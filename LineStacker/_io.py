"""
Cube input/output shared by every LineStacker module.

Cubes can be FITS files or CASA images (read through casatools, CASA 6). Every
cube is returned as an array ordered ``(x, y, stokes, spectral)``, where x is
the longitude (RA) pixel axis and y the latitude (Dec) pixel axis, whatever
the axis order in the file. A :class:`CubeInfo` describes the celestial WCS
and the spectral axis, with a 0-based reference channel.

Arrays and headers are cached per file, so stacking, bootstrapping and the
random-position estimators only read each cube once. A cached cube is read
again when its file changes (for CASA images: any file of the image
directory). Call :func:`clear_cache` to release the memory.
"""

import copy
import os
import re
import tempfile
import threading
import warnings

import numpy as np

# Spectral CTYPE prefixes recognised by the FITS standard (Greisen+2006).
_SPECTRAL_PREFIXES = ('FREQ', 'ENER', 'WAVN', 'VRAD', 'WAVE', 'VOPT', 'ZOPT',
                      'AWAV', 'VELO', 'BETA', 'FELO')
FREQUENCY_TYPES = ('FREQ',)
WAVELENGTH_TYPES = ('WAVE', 'AWAV')
# Units of spectral axes without CUNIT (FITS standard defaults).
_DEFAULT_SPECTRAL_UNITS = {'FREQ': 'Hz', 'ENER': 'J', 'WAVN': '1/m', 'VRAD': 'm/s', 'WAVE': 'm',
                           'VOPT': 'm/s', 'ZOPT': '', 'AWAV': 'm', 'VELO': 'm/s', 'BETA': '',
                           'FELO': 'm/s'}
# Non-linear spectral algorithm codes (Greisen+2006); AIPS frame suffixes
# such as VELO-LSR or FREQ-HEL are linear.
_NONLINEAR_CODES = re.compile(r'^(LOG|TAB|GRI|GRA|[FWVA]2[FWVA])$')

_cache_lock = threading.RLock()
_info_cache = {}
_data_cache = {}


def is_casa_image(path):
    """True if ``path`` is a CASA image (a table directory)."""
    return os.path.isdir(path) and os.path.exists(os.path.join(path, 'table.dat'))


def clear_cache():
    """Forget every cube read so far (frees the memory they use)."""
    with _cache_lock:
        _info_cache.clear()
        _data_cache.clear()


def _signature(path):
    """What changes when an image changes: modification time and size of the
    FITS file, or of every file of a CASA image directory (lock files, which
    CASA touches on every open, excluded)."""
    try:
        if os.path.isdir(path):
            latest, total, count = 0, 0, 0
            for root, _dirs, files in os.walk(path):
                for name in files:
                    if name.endswith('.lock'):
                        continue
                    st = os.stat(os.path.join(root, name))
                    latest = max(latest, st.st_mtime_ns)
                    total += st.st_size
                    count += 1
            return (latest, total, count)
        st = os.stat(path)
        return (st.st_mtime_ns, st.st_size)
    except OSError:
        raise FileNotFoundError('Image not found: {0}'.format(path))


def _key(path, hdul):
    """(absolute path, hdul) and the current signature of the file (None
    when LineStacker.cache_check is off: cached cubes are then trusted)."""
    import LineStacker
    path = os.path.abspath(os.path.expanduser(os.fspath(path)))
    if getattr(LineStacker, 'cache_check', True):
        return (path, hdul), _signature(path)
    if not os.path.exists(path):
        raise FileNotFoundError('Image not found: {0}'.format(path))
    return (path, hdul), None


def _cached(cache, key, sig):
    """Cached value for key if the file has not changed since, else None
    (a changed file's entry is dropped, so old versions do not pile up)."""
    with _cache_lock:
        entry = cache.get(key)
        if entry is None:
            return None
        if sig is None or entry[0] == sig:
            return entry[1]
        del cache[key]
        return None


def _store(cache, key, sig, value):
    if sig is None:                    # cache_check off: record the file as it is now
        sig = _signature(key[0])
    with _cache_lock:
        cache[key] = (sig, value)


def _parse_bunit(text):
    """astropy unit of a BUNIT value ('Jy/beam', 'JY/BEAM', 'mJy/beam', 'K'...),
    or None if empty or not recognised."""
    import astropy.units as u
    if not text:
        return None
    candidates = [text, text.replace('JY', 'Jy').replace('BEAM', 'beam').replace('PIXEL', 'pixel')]
    for cand in candidates:
        for fmt in ('fits', None):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    unit = u.Unit(cand, format=fmt) if fmt else u.Unit(cand)
            except ValueError:
                continue
            if not isinstance(unit, u.UnrecognizedUnit):
                return unit
    return None


def _parse_unit(cunit, base):
    """astropy unit of a spectral axis (FITS default when CUNIT is missing),
    or None if the unit is not recognised."""
    import astropy.units as u
    text = cunit.strip() if cunit else _DEFAULT_SPECTRAL_UNITS.get(base, '')
    if not text:
        return u.dimensionless_unscaled if base in ('ZOPT', 'BETA') else None
    for fmt in ('fits', None):
        try:
            unit = u.Unit(text, format=fmt) if fmt else u.Unit(text)
        except ValueError:
            continue
        if not isinstance(unit, u.UnrecognizedUnit):
            return unit
    return None


class CubeInfo(object):
    """Geometry of one cube: celestial WCS and spectral axis description."""

    def __init__(self, path, header, hdu_index):
        from astropy.wcs import WCS, FITSFixedWarning

        self.path = path
        self.header = header
        self.hdu_index = hdu_index
        naxis = int(header.get('NAXIS', 0))
        if naxis < 2:
            raise ValueError('{0}: HDU {1} holds no image data (NAXIS={2}); '
                             'pass hdul=<index of the data HDU>'.format(path, hdu_index, naxis))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)
            wcs = WCS(header)
        if wcs.wcs.lng < 0 or wcs.wcs.lat < 0:
            raise ValueError('{0}: no celestial (RA/Dec) axes found in the header'.format(path))
        self.wcs = wcs
        self.celestial = wcs.celestial
        # Index (0-based, header order) of every axis we care about.
        self.lng_axis = wcs.wcs.lng
        self.lat_axis = wcs.wcs.lat
        self.spec_axis = -1
        self.stokes_axis = -1
        for i in range(naxis):
            ctype = str(header.get('CTYPE{0}'.format(i + 1), '')).upper()
            if i in (self.lng_axis, self.lat_axis):
                continue
            if ctype.startswith('STOKES'):
                self.stokes_axis = i
            elif ctype.startswith(_SPECTRAL_PREFIXES) and self.spec_axis < 0:
                self.spec_axis = i
        self.naxes = [int(header['NAXIS{0}'.format(i + 1)]) for i in range(naxis)]
        for i, n in enumerate(self.naxes):
            if i not in (self.lng_axis, self.lat_axis, self.spec_axis, self.stokes_axis) and n > 1:
                raise ValueError('{0}: unsupported non-degenerate axis {1} ({2})'.format(
                    path, i + 1, header.get('CTYPE{0}'.format(i + 1))))
        self.nx = self.naxes[self.lng_axis]
        self.ny = self.naxes[self.lat_axis]
        self.nstokes = self.naxes[self.stokes_axis] if self.stokes_axis >= 0 else 1
        self.nchan = self.naxes[self.spec_axis] if self.spec_axis >= 0 else 1

        if self.spec_axis >= 0:
            k = self.spec_axis + 1
            self.spec_ctype = str(header['CTYPE{0}'.format(k)]).upper()
            self.spec_cunit = str(header.get('CUNIT{0}'.format(k), '')).strip()
            self.spec_crval = float(header.get('CRVAL{0}'.format(k), 0.))
            # FITS CRPIX is 1-based; LineStacker indexes arrays from 0.
            self.spec_crpix0 = float(header.get('CRPIX{0}'.format(k), 1.)) - 1.
            if 'CD{0}_{0}'.format(k) in header:
                self.spec_cdelt = float(header['CD{0}_{0}'.format(k)])
            else:
                self.spec_cdelt = (float(header.get('CDELT{0}'.format(k), 1.))
                                   * float(header.get('PC{0}_{0}'.format(k), 1.)))
        else:
            self.spec_ctype, self.spec_cunit = '', ''
            self.spec_crval, self.spec_crpix0, self.spec_cdelt = 0., 0., 1.
        # 'FREQ', 'WAVE', 'VRAD'... and whether values are linear in channel
        self.spec_base = self.spec_ctype[:4]
        code = self.spec_ctype[5:].strip('-') if len(self.spec_ctype) > 5 else ''
        # AIPS 'FELO' is an optical velocity on a frequency-linear grid: not
        # linear in velocity
        self.spec_linear = not _NONLINEAR_CODES.match(code) and self.spec_base != 'FELO'
        self.spec_unit = _parse_unit(self.spec_cunit, self.spec_base) if self.spec_axis >= 0 else None
        self.specsys = str(header['SPECSYS']).strip() if 'SPECSYS' in header else None

        from astropy.wcs.utils import proj_plane_pixel_scales
        cel = self.celestial
        order = [cel.wcs.lng, cel.wcs.lat]
        #: (2, 2) linear transformation from pixel offsets (x along RA, y along
        #: Dec) to projection-plane offsets (degrees, RA then Dec); includes
        #: rotation, skew and the sign of the axes.
        self.cel_matrix = np.array(cel.pixel_scale_matrix)[np.ix_(order, order)]
        self.pixscale_deg = float(np.mean(proj_plane_pixel_scales(cel)))
        #: intensity unit (BUNIT) as written, and as an astropy unit (None if unknown)
        self.bunit = str(header.get('BUNIT', '')).strip()
        self.bunit_unit = _parse_bunit(self.bunit)
        self.beam_major_deg = float(header['BMAJ']) if 'BMAJ' in header else None
        self.beam_minor_deg = float(header['BMIN']) if 'BMIN' in header else None
        self.beam_pa_deg = float(header['BPA']) if 'BPA' in header else None

    @property
    def spectral_kind(self):
        """'frequency', 'wavelength', 'other' or None (no spectral axis)."""
        if self.spec_axis < 0:
            return None
        if self.spec_ctype.startswith(FREQUENCY_TYPES):
            return 'frequency'
        if self.spec_ctype.startswith(WAVELENGTH_TYPES):
            return 'wavelength'
        return 'other'

    def spectral_in(self, unit):
        """This CubeInfo with spectral values (spec_crval, spec_cdelt, and so
        channel_of/value_of) in another unit (an astropy unit). Returns self
        if nothing changes; the cached original is never modified."""
        if unit is None or self.spec_unit is None or self.spec_unit == unit:
            return self
        factor = self.spec_unit.to(unit)
        new = copy.copy(self)
        new.spec_crval = self.spec_crval * factor
        new.spec_cdelt = self.spec_cdelt * factor
        new.spec_unit = unit
        new.spec_cunit = unit.to_string('fits') if unit.to_string('fits') else unit.to_string()
        return new

    def channel_of(self, value):
        """Fractional 0-based channel of a spectral world value (header units)."""
        return self.spec_crpix0 + (value - self.spec_crval) / self.spec_cdelt

    def value_of(self, channel):
        """Spectral world value (header units) of a 0-based channel."""
        return self.spec_crval + (channel - self.spec_crpix0) * self.spec_cdelt

    def world_to_pixel(self, skycoord):
        """0-based (x, y) pixels of an astropy SkyCoord (scalar or array);
        x runs along the RA axis and y along the Dec axis of the file."""
        from astropy.wcs.utils import wcs_to_celestial_frame
        cel = self.celestial
        sc = skycoord.transform_to(wcs_to_celestial_frame(cel))
        lon = np.atleast_1d(sc.spherical.lon.deg)
        lat = np.atleast_1d(sc.spherical.lat.deg)
        world = np.zeros((lon.size, 2))
        world[:, cel.wcs.lng] = lon
        world[:, cel.wcs.lat] = lat
        pix = cel.all_world2pix(world, 0)
        x, y = pix[:, cel.wcs.lng], pix[:, cel.wcs.lat]
        if np.ndim(skycoord.ra) == 0:
            return float(x[0]), float(y[0])
        return x, y

    def pixel_to_world(self, x, y):
        """Inverse of world_to_pixel; returns a SkyCoord in the image frame."""
        from astropy.wcs.utils import wcs_to_celestial_frame
        from astropy.coordinates import SkyCoord
        cel = self.celestial
        x, y = np.atleast_1d(x), np.atleast_1d(y)
        pix = np.zeros((x.size, 2))
        pix[:, cel.wcs.lng], pix[:, cel.wcs.lat] = x, y
        world = cel.all_pix2world(pix, 0)
        return SkyCoord(world[:, cel.wcs.lng], world[:, cel.wcs.lat], unit='deg',
                        frame=wcs_to_celestial_frame(cel))


def _pick_hdu(hdulist, hdul):
    if hdul in (None, 'auto'):
        for i, hdu in enumerate(hdulist):
            if getattr(hdu, 'is_image', False) and hdu.header.get('NAXIS', 0) >= 2:
                return i
        raise ValueError('no image HDU found in {0}'.format(hdulist.filename()))
    hdul = int(hdul)
    if hdul >= len(hdulist):
        raise IndexError('{0} has {1} HDU(s); hdul={2} does not exist. Use hdul="auto" '
                         'or the index of the data HDU.'.format(hdulist.filename(), len(hdulist), hdul))
    return hdul


def _casa_to_fits(path):
    """Export a CASA image to a temporary FITS file and return its path.
    CASA's own exporter writes the WCS, beam and units consistently."""
    from LineStacker._casa_compat import new_tool
    fd, tmp = tempfile.mkstemp(suffix='.fits', prefix='linestacker_')
    os.close(fd)
    ia = new_tool('image')
    try:
        if not ia.open(path):
            raise IOError('CASA could not open {0}'.format(path))
        ia.tofits(outfile=tmp, velocity=False, optical=False, stokeslast=True,
                  overwrite=True, history=False)
    finally:
        ia.done()
    return tmp


def _casa_header(path):
    """FITS header of a CASA image without exporting the whole cube: export
    one spatial plane, then restore the axis lengths."""
    from astropy.io import fits
    from LineStacker._casa_compat import new_tool
    fd, tmp = tempfile.mkstemp(suffix='.fits', prefix='linestacker_')
    os.close(fd)
    ia, rg, cs = new_tool('image'), new_tool('regionmanager'), None
    try:
        if not ia.open(path):
            raise IOError('CASA could not open {0}'.format(path))
        shape = [int(n) for n in ia.shape()]
        cs = ia.coordsys()
        direction = list(cs.findcoordinate('direction')['pixel'])
        stokes = list(cs.findcoordinate('stokes')['pixel']) if 'Stokes' in cs.axiscoordinatetypes() else []
        trc = [shape[i] - 1 if i in direction else 0 for i in range(len(shape))]
        ia.tofits(outfile=tmp, velocity=False, optical=False, stokeslast=True, overwrite=True,
                  history=False, region=rg.box(blc=[0] * len(shape), trc=trc))
        header = fits.getheader(tmp, 0)
    finally:
        if cs is not None:
            cs.done()
        rg.done()
        ia.done()
        os.remove(tmp)
    # tofits(stokeslast=True) moves the Stokes axis to the end
    order = [i for i in range(len(shape)) if i not in stokes] + stokes
    for k, i in enumerate(order):
        header['NAXIS{0}'.format(k + 1)] = shape[i]
    return header


def _read_fits(path, hdul, want_data):
    from astropy.io import fits
    with fits.open(path, memmap=False) as hl:
        idx = _pick_hdu(hl, hdul)
        header = hl[idx].header.copy()
        data = np.array(hl[idx].data) if want_data else None
    return idx, header, data


def _load(path, hdul, want_data):
    """Read header (and data) of a FITS file or CASA image."""
    if is_casa_image(path):
        tmp = _casa_to_fits(path)
        try:
            idx, header, data = _read_fits(tmp, 0, True)
        finally:
            os.remove(tmp)
        return idx, header, data
    if os.path.isdir(path):
        raise IOError('{0} is a directory but not a CASA image'.format(path))
    return _read_fits(path, hdul, want_data)


def _normalise(data, info):
    """Reorder a FITS-ordered array to (x, y, stokes, spectral). Blank
    pixels (NaN in FITS, masked pixels of CASA images) stay NaN."""
    naxis = len(info.naxes)
    # numpy axis of header axis i (0-based) is naxis-1-i.
    np_ax = lambda i: naxis - 1 - i
    order = [np_ax(info.lng_axis), np_ax(info.lat_axis)]
    extra = []
    for name in ('stokes_axis', 'spec_axis'):
        i = getattr(info, name)
        if i >= 0:
            order.append(np_ax(i))
        else:
            extra.append(name)
    leftover = [a for a in range(data.ndim) if a not in order]
    arr = np.transpose(data, order + leftover)
    if leftover:
        arr = arr.reshape(arr.shape[:len(order)])
    if info.stokes_axis < 0:
        arr = arr[:, :, np.newaxis] if info.spec_axis >= 0 else arr[:, :, np.newaxis, np.newaxis]
    if info.spec_axis < 0 and info.stokes_axis >= 0:
        arr = arr[:, :, :, np.newaxis]
    if arr.dtype.kind != 'f':
        arr = arr.astype(np.float32)
    arr = np.ascontiguousarray(arr, dtype=arr.dtype.newbyteorder('='))
    bad = np.isinf(arr)
    if bad.any():
        arr[bad] = np.nan              # infinities are not data either
    return arr


COLLAPSE_METHODS = ('mean', 'median', 'sum')


def _collapse(arr, method):
    """Collapse the spectral axis of an (x, y, stokes, spectral) array,
    ignoring blanked (NaN) channels; keeps a spectral axis of length 1.
    Pixels blank in every channel stay NaN."""
    import warnings
    if method not in COLLAPSE_METHODS:
        raise ValueError('collapse must be one of {0}'.format(COLLAPSE_METHODS))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)       # all-NaN pixels
        func = {'mean': np.nanmean, 'median': np.nanmedian, 'sum': np.nansum}[method]
        out = func(arr, axis=3, keepdims=True)
    out[~np.isfinite(arr).any(axis=3, keepdims=True)] = np.nan
    return np.ascontiguousarray(out, dtype=arr.dtype)


def read_info(path, hdul='auto'):
    """CubeInfo for a FITS file or CASA image (cached)."""
    key, sig = _key(path, hdul)
    info = _cached(_info_cache, key, sig)
    if info is not None:
        return info
    if is_casa_image(key[0]):
        idx, header = 0, _casa_header(key[0])
    else:
        idx, header, _ = _load(key[0], hdul, False)
    info = CubeInfo(key[0], header, idx)
    _store(_info_cache, key, sig, info)
    return info


def n_channels(path, hdul='auto'):
    """Number of spectral channels of an image (1 without a spectral axis),
    without exporting CASA images."""
    key, sig = _key(path, hdul)
    info = _cached(_info_cache, key, sig)
    if info is not None:
        return info.nchan
    if not is_casa_image(key[0]):
        return read_info(path, hdul).nchan
    from LineStacker._casa_compat import new_tool
    ia, cs = new_tool('image'), None
    try:
        if not ia.open(key[0]):
            raise IOError('CASA could not open {0}'.format(key[0]))
        shape = [int(n) for n in ia.shape()]
        cs = ia.coordsys()
        if 'Spectral' not in cs.axiscoordinatetypes():
            return 1
        return shape[int(cs.findcoordinate('spectral')['pixel'][0])]
    finally:
        if cs is not None:
            cs.done()
        ia.done()


def read_cube(path, hdul='auto', collapse=None):
    """Return ``(array, CubeInfo)``, array ordered (x, y, stokes, spectral).

    collapse=None returns the cube as stored. collapse='sum' (or 'mean',
    'median') returns the cube collapsed along its spectral axis (blanked
    channels ignored), shape (x, y, stokes, 1). Collapsed cubes are cached on
    their own, so the full cube does not need to stay in memory. A cube is
    read again if its file changed since it was cached."""
    key, sig = _key(path, hdul)
    ckey = key + (collapse,)
    arr = _cached(_data_cache, ckey, sig)
    info = _cached(_info_cache, key, sig)
    if arr is not None and info is not None:
        return arr, info
    idx, header, data = _load(key[0], hdul, True)
    info = CubeInfo(key[0], header, idx)
    if collapse is None:
        arr = _normalise(data, info)
    else:
        arr = _collapse(_normalise(data, info), collapse)
    del data
    _store(_info_cache, key, sig, info)
    _store(_data_cache, ckey, sig, arr)
    return arr, info


#: Names written as FITS by write_cube (anything else becomes a CASA image).
FITS_SUFFIXES = ('.fits', '.fit', '.fts', '.fits.gz', '.fit.gz', '.fts.gz')


def write_cube(outfile, pixels, header, overwrite=True):
    """Write an (x, y, stokes, spectral) array. A name ending in .fits/.fit/
    .fts (optionally .gz) is written as FITS; any other name (e.g. 'stack.image')
    becomes a CASA image (needs casatools)."""
    from astropy.io import fits
    arr = np.asarray(pixels, dtype=np.float32)
    nx, ny, ns, nc = arr.shape
    lower = str(outfile).lower()
    is_fits = lower.endswith(FITS_SUFFIXES)
    has_spec = 'CTYPE3' in header
    if not has_spec:
        # 2D image (continuum stack of images without a spectral axis)
        out = arr[:, :, 0, 0].T
    elif is_fits:
        # FITS axis order RA, Dec, spectral(, Stokes); numpy order is reversed.
        out = arr[:, :, 0, :].transpose(2, 1, 0) if ns == 1 else arr.transpose(2, 3, 1, 0)
    if is_fits:
        fits.writeto(outfile, out, header, overwrite=overwrite)
        return outfile
    if has_spec:
        # CASA image: same axis order as LineStacker < 2.0 wrote (RA, Dec,
        # Stokes, spectral), so ia.getchunk() returns (x, y, stokes, channel).
        header = header.copy()
        for key in ('CTYPE', 'CUNIT', 'CRPIX', 'CRVAL', 'CDELT'):
            if key + '3' in header:
                header[key + '4'] = header[key + '3']
        header['CTYPE3'], header['CRPIX3'], header['CRVAL3'], header['CDELT3'] = 'STOKES', 1., 1., 1.
        if 'CUNIT3' in header:
            header['CUNIT3'] = ''
        out = arr.transpose(3, 2, 1, 0)
    from LineStacker._casa_compat import new_tool
    import shutil
    fd, tmp = tempfile.mkstemp(suffix='.fits', prefix='linestacker_')
    os.close(fd)
    try:
        fits.writeto(tmp, out, header, overwrite=True)
        if os.path.exists(outfile) and overwrite:
            shutil.rmtree(outfile) if os.path.isdir(outfile) else os.remove(outfile)
        ia = new_tool('image')
        try:
            ia.fromfits(outfile=outfile, infile=tmp, overwrite=overwrite)
        finally:
            ia.done()
    finally:
        os.remove(tmp)
    return outfile

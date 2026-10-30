"""
Stacking in the image domain: spectral lines in data cubes, or continuum.

Inputs are FITS files or CASA images (CASA images need casatools), and can be
mixed. Line stacking extracts a sub-cube around each source, centred on its
line; continuum stacking (``continuum=True``) stacks images. Cubes without a
spectral axis (or with a single channel) are always stacked as continuum;
cubes with a spectral axis are collapsed in memory when continuum stacking is
asked for (with a warning).

What the stacked cube is:

- Stamps are cut on each cube's own pixel grid, centred on the source; cubes
  are not reprojected. Cubes whose pixels differ in size or orientation, or
  whose projection, celestial frame or distortions differ, from the first
  cube's give a warning: reproject them first (e.g. CASA imregrid) to use
  the stacked image spatially.
- The output sky axes are offsets from the stacked position, (0, 0) at the
  stamp centre, with the first cube's pixel size and orientation.
- All cubes must have the same kind of spectral axis (FREQ, WAVE, AWAV, VRAD...);
  units are converted to those of the first cube, in which fEm, lEm and
  observed spectral values are given. Non-linear spectral axes are refused
  for line stacks. Different spectral reference frames (SPECSYS) give a
  warning: LineStacker does not convert frames.
- The beam in the header (BMAJ, BMIN, BPA) is the mean of the input beams, an
  approximation of the stack's effective beam when they differ.
- Blank pixels (NaN) stay NaN and are left out pixel by pixel; stacked
  pixels without any valid data are NaN.
- Pixel units (BUNIT) are compared: scaled units are converted to the first
  cube's (mJy/beam to Jy/beam), recognised but physically different ones
  are refused, unrecognised ones give a warning.
- A continuum stack of collapsed cubes keeps BUNIT, which is then per
  channel ('sum' adds channels); the COLLAPSE keyword says how.
"""
# LineStacker, Python module for stacking of interferometric data.
# Copyright (C) 2014  Lukas Lindroos
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

import os
import warnings

import numpy as np

import LineStacker
from LineStacker import _io


class Prepared(object):
    """Stamps of every source, ready to be combined.

    data      (nsrc, stampsize, stampsize, 1, N_chans) stamps
    valid     (nsrc, N_chans) True where the channel lies inside the cube
    coords    working copies of the pixel coordinates (weights included)
    infos     CubeInfo of each source's cube
    centres   fractional central channel of each source in its cube
    steps     cube channels per stack channel for each source (sign included)
    scales    factor converting each source's pixels to the output BUNIT
    bunit     intensity unit of the stack (None if no cube gives one)
    continuum True for a continuum stack (N_chans == 1)
    collapse  per image: None (used as stored) or how it was collapsed
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def weights(self):
        """(nsrc, N_chans) weight matrix, zero outside the cubes. Weights must
        be finite and >= 0."""
        w = np.empty(self.valid.shape)
        for i, c in enumerate(self.coords):
            w[i] = np.broadcast_to(np.asarray(c.weight, dtype=float), (self.valid.shape[1],))
        check_weights(w, [c.index if c.index is not None else i for i, c in enumerate(self.coords)])
        return w * self.valid


def check_weights(w, labels=None):
    """Raise ValueError if a weight is negative or not finite (rows = sources)."""
    w = np.asarray(w, dtype=float)
    bad = ~np.isfinite(w) | (w < 0)
    if bad.any():
        row = int(np.where(bad.any(axis=-1) if w.ndim > 1 else bad)[0][0])
        value = w[row][bad[row]][0] if w.ndim > 1 else w[row]
        label = labels[row] if labels is not None else row
        raise ValueError('source {0} has weight {1}: weights must be finite and >= 0 '
                         '(set a weight to 0 to leave a source or channel out)'.format(label, value))


def stack(  coords,
            outfile='stackResult.fits',
            stampsize = 32,
            imagenames=None,
            method = 'mean',
            spectralMethod='z',
            weighting = None,
            returnWeights=False,
            maskradius=0,
            fEm = -1,
            lEm = -1,
            N_chans='default',
            plotIt=False,
            regridFromZ=False,
            regridMethod='scaleToMax',
            saveSubCubes=False,
            hdul='auto',
            continuum=False,
            collapse='sum',
            primarybeam=None,
            verbose=None):
    """
        Performs line (or continuum) stacking in the image domain.

        Parameters
        ---------
        coords
            A list of LineStacker.Coord (see LineStacker.readCatalog, LineStacker.makeCoords, LineStacker.readCoords). Set to None to stack the centre of every cube.
        outfile
            Target name for stacked image. Names ending in .fits are written as FITS, other names (e.g. 'stackResult.image') as CASA images (needs casatools). None to skip writing.
        stampsize
            size of target image in pixels
        imagenames
            Cubes or images to extract flux from: FITS files or CASA images (both can be mixed), in the same order as the coordinate files.
            Can be omitted when every coordinate names its cube (readCatalog with an image column, makeCoords(image=...)).
            All cubes must have the same kind of spectral axis (e.g. all FREQ); different units (Hz, GHz) are converted to those of the first cube.
            Stamps are cut on each cube's own pixel grid: cubes with other pixel sizes or orientations are not reprojected (a warning says so).
        method
            'mean' (weighted) or 'median' (over the sources covering each channel; weights are ignored)
        spectralMethod
            What the spectral value of each coordinate (its z, the 3rd column of a coordinate file) is. Possible methods are:\n
            **'z'**: the redshift of the observed line, if used
            **fEm** (emission frequency) or **lEm** (emission wavelength) must be informed as an argument of this Stack function.\n
            **'obsSpec'**: the (observed) central frequency, wavelength or velocity, in the spectral unit of the first cube.\n
            **'channel'**: the channel number (0-based) of the line.
        weighting
            used if method set to **'mean'**, possible values are:\n
            **'sigma2'**: weight times 1/sigma**2, sigma being the standard deviation of the source's stamp (pixels outside **maskradius**, blank pixels excluded)\n
            **'sigma2F'**: the same, computed channel by channel\n
            **'1/A'**, **/!\\\\ only use if lines are visible pre-stacking**: weight divided by the line amplitude (maximum of the central-pixel spectrum; sources with a non-positive amplitude get weight 0)\n
            **'sigma'**: weight times 1/sigma (sigma as for 'sigma2'; if maskradius is 0 the mask radius is twice the beam major axis)\n
            **'pb'**: weight times the squared primary-beam response at each source (needs **primarybeam**); for primary-beam corrected images this is inverse-variance weighting\n
            **None**: the weights of the coordinates (1 by default).\n
            Each weighting multiplies the coordinates' own weights. Weights must be finite and >= 0.
        returnWeights
            If true weights of each source are written in outfile "weightFile.csv" in format: coord.ra, coord.dec, coord.z, coord.weight.\n
            If multiple weights exist for one source (for each spectral channel) all are written.
        maskradius
            allows blanking of centre pixels in weight calculation
        fEm
            rest emission frequency of the line, should be used if spectralMethod = 'z' and spectral axis of image cube is in frequency (in the unit of the first cube, usually Hz).
        lEm
            rest emission wavelength of the line, should be used if spectralMethod = 'z' and spectral axis of image cube is in wavelength (in the unit of the first cube, e.g. Angstrom for MUSE).
        N_chans
            number of channels of the resulting stack, default is number of channels in the first image.
        plotIt
            direct plot option
        regridFromZ
            if set to True each source's spectrum is resampled (linear interpolation) so that all sources share the same rest-frame channel width, correcting the redshift dependence of observed line widths. Requires spectralMethod='z'.
        regridMethod
            Used if regridFromZ is True. Can be set either to 'scaleToMin' or 'scaleToMax'.
            The common rest-frame channel is the one of the source with the smallest ('scaleToMin') or highest ('scaleToMax') redshift.
        saveSubCubes
            If set to True (or to str) sub cubes will be saved as a numpy file. (outSubCubes.npy if set to True, user defined is set to str)
        hdul
            index of the FITS HDU holding the cube; 'auto' (default) picks the first HDU containing an image.
        continuum
            If True, stack continuum images instead of lines (spectralMethod, fEm, lEm, N_chans and regridFromZ are then not used).
            Continuum stacking is also used automatically when no image has a spectral axis (or only one channel),
            and when the coordinates carry no spectral value (readCoords(..., lineON=False)).
            Cubes with a spectral axis are then collapsed in memory (with a warning), see **collapse**.
        collapse
            How cubes are collapsed for continuum stacking: 'sum' (default), 'mean' or 'median' over channels. Blanked channels are ignored.
        primarybeam
            Used by weighting='pb': the primary-beam image of each cube (FITS or CASA, e.g. tclean's .pb), as one path for a single cube,
            a list in the order of imagenames or a dict {cube: pb image}; the response is read at each source's position (and line channel).
            Or a function f(dx, dy, nu) returning the response at offsets dx, dy (radians, east and north on the sky) from the image reference pixel
            and frequency nu (Hz, None without a frequency axis), e.g. **GaussianPrimaryBeam(dish=12.)**. Entries of a list or dict can also be functions.
        verbose
            If True, print the path of every file written (stacked cube, weights, sub-cubes).
            Default: LineStacker.verbose (True unless changed).

        Returns
        -------
        list
            [central pixel at the central channel, spectrum of the central pixel,
            stacked cube (x, y, stokes, channel), number of sources missing from
            each channel]. For a continuum stack the spectrum has one channel.
    """
    prep = _prepare(coords, imagenames, stampsize=stampsize, N_chans=N_chans,
                    spectralMethod=spectralMethod, fEm=fEm, lEm=lEm, hdul=hdul,
                    regridFromZ=regridFromZ, regridMethod=regridMethod,
                    method=method, weighting=weighting, maskradius=maskradius,
                    continuum=continuum, collapse=collapse, primarybeam=primarybeam)
    nchan = prep.valid.shape[1]

    if returnWeights:
        _write_weights('weightFile.csv', coords, prep)
        LineStacker._info('weights written to ' + os.path.abspath('weightFile.csv'), verbose)

    stacked_im = _combine(prep.data, prep.weights(), prep.valid, method)

    s = prep.stampsize
    what = '{0} {1} stack of {2} source{3} ({4}x{4} pixels{5})'.format(
        method, 'continuum' if prep.continuum else 'line', len(prep.coords),
        's' if len(prep.coords) != 1 else '', s, '' if prep.continuum else ', {0} channels'.format(nchan))
    if outfile is not None:
        _write_stacked_image(outfile, stacked_im, prep, spectralMethod=spectralMethod,
                             fEm=fEm, lEm=lEm, method=method)
        LineStacker._info('{0} written to {1}'.format(what, os.path.abspath(outfile)), verbose)
    else:
        LineStacker._info(what + ' done (outfile=None: not written to disk)', verbose)

    if saveSubCubes is not False:
        base = saveSubCubes if isinstance(saveSubCubes, str) else None
        npy = (base or 'outSubCubes') + '.npy'
        txt = (base + '_outCoordsFile.txt') if base else 'outCoordsFile.txt'
        np.save(npy, prep.data)
        with open(txt, 'w') as fh:
            for coord in prep.coords:
                fh.write(str(coord.ra)+', '+str(coord.dec)+', '+str(coord.spec)+'\n')
        LineStacker._info('sub-cubes written to {0}, their pixel positions to {1}'.format(
            os.path.abspath(npy), os.path.abspath(txt)), verbose)

    centre = int(prep.data.shape[1] / 2)
    missing = (~prep.valid).sum(axis=0).astype(float)      # sources without each channel
    if plotIt and prep.continuum:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        im = ax.imshow(stacked_im[:, :, 0, 0].T, origin='lower')
        fig.colorbar(im, ax=ax)
        ax.set_title('Continuum stack of {0} sources'.format(len(prep.coords)))
        fig.show()
    elif plotIt:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(11, 4))
        axes[0].step(range(nchan), stacked_im[centre, centre, 0, :], color='red', where='mid')
        axes[0].set_title('Stack flux in the central pixel')
        axes[0].set_xlabel('Channel number')
        axes[1].step(range(nchan), len(prep.coords) - missing, color='blue', where='mid')
        axes[1].set_title('Number of sources per channel')
        axes[1].set_xlabel('Channel number')
        fig.tight_layout()
        fig.show()

    return [stacked_im[centre, centre, 0, int(nchan/2)], stacked_im[centre, centre, 0, :],
            stacked_im, missing]


def _centre_coords(imagenames, hdul):
    """One pixel coordinate at the centre (space and spectrum) of each cube."""
    out = []
    for k, name in enumerate(imagenames):
        inf = _io.read_info(name, hdul)
        out.append(LineStacker.Coord(inf.nx // 2, inf.ny // 2, spec=inf.nchan // 2, image=name,
                                     imagenumber=k, type='pixel'))
    return out


def _centre_channel(coord, info, spectralMethod, fEm, lEm):
    """Fractional 0-based channel of a source's line in its cube."""
    if spectralMethod == 'channel':
        return float(coord.spec)
    if spectralMethod == 'obsSpec':
        return info.channel_of(float(coord.spec))
    if spectralMethod == 'z':
        kind = info.spectral_kind
        z = float(coord.spec)
        if kind == 'frequency':
            if fEm == -1:
                raise ValueError('fEm should be specified if spectralMethod is z and the spectral axis '
                                 'of {0} is in frequency'.format(info.path))
            return info.channel_of(fEm / (1. + z))
        if kind == 'wavelength':
            if lEm == -1:
                raise ValueError('lEm should be specified if spectralMethod is z and the spectral axis '
                                 'of {0} is in wavelength'.format(info.path))
            return info.channel_of(lEm * (1. + z))
        raise ValueError("spectralMethod='z' needs a frequency or wavelength axis; {0} has {1!r}. "
                         "Use spectralMethod='obsSpec' or 'channel'.".format(info.path, info.spec_ctype))
    raise ValueError("spectralMethod must be 'z', 'obsSpec' or 'channel'")


def _rest_channel_width(info, z):
    """Rest-frame width of one channel of a cube for a source at redshift z."""
    if info.spectral_kind == 'frequency':
        return abs(info.spec_cdelt) * (1. + z)
    if info.spectral_kind == 'wavelength':
        return abs(info.spec_cdelt) / (1. + z)
    raise ValueError('regridFromZ needs a frequency or wavelength axis ({0})'.format(info.path))


def _extract(cube, info, x, y, centre, step, stampsize, nchan_out):
    """Stamp of stampsize x stampsize pixels around (x, y), sampled at channels
    centre + o*step (o = -nchan_out//2 ...). Integer positions are copied,
    fractional ones linearly interpolated. Returns (stamp, inside_cube)."""
    s = stampsize
    blcx = int(np.floor(x - int(s / 2.) + 0.5))
    blcy = int(np.floor(y - int(s / 2.) + 0.5))
    if blcx < 0 or blcy < 0 or blcx + s > info.nx or blcy + s > info.ny:
        raise Exception('the source located at '+str(x)+' '+str(y)+' on image '+str(info.path) +
                        ' is too close to the edge, trying to stack out of boundaries, try stacking with smaller stampsize')
    sub = cube[blcx:blcx + s, blcy:blcy + s, 0, :]
    nchan = sub.shape[-1]
    pos = centre + (np.arange(nchan_out) - int(nchan_out / 2.)) * step
    inside = (pos > -1e-6) & (pos < nchan - 1 + 1e-6)
    pos = np.clip(pos, 0, nchan - 1)
    j0 = np.floor(pos + 1e-9).astype(int)
    j1 = np.minimum(j0 + 1, nchan - 1)
    t = np.clip(pos - j0, 0., 1.)
    out = np.full((s, s, 1, nchan_out), np.nan)
    out[:, :, 0, :] = _interp_channels(sub[..., j0], sub[..., j1], t, inside)
    return out, inside


def _interp_channels(a, b, t, inside):
    """a*(1-t) + b*t along the last axis, NaN outside the cube. Samples that
    fall on a channel (t = 0 or 1) copy it, so a blank neighbour does not
    blank them (NaN*0 is NaN)."""
    out = np.full(np.broadcast(a, b).shape, np.nan)
    left = inside & (t == 0.)
    right = inside & (t == 1.)
    mid = inside & ~(left | right)
    out[..., left] = a[..., left]
    out[..., right] = b[..., right]
    out[..., mid] = a[..., mid] * (1. - t[mid]) + b[..., mid] * t[mid]
    return out


def _extract_many(cube, info, xs, ys, centre, step, stampsize, nchan_out):
    """_extract for many positions of one source at once (same spectral
    sampling). Returns stamps (k, s, s, 1, nchan_out); positions must be
    inside the image."""
    s = stampsize
    bx = np.floor(np.asarray(xs) - int(s / 2.) + 0.5).astype(int)
    by = np.floor(np.asarray(ys) - int(s / 2.) + 0.5).astype(int)
    nchan = cube.shape[3]
    pos = centre + (np.arange(nchan_out) - int(nchan_out / 2.)) * step
    inside = (pos > -1e-6) & (pos < nchan - 1 + 1e-6)
    pos = np.clip(pos, 0, nchan - 1)
    j0 = np.floor(pos + 1e-9).astype(int)
    j1 = np.minimum(j0 + 1, nchan - 1)
    t = np.clip(pos - j0, 0., 1.)
    off = np.arange(s)
    ix = (bx[:, None] + off[None, :])[:, :, None]
    iy = (by[:, None] + off[None, :])[:, None, :]
    plane = cube[:, :, 0, :]
    stamps = plane[ix, iy]
    out = _interp_channels(stamps[..., j0], stamps[..., j1], t, inside)
    return out[:, :, :, None, :]


def _continuum_mode(coords, imagenames, nchans, continuum):
    """Decide between line and continuum stacking; returns (continuum, reason)."""
    if continuum:
        return True, 'continuum=True'
    if all(n <= 1 for n in nchans):
        return True, None                     # nothing spectral to stack on
    if coords is not None and len(coords) and all(isinstance(c.spec, str) for c in coords):
        return True, 'the coordinates have no spectral value'
    flat = [name for name, n in zip(imagenames, nchans) if n <= 1]
    if flat:
        raise ValueError('line stacking needs a spectral axis, but {0} image(s) have none ({1}{2}). '
                         'Use continuum=True to stack everything as continuum.'.format(
                             len(flat), ', '.join(flat[:3]), '...' if len(flat) > 3 else ''))
    return False, None


def _cubes(prep, hdul='auto'):
    """Arrays the stamps of a Prepared stack were cut from (cached)."""
    return [_io.read_cube(name, hdul, collapse=col)[0] for name, col in zip(prep.imagenames, prep.collapse)]


def _prepare(coords, imagenames, stampsize=32, N_chans='default', spectralMethod='z',
             fEm=-1, lEm=-1, hdul='auto', regridFromZ=False, regridMethod='scaleToMax',
             method='mean', weighting=None, maskradius=0, continuum=False, collapse='sum',
             primarybeam=None):
    """Read cubes, extract every stamp and compute weights. The user's
    coordinates are never modified."""
    imagenames = LineStacker._cube_list(coords, imagenames)
    nchans = [_io.n_channels(name, hdul) for name in imagenames]
    continuum, reason = _continuum_mode(coords, imagenames, nchans, continuum)
    if continuum and regridFromZ:
        raise ValueError('regridFromZ does not apply to continuum stacking')
    if coords is None:
        coords = _centre_coords(imagenames, hdul)
        spectralMethod = 'channel'
    if len(coords) == 0:
        raise ValueError('no coordinates to stack')
    if regridFromZ and spectralMethod != 'z':
        raise ValueError("regridFromZ needs spectralMethod='z'")

    if continuum:
        if collapse not in _io.COLLAPSE_METHODS:
            raise ValueError('collapse must be one of {0}'.format(_io.COLLAPSE_METHODS))
        collapsing = [name for name, n in zip(imagenames, nchans) if n > 1]
        if collapsing:
            warnings.warn('continuum stacking ({0}): collapsing {1} cube(s) with a spectral axis '
                          '({2} over channels) in memory'.format(reason or 'continuum=True', len(collapsing), collapse),
                          stacklevel=3)
        collapses = [collapse if n > 1 else None for n in nchans]
    else:
        collapses = [None] * len(imagenames)
    cubes, infos = [], []
    for name, col in zip(imagenames, collapses):
        cube, inf = _io.read_cube(name, hdul, collapse=col)
        cubes.append(cube)
        infos.append(inf)
    infos = _common_spectral_axis(infos, continuum, spectralMethod)
    _check_grids(infos)
    image_scales, bunit = _intensity_scales(infos)
    work = LineStacker.getPixelCoords(coords, imagenames, hdul=hdul)
    if continuum:
        N_chans = 1
    elif N_chans == 'default':
        N_chans = infos[0].nchan
    N_chans = int(N_chans)

    src_infos = [infos[c.imagenumber] for c in work]
    if continuum:
        centres = np.zeros(len(work))
    else:
        centres = np.array([_centre_channel(c, inf, spectralMethod, fEm, lEm)
                            for c, inf in zip(work, src_infos)])
    # Cubes whose spectral axis runs the other way are read backwards so all
    # spectra share the orientation of the first source's cube.
    ref_sign = np.sign(src_infos[0].spec_cdelt) or 1.
    orient = np.array([(np.sign(inf.spec_cdelt) or 1.) * ref_sign for inf in src_infos])
    if continuum:
        steps = np.ones(len(work))
    elif regridFromZ:
        zs = np.array([float(c.spec) for c in work])
        ref = int(np.argmin(zs)) if regridMethod == 'scaleToMin' else int(np.argmax(zs))
        target = _rest_channel_width(src_infos[ref], zs[ref])
        steps = orient * np.array([target / _rest_channel_width(inf, z) for inf, z in zip(src_infos, zs)])
    else:
        # Classic behaviour: nearest channel, one cube channel per stack channel.
        centres = np.array([LineStacker.round_half_up(cn) for cn in centres], dtype=float)
        steps = orient

    n = len(work)
    stamps = np.zeros((n, stampsize, stampsize, 1, N_chans))
    inside = np.zeros((n, N_chans), dtype=bool)
    for i, c in enumerate(work):
        stamps[i], inside[i] = _extract(cubes[c.imagenumber], src_infos[i], c.ra, c.dec,
                                        centres[i], steps[i], stampsize, N_chans)
    scales = np.array([image_scales[c.imagenumber] for c in work])
    if np.any(scales != 1):
        stamps *= scales[:, None, None, None, None]

    prep = Prepared(data=stamps, valid=inside, coords=work, infos=src_infos, centres=centres,
                    steps=steps, stampsize=stampsize, regridFromZ=regridFromZ,
                    imagenames=imagenames, spectralMethod=spectralMethod,
                    continuum=continuum, collapse=collapses, scales=scales, bunit=bunit)
    if method == 'mean' and weighting is not None:
        _apply_weighting(prep, weighting, maskradius, primarybeam, hdul)
    return prep


def _common_spectral_axis(infos, continuum, spectralMethod):
    """Check that the cubes' spectral axes can be stacked together and
    express them all in the unit of the first cube (copies; the cached
    CubeInfo objects are not changed)."""
    spec = [i for i in infos if i.spec_axis >= 0 and i.nchan > 1]
    if not spec:
        return list(infos)
    ref = spec[0]
    bases = sorted(set(i.spec_base for i in spec))
    if len(bases) > 1:
        if continuum:
            return list(infos)        # the stack header then has no spectral axis
        raise ValueError('the cubes have different kinds of spectral axis ({0}); LineStacker stacks channels '
                         'of one kind only (e.g. all FREQ, or all WAVE). Convert the cubes first.'.format(
                             ', '.join('{0} ({1})'.format(b, os.path.basename(next(
                                 i.path for i in spec if i.spec_base == b))) for b in bases)))
    if not continuum:
        bad = [i for i in spec if not i.spec_linear]
        if bad:
            raise ValueError('{0} has a non-linear spectral axis ({1}): LineStacker needs channels evenly spaced '
                             'in {2} (also to write the spectral axis of the stack). Regrid the cube to a linear '
                             'axis first (e.g. CASA imregrid, or tclean with a linear spectral grid).'.format(
                                 bad[0].path, bad[0].spec_ctype, bad[0].spec_base))
    out = []
    for inf in infos:
        if inf.spec_axis >= 0 and inf.nchan > 1 and inf.spec_unit != ref.spec_unit:
            if inf.spec_unit is None or ref.spec_unit is None:
                if inf.spec_cunit != ref.spec_cunit:
                    raise ValueError('cannot compare the spectral units of {0} ({1!r}) and {2} ({3!r})'.format(
                        ref.path, ref.spec_cunit, inf.path, inf.spec_cunit))
            else:
                inf = inf.spectral_in(ref.spec_unit)
        out.append(inf)
    frames = sorted(set(i.specsys or 'not given' for i in spec))
    if len(frames) > 1 and not continuum:
        warnings.warn('the cubes use different spectral reference frames (SPECSYS: {0}): the same frequency is a '
                      'different velocity in each (up to tens of km/s). Regrid the cubes to one frame (e.g. LSRK '
                      'with CASA imregrid or tclean outframe), or check the frame of the cubes without SPECSYS, '
                      'for an accurate line stack.'.format(', '.join(frames)), stacklevel=4)
    return out


def _intensity_scales(infos):
    """Factor converting each cube's pixel values to the intensity unit
    (BUNIT) of the first cube that has one. Compatible units are converted
    (mJy/beam to Jy/beam); physically different units (K and Jy/beam,
    Jy/pixel and Jy/beam) are an error; cubes without BUNIT, or with units
    astropy cannot read, give a warning."""
    have = [i for i in infos if i.bunit]
    if not have:
        return np.ones(len(infos)), None
    ref = have[0]
    scales, missing, unknown = np.ones(len(infos)), [], []
    for k, inf in enumerate(infos):
        if not inf.bunit:
            missing.append(os.path.basename(inf.path.rstrip('/')))
            continue
        if inf.bunit == ref.bunit:
            continue
        if inf.bunit_unit is None or ref.bunit_unit is None:
            if inf.bunit.replace(' ', '').lower() != ref.bunit.replace(' ', '').lower():
                unknown.append('{0} ({1!r})'.format(os.path.basename(inf.path.rstrip('/')), inf.bunit))
            continue
        try:
            scales[k] = inf.bunit_unit.to(ref.bunit_unit)
        except Exception:
            raise ValueError('{0} is in {1!r} but {2} in {3!r}: these intensity units cannot be stacked '
                             'together. Convert the cubes to one unit first.'.format(
                                 inf.path, inf.bunit, ref.path, ref.bunit))
    if missing:
        warnings.warn('{0} cube(s) have no BUNIT ({1}{2}): their pixels are taken to be in {3!r}, the unit of '
                      '{4}'.format(len(missing), ', '.join(missing[:3]), ', ...' if len(missing) > 3 else '',
                                   ref.bunit, os.path.basename(ref.path.rstrip('/'))), stacklevel=4)
    if unknown:
        warnings.warn('intensity units that cannot be compared with {0!r}: {1}; their pixels are stacked as they '
                      'are'.format(ref.bunit, ', '.join(unknown[:3])), stacklevel=4)
    return scales, ref.bunit


def _grid_difference(m, ref):
    """How the pixel grid m (CubeInfo.cel_matrix) differs from ref, as text
    ('' if it is the same)."""
    s_ref = np.sqrt(abs(np.linalg.det(ref)))
    if np.allclose(m, ref, rtol=0, atol=1e-4 * s_ref):
        return ''
    scale = np.sqrt(abs(np.linalg.det(m)))
    parts = []
    if abs(scale / s_ref - 1) > 1e-4:
        parts.append('{0:.4g}" pixels instead of {1:.4g}"'.format(scale * 3600, s_ref * 3600))
    if np.sign(np.linalg.det(m)) != np.sign(np.linalg.det(ref)):
        parts.append('mirrored')
    angle = np.degrees(np.arctan2(m[1, 0], m[0, 0]) - np.arctan2(ref[1, 0], ref[0, 0]))
    angle = (angle + 180.) % 360. - 180.
    if abs(angle) > 0.01:
        parts.append('rotated by {0:.2f} deg'.format(angle))
    return ', '.join(parts) or 'different pixel grid'


def _sky_geometry(inf):
    """(projection, celestial frame, distortion) of a cube, as text."""
    from astropy.wcs.utils import wcs_to_celestial_frame
    cel = inf.celestial
    projection = cel.wcs.ctype[cel.wcs.lng][4:].strip('-') or 'none'
    try:
        fr = wcs_to_celestial_frame(cel)
        frame = fr.name + (' {0}'.format(fr.equinox) if getattr(fr, 'equinox', None) is not None else '')
    except Exception:
        frame = 'unknown'
    return projection, frame, _distortion_signature(cel)


def _distortion_signature(cel):
    """The distortion terms of a celestial WCS: SIP orders and coefficient
    arrays, and the lookup tables (cpdis, det2im) with their reference
    pixels and steps. None where a term is absent."""
    sip = getattr(cel, 'sip', None)
    sig = {'sip': None if sip is None else
           tuple(None if a is None else np.asarray(a, dtype=float) for a in (sip.a, sip.b, sip.ap, sip.bp))}
    for name in ('cpdis1', 'cpdis2', 'det2im1', 'det2im2'):
        t = getattr(cel, name, None)
        sig[name] = None if t is None else (np.asarray(t.data, dtype=float), tuple(t.crpix), tuple(t.crval),
                                            tuple(t.cdelt))
    return sig


def _same_distortion(a, b):
    """True if two distortion signatures describe the same distortions."""
    def same(x, y):
        if x is None or y is None:
            return x is y
        if isinstance(x, tuple):
            return len(x) == len(y) and all(same(u, v) for u, v in zip(x, y))
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        return x.shape == y.shape and np.allclose(x, y, rtol=1e-9, atol=1e-12)
    return a.keys() == b.keys() and all(same(a[k], b[k]) for k in a)


def _check_grids(infos):
    """Warn when the cubes do not share one pixel grid (scale, rotation,
    orientation, projection, celestial frame, distortions): stamps are cut
    on each cube's own grid, not reprojected."""
    ref = infos[0]
    ref_geo = _sky_geometry(ref)
    issues = []
    for inf in infos[1:]:
        what = _grid_difference(inf.cel_matrix, ref.cel_matrix)
        geo = _sky_geometry(inf)
        extra = []
        if geo[0] != ref_geo[0]:
            extra.append('{0} projection instead of {1}'.format(geo[0], ref_geo[0]))
        if geo[1] != ref_geo[1]:
            extra.append('frame {0} instead of {1}'.format(geo[1], ref_geo[1]))
        if not _same_distortion(geo[2], ref_geo[2]):
            has, ref_has = (any(v is not None for v in g.values()) for g in (geo[2], ref_geo[2]))
            extra.append('different distortion terms' if has and ref_has else
                         'distortion terms' if has else 'no distortion terms')
        what = ', '.join([w for w in [what] + extra if w])
        if what:
            issues.append('{0}: {1}'.format(os.path.basename(inf.path.rstrip('/')), what))
    if issues:
        warnings.warn('{0} cube(s) do not share the pixel grid of {1} ({2}{3}). Every stamp is centred on its '
                      'source but cut on its own cube\'s grid, without reprojection, so the stacked image mixes '
                      'grids. Reproject the cubes to a common grid (e.g. CASA imregrid) to use the stacked image '
                      'spatially.'.format(len(issues), os.path.basename(ref.path.rstrip('/')), '; '.join(issues[:3]),
                                          '; ...' if len(issues) > 3 else ''), stacklevel=4)


def _mask(stampsize, maskradius):
    X = np.arange(0, stampsize) - int(stampsize / 2)
    X, Y = np.meshgrid(X, X)
    return (X**2 + Y**2 > maskradius**2).T.astype(float)


def _apply_weighting(prep, weighting, maskradius=0, primarybeam=None, hdul='auto'):
    if weighting == 'sigma2':
        _sigma2_weights(prep.data, prep.coords, maskradius)
    elif weighting == 'sigma2F':
        _sigma2_weights_spectral(prep.data, prep.coords, maskradius, prep.valid)
    elif weighting == '1/A':
        _amp_weights(prep.data, prep.coords, fit=False)
    elif weighting == 'sigma':
        _sigma_weights(prep, maskradius)
    elif weighting == 'pb':
        if primarybeam is None:
            raise ValueError("weighting='pb' needs primarybeam: the primary-beam image(s) of the cubes, "
                             "or a function such as LineStacker.image.GaussianPrimaryBeam(dish=12.)")
        _pb_weights(prep, primarybeam, hdul)
    else:
        raise ValueError("weighting must be None, 'sigma2', 'sigma2F', '1/A', 'sigma' or 'pb'")


def _sigma_weights(prep, maskradius=0):
    """weight *= 1/sigma, sigma being the standard deviation of the unmasked,
    non-blank stamp pixels. Without maskradius the mask is 2 beams wide."""
    s = prep.stampsize
    for i, coord in enumerate(prep.coords):
        radius = maskradius
        inf = prep.infos[i]
        if not radius:
            radius = 2 * abs(inf.beam_major_deg / inf.pixscale_deg) if inf.beam_major_deg else 10.
        sigma = _stamp_sigma(prep.data[i], _mask(s, min(radius, s / 2.)).astype(bool))
        if sigma is None:
            coord.weight = coord.weight * 0.
            warnings.warn('source {0}: fewer than 2 valid pixels, or zero noise, outside the mask of its stamp; '
                          'weight set to 0'.format(coord.index if coord.index is not None else i), stacklevel=4)
        else:
            coord.weight = coord.weight * 1. / sigma


class GaussianPrimaryBeam(object):
    """
        Gaussian primary beam (power pattern) of a dish, for weighting='pb':
        response exp(-4 ln2 r**2 / FWHM**2) with FWHM = factor * c / (nu * dish),
        r being the offset from the image reference pixel (the pointing centre
        of single-pointing images). factor=1.13 is the ALMA 12 m value, as in
        LineStacker.uv.

        Parameters
        ---------
        dish
            dish diameter in metres (default 12)
        factor
            FWHM in units of lambda/D (default 1.13)
        nu0
            frequency (Hz) to use for images without a frequency axis
    """

    def __init__(self, dish=12., factor=1.13, nu0=None):
        self.dish, self.factor, self.nu0 = float(dish), float(factor), nu0

    def fwhm(self, nu):
        """FWHM in radians at frequency nu (Hz)."""
        return self.factor * 299792458. / (nu * self.dish)

    def __call__(self, dx, dy, nu=None):
        nu = self.nu0 if nu is None else nu
        if nu is None:
            raise ValueError('GaussianPrimaryBeam needs a frequency: the cube has no frequency axis, give nu0 (Hz)')
        return np.exp(-4. * np.log(2.) * (dx**2 + dy**2) / self.fwhm(nu)**2)

    def __repr__(self):
        return 'GaussianPrimaryBeam(dish={0}, factor={1}, nu0={2})'.format(self.dish, self.factor, self.nu0)


def _pb_models(primarybeam, imagenames):
    """One primary-beam model (function or image path) per cube."""
    if callable(primarybeam):
        return [primarybeam] * len(imagenames)
    if isinstance(primarybeam, (str, os.PathLike)):
        if len(imagenames) != 1:
            raise ValueError('one primary-beam image was given for {0} cubes: give one per cube (a list in the '
                             'order of imagenames, or a dict {{cube: pb image}}), or a function'.format(len(imagenames)))
        return [os.fspath(primarybeam)]
    if isinstance(primarybeam, dict):
        norm = {os.path.normpath(os.path.abspath(os.fspath(k))): v for k, v in primarybeam.items()}
        out = []
        for name in imagenames:
            key = os.path.normpath(os.path.abspath(name))
            if key not in norm:
                raise ValueError('primarybeam has no entry for the cube {0}'.format(name))
            out.append(norm[key])
        return out
    seq = list(primarybeam)
    if len(seq) != len(imagenames):
        raise ValueError('primarybeam has {0} entries for {1} cubes'.format(len(seq), len(imagenames)))
    return seq


def _pb_image_response(pbname, inf, coord, value, hdul):
    """Primary-beam image value at a source (nearest pixel; channel of the
    spectral value `value`, or the mean over channels if value is None)."""
    pb, pinf = _io.read_cube(os.fspath(pbname), hdul)
    sky = inf.pixel_to_world(coord.ra, coord.dec)
    px, py = pinf.world_to_pixel(sky)
    ix, iy = LineStacker.round_half_up(float(np.atleast_1d(px)[0])), LineStacker.round_half_up(float(np.atleast_1d(py)[0]))
    if not (0 <= ix < pinf.nx and 0 <= iy < pinf.ny):
        return None
    spectrum = pb[ix, iy, 0, :]
    if pinf.nchan > 1 and value is not None:
        if pinf.spec_base != inf.spec_base or (pinf.spec_unit is None) != (inf.spec_unit is None):
            raise ValueError('the primary-beam image {0} has a {1} spectral axis but the cube {2} a {3} one: '
                             'give a primary beam on the same kind of axis, or a single-plane one'.format(
                                 pbname, pinf.spec_base or 'missing', inf.path, inf.spec_base or 'missing'))
        pin = pinf.spectral_in(inf.spec_unit) if inf.spec_unit is not None else pinf
        chan = pin.channel_of(value)
        if not (-0.5 <= chan <= pinf.nchan - 0.5):
            warnings.warn('the line of a source in {0} is outside the channels of its primary-beam image {1}; '
                          'the nearest primary-beam channel is used'.format(inf.path, pbname), stacklevel=5)
        ch = int(np.clip(LineStacker.round_half_up(chan), 0, pinf.nchan - 1))
        return float(spectrum[ch])
    # continuum, or a single-plane primary beam: mean of its non-blank planes
    good = spectrum[np.isfinite(spectrum)]
    return float(np.mean(good)) if good.size else np.nan     # blank: weight 0, with a warning


def _pb_weights(prep, primarybeam, hdul='auto'):
    """weight *= response**2 of the primary beam at each source: read from the
    cube's primary-beam image, or primarybeam(dx, dy, nu) with (dx, dy) the
    spherical offsets (radians, east and north on the sky) of the source from
    the image reference pixel. Sources outside their primary-beam image, or
    where the response is blank (NaN) or negative, get weight 0."""
    models = _pb_models(primarybeam, prep.imagenames)
    outside, blank = [], []
    for i, coord in enumerate(prep.coords):
        inf = prep.infos[i]
        model = models[coord.imagenumber]
        if prep.continuum or inf.spec_axis < 0:
            value = None
            centre = (inf.nchan - 1) / 2.
        else:
            centre = prep.centres[i]
            value = inf.value_of(centre)
        if callable(model):
            # spherical offsets of the source from the reference pixel's sky
            # position (exact for any projection, distortions included)
            cel = inf.celestial
            src = inf.pixel_to_world(coord.ra, coord.dec)
            ref = inf.pixel_to_world(cel.wcs.crpix[cel.wcs.lng] - 1., cel.wcs.crpix[cel.wcs.lat] - 1.)
            dlon, dlat = ref.spherical_offsets_to(src)
            dx, dy = float(dlon.rad[0]), float(dlat.rad[0])
            nu = None
            if inf.spectral_kind == 'frequency':
                nu = inf.value_of(centre) * (inf.spec_unit.to('Hz') if inf.spec_unit is not None else 1.)
            response = float(model(dx, dy, nu))
        else:
            response = _pb_image_response(model, inf, coord, value, hdul)
            if response is None:
                outside.append(coord.index if coord.index is not None else i)
                response = 0.
        if not np.isfinite(response) or response < 0:
            blank.append(coord.index if coord.index is not None else i)
            response = 0.
        coord.weight = coord.weight * response**2
    if outside:
        warnings.warn('{0} source(s) ({1}) are outside their primary-beam image: weight set to 0'.format(
            len(outside), ', '.join(str(k) for k in outside[:10])), stacklevel=4)
    if blank:
        warnings.warn('{0} source(s) ({1}) have a blank (NaN) or negative primary-beam response, e.g. beyond '
                      'the primary-beam cutoff: weight set to 0'.format(
                          len(blank), ', '.join(str(k) for k in blank[:10])), stacklevel=4)


def _stamp_sigma(stamp, keep):
    """Standard deviation of the valid pixels of a (x, y, stokes, chan)
    stamp: inside `keep` (x, y) and not blank (NaN, as everything outside
    the cube). None with fewer than 2 such pixels or zero noise."""
    good = keep[:, :, None, None] & np.isfinite(stamp)
    if good.sum() < 2:
        return None
    sigma = float(np.std(stamp[good]))
    return sigma if np.isfinite(sigma) and sigma > 0 else None


def _sigma2_weights(stamps, coords, maskradius=0.):
    s = stamps.shape[1]
    if maskradius >= s:
        raise Exception('the mask radius is bigger than the stamp size')
    if s <= 4:
        warnings.warn('stamp size is very small, sigma2 weights are probably not meaningful', stacklevel=4)
    keep = _mask(s, maskradius).astype(bool) if maskradius != 0 else np.ones((s, s), bool)
    for i, coord in enumerate(coords):
        sigma = _stamp_sigma(stamps[i], keep)
        if sigma is None:
            coord.weight = coord.weight * 0.
            warnings.warn('source {0}: fewer than 2 valid pixels, or zero noise, outside the mask of its stamp '
                          '(empty map?); weight set to 0'.format(coord.index if coord.index is not None else i),
                          stacklevel=4)
        else:
            coord.weight = coord.weight * 1. / sigma**2


def _sigma2_weights_spectral(stamps, coords, maskradius=0., inside=None):
    """weight *= 1/sigma_k**2 in every channel k, sigma_k being the standard
    deviation of the stamp's valid pixels in that channel: outside the mask,
    not blank (NaN, as everything outside the cube)."""
    s, nchan = stamps.shape[1], stamps.shape[-1]
    if maskradius >= s:
        raise Exception('the mask radius is bigger than the stamp size')
    if s <= 4:
        warnings.warn('stamp size is very small, sigma2F weights may be biased', stacklevel=4)
    keep = _mask(s, maskradius).astype(bool) if maskradius != 0 else np.ones((s, s), bool)
    starved = []
    for i, coord in enumerate(coords):
        st = stamps[i, :, :, 0, :]
        good = keep[:, :, None] & np.isfinite(st)
        n = good.sum(axis=(0, 1))
        with np.errstate(invalid='ignore', divide='ignore'):
            mean = np.where(good, st, 0.).sum(axis=(0, 1)) / n
            sigmaF = np.sqrt(np.where(good, (st - mean)**2, 0.).sum(axis=(0, 1)) / n)
        base = np.broadcast_to(np.asarray(coord.weight, dtype=float), (nchan,))
        use = (base != 0) & np.isfinite(base) & (n >= 2) & np.isfinite(sigmaF) & (sigmaF > 0)
        if inside is not None:
            use &= inside[i]
            if np.any(inside[i] & (base != 0) & (n < 2)):
                starved.append(coord.index if coord.index is not None else i)
        w = np.zeros(nchan)
        w[use] = base[use] / sigmaF[use]**2
        coord.weight = w
    if starved:
        warnings.warn('sigma2F: {0} source(s) ({1}) have channels with fewer than 2 valid pixels outside the mask; '
                      'those channels get weight 0'.format(len(starved), ', '.join(str(k) for k in starved[:10])),
                      stacklevel=4)


def _amp_weights(stamps, coords, fit=False):
    """weight /= line amplitude (central-pixel spectrum); a non-positive or
    undefined amplitude gives weight 0."""
    c = int(stamps.shape[1] / 2.)
    bad = []
    for i, coord in enumerate(coords):
        spectrum = stamps[i, c, c, 0, :]
        if fit:
            import LineStacker.tools.fit as fitTool
            try:
                amp = fitTool.GaussFit(fctToFit=spectrum, returnInfos=True)[1][0]
            except RuntimeError:
                amp = np.nan
        else:
            amp = np.nanmax(spectrum) if np.any(np.isfinite(spectrum)) else np.nan
        if not np.isfinite(amp) or amp <= 0:
            bad.append(coord.index if coord.index is not None else i)
            coord.weight = coord.weight * 0.
        else:
            coord.weight = coord.weight * 1. / amp
    if bad:
        warnings.warn("weighting='1/A': {0} source(s) ({1}) have no positive line amplitude; weight set to 0".format(
            len(bad), ', '.join(str(k) for k in bad[:10])), stacklevel=4)


def _combine(stamps, weights, inside, method):
    """Combine (nsrc, x, y, stokes, chan) stamps over sources, pixel by
    pixel. Blank (NaN) pixels are left out. Mean: weighted, NaN where no
    source has a valid pixel with weight > 0. Median: over the sources whose
    cube covers the channel, NaN where none does."""
    if method == 'mean':
        finite = np.isfinite(stamps)
        if finite.all():
            num = np.einsum('nabcf,nf->abcf', stamps, weights)
            den = weights.sum(axis=0)
        else:
            num = np.einsum('nabcf,nf->abcf', np.where(finite, stamps, 0.), weights)
            den = np.einsum('nabcf,nf->abcf', finite.astype(float), weights)
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(den > 0, num / den, np.nan)
    if method == 'median':
        tmp = np.where(inside[:, None, None, None, :], stamps, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            return np.nanmedian(tmp, axis=0)
    raise ValueError("method must be 'mean' or 'median'")


def _write_weights(path, user_coords, prep):
    with open(path, 'w') as fh:
        for c in prep.coords:
            orig = user_coords[c.index] if user_coords is not None else c
            ra = orig.ra.deg if hasattr(orig.ra, 'deg') else orig.ra
            dec = orig.dec.deg if hasattr(orig.dec, 'deg') else orig.dec
            w = np.atleast_1d(np.asarray(c.weight, dtype=float))
            ws = str(w[0]) if np.all(w == w[0]) else ', '.join(str(v) for v in w)
            fh.write('{0}, {1}, {2}, {3}\n'.format(ra, dec, c.spec, ws))


def _output_header(prep, nchan, spectralMethod, fEm, lEm, method):
    """FITS header of the stacked cube: spatial offsets around (0, 0) at the
    stamp centre, spectral axis centred on the line."""
    from astropy.io import fits
    infos = prep.infos
    t = infos[0]
    cel = t.celestial
    lng, lat = cel.wcs.lng, cel.wcs.lat
    h = fits.Header()
    s = prep.stampsize
    h['CTYPE1'], h['CTYPE2'] = cel.wcs.ctype[lng], cel.wcs.ctype[lat]
    h['CRVAL1'], h['CRVAL2'] = 0., 0.
    h['CRPIX1'] = h['CRPIX2'] = int(s / 2) + 1
    cdelt, pc = _cdelt_pc(t.cel_matrix)
    h['CDELT1'], h['CDELT2'] = float(cdelt[0]), float(cdelt[1])
    if not np.allclose(pc, np.eye(2), rtol=0, atol=1e-12):
        # rotated or skewed grid (stamps keep the first cube's orientation)
        for (a, b) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            h['PC{0}_{1}'.format(a + 1, b + 1)] = float(pc[a, b])
    h['CUNIT1'] = h['CUNIT2'] = 'deg'
    for key in ('RADESYS', 'EQUINOX'):
        if key in t.header:
            h[key] = t.header[key]

    if prep.continuum:
        _continuum_header(h, prep)
        _finish_header(h, prep, method, 'continuum ({0})'.format(
            ', '.join(sorted(set(c for c in prep.collapse if c))) + ' collapse' if any(prep.collapse) else 'images'))
        return h

    zs = np.array([float(c.spec) for c in prep.coords]) if spectralMethod != 'channel' else None
    if prep.regridFromZ:
        cdelts = np.array([_rest_channel_width(inf, z) * np.sign(infos[0].spec_cdelt)
                           for inf, z in zip(infos, zs)]) * np.abs(prep.steps)
    else:
        cdelts = np.array([abs(inf.spec_cdelt) * np.sign(infos[0].spec_cdelt) for inf in infos])
    if spectralMethod == 'z':
        kind = t.spectral_kind
        if kind == 'frequency':
            crval = fEm
            cdelt = np.mean(cdelts if prep.regridFromZ else cdelts * (1. + zs))
            h['RESTFRQ'] = fEm
        else:
            crval = lEm
            cdelt = np.mean(cdelts if prep.regridFromZ else cdelts / (1. + zs))
            h['RESTWAV'] = lEm
    elif spectralMethod == 'obsSpec':
        crval, cdelt = float(np.mean(zs)), float(np.mean(cdelts))
    else:
        crval = float(np.mean([inf.value_of(cn) for inf, cn in zip(infos, prep.centres)]))
        cdelt = float(np.mean(cdelts))
    h['CTYPE3'] = t.spec_ctype or 'FREQ'
    if t.spec_cunit:
        h['CUNIT3'] = t.spec_cunit
    h['CRPIX3'] = int(nchan / 2) + 1
    h['CRVAL3'] = float(crval)
    h['CDELT3'] = float(cdelt)
    if 'SPECSYS' in t.header:
        h['SPECSYS'] = 'SOURCE' if spectralMethod == 'z' else t.header['SPECSYS']
    _finish_header(h, prep, method, 'line, spectralMethod={0}{1}'.format(
        spectralMethod, ', regridFromZ' if prep.regridFromZ else ''))
    return h


def _cdelt_pc(m):
    """FITS CDELT and PC of a (2, 2) pixel-to-sky matrix (CubeInfo.cel_matrix):
    CDELT_i is the length of row i, with the sign of its diagonal element, so
    an unrotated grid gives PC = identity."""
    norms = np.sqrt((m**2).sum(axis=1))
    signs = np.where(np.diag(m) != 0, np.sign(np.diag(m)), 1.)
    cdelt = norms * signs
    return cdelt, m / cdelt[:, None]


def _continuum_header(h, prep):
    """Spectral axis of a continuum stack: one channel covering the bands
    of the stacked images (none if no image has a spectral axis)."""
    spec = [i for i in prep.infos if i.spec_axis >= 0]
    if not spec or len(set(i.spec_base for i in spec)) > 1 or not all(i.spec_linear for i in spec):
        return                    # no spectral axis, kinds that cannot be averaged, or non-linear axes
    t = spec[0]
    h['CTYPE3'] = t.spec_ctype
    if t.spec_cunit:
        h['CUNIT3'] = t.spec_cunit
    h['CRPIX3'] = 1.
    h['CRVAL3'] = float(np.mean([i.value_of((i.nchan - 1) / 2.) for i in spec]))
    h['CDELT3'] = float(np.mean([i.nchan * abs(i.spec_cdelt) for i in spec])) * (np.sign(t.spec_cdelt) or 1.)
    for key in ('SPECSYS', 'RESTFRQ', 'RESTWAV'):
        if key in t.header:
            h[key] = t.header[key]
    methods = sorted(set(c for c in prep.collapse if c))
    if methods:
        h['COLLAPSE'] = (', '.join(methods), 'channels combined (BUNIT is per channel)')


def _finish_header(h, prep, method, what):
    infos = prep.infos
    t = infos[0]
    beams = [(i.beam_major_deg, i.beam_minor_deg, i.beam_pa_deg) for i in infos if i.beam_major_deg]
    if beams:
        b = np.array([[np.nan if v is None else v for v in beam] for beam in beams], dtype=float)
        h['BMAJ'] = (float(np.mean(b[:, 0])), 'mean of the stacked images')
        if np.any(np.isfinite(b[:, 1])):
            h['BMIN'] = (float(np.nanmean(b[:, 1])), 'mean of the stacked images')
        pa = b[np.isfinite(b[:, 2]), 2]
        if pa.size:
            h['BPA'] = (_axial_mean_deg(pa), 'axial mean of the stacked images')
    bunit = getattr(prep, 'bunit', None) or t.header.get('BUNIT')
    if bunit:
        h['BUNIT'] = bunit
    h['OBJECT'] = 'LineStacker stack'
    h['NSTACK'] = (len(prep.coords), 'number of stacked sources')
    h['HISTORY'] = 'LineStacker {0} {1} stack, {2}'.format(LineStacker.__version__, method, what)


def _axial_mean_deg(angles):
    """Mean of position angles (degrees) defined modulo 180: 179 and 1 give 0."""
    a = np.deg2rad(np.asarray(angles, dtype=float)) * 2.
    return float(np.degrees(np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a)))) / 2.)


def _write_stacked_image(imagename, pixels, prep, spectralMethod='z', method='mean', fEm=-1, lEm=-1):
    header = _output_header(prep, pixels.shape[-1], spectralMethod, fEm, lEm, method)
    _io.write_cube(imagename, pixels, header)


def readCube(imagename, hdul='auto', collapse=None):
    """
        Reads a cube or image (FITS file or CASA image) as a numpy array
        ordered (x, y, stokes, channel), x along RA and y along Dec, as
        LineStacker uses it. Blank pixels are NaN.

        Parameters
        ---------
        imagename
            path of the FITS file or CASA image
        hdul
            FITS HDU holding the data ('auto': first HDU with an image)
        collapse
            None (default) for the full cube, or 'sum', 'mean', 'median' to
            collapse the spectral axis (blank channels ignored).
    """
    return _io.read_cube(imagename, hdul, collapse=collapse)[0]


def clear_cache():
    """
        Forget the cubes read so far: LineStacker reads each cube once and keeps
        it in memory for the following stacks and analysis tools (a cube is read
        again when its file changes). This frees that memory.
    """
    _io.clear_cache()


def getFlux(imagename, hdul='auto'):
    """Pixel value at the reference pixel of an image, e.g. the centre of a
    stacked image (first Stokes plane and channel)."""
    cube, info = _io.read_cube(imagename, hdul)
    cel = info.celestial
    x = int(round(cel.wcs.crpix[cel.wcs.lng] - 1))
    y = int(round(cel.wcs.crpix[cel.wcs.lat] - 1))
    return float(cube[x, y, 0, 0])

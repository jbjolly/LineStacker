"""Synthetic cubes with known content for the LineStacker tests.

Run from the repository root with:  python -m pytest tests
Tests that need CASA (casatools) are skipped when it is not installed.
"""
import os
import sys
import warnings

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    # through LineStacker, so CASA's log goes where LineStacker.casa_logfile
    # says (not a new casa-<date>.log in the folder pytest runs from)
    from LineStacker._casa_compat import new_tool as _new_tool
    _new_tool('image').done()
    HAVE_CASA = True
except Exception:
    HAVE_CASA = False

needs_casa = pytest.mark.skipif(not HAVE_CASA, reason='casatools not installed')

FEM = 1897420620253.1646        # [CII] rest frequency used by the examples (Hz)


def wcs_header(nx, ny, nchan, crval3, cdelt3, crpix3, kind='freq', stokes=True, cd=False, spectral=True):
    from astropy.io import fits
    h = fits.Header()
    h['CTYPE1'], h['CTYPE2'] = 'RA---SIN', 'DEC--SIN'
    h['CRVAL1'], h['CRVAL2'] = 53.1, -27.8
    h['CRPIX1'], h['CRPIX2'] = nx / 2. + 0.5, ny / 2. + 0.5
    if cd:
        h['CD1_1'], h['CD2_2'] = -0.25 / 3600, 0.25 / 3600
    else:
        h['CDELT1'], h['CDELT2'] = -0.25 / 3600, 0.25 / 3600
    h['CUNIT1'] = h['CUNIT2'] = 'deg'
    h['RADESYS'] = 'ICRS'
    h['BMAJ'], h['BMIN'], h['BPA'] = 1.0 / 3600, 0.8 / 3600, 30.
    h['BUNIT'] = 'Jy/beam'
    if not spectral:
        return h
    if kind == 'freq':
        h['CTYPE3'], h['CUNIT3'] = 'FREQ', 'Hz'
    else:
        h['CTYPE3'], h['CUNIT3'] = 'AWAV', 'Angstrom'
    h['CRVAL3'], h['CRPIX3'] = crval3, crpix3
    if cd:
        h['CD3_3'] = cdelt3
    else:
        h['CDELT3'] = cdelt3
    if stokes:
        h['CTYPE4'], h['CRVAL4'], h['CDELT4'], h['CRPIX4'] = 'STOKES', 1., 1., 1.
    return h


class Cube(object):
    """One synthetic cube: array (x, y, chan) plus how to write it.
    nchan=None makes a 2D image (no spectral axis)."""

    def __init__(self, nx=40, ny=40, nchan=60, crval3=232e9, cdelt3=-7.6e7, crpix3=1.,
                 kind='freq', stokes=True, cd=False):
        self.spectral = nchan is not None
        nchan = nchan or 1
        self.nx, self.ny, self.nchan = nx, ny, nchan
        self.header = wcs_header(nx, ny, nchan, crval3, cdelt3, crpix3, kind, stokes and self.spectral,
                                 cd, spectral=self.spectral)
        self.crval3, self.cdelt3, self.crpix0 = crval3, cdelt3, crpix3 - 1.
        self.kind, self.stokes = kind, stokes
        self.arr = np.zeros((nx, ny, nchan))

    def channel_of(self, value):
        return self.crpix0 + (value - self.crval3) / self.cdelt3

    def add_line(self, x, y, chan, amp, sig_chan, sig_pix=1.5):
        X, Y = np.meshgrid(np.arange(self.nx), np.arange(self.ny), indexing='ij')
        spat = np.exp(-((X - x)**2 + (Y - y)**2) / (2 * sig_pix**2))
        spec = np.exp(-(np.arange(self.nchan) - chan)**2 / (2 * sig_chan**2))
        self.arr += amp * spat[:, :, None] * spec[None, None, :]

    def add_noise(self, sigma, rng):
        self.arr += rng.normal(0, sigma, self.arr.shape)

    def write(self, path, hdu=0, casa=False):
        from astropy.io import fits
        a = self.arr.transpose(2, 1, 0).astype('>f4')       # FITS (chan, y, x)
        if not self.spectral:
            a = a[0]
        elif self.stokes:
            a = a[None]
        if hdu == 0:
            hl = fits.HDUList([fits.PrimaryHDU(a, self.header)])
        else:
            hl = fits.HDUList([fits.PrimaryHDU()] + [fits.ImageHDU(a, self.header)] * hdu)
        fitspath = path if not casa else path + '.tmp.fits'
        hl.writeto(fitspath, overwrite=True)
        if casa:
            from casatools import image
            ia = image()
            ia.fromfits(outfile=path, infile=fitspath, overwrite=True)
            ia.done()
            os.remove(fitspath)
        return path

    def sky(self, x, y):
        from astropy.wcs import WCS
        from astropy.coordinates import SkyCoord
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            w = WCS(self.header).celestial
        ra, dec = w.all_pix2world([[x, y]], 0)[0]
        return SkyCoord(ra, dec, unit='deg', frame='icrs')


def casa_style(sky):
    """'HH:MM:SS.ss, -DD.MM.SS.sss' as written by CASA (dots in Dec)."""
    ra = sky.ra.to_string(unit='hourangle', sep=':', precision=4, pad=True)
    dec = sky.dec.to_string(unit='deg', sep='.', precision=4, alwayssign=True, pad=True)
    return ra, dec


def reference_stack(cubes, pix, centres, nchan_out, stampsize, weights, method):
    """Straightforward re-implementation of the LineStacker stacking rule,
    loop by loop, to check the vectorised code against."""
    n = len(cubes)
    s = stampsize
    stamps = np.zeros((n, s, s, nchan_out))
    ok = np.zeros((n, nchan_out), bool)
    for i, (cube, (x, y), c) in enumerate(zip(cubes, pix, centres)):
        bx = int(np.floor(x - s // 2 + 0.5))
        by = int(np.floor(y - s // 2 + 0.5))
        for k in range(nchan_out):
            ch = c + k - nchan_out // 2
            if 0 <= ch < cube.shape[2]:
                stamps[i, :, :, k] = cube[bx:bx + s, by:by + s, ch]
                ok[i, k] = True
    out = np.zeros((s, s, nchan_out))
    for k in range(nchan_out):
        use = np.where(ok[:, k])[0]
        if len(use) == 0:
            continue
        if method == 'mean':
            w = np.array([np.broadcast_to(weights[i], (nchan_out,))[k] for i in use])
            if w.sum() > 0:
                out[:, :, k] = np.tensordot(w, stamps[use, :, :, k], axes=(0, 0)) / w.sum()
        else:
            out[:, :, k] = np.median(stamps[use, :, :, k], axis=0)
    return out, ok


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from LineStacker import _io
    _io.clear_cache()
    return tmp_path

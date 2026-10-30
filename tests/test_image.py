import warnings

import numpy as np
import pytest

from conftest import (Cube, FEM, casa_style, needs_casa, reference_stack)

import LineStacker
import LineStacker.image as li
import LineStacker.tools.fit as fitTools


def build(n=8, seed=1, flip=(), edge=None, noise=0., casa=False, sig_chan=3.,
          integer=True, zrange=None, rest_sigma_hz=None, nchan=60):
    """n cubes with one line each; returns coordfiles, imagenames, cubes, truth."""
    rng = np.random.default_rng(seed)
    coordfiles, imagenames, cubes, truth = [], [], [], []
    for i in range(n):
        if zrange is not None:
            # each cube observed around its own line, as in real data
            z = rng.uniform(*zrange)
            chan = rng.uniform(nchan / 2 - 8, nchan / 2 + 8)
            cube = Cube(nchan=nchan, crval3=FEM / (1 + z) + chan * 7.6e7)
        elif i in flip:
            # same frequency grid, stored in increasing order
            cube = Cube(nchan=nchan, crval3=232e9 - (nchan - 1) * 7.6e7, cdelt3=7.6e7, crpix3=1.)
        else:
            cube = Cube(nchan=nchan)
        x, y = (float(rng.integers(12, 28)), float(rng.integers(12, 28))) if integer else \
            (rng.uniform(12, 28), rng.uniform(12, 28))
        if zrange is not None:
            pass
        else:
            chan_int = int(rng.integers(15, nchan - 15)) if edge is None or i not in edge else edge[i]
            # flipped cubes: same physical frequency -> mirrored channel
            f = 232e9 - chan_int * 7.6e7
            chan = cube.channel_of(f)
            z = FEM / f - 1.
        amp = rng.uniform(1., 3.)
        sc = sig_chan if rest_sigma_hz is None else rest_sigma_hz / (1 + z) / abs(cube.cdelt3)
        cube.add_line(x, y, chan, amp, sc)
        if noise:
            cube.add_noise(noise, rng)
        name = 'cube_%02d' % i + ('.image' if casa else '.fits')
        cube.write(name, casa=casa)
        ra, dec = casa_style(cube.sky(x, y))
        cf = 'cube_%02d_coord.txt' % i
        with open(cf, 'w') as fh:
            fh.write('%s, %s, %r\n' % (ra, dec, z))
        coordfiles.append(cf)
        imagenames.append(name)
        cubes.append(cube)
        truth.append(dict(x=x, y=y, chan=chan, amp=amp, z=z, sig=sc))
    return coordfiles, imagenames, cubes, truth


def test_exact_line_recovery(workdir):
    cf, names, cubes, truth = build()
    coords = LineStacker.readCoords(cf)
    res = li.stack(coords, imagenames=names, fEm=FEM, stampsize=12, N_chans=21)
    o = np.arange(21) - 10
    expected = np.mean([t['amp'] for t in truth]) * np.exp(-o**2 / (2 * 3.**2))
    np.testing.assert_allclose(res[1], expected, rtol=2e-6, atol=1e-7)
    assert np.argmax(res[1]) == 10
    np.testing.assert_allclose(res[0], np.mean([t['amp'] for t in truth]), rtol=2e-6)
    ref, ok = reference_stack([c.arr.astype('f4') for c in cubes], [(t['x'], t['y']) for t in truth],
                              [int(round(t['chan'])) for t in truth], 21, 12, [1.] * len(truth), 'mean')
    np.testing.assert_allclose(res[2][:, :, 0, :], ref, rtol=1e-6, atol=1e-7)
    # user's coordinates are not modified
    assert coords[0].type == 'physical' and coords[0].weight == 1.


def test_reversed_spectral_axis_gives_same_stack(workdir):
    cf, names, _, _ = build(seed=4)
    a = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=12, N_chans=21, outfile=None)
    cf, names, _, _ = build(seed=4, flip=(1, 3, 5))
    b = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=12, N_chans=21, outfile=None)
    np.testing.assert_allclose(a[2], b[2], rtol=1e-6, atol=1e-7)


@needs_casa
def test_casa_images_equal_fits(workdir):
    cf, names, _, _ = build(seed=7, noise=0.1)
    a = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=12, N_chans=21, outfile=None)
    cf, names, _, _ = build(seed=7, noise=0.1, casa=True)
    b = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=12, N_chans=21, outfile=None)
    np.testing.assert_array_equal(a[2], b[2])


@pytest.mark.parametrize('method', ['mean', 'median'])
def test_spectral_edges(workdir, method):
    cf, names, cubes, truth = build(n=6, seed=3, edge={0: 3, 1: 56, 2: 1}, noise=0.05)
    res = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=10,
                   N_chans=21, method=method, outfile=None)
    ref, ok = reference_stack([c.arr.astype('f4') for c in cubes], [(t['x'], t['y']) for t in truth],
                              [int(round(t['chan'])) for t in truth], 21, 10, [1.] * 6, method)
    np.testing.assert_allclose(res[2][:, :, 0, :], ref, rtol=1e-6, atol=1e-6)
    np.testing.assert_array_equal(res[3], (~ok).sum(0))


def test_sigma2_weights_match_reference_and_do_not_compound(workdir):
    cf, names, cubes, truth = build(n=6, seed=5, noise=0.3)
    coords = LineStacker.readCoords(cf)
    kw = dict(imagenames=names, fEm=FEM, stampsize=12, N_chans=21, weighting='sigma2', maskradius=3, outfile=None)
    r1 = li.stack(coords, **kw)
    r2 = li.stack(coords, **kw)
    np.testing.assert_array_equal(r1[2], r2[2])
    # independent weights
    s = 12
    X = np.arange(s) - s // 2
    X, Y = np.meshgrid(X, X)
    mask = (X**2 + Y**2 > 9).astype(float)
    ref_w, stamps = [], []
    for c, t in zip(cubes, truth):
        ch = int(round(t['chan']))
        bx, by = int(t['x']) - 6, int(t['y']) - 6
        st = c.arr.astype('f4')[bx:bx + s, by:by + s, ch - 10:ch + 11].astype(float)
        tmp = st * mask[:, :, None]
        ref_w.append(1. / np.std(tmp[np.nonzero(tmp)])**2)
        stamps.append(st)
    ref = np.tensordot(ref_w, stamps, axes=(0, 0)) / np.sum(ref_w)
    np.testing.assert_allclose(r1[2][:, :, 0, :], ref, rtol=1e-6)


def test_regrid_from_z_equalises_rest_frame_widths(workdir):
    rest_sigma = 4 * 7.6e7 * 6.0           # Hz, rest frame
    cf, names, _, truth = build(n=10, seed=11, zrange=(4.0, 7.0), rest_sigma_hz=rest_sigma, integer=True)
    coords = LineStacker.readCoords(cf)
    from astropy.io import fits
    zs = np.array([t['z'] for t in truth])
    for method, zref in (('scaleToMax', zs.max()), ('scaleToMin', zs.min())):
        res = li.stack(coords, imagenames=names, fEm=FEM, stampsize=10, N_chans=41,
                       regridFromZ=True, regridMethod=method, outfile='regrid.fits')
        h = fits.getheader('regrid.fits')
        # every source now has the rest-frame channel width of the reference source
        np.testing.assert_allclose(abs(h['CDELT3']), 7.6e7 * (1 + zref), rtol=1e-9)
        prep = li._prepare(coords, names, stampsize=10, N_chans=41, fEm=FEM,
                           regridFromZ=True, regridMethod=method)
        widths = []
        for i in range(len(truth)):
            p = fitTools.GaussFit(fctToFit=prep.data[i, 5, 5, 0, :], returnInfos=True)[1]
            assert abs(p[1] - 20) < 0.05          # fractional centres are honoured
            widths.append(abs(p[2]) * abs(h['CDELT3']))
        # linear interpolation broadens a little when under-sampling (scaleToMax)
        np.testing.assert_allclose(widths, rest_sigma, rtol=0.03 if method == 'scaleToMax' else 0.01)
        p = fitTools.GaussFit(fctToFit=res[1], returnInfos=True)[1]
        assert abs(abs(p[2]) * abs(h['CDELT3']) / rest_sigma - 1) < 0.03
    # without regridding, observed widths differ by (1+z) and the stack is broader
    res = li.stack(coords, imagenames=names, fEm=FEM, stampsize=10, N_chans=41, outfile='raw.fits')
    h = fits.getheader('raw.fits')
    p = fitTools.GaussFit(fctToFit=res[1], returnInfos=True)[1]
    spread = np.std([rest_sigma / (1 + z) / 7.6e7 for z in zs])
    assert spread > 0.3


def test_output_header_wcs(workdir):
    cf, names, _, truth = build(n=5, seed=2)
    res = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=12, N_chans=21,
                   outfile='out.fits')
    from astropy.io import fits
    from astropy.wcs import WCS
    d, h = fits.getdata('out.fits'), fits.getheader('out.fits')
    assert d.shape == (21, 12, 12)
    np.testing.assert_allclose(d.transpose(2, 1, 0), res[2][:, :, 0, :], rtol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        w = WCS(h)
    peak = int(np.argmax(res[1]))
    f_peak = w.spectral.pixel_to_world(peak).to_value('Hz')
    assert abs(f_peak - FEM) < 1.
    assert h['CRPIX1'] == 7 and h['NSTACK'] == 5
    zs = np.array([t['z'] for t in truth])
    np.testing.assert_allclose(h['CDELT3'], np.mean(-7.6e7 * (1 + zs)), rtol=1e-9)


@needs_casa
def test_casa_output_image(workdir):
    cf, names, _, _ = build(n=4, seed=8)
    res = li.stack(LineStacker.readCoords(cf), imagenames=names, fEm=FEM, stampsize=8, N_chans=15,
                   outfile='stackResult.image')
    from casatools import image
    ia = image()
    ia.open('stackResult.image')
    chunk = ia.getchunk()
    names_ax = ia.coordsys().names()
    ia.done()
    assert chunk.shape == (8, 8, 1, 15)
    assert names_ax == ['Right Ascension', 'Declination', 'Stokes', 'Frequency']
    np.testing.assert_allclose(chunk, res[2], rtol=1e-6)


def test_wavelength_cube_in_extension(workdir):
    lem = 6564.61
    rng = np.random.default_rng(0)
    cf, names, amps = [], [], []
    for i in range(4):
        cube = Cube(nx=30, ny=30, nchan=80, crval3=8000., cdelt3=1.25, crpix3=1., kind='wave',
                    stokes=False, cd=True)
        chan = int(rng.integers(20, 60))
        lam = cube.crval3 + chan * 1.25
        z = lam / lem - 1
        a = rng.uniform(1, 2)
        cube.add_line(15, 15, chan, a, 2.)
        cube.write('muse_%d.fits' % i, hdu=1)
        ra, dec = casa_style(cube.sky(15, 15))
        with open('muse_%d_coord.txt' % i, 'w') as fh:
            fh.write('%s %s %r\n' % (ra, dec, z))      # space separated
        cf.append('muse_%d_coord.txt' % i)
        names.append('muse_%d.fits' % i)
        amps.append(a)
    res = li.stack(LineStacker.readCoords(cf), imagenames=names, lEm=lem, stampsize=8, N_chans=11,
                   outfile=None)
    assert np.argmax(res[1]) == 5
    np.testing.assert_allclose(res[0], np.mean(amps), rtol=1e-6)


def test_obsspec_and_channel_methods(workdir):
    cf, names, cubes, truth = build(n=4, seed=9)
    coords = LineStacker.readCoords(cf)
    a = li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=11, outfile=None)
    for c, t, cube in zip(coords, truth, cubes):
        c.spec = cube.crval3 + (t['chan'] - cube.crpix0) * cube.cdelt3
    b = li.stack(coords, imagenames=names, spectralMethod='obsSpec', stampsize=8, N_chans=11, outfile=None)
    for c, t in zip(coords, truth):
        c.spec = t['chan']
    c_ = li.stack(coords, imagenames=names, spectralMethod='channel', stampsize=8, N_chans=11, outfile=None)
    np.testing.assert_allclose(a[2], b[2])
    np.testing.assert_allclose(a[2], c_[2])


def test_coords_none_stacks_cube_centres(workdir):
    cube = Cube(nx=20, ny=20, nchan=21)
    cube.add_line(10, 10, 10, 2., 2.)
    cube.write('c.fits')
    res = li.stack(None, imagenames=['c.fits', 'c.fits'], stampsize=6, N_chans=9, outfile=None)
    np.testing.assert_allclose(res[0], 2., rtol=1e-6)

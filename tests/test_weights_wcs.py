"""Weighting, spectral-axis and sky-grid handling of LineStacker.image
(fixes from the 2.0 code review)."""
import warnings

import numpy as np
import pytest

from conftest import Cube, FEM

import LineStacker
import LineStacker.image as li
import LineStacker.analysisTools as A
from LineStacker import _io


def line_cubes(n=4, seed=3, noise=0.05, edit=None, nchan=50, nx=40):
    """n FITS cubes with one line each (same frequency grid unless `edit`
    changes a header); returns coords (from the sky positions), names, truth."""
    rng = np.random.default_rng(seed)
    names, truth, coords = [], [], []
    for i in range(n):
        cube = Cube(nx=nx, ny=nx, nchan=nchan)
        if edit is not None:
            edit(i, cube)
        ch = int(rng.integers(15, nchan - 15))
        f = 232e9 - ch * 7.6e7
        z = FEM / f - 1.
        x, y = float(rng.integers(14, nx - 14)), float(rng.integers(14, nx - 14))
        cube.add_line(x, y, ch, rng.uniform(1., 3.), 2.5)
        if noise:
            cube.add_noise(noise, rng)
        name = 'w%02d.fits' % i
        cube.write(name)
        sky = cube.sky(x, y)
        coords.append(LineStacker.Coord(sky.ra, sky.dec, spec=z, image=name))
        names.append(name)
        truth.append(dict(x=x, y=y, z=z, chan=ch, cube=cube))
    return coords, names, truth


def ghz(i, cube):
    """Same frequencies, header in GHz for odd cubes."""
    if i % 2:
        h = cube.header
        h['CRVAL3'], h['CDELT3'], h['CUNIT3'] = h['CRVAL3'] / 1e9, h['CDELT3'] / 1e9, 'GHz'


def no_grid_warning(record):
    return not [w for w in record if 'pixel grid' in str(w.message)]


# -- weights -----------------------------------------------------------------

def test_sigma2F_keeps_the_source_weights(workdir):
    coords, names, _ = line_cubes()
    p1 = li._prepare(coords, names, stampsize=16, N_chans=11, fEm=FEM, weighting='sigma2F')
    for c in coords:
        c.weight = 2.5
    p2 = li._prepare(coords, names, stampsize=16, N_chans=11, fEm=FEM, weighting='sigma2F')
    np.testing.assert_allclose(p2.weights(), 2.5 * p1.weights(), rtol=1e-12)
    # and a weight per source changes only that source
    coords[1].weight = 0.
    p3 = li._prepare(coords, names, stampsize=16, N_chans=11, fEm=FEM, weighting='sigma2F')
    assert np.all(p3.weights()[1] == 0) and np.all(p3.weights()[0] > 0)


def test_sigma2F_uses_only_valid_pixels():
    rng = np.random.default_rng(0)
    s, nchan = 12, 5
    stamps = rng.normal(0., 1., (1, s, s, 1, nchan))
    stamps[0, :3, :, 0, :] = np.nan               # blank pixels
    stamps[0, 5, 0, 0, :] = 0.                     # a genuine zero counts
    coord = LineStacker.Coord(0, 0, spec=0, weight=3.)
    li._sigma2_weights_spectral(stamps, [coord], maskradius=3, inside=np.ones((1, nchan), bool))
    keep = li._mask(s, 3).astype(bool)
    for k in range(nchan):
        v = stamps[0, :, :, 0, k][keep & np.isfinite(stamps[0, :, :, 0, k])]
        assert 0. in v
        assert coord.weight[k] == pytest.approx(3. / np.std(v)**2, rel=1e-12)


def test_sigma2F_channels_without_pixels_get_zero_weight():
    s, nchan = 8, 4
    stamps = np.random.default_rng(1).normal(0., 1., (1, s, s, 1, nchan))
    stamps[0, :, :, 0, 2] = np.nan                # blank channel
    coord = LineStacker.Coord(0, 0, spec=0)
    with pytest.warns(UserWarning, match='fewer than 2 valid pixels'):
        li._sigma2_weights_spectral(stamps, [coord], inside=np.ones((1, nchan), bool))
    assert coord.weight[2] == 0 and np.all(coord.weight[[0, 1, 3]] > 0)


def test_negative_or_nan_weights_are_rejected(workdir):
    coords, names, _ = line_cubes(noise=0.)
    coords[2].weight = -1.
    with pytest.raises(ValueError, match='source 2 has weight -1'):
        li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None)
    coords[2].weight = np.nan
    with pytest.raises(ValueError, match='finite'):
        li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None)
    with pytest.raises(ValueError, match='weight'):
        LineStacker.makeCoords([1.], [2.], weight=[-1.])


def test_amplitude_weights_need_a_positive_line(workdir):
    coords, names, truth = line_cubes(noise=0.)
    c = truth[1]['cube']
    c.arr *= -1.                                   # an absorption "line" only: max <= 0
    c.write(names[1])
    _io.clear_cache()
    with pytest.warns(UserWarning, match='no positive line amplitude'):
        prep = li._prepare(coords, names, stampsize=8, N_chans=21, fEm=FEM, weighting='1/A')
    w = prep.weights()
    assert np.all(w[1] == 0) and np.all(np.isfinite(w)) and np.all(w[0] > 0)


# -- spectral axes -------------------------------------------------------------

def test_hz_and_ghz_cubes_stack_identically(workdir):
    coords, names, _ = line_cubes(seed=5)
    ref = li.stack(coords, imagenames=names, fEm=FEM, stampsize=10, N_chans=15, outfile='hz.fits')
    _io.clear_cache()
    coords2, names2, _ = line_cubes(seed=5, edit=ghz)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        new = li.stack(coords2, imagenames=names2, fEm=FEM, stampsize=10, N_chans=15, outfile='mixed.fits')
    assert no_grid_warning(rec)
    np.testing.assert_allclose(new[2], ref[2], rtol=1e-6, atol=1e-7)
    from astropy.io import fits
    assert fits.getheader('mixed.fits')['CDELT3'] == pytest.approx(fits.getheader('hz.fits')['CDELT3'])
    # with a GHz cube first, fEm is in GHz
    _io.clear_cache()
    first = [names2[1]] + [n for k, n in enumerate(names2) if k != 1]
    res = li.stack(coords2, imagenames=first, fEm=FEM / 1e9, stampsize=10, N_chans=15, outfile=None)
    np.testing.assert_allclose(res[1], ref[1], rtol=1e-6, atol=1e-7)


def test_frequency_and_wavelength_cubes_are_not_mixed(workdir):
    def wave(i, cube):
        if i == 2:
            cube.header['CTYPE3'], cube.header['CUNIT3'] = 'AWAV', 'Angstrom'
            cube.header['CRVAL3'], cube.header['CDELT3'] = 5000., 1.25
    coords, names, _ = line_cubes(edit=wave)
    with pytest.raises(ValueError, match='different kinds of spectral axis'):
        li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None)


def test_nonlinear_spectral_axis_is_refused(workdir):
    cube = Cube(nchan=30, kind='wave', crval3=5000e-10, cdelt3=1.25e-10)
    cube.header['CTYPE3'], cube.header['CUNIT3'] = 'WAVE-LOG', 'm'
    cube.write('log.fits')
    with pytest.raises(ValueError, match='non-linear spectral axis'):
        li.stack([LineStacker.Coord(20, 20, spec=0.1, image='log.fits', type='pixel')], imagenames=['log.fits'],
                 lEm=4500e-10, stampsize=8, N_chans=5, outfile=None)


def test_velocity_axes_in_m_s_and_km_s(workdir):
    def vrad(i, cube):
        h = cube.header
        h['CTYPE3'] = 'VRAD'
        h['CRPIX3'] = 1.
        if i % 2:
            h['CRVAL3'], h['CDELT3'], h['CUNIT3'] = -100000., 10000., 'm/s'
        else:
            h['CRVAL3'], h['CDELT3'], h['CUNIT3'] = -100., 10., 'km/s'
    coords, names, truth = line_cubes(edit=vrad, noise=0.)
    for c, t in zip(coords, truth):
        c.spec = -100. + 10. * t['chan']            # km/s: the unit of the first cube
    res = li.stack(coords, imagenames=names, spectralMethod='obsSpec', stampsize=8, N_chans=11, outfile='v.fits')
    assert np.argmax(res[1]) == 5
    from astropy.io import fits
    h = fits.getheader('v.fits')
    assert h['CUNIT3'].strip() in ('km/s', 'km s-1') and h['CDELT3'] == pytest.approx(10.)
    with pytest.raises(ValueError, match='frequency or wavelength'):
        li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None)


def test_different_spectral_frames_warn(workdir):
    def frames(i, cube):
        cube.header['SPECSYS'] = 'LSRK' if i else 'TOPOCENT'
    coords, names, _ = line_cubes(edit=frames)
    with pytest.warns(UserWarning, match='spectral reference frames'):
        li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None)


# -- sky grids -----------------------------------------------------------------

def rotate(angle_deg):
    def edit(i, cube):
        a = np.deg2rad(angle_deg)
        h = cube.header
        h['PC1_1'], h['PC1_2'], h['PC2_1'], h['PC2_2'] = np.cos(a), -np.sin(a), np.sin(a), np.cos(a)
    return edit


def test_rotated_grid_is_kept_in_the_output(workdir):
    from astropy.io import fits
    coords, names, _ = line_cubes(edit=rotate(30.), noise=0.)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        res = li.stack(coords, imagenames=names, fEm=FEM, stampsize=12, N_chans=11, outfile='rot.fits')
    assert no_grid_warning(rec)
    assert np.argmax(res[1]) == 5
    h = fits.getheader('rot.fits')
    a = np.deg2rad(30.)
    np.testing.assert_allclose([h['PC1_1'], h['PC1_2'], h['PC2_1'], h['PC2_2']],
                               [np.cos(a), -np.sin(a), np.sin(a), np.cos(a)], atol=1e-12)
    assert h['CDELT1'] == pytest.approx(-0.25 / 3600) and h['CDELT2'] == pytest.approx(0.25 / 3600)
    # an unrotated stack has no PC keywords, as before
    _io.clear_cache()
    coords, names, _ = line_cubes(noise=0.)
    li.stack(coords, imagenames=names, fEm=FEM, stampsize=12, N_chans=11, outfile='plain.fits')
    assert 'PC1_2' not in fits.getheader('plain.fits')


@pytest.mark.parametrize('edit,match', [(lambda i, c: c.header.update(CDELT1=-0.3 / 3600, CDELT2=0.3 / 3600)
                                         if i == 1 else None, 'pixels instead of'),
                                        (lambda i, c: rotate(20.)(i, c) if i == 2 else None, 'rotated by'),
                                        (lambda i, c: c.header.update(CDELT1=0.25 / 3600) if i == 3 else None,
                                         'mirrored')])
def test_different_grids_warn(workdir, edit, match):
    coords, names, _ = line_cubes(edit=edit)
    with pytest.warns(UserWarning, match=match):
        li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None)


def test_primary_beam_offsets_on_a_rotated_grid(workdir):
    coords, names, truth = line_cubes(edit=rotate(35.), noise=0.)
    pb = li.GaussianPrimaryBeam(dish=50.)      # FWHM ~5", offsets up to ~2"
    prep = li._prepare(coords, names, stampsize=8, N_chans=5, fEm=FEM, weighting='pb', primarybeam=pb)
    for c, t, wrow in zip(prep.coords, truth, prep.weights()):
        cube = t['cube']
        from astropy.wcs import WCS
        centre = WCS(cube.header).celestial.all_pix2world([[cube.header['CRPIX1'] - 1, cube.header['CRPIX2'] - 1]], 0)[0]
        from astropy.coordinates import SkyCoord
        sep = cube.sky(t['x'], t['y']).separation(SkyCoord(*centre, unit='deg')).rad
        nu = FEM / (1 + t['z'])
        expected = np.exp(-4 * np.log(2) * sep**2 / pb.fwhm(nu)**2)**2
        assert wrow[2] == pytest.approx(expected, rel=1e-6)


def test_primary_beam_images(workdir):
    from astropy.io import fits
    coords, names, truth = line_cubes(n=2, noise=0.)
    pbs = []
    for i, t in enumerate(truth):
        h = t['cube'].header.copy()
        X, Y = np.meshgrid(np.arange(40), np.arange(40), indexing='ij')
        pbmap = 0.5 + 0.01 * X + 0.001 * Y + 0.1 * i          # (x, y)
        arr = np.repeat(pbmap.T[None, None], 50, axis=1)        # (stokes, chan, y, x)
        fits.writeto('pb%d.fits' % i, arr.astype('f4'), h, overwrite=True)
        pbs.append('pb%d.fits' % i)
    expected = [(0.5 + 0.01 * t['x'] + 0.001 * t['y'] + 0.1 * i)**2 for i, t in enumerate(truth)]
    for primarybeam in (pbs, dict(zip(names, pbs))):
        prep = li._prepare(coords, names, stampsize=8, N_chans=5, fEm=FEM, weighting='pb', primarybeam=primarybeam)
        np.testing.assert_allclose(prep.weights()[:, 2], expected, rtol=1e-5)
    with pytest.raises(ValueError, match='one per cube'):
        li._prepare(coords, names, stampsize=8, N_chans=5, fEm=FEM, weighting='pb', primarybeam=pbs[0])


def test_beam_position_angles_are_averaged_as_axes(workdir):
    from astropy.io import fits

    def pa(i, cube):
        cube.header['BPA'] = (179., 1., 178., 2.)[i]
    coords, names, _ = line_cubes(edit=pa)
    li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile='pa.fits')
    bpa = fits.getheader('pa.fits')['BPA']
    assert abs(((bpa + 90.) % 180.) - 90.) < 1e-9
    assert li._axial_mean_deg([30., 30.]) == pytest.approx(30.)
    assert li._axial_mean_deg([-60.]) == pytest.approx(-60.)


# -- regridFromZ (cube writer) ---------------------------------------------------

def test_regridded_cubes_follow_the_documented_contract(workdir):
    from astropy.io import fits
    from astropy.wcs import WCS
    rng = np.random.default_rng(2)
    names, coords, sky = [], [], []
    for i in range(4):
        z = 5.6 + 0.3 * i
        cube = Cube(nx=40, ny=40, nchan=80, crval3=FEM / (1 + z) + 40 * 7.6e7)
        chan = cube.channel_of(FEM / (1 + z))
        x, y = float(rng.integers(14, 26)), float(rng.integers(14, 26))
        cube.add_line(x, y, chan, 2., 3e9 / (1 + z) / 7.6e7)     # same rest-frame width
        cube.write('g%d.fits' % i)
        s = cube.sky(x, y)
        sky.append(s)
        coords.append(LineStacker.Coord(s.ra, s.dec, spec=z, image='g%d.fits' % i))
        names.append('g%d.fits' % i)
    new_names, new_coords = A.regridFromZ(coords, names, stampsize=12, fEm=FEM, regridMethod='scaleToMax',
                                          outdir='regridded')
    zmax = 5.6 + 0.9
    for name, c, s in zip(new_names, new_coords, sky):
        h = fits.getheader(name)
        if name.startswith('regridded'):
            assert c.spec == 0
            # rest-frame channels: same width as the reference source's
            assert abs(h['CDELT3']) == pytest.approx(7.6e7 * (1 + zmax), rel=1e-9)
            # the stamp keeps the sky position of its source
            w = WCS(h).celestial
            ra, dec = w.all_pix2world([[c.ra, c.dec]], 0)[0]
            assert abs(ra - s.ra.deg) < 1e-7 and abs(dec - s.dec.deg) < 1e-7
        else:
            assert c.spec == pytest.approx(zmax)
    res = li.stack(new_coords, imagenames=new_names, fEm=FEM, stampsize=8, N_chans=21, outfile='rz.fits')
    assert np.argmax(res[1]) == 10
    h = fits.getheader('rz.fits')
    assert abs(h['CDELT3']) == pytest.approx(7.6e7 * (1 + zmax), rel=1e-9)

"""Fixes from the second code review: non-linear spectral axes, intensity
units, sky geometry, regridFromZ pixels, 1D inputs, fits, validation,
caching, parallel environment and CASA's log."""
import os
import subprocess
import sys
import threading
import warnings

import numpy as np
import pytest

from conftest import Cube, FEM, needs_casa
from test_weights_wcs import line_cubes

import LineStacker
import LineStacker.OneD_Stacker as OneD
import LineStacker.analysisTools as A
import LineStacker.image as li
import LineStacker.tools.fit as fitTools
from LineStacker import _io
from LineStacker.tools import parallel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def stack(coords, names, **kw):
    kw.setdefault('fEm', FEM)
    kw.setdefault('stampsize', 8)
    kw.setdefault('N_chans', 5)
    kw.setdefault('outfile', None)
    return li.stack(coords, imagenames=names, **kw)


# -- spectral axes ---------------------------------------------------------------

def test_nonlinear_axes_are_refused_in_every_line_mode(workdir):
    from astropy.io import fits
    cube = Cube(nchan=30, kind='wave', crval3=5000e-10, cdelt3=1.25e-10)
    cube.header['CTYPE3'], cube.header['CUNIT3'] = 'WAVE-LOG', 'm'
    cube.write('log.fits')
    c = [LineStacker.Coord(20, 20, spec=10, image='log.fits', type='pixel')]
    with pytest.raises(ValueError, match='non-linear spectral axis'):
        li.stack(c, imagenames=['log.fits'], spectralMethod='channel', stampsize=8, N_chans=5, outfile=None)
    # continuum: the stack is fine, but no spectral axis is claimed
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        li.stack(c, imagenames=['log.fits'], continuum=True, stampsize=8, outfile='cont.fits')
    assert 'CTYPE3' not in fits.getheader('cont.fits')


def test_missing_specsys_warns(workdir):
    def frames(i, cube):
        if i:
            cube.header['SPECSYS'] = 'LSRK'
    coords, names, _ = line_cubes(edit=frames)
    with pytest.warns(UserWarning, match='not given'):
        stack(coords, names)


# -- intensity units ---------------------------------------------------------------

def test_intensity_units_scale_exactly(workdir):
    coords, names, truth = line_cubes(seed=9, noise=0.)
    ref = stack(coords, names)[2]
    # write cube 1 again in mJy/beam: same physical values
    c = truth[1]['cube']
    c.arr *= 1000.
    c.header['BUNIT'] = 'mJy/beam'
    c.write(names[1])
    _io.clear_cache()
    np.testing.assert_allclose(stack(coords, names, outfile='u.fits')[2], ref, rtol=1e-6, atol=1e-9)
    from astropy.io import fits
    assert fits.getheader('u.fits')['BUNIT'] == 'Jy/beam'      # the first cube's unit


@pytest.mark.parametrize('unit', ['K', 'Jy/pixel'])
def test_incompatible_intensity_units_are_refused(workdir, unit):
    def other(i, cube):
        if i == 2:
            cube.header['BUNIT'] = unit
    coords, names, _ = line_cubes(edit=other)
    with pytest.raises(ValueError, match='intensity units cannot be stacked'):
        stack(coords, names)


def test_missing_bunit_warns(workdir):
    def nounit(i, cube):
        if i == 1:
            del cube.header['BUNIT']
    coords, names, _ = line_cubes(edit=nounit)
    with pytest.warns(UserWarning, match='no BUNIT'):
        stack(coords, names)


# -- sky geometry ------------------------------------------------------------------

def test_projection_and_frame_differences_warn(workdir):
    def tan(i, cube):
        if i == 1:
            cube.header['CTYPE1'], cube.header['CTYPE2'] = 'RA---TAN', 'DEC--TAN'
        if i == 2:
            cube.header['RADESYS'] = 'FK5'
            cube.header['EQUINOX'] = 2000.
    coords, names, _ = line_cubes(edit=tan)
    with pytest.warns(UserWarning) as rec:
        stack(coords, names)
    text = ' '.join(str(w.message) for w in rec)
    assert 'TAN projection instead of SIN' in text and 'frame fk5' in text


# -- regridFromZ ----------------------------------------------------------------------

def test_regridded_coordinates_keep_fractional_pixels(workdir):
    from astropy.io import fits
    from astropy.wcs import WCS
    names, coords, sky = [], [], []
    for i in range(3):
        z = 5.6 + 0.4 * i
        cube = Cube(nx=40, ny=40, nchan=80, crval3=FEM / (1 + z) + 40 * 7.6e7)
        x, y = 20.3 + i * 0.2, 19.6 - i * 0.3                 # fractional positions
        cube.add_line(x, y, cube.channel_of(FEM / (1 + z)), 2., 3.)
        cube.write('f%d.fits' % i)
        s = cube.sky(x, y)
        sky.append(s)
        coords.append(LineStacker.Coord(s.ra, s.dec, spec=z, image='f%d.fits' % i))
        names.append('f%d.fits' % i)
    new_names, new_coords = A.regridFromZ(coords, names, stampsize=12, fEm=FEM, regridMethod='scaleToMax',
                                          outdir='rz', verbose=False)
    for name, c, s in zip(new_names, new_coords, sky):
        if not name.startswith('rz'):
            continue
        assert c.ra != int(c.ra) or c.dec != int(c.dec)
        ra, dec = WCS(fits.getheader(name)).celestial.all_pix2world([[c.ra, c.dec]], 0)[0]
        assert abs(ra - s.ra.deg) < 1e-8 and abs(dec - s.dec.deg) < 1e-8


# -- 1D inputs and fits -----------------------------------------------------------------

def test_infinite_fluxes_are_left_out():
    a = OneD.Spectrum(amp=[1., 1., np.inf, 1., 1.], centerIndex=2)
    b = OneD.Spectrum(amp=[3., 3., 3., 3., -np.inf], centerIndex=2)
    for method, expected in (('mean', [2., 2., 3., 2., 1.]), ('median', [2., 2., 3., 2., 1.])):
        st, n = OneD.Stack([a, b], method=method)
        np.testing.assert_allclose(st, expected)
    np.testing.assert_allclose(OneD.Stack([a, b])[1], [2, 2, 1, 2, 1])


def test_spectrum_shape_checks():
    with pytest.raises(ValueError, match='at least 2'):
        OneD.Spectrum(amp=[1.], centerIndex=0)
    with pytest.raises(ValueError, match='5 fluxes but 4 velocities'):
        OneD.Spectrum(amp=np.ones(5), velocities=np.arange(4.), centerIndex=2)
    with pytest.raises(ValueError, match='frequencies'):
        OneD.Spectrum(amp=np.ones(5), frequencies=np.arange(6.), velOrFreq='freq', centerIndex=2)


def test_fits_raise_fiterror_on_short_or_nan_axes():
    for bad in ([], [1.], [1., 2.]):
        with pytest.raises(fitTools.FitError):
            fitTools.DoublePeakFit(fctToFit=bad)
    # leading NaN channels (and axis values) do not spoil the default width
    x = np.arange(40.)
    y = 2. * np.exp(-(x - 22.)**2 / (2 * 3.**2))
    f = x.copy()
    y[:2], f[:2] = np.nan, np.nan
    _, p = fitTools.GaussFit(fctToFit=y, fullFreq=f, returnInfos=True)
    np.testing.assert_allclose([p[0], p[1], abs(p[2])], [2., 22., 3.], rtol=1e-6)


# -- validation ---------------------------------------------------------------------------

def test_subsample_sizes_must_be_integers():
    sp = [OneD.Spectrum(amp=np.random.default_rng(k).normal(size=20), centerIndex=10) for k in range(6)]
    with pytest.raises(ValueError, match='minSize must be an integer'):
        A.subsample_OneD(sp, nRandom=5, minSize=2.7, verbose=False)
    with pytest.raises(ValueError, match='maxSize must be an integer'):
        A.subsample_OneD(sp, nRandom=5, maxSize=4.5, verbose=False)


def test_randomize_coords_accepts_one_image_path(workdir):
    Cube(nx=60, ny=60, nchan=3).write('one.fits')
    c = [LineStacker.Coord(30., 30., spec=1, image='one.fits', type='pixel')]
    out = LineStacker.randomizeCoords(c, beam=np.deg2rad(1. / 3600), maxBeamRange=3, imagenames='one.fits', rng=1)
    assert len(out) == 1 and out[0].ra != 30.
    out = A.randomizeCoords(c, lowerLimit=np.deg2rad(1. / 3600), upperLimit=np.deg2rad(2. / 3600),
                            imagenames='one.fits', seed=1)
    r = np.hypot(out[0].ra - 30., out[0].dec - 30.) * 0.25
    assert 1. - 1e-9 <= r <= 2. + 1e-9
    with pytest.raises(ValueError, match='beam'):
        LineStacker.randomizeCoords(c, beam=0., imagenames='one.fits')
    with pytest.raises(ValueError, match='maxBeamRange'):
        LineStacker.randomizeCoords(c, beam=1e-5, maxBeamRange=0.5, imagenames='one.fits')
    with pytest.raises(ValueError, match='lowerLimit <= upperLimit'):
        A.randomizeCoords(c, lowerLimit=2e-5, upperLimit=1e-5, imagenames='one.fits')
    with pytest.raises(ValueError, match='beam'):
        A.randomizeCoords(c, imagenames='one.fits')


def test_congrid_minusone_singleton():
    with pytest.raises(ValueError, match='minusone'):
        A.congrid(np.ones((4, 4)), [1, 4], minusone=True)
    assert A.congrid(np.ones((4, 4)), [2, 3], minusone=True).shape == (2, 3)


@pytest.mark.parametrize('bad', [None, np.nan, 'nan', 'None'])
def test_catalogue_missing_image_names(workdir, bad):
    from astropy.table import Table
    t = Table({'ra': [53.1, 53.2], 'dec': [-27.8, -27.8], 'z': [6.5, 6.6],
               'image': np.array(['a.fits', bad], dtype=object)})
    with pytest.raises(ValueError, match='empty in rows'):
        LineStacker.readCatalog(t, verbose=False)


def test_primary_beam_on_another_spectral_axis_is_refused(workdir):
    from astropy.io import fits
    coords, names, truth = line_cubes(n=1, noise=0.)
    h = truth[0]['cube'].header.copy()
    h['CTYPE3'], h['CUNIT3'], h['CRVAL3'], h['CDELT3'] = 'WAVE', 'm', 1e-3, 1e-6
    fits.writeto('pbw.fits', np.ones((1, 50, 40, 40), 'f4'), h, overwrite=True)
    with pytest.raises(ValueError, match='primary-beam image'):
        li._prepare(coords, names, stampsize=8, N_chans=5, fEm=FEM, weighting='pb', primarybeam='pbw.fits')


# -- caching -------------------------------------------------------------------------------

def test_cache_check_can_be_turned_off(workdir, monkeypatch):
    import time
    cube = Cube(nchan=4)
    cube.arr[:] = 1.
    cube.write('k.fits')
    assert _io.read_cube('k.fits')[0].max() == 1.
    monkeypatch.setattr(LineStacker, 'cache_check', False)
    time.sleep(0.01)
    cube.arr[:] = 2.
    cube.write('k.fits')
    assert _io.read_cube('k.fits')[0].max() == 1.          # trusted until cleared
    li.clear_cache()
    assert _io.read_cube('k.fits')[0].max() == 2.
    monkeypatch.setattr(LineStacker, 'cache_check', True)
    time.sleep(0.01)
    cube.arr[:] = 3.
    cube.write('k.fits')
    assert _io.read_cube('k.fits')[0].max() == 3.


# -- parallel environment ---------------------------------------------------------------------

def _square(sl):
    return [k * k for k in range(sl.start, sl.stop)]


def test_concurrent_process_jobs_restore_the_environment(monkeypatch):
    monkeypatch.setenv('OMP_NUM_THREADS', '7')
    results, errors = [], []

    def job():
        try:
            out = parallel.map_chunks(_square, [slice(0, 3), slice(3, 6)], n_jobs=2, backend='process',
                                      progress=False)
            results.append([v for chunk in out for v in chunk])
        except Exception as err:
            errors.append(err)
    threads = [threading.Thread(target=job) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors and results == [[k * k for k in range(6)]] * 2
    assert os.environ['OMP_NUM_THREADS'] == '7'


# -- CASA's log ----------------------------------------------------------------------------------

def _run(code, cwd):
    out = subprocess.run([sys.executable, '-c', 'import sys; sys.path.insert(0, {0!r}); '.format(ROOT) + code],
                         capture_output=True, text=True, cwd=cwd, timeout=300)
    assert out.returncode == 0, out.stderr
    return out


def test_fits_use_starts_no_casa(workdir):
    Cube(nchan=6).write('n.fits')
    os.makedirs('run')
    out = _run('import LineStacker; a = "casatools" in sys.modules; '
               'from LineStacker import _io; _io.read_cube("../n.fits"); '
               'print(a, "casatools" in sys.modules)', os.path.join(str(workdir), 'run'))
    assert out.stdout.split() == ['False', 'False']
    assert not os.listdir('run')


@needs_casa
@pytest.mark.parametrize('setting', ['auto', 'path', 'none', 'casa'])
def test_casa_log_settings(workdir, setting):
    from LineStacker import _casa_compat
    Cube(nchan=4).write('l.image', casa=True)
    os.makedirs('run')
    run = os.path.join(str(workdir), 'run')
    log = os.path.join(str(workdir), 'my.log')
    pre = {'auto': '', 'path': 'LineStacker.casa_logfile = {0!r}; '.format(log),
           'none': 'LineStacker.casa_logfile = None; ', 'casa': "LineStacker.casa_logfile = 'casa'; "}[setting]
    default = _casa_compat.default_logfile()
    before = os.path.getsize(default) if os.path.exists(default) else 0
    out = _run('import LineStacker; {0}from LineStacker import _io; '
               'print("casatools" in sys.modules); _io.clear_cache(); print(_io.read_cube("../l.image")[0].shape); '
               'print("casatools" in sys.modules)'.format(pre), run)
    lines = out.stdout.strip().splitlines()
    assert lines[0] == 'False' and lines[-1] == 'True'          # CASA started only for the CASA image
    logs_here = [f for f in os.listdir(run) if f.endswith('.log')]
    if setting == 'casa':
        assert logs_here                                          # CASA's own default
        return
    assert not logs_here
    # CASA's tool messages go to the log, not to the terminal
    assert 'INFO' not in out.stdout + out.stderr
    if setting == 'path':
        assert '::casa' in open(log).read()
    elif setting == 'auto':
        assert os.path.getsize(default) > before


# -- blank pixels are NaN ----------------------------------------------------------------------

def test_blank_pixels_stay_nan(workdir):
    coords, names, truth = line_cubes(n=3, noise=0.05)
    c = truth[0]['cube']
    x, y = int(truth[0]['x']), int(truth[0]['y'])
    c.arr[x, y, :] = np.nan                          # blank spectrum at source 0
    c.arr[:5, :, :] = np.nan                         # blank edge
    c.write(names[0])
    _io.clear_cache()
    cube = li.readCube(names[0])
    assert np.isnan(cube[x, y, 0, :]).all() and np.isnan(cube[:5]).all()
    # the mean leaves blank pixels out: the centre is the mean of the 2 others
    res = stack(coords, names, N_chans=11, outfile='nan.fits')
    pix = LineStacker.getPixelCoords(coords, names)
    others = []
    for k in (1, 2):
        arr = truth[k]['cube'].arr
        cc = int(round(truth[k]['chan']))
        others.append(arr[int(round(pix[k].ra)), int(round(pix[k].dec)), cc - 5:cc + 6])
    np.testing.assert_allclose(res[1], np.mean(others, axis=0), rtol=1e-5, atol=1e-6)
    # a stack channel no source covers is NaN, not 0
    res = stack(coords, names, N_chans=200, outfile=None)
    assert np.isnan(res[1][0]) and np.isnan(res[1][-1]) and np.isfinite(res[1][100])


def test_collapse_keeps_fully_blank_pixels_nan(workdir):
    cube = Cube(nchan=6)
    cube.arr[:] = 1.
    cube.arr[3, 4, :] = np.nan
    cube.arr[5, 5, :2] = np.nan
    cube.write('c.fits')
    col = li.readCube('c.fits', collapse='sum')
    assert np.isnan(col[3, 4, 0, 0]) and col[5, 5, 0, 0] == 4. and col[0, 0, 0, 0] == 6.


def test_zero_pixels_are_data(workdir):
    s = 10
    stamp = np.random.default_rng(2).normal(0, 1, (s, s, 1, 3))
    stamp[2:4] = 0.                                     # real zeros
    keep = np.ones((s, s), bool)
    assert li._stamp_sigma(stamp, keep) == pytest.approx(np.std(stamp))
    stamp[2:4] = np.nan
    assert li._stamp_sigma(stamp, keep) == pytest.approx(np.nanstd(stamp))


# -- third review: distortions, PB offsets, small checks -----------------------------------------

def sip_images(coeffs, seed=4):
    """2D images (SIP needs a 2D WCS) with a source each; coords, names, truth."""
    rng = np.random.default_rng(seed)
    coords, names, truth = [], [], []
    for i, coeff in enumerate(coeffs):
        img = Cube(nchan=None)
        h = img.header
        h['CTYPE1'], h['CTYPE2'] = 'RA---TAN-SIP', 'DEC--TAN-SIP'
        h['A_ORDER'] = h['B_ORDER'] = 2
        h['A_0_2'], h['B_2_0'] = coeff, coeff
        x, y = float(rng.integers(12, 28)), float(rng.integers(12, 28))
        img.add_line(x, y, 0, 1., 1.)
        img.write('s%d.fits' % i)
        sky = img.sky(x, y)
        coords.append(LineStacker.Coord(sky.ra, sky.dec, image='s%d.fits' % i))
        names.append('s%d.fits' % i)
        truth.append(dict(x=x, y=y, img=img))
    return coords, names, truth


def test_different_sip_distortions_warn(workdir):
    coords, names, _ = sip_images([1e-6, 2e-6, 1e-6])
    with pytest.warns(UserWarning, match='different distortion terms'):
        li.stack(coords, imagenames=names, stampsize=8, outfile=None)
    _io.clear_cache()
    coords, names, _ = sip_images([1e-6, 1e-6, 1e-6])
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        li.stack(coords, imagenames=names, stampsize=8, outfile=None)
    assert not [w for w in rec if 'pixel grid' in str(w.message)]


def test_gaussian_primary_beam_uses_sky_offsets(workdir):
    coords, names, truth = sip_images([3e-5, 3e-5])
    pb = li.GaussianPrimaryBeam(dish=50., nu0=230e9)
    prep = li._prepare(coords, names, stampsize=8, weighting='pb', primarybeam=pb)
    for t, w in zip(truth, prep.weights()):
        h = t['img'].header
        centre = t['img'].sky(h['CRPIX1'] - 1, h['CRPIX2'] - 1)
        sep = t['img'].sky(t['x'], t['y']).separation(centre).rad
        assert w[0] == pytest.approx(np.exp(-4 * np.log(2) * sep**2 / pb.fwhm(230e9)**2)**2, rel=1e-6)


def test_integer_and_range_checks(workdir):
    Cube(nx=30, ny=30, nchan=3).write('r.fits')
    with pytest.raises(ValueError, match='ncoords must be an integer'):
        LineStacker.randomCoords('r.fits', ncoords=2.5)
    c = [LineStacker.Coord(15., 15., spec=1, image='r.fits', type='pixel')]
    with pytest.raises(ValueError, match='maxTries must be an integer'):
        LineStacker.randomizeCoords(c, beam=1e-5, imagenames='r.fits', maxTries=2.5)
    with pytest.raises(ValueError, match='maxBeamRange'):
        LineStacker.randomizeCoords(c, beam=1e-5, imagenames='r.fits', maxBeamRange=np.inf)
    with pytest.raises(ValueError, match='nRandom must be an integer >= 1'):
        A.noise_estimator('r.fits', nRandom=0, plotIt=False, verbose=False)
    coords, names, _ = line_cubes(n=2)
    with pytest.raises(ValueError, match='nRandom must be an integer >= 1'):
        A.stack_estimator(coords, nRandom=0, imagenames=names, fEm=FEM, lowerLimit=1e-5, upperLimit=2e-5)


def test_default_fit_width_ignores_gaps():
    x = np.arange(60.)
    y = 3. * np.exp(-(x - 30.)**2 / (2 * 4.**2))
    y[1] = np.nan                                       # gap between the first valid samples
    _, p = fitTools.GaussFit(fctToFit=y, returnInfos=True)
    np.testing.assert_allclose([p[0], p[1], abs(p[2])], [3., 30., 4.], rtol=1e-6)


def test_spectrum_axis_must_be_usable():
    with pytest.raises(ValueError, match='finite'):
        OneD.Spectrum(amp=np.ones(4), velocities=[0., np.nan, 2., 3.], centerIndex=1)
    with pytest.raises(ValueError, match='increase or decrease'):
        OneD.Spectrum(amp=np.ones(4), velocities=[0., 1., 1., 2.], centerIndex=1)
    OneD.Spectrum(amp=np.ones(4), velocities=[3., 2., 1., 0.], centerIndex=1)       # decreasing is fine

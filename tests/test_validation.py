"""1D velocities, fitting, input validation, coordinate frames, parallel
state and CASA I/O (fixes from the 2.0 code review)."""
import os
import subprocess
import sys
import threading

import numpy as np
import pytest

from conftest import Cube, needs_casa

import LineStacker
import LineStacker.OneD_Stacker as OneD
import LineStacker.analysisTools as A
import LineStacker.tools.fit as fitTools
from LineStacker import _io
from LineStacker.tools import parallel

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
C_KMS = 299792.458


def gauss_spectra(n=6, L=40, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        x = np.arange(L, dtype=float)
        amp = 2. * np.exp(-(x - 20.) ** 2 / 8.) + rng.normal(0, 0.1, L)
        sp = OneD.Spectrum(amp=amp, centerIndex=20, name=str(i))
        sp.z = 1. + 0.1 * i                          # for the regridding checks
        out.append(sp)
    return out


# -- 1D velocities -------------------------------------------------------------

@pytest.mark.parametrize('convention,vline', [('radio', lambda z: C_KMS * z / (1 + z)),
                                              ('optical', lambda z: C_KMS * z),
                                              ('relativistic', lambda z: C_KMS * ((1 + z)**2 - 1) / ((1 + z)**2 + 1))])
def test_velocity_spectra_centred_from_z(convention, vline):
    z, f0 = 0.02, 230.538e9
    v = np.arange(-2000., 12000., 10.)                  # km/s relative to the rest frequency
    amp = np.exp(-(v - vline(z))**2 / (2 * 50.**2))
    sp = OneD.Spectrum(amp=amp, velocities=v, z=z, fEmLine=f0, velConvention=convention)
    assert sp.centerIndex == int(np.argmin(np.abs(v - vline(z))))
    # the frequencies of the channels follow the same convention
    # (to within half a 10 km/s channel)
    assert abs(sp.frequencies[sp.centerIndex] - f0 / (1 + z)) <= f0 * 5.01 / C_KMS
    assert sp.velToFreq(0., 0., f0) == pytest.approx(f0)
    assert sp.freqToVel(sp.velToFreq(1234., 0., f0), 0., f0) == pytest.approx(1234.)


def test_unknown_velocity_convention():
    with pytest.raises(ValueError, match='velConvention'):
        OneD.Spectrum(amp=np.ones(5), velConvention='lsr')


def test_nan_channels_carry_no_weight():
    a = OneD.Spectrum(amp=[1., 1., np.nan, 1., 1.], centerIndex=2)
    b = OneD.Spectrum(amp=[3., 3., 3., 3., 3.], centerIndex=2)
    st, n = OneD.Stack([a, b])
    np.testing.assert_allclose(st, [2., 2., 3., 2., 2.])
    np.testing.assert_allclose(n, [2, 2, 1, 2, 2])


def test_oned_weights_must_be_valid():
    good = OneD.Spectrum(amp=np.ones(5), centerIndex=2)
    for w in (-1., np.nan):
        bad = OneD.Spectrum(amp=np.ones(5), centerIndex=2, weights=w)
        with pytest.raises(ValueError, match='weights must be finite and >= 0'):
            OneD.Stack([good, bad])
    with pytest.warns(UserWarning, match='1/A'):
        sp = OneD.Spectrum(amp=-np.ones(5), centerIndex=2, weights='1/A')
    assert sp.weights == 0


def test_empty_and_invalid_inputs():
    with pytest.raises(ValueError, match='no spectra'):
        OneD.Stack([])
    with pytest.raises(ValueError, match='no spectra'):
        A.bootstraping_OneD([])
    sp = gauss_spectra()
    with pytest.raises(ValueError, match='nRandom'):
        A.bootstraping_OneD(sp, nRandom=-1, verbose=False)
    with pytest.raises(ValueError, match='nRandom'):
        A.subsample_OneD(sp, nRandom=2.5, verbose=False)
    for kw in (dict(minSize=0), dict(minSize=4, maxSize=3), dict(maxSize=99)):
        with pytest.raises(ValueError, match='minSize'):
            A.subsample_OneD(sp, nRandom=5, verbose=False, **kw)
    with pytest.raises(ValueError, match='regridMethod'):
        A.regridFromZ1D(sp, regridMethod='scaleToMiddle', verbose=False)
    with pytest.raises(ValueError, match='regridMethod'):
        OneD.Stack(sp, regridFromZ=True, regridMethod='nope', verbose=False)
    with pytest.raises(ValueError, match='z of every spectrum'):
        A.regridFromZ1D([OneD.Spectrum(amp=np.ones(5), centerIndex=2)], verbose=False)
    with pytest.raises(ValueError, match='congrid'):
        A.congrid(np.ones((4, 4)), [0, 4])


# -- fitting -----------------------------------------------------------------

def test_gauss_fit_ignores_nan_and_reports_failures():
    x = np.arange(60.)
    y = 3. * np.exp(-(x - 31.)**2 / (2 * 4.**2))
    y[[5, 30, 44]] = np.nan
    model, p = fitTools.GaussFit(fctToFit=y, returnInfos=True)
    np.testing.assert_allclose([p[0], p[1], abs(p[2])], [3., 31., 4.], rtol=1e-6)
    assert model.shape == y.shape and np.all(np.isfinite(model))
    assert issubclass(fitTools.FitError, RuntimeError)
    for bad in (np.ones(20), [1., 2.], np.full(10, np.nan)):
        with pytest.raises(fitTools.FitError):
            fitTools.GaussFit(fctToFit=bad)
    # callers turn failures into NaN / 0 instead of crashing
    assert A.maximizeSNR(np.ones(30)) == 0


def test_double_peak_fit_on_an_offset_axis():
    f = np.linspace(100., 220., 121)
    truth = fitTools.DoubleGaussPeak(f, amp=2., dist=30., sigma=6.)
    peaks = f[[i for i in range(1, 120) if truth[i] > truth[i - 1] and truth[i] > truth[i + 1]]]
    np.testing.assert_allclose(peaks, [145., 175.], atol=1.5)          # around the middle, 160
    _, p = fitTools.DoublePeakFit(fctToFit=truth, fullFreq=f, returnInfos=True)
    np.testing.assert_allclose([p[0], abs(p[1]), abs(p[2])], [2., 30., 6.], rtol=1e-4)


# -- coordinates and catalogues ------------------------------------------------

def test_mixed_frames_are_converted(workdir):
    from astropy.coordinates import SkyCoord
    cube = Cube(nchan=5)
    cube.write('f.fits')
    icrs = [cube.sky(10, 12), cube.sky(25, 30)]
    fk5 = icrs[1].transform_to('fk5')
    coords = [LineStacker.Coord(icrs[0].ra, icrs[0].dec, spec=0, image='f.fits'),
              LineStacker.Coord(fk5.ra, fk5.dec, spec=0, image='f.fits', frame='fk5')]
    pix = LineStacker.getPixelCoords(coords, ['f.fits'])
    np.testing.assert_allclose([pix[0].ra, pix[0].dec, pix[1].ra, pix[1].dec], [10, 12, 25, 30], atol=1e-6)
    # first coordinate in FK5: the ICRS one is converted too
    pix = LineStacker.getPixelCoords(coords[::-1], ['f.fits'])
    np.testing.assert_allclose([pix[0].ra, pix[0].dec, pix[1].ra, pix[1].dec], [25, 30, 10, 12], atol=1e-6)
    with pytest.raises(ValueError, match='RA/Dec'):
        LineStacker.makeCoords([1.], [2.], frame='galactic')
    gal = SkyCoord(l=[10.], b=[20.], unit='deg', frame='galactic')
    c = LineStacker.makeCoords(gal)
    assert c[0].frame == 'icrs' and c[0].ra.deg == pytest.approx(gal.icrs.ra.deg[0])


def test_catalogue_blanks_and_bad_values(workdir):
    with open('blank.csv', 'w') as fh:
        fh.write('ra,dec,z,image\n53.1,-27.8,6.5,a.fits\n53.2,-27.8,6.6,  \n')
    with pytest.raises(ValueError, match='image'):
        LineStacker.readCatalog('blank.csv', verbose=False)
    for col, val, match in (('z', 'nan', 'non-finite'), ('weight', '-2', 'negative'), ('weight', 'inf', 'non-finite')):
        with open('bad.csv', 'w') as fh:
            fh.write('ra,dec,z,weight\n53.1,-27.8,6.5,1\n53.2,-27.8,{0},{1}\n'.format(
                val if col == 'z' else 6.6, val if col == 'weight' else 1))
        with pytest.raises(ValueError, match=match):
            LineStacker.readCatalog('bad.csv', verbose=False)
    with pytest.raises(ValueError, match='z is'):
        LineStacker.makeCoords([1., 2.], [3., 4.], z=[0.1, np.nan], verbose=False)
    with pytest.raises(ValueError, match='empty'):
        LineStacker.makeCoords([1.], [3.], z=[0.1], image=['  '], verbose=False)


def test_catalogue_read_errors_are_kept(workdir):
    with open('broken.fits', 'wb') as fh:
        fh.write(b'not a fits file')
    with pytest.raises(Exception) as err:
        LineStacker.readCatalog('broken.fits', verbose=False)
    assert 'as a text table' not in str(err.value)        # FITS: astropy's own error, no text fallback
    with open('broken.txt', 'w') as fh:
        fh.write('ra dec\n1 2 3 4\n"\n')
    with pytest.raises(ValueError, match='as a text table'):
        LineStacker.readCatalog('broken.txt', verbose=False)


def test_random_coords_validation(workdir):
    Cube(nx=30, ny=20, nchan=4).write('r.fits')
    out = LineStacker.randomCoords('r.fits', ncoords=5, margin=3, rng=1)
    assert len(out) == 5 and all(o.image == 'r.fits' for o in out)
    assert all(3 <= o.ra <= 26 and 3 <= o.dec <= 16 for o in out)
    with pytest.raises(ValueError, match='at least one image'):
        LineStacker.randomCoords([], ncoords=3)
    with pytest.raises(ValueError, match='margin'):
        LineStacker.randomCoords(['r.fits'], ncoords=3, margin=10)
    with pytest.raises(ValueError, match='ncoords'):
        LineStacker.randomCoords(['r.fits'], ncoords=-1)


def test_derived_names_of_compressed_fits():
    assert A._derived_name('d/cube.fits.gz', '_X') == 'd/cube_X.fits'
    assert A._derived_name('cube.FTS.gz', '_X') == 'cube_X.fits'
    assert A._derived_name('cube.fits', '_X', outdir='o') == os.path.join('o', 'cube_X.fits')
    assert A._derived_name('cube', '_X') == 'cube_X.fits'


def test_rebin_width_validation(workdir):
    c = Cube(nchan=20)
    c.write('a.fits')
    c.write('b.fits')
    coords = LineStacker.makeCoords([20, 20], [20, 20], z=[0, 0], image=['a.fits', 'b.fits'], type='pixel',
                                    verbose=False)
    for widths in ([1., 0.], [1.], [1., np.nan]):
        with pytest.raises(ValueError):
            A.rebin_CubesSpectra(coords, widths=widths, verbose=False)


def test_noise_estimator_with_cancelling_pixels(workdir):
    cube = Cube(nx=10, ny=10, nchan=2)
    cube.arr[:] = np.where(np.indices((10, 10, 2)).sum(axis=0) % 2, 1., -1.)   # sums to 0
    cube.write('n.fits')
    sigma = A.noise_estimator('n.fits', nRandom=200, maskradius=2, plotIt=False, seed=1, verbose=False)
    assert sigma == pytest.approx(1., rel=0.1)


# -- parallel ------------------------------------------------------------------

def _scaled(sl):
    return [parallel.get_shared('factor') * k for k in range(sl.start, sl.stop)]


def test_concurrent_thread_jobs_keep_their_inputs():
    results, errors = {}, []

    def job(factor):
        try:
            for _ in range(20):
                out = parallel.map_chunks(_scaled, [slice(i, i + 5) for i in range(0, 50, 5)], n_jobs=3,
                                          backend='thread', shared={'factor': factor}, progress=False)
                got = [v for chunk in out for v in chunk]
                assert got == [factor * k for k in range(50)]
            results[factor] = True
        except Exception as err:                     # reported below
            errors.append(err)
    threads = [threading.Thread(target=job, args=(f,)) for f in (1, 7, 13)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not errors and sorted(results) == [1, 7, 13]
    with pytest.raises(RuntimeError, match='map_chunks'):
        parallel.get_shared('factor')


def test_serial_and_parallel_backends_agree():
    sp = gauss_spectra(n=10)
    ref = A.bootstraping_OneD(sp, nRandom=60, seed=3, save='ampAndWidth', verbose=False)
    for kw in (dict(n_jobs=2, backend='thread'), dict(n_jobs=2, backend='process')):
        out = A.bootstraping_OneD(sp, nRandom=60, seed=3, save='ampAndWidth', verbose=False, **kw)
        np.testing.assert_allclose(np.array(out[1:], float), np.array(ref[1:], float), rtol=1e-12)


# -- CASA I/O ------------------------------------------------------------------

def test_fits_only_use_never_imports_casa(workdir):
    Cube(nchan=8).write('only.fits')
    code = ('import sys; sys.path.insert(0, {0!r}); import LineStacker, LineStacker.image as li; '
            'c = LineStacker.makeCoords([20], [20], z=[0], image=["only.fits"], type="pixel", verbose=False); '
            'li.stack(c, spectralMethod="channel", stampsize=4, N_chans=3, outfile="o.fits", verbose=False); '
            'print("casatools" in sys.modules, "casatasks" in sys.modules)').format(ROOT)
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, cwd=str(workdir), timeout=300)
    assert out.stdout.split() == ['False', 'False'], out.stderr
    assert not [f for f in os.listdir(workdir) if f.startswith('casa-') and f.endswith('.log')]


@needs_casa
def test_modified_casa_image_is_read_again(workdir):
    import time
    cube = Cube(nchan=6)
    cube.arr[:] = 1.
    cube.write('m.image', casa=True)
    a, _ = _io.read_cube('m.image')
    assert float(a.max()) == 1.
    time.sleep(0.01)
    from LineStacker._casa_compat import new_tool
    ia = new_tool('image')
    ia.open('m.image')
    ia.putchunk(ia.getchunk() * 3.)
    ia.done()
    b, _ = _io.read_cube('m.image')
    assert float(b.max()) == 3.
    # the stale entry is replaced, not kept next to the new one
    assert len([k for k in _io._data_cache if k[0].endswith('m.image')]) == 1


@needs_casa
def test_repeated_casa_header_reads(workdir):
    Cube(nchan=6).write('h.image', casa=True)
    for _ in range(30):
        _io.clear_cache()
        assert _io.read_info('h.image').nchan == 6
        assert _io.n_channels('h.image') == 6


@needs_casa
@pytest.mark.parametrize('setting', ['path', 'auto'])
def test_casa_log_stays_out_of_the_working_directory(workdir, setting):
    Cube(nchan=4).write('l.image', casa=True)
    os.makedirs('run')
    log = os.path.join(str(workdir), 'my.log')
    pre = 'LineStacker.casa_logfile = {0!r}; '.format(log) if setting == 'path' else ''
    code = ('import sys; sys.path.insert(0, {0!r}); import LineStacker; {1}'
            'from LineStacker import _io, _casa_compat; print(_io.read_info("../l.image").nchan); '
            'print(_casa_compat.default_logfile())').format(ROOT, pre)
    out = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True,
                         cwd=os.path.join(str(workdir), 'run'), timeout=300)
    lines = out.stdout.strip().splitlines()
    assert lines and lines[0] == '4', out.stderr
    assert not [f for f in os.listdir('run') if f.endswith('.log')]
    assert os.path.exists(log if setting == 'path' else lines[1])

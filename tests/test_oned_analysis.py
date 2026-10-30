import os
import warnings

import numpy as np
import pytest

from conftest import Cube, FEM, casa_style
import _reference_oned

import LineStacker
import LineStacker.OneD_Stacker as OneD
import LineStacker.analysisTools as A
import LineStacker.image as li

EXAMPLE = os.environ.get('LINESTACKER_EXAMPLES',
                         os.path.join(os.path.dirname(__file__), '..', 'LineStacker', 'Example'))
warnings.filterwarnings('ignore', category=RuntimeWarning)


def spectra(n=12, seed=0, ragged=True, perchan=False, L=60):
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n):
        Li = int(rng.integers(L - 20, L + 20)) if ragged else L
        v = np.arange(Li) * 10. - 5. * Li
        ci = int(rng.integers(10, Li - 10))
        a = rng.normal(0, 0.5, Li) + rng.uniform(1, 4) * np.exp(-(np.arange(Li) - ci)**2 / 8.)
        w = rng.uniform(0.5, 2, Li) if perchan else float(rng.uniform(0.5, 2))
        out.append((v, a, w, ci))
    return out


def as_spectra(data, listweights=False):
    return [OneD.Spectrum(spectrum=np.array([v, a]), centerIndex=ci, name='s%d' % k,
                       weights=(list(w) if (listweights and np.ndim(w)) else w))
            for k, (v, a, w, ci) in enumerate(data)]


@pytest.mark.parametrize('method', ['mean', 'median'])
@pytest.mark.parametrize('chans', ['full', 25, 24])
def test_stack_matches_reference(method, chans):
    for seed in range(5):
        data = spectra(seed=seed, perchan=(chans == 'full' and seed % 2 == 1))
        new = OneD.Stack(as_spectra(data), chansStack=chans, method=method)
        ref = _reference_oned.Stack(as_spectra(data, listweights=True), chansStack=chans, method=method)
        np.testing.assert_allclose(new[0], ref[0], rtol=1e-12, atol=1e-12)
        np.testing.assert_array_equal(new[1], ref[1])


def test_per_channel_weights_in_fixed_width_mode():
    """The weight of each spectrum must follow its own line centre."""
    data = spectra(n=4, seed=3, perchan=True, ragged=False)
    ims = as_spectra(data)
    st, _ = OneD.Stack(ims, chansStack=11, method='mean')
    expected = []
    for j in range(11):
        vals = [im.amp[im.centerIndex - 5 + j] for im in ims]
        ws = [im.weights[im.centerIndex - 5 + j] for im in ims]
        expected.append(np.average(vals, weights=ws))
    np.testing.assert_allclose(st, expected, rtol=1e-12)


@pytest.mark.parametrize('method', ['mean', 'median'])
@pytest.mark.parametrize('save', ['amp', 'all'])
def test_bootstrap_oned_matches_bruteforce(method, save):
    data = spectra(n=9, seed=7)
    ims = as_spectra(data)
    res = A.bootstraping_OneD(ims, nRandom=40, method=method, save=save, seed=123)
    idx = np.random.default_rng(123).integers(0, 9, size=(40, 9))
    assert len(res) == 41
    for r in range(40):
        st = _reference_oned.Stack([ims[k] for k in idx[r]], method=method)[0]
        if save == 'amp':
            assert res[r + 1] == pytest.approx(np.max(st), rel=1e-12, abs=1e-12)
        else:
            np.testing.assert_allclose(res[r + 1], st, rtol=1e-12, atol=1e-12)


def test_bootstrap_oned_parallel_identical():
    ims = as_spectra(spectra(n=15, seed=2))
    base = A.bootstraping_OneD(ims, nRandom=300, method='median', save='amp', seed=5)
    for backend in ('thread', 'process'):
        par = A.bootstraping_OneD(ims, nRandom=300, method='median', save='amp', seed=5,
                                  n_jobs=3, backend=backend)
        assert par == base
    fit1 = A.bootstraping_OneD(ims, nRandom=30, save='ampAndWidth', seed=1)
    fit3 = A.bootstraping_OneD(ims, nRandom=30, save='ampAndWidth', seed=1, n_jobs=3)
    np.testing.assert_allclose(fit1, fit3)


def test_bootstrap_respects_global_seed():
    ims = as_spectra(spectra(n=8, seed=4))
    np.random.seed(42)
    a = A.bootstraping_OneD(ims, nRandom=20, save='amp')
    np.random.seed(42)
    b = A.bootstraping_OneD(ims, nRandom=20, save='amp')
    assert a == b


def _ampgrade(s):
    return float(np.max(s))


@pytest.mark.parametrize('method', ['mean', 'median'])
def test_subsample_matches_bruteforce(method):
    data = spectra(n=10, seed=9)
    ims = as_spectra(data)
    grades = A.subsample_OneD(ims, nRandom=60, maxTest=_ampgrade, seed=8, method=method)
    rng = np.random.default_rng(8)
    sizes = rng.integers(3, 11, size=60)
    perm = np.argsort(rng.random((60, 10)), axis=1)
    ref = {im.name: 0. for im in ims}
    for r in range(60):
        sub = [ims[k] for k in perm[r, :sizes[r]]]
        g = _ampgrade(_reference_oned.Stack(sub, method=method)[0])
        for im in sub:
            ref[im.name] += g
    for k in ref:
        assert grades[k] == pytest.approx(ref[k], rel=1e-12)
    par = A.subsample_OneD(ims, nRandom=60, maxTest=_ampgrade, seed=8, method=method, n_jobs=2)
    assert par == grades


def test_subsample_with_lambda_falls_back_to_threads():
    ims = as_spectra(spectra(n=6, seed=1))
    with pytest.warns(UserWarning):
        g = A.subsample_OneD(ims, nRandom=20, maxTest=lambda s: float(np.max(s)), seed=3, n_jobs=2)
    assert g == A.subsample_OneD(ims, nRandom=20, maxTest=_ampgrade, seed=3)


def test_congrid_and_regrid1d_numpy2():
    a = np.sin(np.linspace(0, 3, 40))
    for m in ('linear', 'nearest', 'spline'):
        out = A.congrid(a, [30], method=m)
        assert out.shape == (30,)
    assert A.congrid(np.arange(12.).reshape(3, 4), [6, 8], method='neighbour').shape == (6, 8)
    ims = []
    for z in (1.0, 1.5, 2.0):
        v = np.arange(80) * 5.
        a = np.exp(-(np.arange(80) - 40)**2 / (2 * (6 / (1 + z))**2))
        ims.append(OneD.Spectrum(spectrum=np.array([v, a]), z=z, fEmLine=1e10, centerIndex=40,
                              frequencies=np.arange(80) * 1e6 + 1e10 / (1 + z) - 40e6, velOrFreq='freq'))
    st = OneD.Stack(ims, regridFromZ=True, regridMethod='scaleToMin', center=40)
    # channels beyond every regridded spectrum have no data: NaN, the others finite
    assert np.array_equal(np.isfinite(st[0]), st[1] > 0) and np.isfinite(st[0]).sum() > 70


def make_cube_set(n, seed, noise=0.2, nchan=40):
    rng = np.random.default_rng(seed)
    cf, names = [], []
    for i in range(n):
        cube = Cube(nx=48, ny=48, nchan=nchan)
        ch = int(rng.integers(12, nchan - 12))
        f = 232e9 - ch * 7.6e7
        z = FEM / f - 1
        x, y = float(rng.integers(20, 28)), float(rng.integers(20, 28))
        cube.add_line(x, y, ch, rng.uniform(1, 3), 2.5)
        cube.add_noise(noise, rng)
        cube.write('c%02d.fits' % i)
        ra, dec = casa_style(cube.sky(x, y))
        with open('c%02d_coord.txt' % i, 'w') as fh:
            fh.write('%s, %s, %r\n' % (ra, dec, z))
        cf.append('c%02d_coord.txt' % i)
        names.append('c%02d.fits' % i)
    return LineStacker.readCoords(cf), names


@pytest.mark.parametrize('method,save', [('mean', 'amp'), ('mean', 'ampAndWidth'), ('median', 'amp'), ('mean', 'all')])
def test_bootstrap_cube_matches_bruteforce(workdir, method, save):
    coords, names = make_cube_set(7, 1)
    kw = dict(imagenames=names, fEm=FEM, stampsize=8, N_chans=15, method=method, save=save, nRandom=12)
    res = A.bootStraping_Cube(coords, seed=4, **kw)
    prep = li._prepare(coords, names, stampsize=8, N_chans=15, fEm=FEM, method=method)
    W = prep.weights()
    idx = np.random.default_rng(4).integers(0, 7, size=(12, 7))
    import LineStacker.tools.fit as fitTools
    for r in range(12):
        cube = li._combine(prep.data[idx[r]], W[idx[r]], prep.valid[idx[r]], method)
        spec = cube.sum(axis=(0, 1, 2))
        if save == 'amp':
            assert res[r + 1] == pytest.approx(np.max(spec), rel=1e-10)
        elif save == 'all':
            np.testing.assert_allclose(res[r + 1], cube, rtol=1e-10)
        else:
            p = fitTools.GaussFit(fctToFit=spec, returnInfos=True)[1]
            np.testing.assert_allclose(res[r + 1], [p[0], p[2]], rtol=1e-6)
    par = A.bootStraping_Cube(coords, seed=4, n_jobs=2, **kw)
    if save == 'all':
        for a, b in zip(res, par):
            np.testing.assert_array_equal(a, b)
    else:
        np.testing.assert_allclose(np.array(par, dtype=float), np.array(res, dtype=float))


def test_stack_estimator(workdir):
    coords, names = make_cube_set(6, 2, noise=1.0)
    out = A.stack_estimator(coords, nRandom=40, imagenames=names, fEm=FEM, N_chans=11, seed=3)
    out2 = A.stack_estimator(coords, nRandom=40, imagenames=names, fEm=FEM, N_chans=11, seed=3, n_jobs=3)
    assert len(out[2]) == 40 and out[2][0].shape == (11,)
    np.testing.assert_array_equal(np.array(out[2]), np.array(out2[2]))
    # pure noise of sigma 1 averaged over 6 sources
    assert abs(out[1]) < 0.2 and 0.25 < out[0] < 0.6
    allst = A.stack_estimator(coords, nRandom=3, imagenames=names, fEm=FEM, N_chans=11, stampsize=5, save='all', seed=1)
    assert allst[2][0].shape == (5, 5, 1, 11)


def test_noise_estimator(workdir):
    coords, names = make_cube_set(3, 3, noise=0.5)
    s = A.noise_estimator(names, nRandom=20000, maskradius=6, plotIt=False, seed=1)
    assert abs(s - 0.5) < 0.03
    per = A.noise_estimator(names, nRandom=2000, continuum=False, plotIt=False, seed=1)
    assert per.shape == (40,) and np.all(np.abs(per - 0.5) < 0.1)


def test_rebin_cubes_spectra(workdir):
    import LineStacker.tools.fit as fitTools
    from LineStacker import _io
    names, coords, freqs = [], [], []
    for i, sig in enumerate([2., 3., 4.]):
        cube = Cube(nx=30, ny=30, nchan=60, crpix3=7.)        # reference pixel not at 0
        cube.add_line(15, 15, 30, 2., sig)
        cube.write('r%d.fits' % i)
        names.append('r%d.fits' % i)
        coords.append(LineStacker.Coord(15., 15., spec=0., image='r%d.fits' % i, type='pixel'))
        freqs.append(cube.crval3 + (30 - cube.crpix0) * cube.cdelt3)
    out = A.rebin_CubesSpectra(coords, names)
    assert out == ['r0_SpectralRebinned.fits', 'r1_SpectralRebinned.fits', 'r2_SpectralRebinned.fits']
    for name, f in zip(out, freqs):
        cube, info = _io.read_cube(name)
        p = fitTools.GaussFit(fctToFit=cube[15, 15, 0, :], returnInfos=True)[1]
        assert abs(abs(p[2]) - 2.) < 0.1
        # the line keeps its frequency in the new WCS
        assert abs(info.value_of(p[1]) - f) < 0.05 * abs(info.spec_cdelt)


def test_regrid_cube_writer(workdir):
    coords, names = make_cube_set(4, 5, noise=0.)
    new_names, new_coords = A.regridFromZ(coords, names, stampsize=10, fEm=FEM, regridMethod='scaleToMax')
    assert len(new_names) == 4
    res = li.stack(new_coords, imagenames=new_names, fEm=FEM, stampsize=8, N_chans=11, outfile=None)
    assert np.argmax(res[1]) == 5


@pytest.mark.skipif(not os.path.isdir(os.path.join(EXAMPLE, 'oneD', 'dataExtended')), reason='example data not found')
def test_oned_example_data():
    d = os.path.join(EXAMPLE, 'oneD', 'dataExtended')
    ims = [OneD.Spectrum(spectrum=np.loadtxt(os.path.join(d, 'spectra_%d' % i)),
                      centralVelocity=np.loadtxt(os.path.join(d, 'central_velocity_%d' % i)),
                      name='spectrum_%d' % i) for i in range(50)]
    for m in ('mean', 'median'):
        new = OneD.Stack(ims, method=m)
        ref = _reference_oned.Stack(ims, method=m)
        np.testing.assert_allclose(new[0], ref[0], rtol=1e-12, atol=1e-14)


def test_maximize_snr_ignores_lines_outside_the_spectrum(monkeypatch):
    import LineStacker.tools.fit as fitTools
    rng = np.random.default_rng(0)
    spec = rng.normal(0, 1, 100) + 5 * np.exp(-(np.arange(100) - 50.)**2 / 8.)
    good = A.maximizeSNR(spec)
    assert 3 < good < 20
    # fit locked onto the edge, centred outside the spectrum: no detection
    monkeypatch.setattr(fitTools, 'GaussFit', lambda **k: (None, np.array([600., -22.7, 6.6])))
    assert A.maximizeSNR(spec) == 0
    # negative sigma is the same Gaussian
    monkeypatch.setattr(fitTools, 'GaussFit', lambda **k: (None, np.array([5., 50., -2.])))
    neg = A.maximizeSNR(spec)
    monkeypatch.setattr(fitTools, 'GaussFit', lambda **k: (None, np.array([5., 50., 2.])))
    assert neg == A.maximizeSNR(spec)


@pytest.mark.skipif(not os.path.isdir(os.path.join(EXAMPLE, 'oneD', 'dataExtended')), reason='example data not found')
def test_subsampling_finds_the_bright_spectra():
    """Documented result: spectra 40-49 (line amplitude 10) get the highest grades."""
    d = os.path.join(EXAMPLE, 'oneD', 'dataExtended')
    spectra = [OneD.Spectrum(spectrum=np.loadtxt(os.path.join(d, 'spectra_%d' % i)),
                             centralVelocity=np.loadtxt(os.path.join(d, 'central_velocity_%d' % i)),
                             name='spectrum_%d' % i) for i in range(50)]
    g = A.subsample_OneD(spectra, nRandom=10000, method='median', seed=1, n_jobs=-1)
    v = np.array([g['spectrum_%d' % i] for i in range(50)])
    assert len(set(np.argsort(-v)[:10]) & set(range(40, 50))) >= 9
    assert v[40:].min() > np.median(v[:40])

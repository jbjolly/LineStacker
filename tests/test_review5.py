"""Fixes from the fifth code review and the release review: local spline
resampling, rebinned and regridded spectra (axes, centres, widths), evenly
spaced spectral axes, decreasing axes, channels without data in 1D stacks,
line centres outside a spectrum, regridded cube labels, and smaller checks."""
import warnings

import numpy as np
import pytest

from conftest import Cube, FEM
from test_weights_wcs import line_cubes

import LineStacker
import LineStacker.OneD_Stacker as OneD
import LineStacker.analysisTools as A
import LineStacker.image as li
import LineStacker.tools.fit as fitTools


def _gauss(n, centre, sigma, amp=1.):
    x = np.arange(n, dtype=float)
    return amp * np.exp(-0.5 * ((x - centre) / sigma) ** 2)


# -- spline resampling ------------------------------------------------------------

def test_spline_values_never_depend_on_data_beyond_a_blank():
    rng = np.random.default_rng(0)
    a = rng.normal(size=60)
    a[30] = np.nan
    pos = np.linspace(0, 59, 101)
    one = A._spline_axis(a, 0, pos)
    b = a.copy()
    b[31:] = rng.normal(size=29) * 100.               # change everything after the blank
    two = A._spline_axis(b, 0, pos)
    left = pos < 30
    np.testing.assert_array_equal(one[left], two[left])    # bit-identical
    assert np.isnan(one[(pos > 29) & (pos < 31)]).all()
    # far from the blank the values are those of a spline through all the data
    full = A._spline_axis(np.where(np.isnan(a), 0., a), 0, pos)
    far = pos < 15
    np.testing.assert_allclose(one[far], full[far], atol=1e-4)


def test_spline_is_nan_beyond_the_last_sample():
    a = np.linspace(0., 9., 10)
    out = A._congrid_nan(a, [13], method='spline')        # positions up to 9.23
    assert np.isnan(out[-1])
    assert A.congrid(a, [13], method='spline')[-1] == 0.   # plain congrid gives 0 there
    np.testing.assert_allclose(out[:-1], A.congrid(a, [13], method='spline')[:-1])


# -- rebinned and regridded spectra ------------------------------------------------

def test_regridded_spectra_keep_even_axes_and_their_line():
    v = np.linspace(-500., 500., 101)
    out = []
    for z in (1.0, 1.4):
        line = np.exp(-0.5 * (v / (60. / (1 + z))) ** 2)
        out.append(OneD.Spectrum(spectrum=np.array([v, line]), z=z, fEmLine=230e9))
    # z sets the centre: put it back on the line for this test
    for sp in out:
        sp.centerIndex = 50
    new = A.regridFromZ1D(out, regridMethod='scaleToMin', verbose=False)[1]
    ratio = 2.4 / 2.0
    assert len(new.amp) == int(101 * ratio)
    steps = np.diff(new.velocities)
    np.testing.assert_allclose(steps, 10. / ratio, rtol=1e-12)           # even, the right width
    assert new.velBin == pytest.approx(10. / ratio)
    assert new.velocities[-1] != 0.
    assert new.centerIndex == LineStacker.round_half_up(50 * ratio)
    assert abs(int(np.nanargmax(new.amp)) - new.centerIndex) <= 1
    assert new.velocities[new.centerIndex] == pytest.approx(0., abs=10. / ratio)
    np.testing.assert_allclose(new.spectrum[0], new.velocities)          # [2, N] array follows
    np.testing.assert_allclose(new.spectrum[1], new.amp)


def test_rebinned_spectra_keep_their_line_centre():
    a = OneD.Spectrum(amp=_gauss(60, 40., 6.), centerIndex=40, name='wide')
    b = OneD.Spectrum(amp=_gauss(60, 40., 3.), centerIndex=40, name='narrow')
    new = A.rebin_OneD([a, b], verbose=False)
    assert len(new[0].amp) == 30
    assert abs(new[0].centerIndex - 20) <= 1
    assert abs(int(np.argmax(new[0].amp)) - new[0].centerIndex) <= 1
    assert new[1].centerIndex == 40 and len(new[1].amp) == 60
    st, n = OneD.Stack(new, chansStack=21)
    assert np.argmax(st) == 10


@pytest.mark.filterwarnings('ignore:Covariance of the parameters')   # noise-free lines
def test_rebinned_widths_are_in_channels():
    # same line in km/s, but the second spectrum has channels twice as wide
    v1, v2 = np.arange(-300., 301., 5.), np.arange(-600., 601., 10.)
    s1 = OneD.Spectrum(amp=np.exp(-0.5 * (v1 / 40.) ** 2), velocities=v1, centralVelocity=0., fit=True)
    s2 = OneD.Spectrum(amp=np.exp(-0.5 * (v2 / 40.) ** 2), velocities=v2, centralVelocity=0., fit=True)
    new = A.rebin_OneD([s1, s2], verbose=False)
    # s2's line spans half as many channels: it is the narrowest, s1 is halved
    assert len(new[1].amp) == len(v2)
    assert abs(len(new[0].amp) - len(v1) / 2.) <= 1
    assert new[0].velBin == pytest.approx(10., rel=0.02)


def test_per_channel_weights_follow_the_resampling():
    w = np.linspace(1., 2., 41)
    sp = OneD.Spectrum(amp=_gauss(41, 20., 6.), centerIndex=20, weights=w, z=1., fEmLine=1e9,
                       frequencies=np.linspace(4.9e8, 5.1e8, 41), velOrFreq='freq')
    sp.centerIndex = 20
    ref = OneD.Spectrum(amp=_gauss(41, 20., 3.), centerIndex=20, z=0.5, fEmLine=1e9,
                        frequencies=np.linspace(6.5e8, 6.8e8, 41), velOrFreq='freq')
    ref.centerIndex = 20
    new = A.regridFromZ1D([sp, ref], verbose=False)[0]
    assert len(new.weights) == len(new.amp)
    assert np.all(new.weights[np.isfinite(new.weights)] >= 1.)


# -- spectral axes ----------------------------------------------------------------

def test_uneven_axes_are_refused():
    v = np.arange(50.) * 10.
    v[25:] += 20.                                     # two channels missing
    with pytest.raises(ValueError, match='evenly spaced'):
        OneD.Spectrum(amp=np.ones(50), velocities=v, centerIndex=10)
    f = np.geomspace(1e9, 1.3e9, 200)                 # log-spaced
    with pytest.raises(ValueError, match='evenly spaced'):
        OneD.Spectrum(amp=np.ones(200), frequencies=f, velOrFreq='freq', centerIndex=10)


def test_rounded_axes_are_accepted_with_their_mean_width():
    v = np.round(np.arange(200) * 10.05 - 1000., 1)  # values written with one decimal
    sp = OneD.Spectrum(amp=np.ones(200), velocities=v, centralVelocity=v[150])
    assert sp.velBin == pytest.approx((v[-1] - v[0]) / 199.)
    assert sp.centerIndex == 150
    assert isinstance(sp.velocities, np.ndarray)


def test_sigma2_with_a_fit_on_decreasing_axes():
    rng = np.random.default_rng(3)
    v = np.arange(-400., 401., 10.)
    amp = np.exp(-0.5 * (v / 30.) ** 2) * 5 + rng.normal(0., 0.2, len(v))
    up = OneD.Spectrum(amp=amp, velocities=v, centralVelocity=0., fit=True, weights='sigma2')
    down = OneD.Spectrum(amp=amp[::-1], velocities=v[::-1], centralVelocity=0., fit=True, weights='sigma2')
    # the line is left out of the noise either way (it was not for decreasing
    # axes: the window was reversed and the noise included the line)
    off = np.abs(v) > 30. * 2.35 + 10.
    for sp in (up, down):
        assert sp.weights == pytest.approx(1. / np.std(amp[off]) ** 2, rel=0.1)
        assert sp.weights > 5. / np.var(amp)
    assert up.weights == pytest.approx(down.weights, rel=0.05)


# -- 1D stacks ----------------------------------------------------------------------

def test_channels_without_data_are_nan_in_1d_stacks():
    a = OneD.Spectrum(amp=np.ones(11), centerIndex=5)
    b = OneD.Spectrum(amp=np.full(11, 3.), centerIndex=5)
    for method in ('mean', 'median'):
        st, n = OneD.Stack([a, b], chansStack=21, method=method)
        assert np.isnan(st[:5]).all() and np.isnan(st[16:]).all()
        np.testing.assert_allclose(st[5:16], 2.)
        assert np.array_equal(n[:5], np.zeros(5))
    # spectra whose weight is 0 leave their channels without data
    c = OneD.Spectrum(amp=np.ones(11), centerIndex=5, weights=0.)
    st, _ = OneD.Stack([c])
    assert np.isnan(st).all()
    boot = A.bootstraping_OneD([a, b], nRandom=5, save='amp', chansStack=21, seed=1, verbose=False)
    assert np.all(np.isfinite(boot))
    assert A.maximizeAmp(np.array([np.nan, 1., 2.])) == 2.
    assert A.maximizeAmp(np.full(3, np.nan)) == 0.


def test_line_centre_outside_the_spectrum_is_an_error():
    v = np.linspace(-500., 500., 41)
    sp = OneD.Spectrum(amp=np.ones(41), velocities=v, z=2., fEmLine=230e9)   # line at 200 000 km/s
    with pytest.raises(ValueError, match='outside the spectrum'):
        OneD.Stack([sp])
    # another centring can still be used
    st, n = OneD.Stack([sp], center='zero_vel')
    assert len(st) == 41
    with pytest.raises(ValueError, match='zero velocity'):
        OneD.Stack([OneD.Spectrum(amp=np.ones(5), velocities=np.arange(100., 105.))], center='zero_vel')
    with pytest.raises(ValueError, match='needs velocities'):
        OneD.Stack([OneD.Spectrum(amp=np.ones(5), frequencies=np.arange(5.), velOrFreq='freq')],
                   center='zero_vel')


def test_fit_centring_of_frequency_spectra_uses_frequencies():
    f = np.linspace(100e9, 101e9, 101)
    amp = np.exp(-0.5 * ((f - 100.6e9) / 3e7) ** 2)
    sp = OneD.Spectrum(amp=amp, frequencies=f, velOrFreq='freq', centerIndex=0)
    sp.velocities = np.arange(101.)                 # unrelated velocities are not used
    OneD._set_centers([sp], 'fit')
    assert sp.centerIndex == 60


# -- regridded cubes ----------------------------------------------------------------

def test_regridded_cube_labels_put_the_line_at_fem(workdir):
    from astropy.io import fits
    names, coords, chans = [], [], []
    for i, z in enumerate((5.6, 6.3)):
        cube = Cube(nx=30, ny=30, nchan=80, crval3=FEM / (1 + z) + 37.3 * 7.6e7)
        chan = cube.channel_of(FEM / (1 + z))            # fractional channel
        cube.add_line(15., 15., chan, 2., 2.)
        cube.write('c%d.fits' % i)
        s = cube.sky(15., 15.)
        coords.append(LineStacker.Coord(s.ra, s.dec, spec=z, image='c%d.fits' % i))
        names.append('c%d.fits' % i)
    new_names, new_coords = A.regridFromZ(coords, names, stampsize=8, fEm=FEM, regridMethod='scaleToMax',
                                          outdir='rz', verbose=False)
    h = fits.getheader(new_names[0])
    data = li.readCube(new_names[0])[4, 4, 0, :]
    ok = np.isfinite(data)
    k = np.arange(len(data))[ok]
    # line centre from the moments of the regridded spectrum
    w = np.clip(data[ok], 0, None) ** 2
    centre = np.sum(k * w) / np.sum(w)
    rest = h['CRVAL3'] + (centre + 1 - h['CRPIX3']) * h['CDELT3']
    assert abs(rest - FEM) < 0.1 * abs(h['CDELT3'])
    assert abs(h['CDELT3']) == pytest.approx(7.6e7 * 7.3, rel=1e-12)


# -- smaller checks -------------------------------------------------------------------

def test_felo_axes_are_non_linear(workdir):
    cube = Cube(nchan=30)
    cube.header['CTYPE3'], cube.header['CUNIT3'] = 'FELO-HEL', 'm/s'
    cube.header['CRVAL3'], cube.header['CDELT3'] = 1e6, 1e4
    cube.header['RESTFRQ'] = 1.420405752e9
    cube.write('felo.fits')
    from LineStacker import _io
    assert not _io.read_info('felo.fits').spec_linear


def test_region_outside_the_cube_is_an_error(workdir):
    coords, names, _ = line_cubes(n=2, noise=0.)
    with pytest.raises(ValueError, match='not inside'):
        A.rebin_CubesSpectra(coords, imagenames=names, regionSize=200, verbose=False)


def test_blank_primary_beam_response_warns(workdir):
    coords, names, _ = line_cubes(n=2, noise=0.05)

    def pb(dx, dy, nu):
        return np.nan
    with pytest.warns(UserWarning, match='blank'):
        res = li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=5, outfile=None,
                       weighting='pb', primarybeam=pb, verbose=False)
    assert np.isnan(res[2]).all()                 # every weight is 0: no data


def test_congrid_has_its_docstring():
    assert A.congrid.__doc__ and 'resampling' in A.congrid.__doc__
    for f in (A.maximizeAmp, A.maximizeOutflow, A.RebinToRest):
        assert f.__doc__


# -- final review --------------------------------------------------------------------

def test_blank_single_plane_primary_beam_warns(workdir):
    coords, names, truth = line_cubes(n=2, noise=0.05)
    pbs = []
    for k, t in enumerate(truth):
        pb = Cube(nx=40, ny=40, nchan=None)
        pb.arr[:] = 0.8
        if k == 0:
            pb.arr[int(t['x']) - 2:int(t['x']) + 3, int(t['y']) - 2:int(t['y']) + 3] = np.nan
        pb.write('pb%d.fits' % k)
        pbs.append('pb%d.fits' % k)
    with pytest.warns(UserWarning, match='blank'):
        prep = li._prepare(coords, names, stampsize=8, N_chans=5, fEm=FEM, weighting='pb', primarybeam=pbs)
    w = [float(np.max(c.weight)) for c in prep.coords]
    assert w[0] == 0. and w[1] == pytest.approx(0.64)


@pytest.mark.parametrize('size', [4, 5])
def test_region_spectra_have_the_requested_size(workdir, monkeypatch, size):
    coords, names, truth = line_cubes(n=2, noise=0.)
    li.clear_cache()
    seen = []
    fit = fitTools.GaussFit

    def spy(fctToFit=None, **kw):
        seen.append(np.array(fctToFit, dtype=float))
        return fit(fctToFit=fctToFit, **kw)
    monkeypatch.setattr(fitTools, 'GaussFit', spy)
    A.rebin_CubesSpectra(coords, imagenames=names, regionSize=size, verbose=False)
    t = truth[0]
    cube = li.readCube(names[0])
    bx = int(np.floor(round(t['x']) - size // 2 + 0.5))
    by = int(np.floor(round(t['y']) - size // 2 + 0.5))
    np.testing.assert_allclose(seen[0], cube[bx:bx + size, by:by + size, 0, :].sum(axis=(0, 1)))


def test_congrid_nan_checks_its_dimensions():
    for bad in ([0], [2.5], [3, 3]):
        with pytest.raises(ValueError, match='dimensions'):
            A._congrid_nan(np.arange(10.), bad, method='spline')
    with pytest.raises(ValueError, match='method'):
        A._congrid_nan(np.arange(10.), [5], method='cubic')

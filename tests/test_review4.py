"""Fixes from the fourth code review: blank (NaN) channels next to exact
samples, NaN-aware rebinning and regridding, the SNR grade's noise, the
order of Spectrum's axis checks and distortion signatures."""
import types

import numpy as np
import pytest

from conftest import Cube
from test_weights_wcs import line_cubes

import LineStacker
import LineStacker.OneD_Stacker as OneD
import LineStacker.analysisTools as A
import LineStacker.image as li
import LineStacker.tools.fit as fitTools


def _gauss(n, centre, sigma, amp=1.):
    x = np.arange(n, dtype=float)
    return amp * np.exp(-0.5 * ((x - centre) / sigma) ** 2)


# -- exact channels next to a blank one ------------------------------------------

def _blank_middle_cube(nchan=3):
    cube = np.arange(8 * 8 * nchan, dtype=float).reshape(8, 8, 1, nchan)
    cube[..., 1] = np.nan
    return cube, types.SimpleNamespace(nx=8, ny=8, path='test')


def test_extract_exact_channel_ignores_blank_neighbour():
    cube, info = _blank_middle_cube()
    # samples at channels -2 (outside), 0 and 2; channel 1 is blank
    stamp, inside = li._extract(cube, info, 4, 4, 0., 2., 4, 3)
    assert list(inside) == [False, True, True]
    np.testing.assert_array_equal(stamp[:, :, 0, 1], cube[2:6, 2:6, 0, 0])
    np.testing.assert_array_equal(stamp[:, :, 0, 2], cube[2:6, 2:6, 0, 2])
    assert np.isnan(stamp[:, :, 0, 0]).all()
    # a fractional sample uses the blank channel: blank
    stamp, _ = li._extract(cube, info, 4, 4, 0.5, 1., 4, 1)
    assert np.isnan(stamp).all()


def test_extract_many_exact_channel_ignores_blank_neighbour():
    cube, info = _blank_middle_cube()
    stamps = li._extract_many(cube, info, [3, 4], [4, 5], 0., 2., 4, 3)
    assert stamps.shape == (2, 4, 4, 1, 3)
    np.testing.assert_array_equal(stamps[0, :, :, 0, 1], cube[1:5, 2:6, 0, 0])
    np.testing.assert_array_equal(stamps[1, :, :, 0, 2], cube[2:6, 3:7, 0, 2])
    assert np.isnan(stamps[..., 0]).all()
    for x, y, s in zip([3, 4], [4, 5], stamps):
        one, _ = li._extract(cube, info, x, y, 0., 2., 4, 3)
        np.testing.assert_array_equal(s, one)


def test_blank_channel_does_not_blank_exact_stack_channels(workdir):
    coords, names, truth = line_cubes(n=2, noise=0.)
    for t, name in zip(truth, names):
        cube = t['cube']
        cube.arr[..., t['chan'] + 1] = np.nan     # the channel after the line
        cube.write(name)
    li.clear_cache()
    pixel = [LineStacker.Coord(t['x'], t['y'], spec=t['chan'], image=n, type='pixel')
             for t, n in zip(truth, names)]
    peak, spec, stacked, missing = li.stack(pixel, imagenames=names, spectralMethod='channel',
                                            stampsize=8, N_chans=5, outfile=None)
    assert np.isfinite(spec[[0, 1, 2, 4]]).all()   # only the blank channel is blank
    assert np.isnan(spec[3])
    assert np.isfinite(peak) and peak > 0.5


# -- rebinning and regridding ------------------------------------------------------

def test_congrid_nan_linear_blanks_only_what_it_uses():
    a = np.arange(1., 9.)
    a[2] = np.nan
    out = A._congrid_nan(a, [4])                   # samples 0, 2, 4, 6
    assert np.isnan(out[1]) and np.isfinite(out[[0, 2, 3]]).all()
    np.testing.assert_allclose(out[[0, 2, 3]], [1., 5., 7.])
    a[2] = 3.
    np.testing.assert_allclose(A._congrid_nan(a, [4]), A.congrid(a, [4]))


def test_congrid_nan_spline_keeps_blanks_local():
    a = _gauss(40, 20., 4.)
    a[10] = np.nan
    out = A._congrid_nan(a, [30], method='spline')
    bad = np.flatnonzero(np.isnan(out))
    assert 0 < len(bad) <= 4
    assert np.all(np.abs(bad * 40. / 30. - 10.) < 2.)
    # a plain spline spreads the blank everywhere
    assert np.isnan(A.congrid(a, [30], method='spline')).sum() > 4
    # with no blanks it is congrid
    b = _gauss(40, 20., 4.)
    np.testing.assert_allclose(A._congrid_nan(b, [30], method='spline'), A.congrid(b, [30], method='spline'))


def test_congrid_nan_keeps_blanks_on_unresampled_axes():
    a = np.ones((3, 3, 20))
    a[1, 1, 5] = np.nan
    out = A._congrid_nan(a, [3, 3, 30], method='spline')   # upsampled spectrally
    # the last new channel lies beyond the last old one (19.33): no data
    assert np.isnan(out[..., -1]).all()
    body = out[..., :-1]
    assert np.isfinite(body[[0, 2]]).all() and np.isfinite(body[:, [0, 2]]).all()
    # pixel (1, 1): only the new channels between old channels 4 and 6
    bad = np.flatnonzero(np.isnan(body[1, 1]))
    assert len(bad) and np.all(np.abs(bad * 20. / 30. - 5.) < 1.)


def test_rebin_and_regrid_1d_keep_blanks_local():
    wide = _gauss(41, 20., 6.)
    wide[6] = np.nan                                  # sampled (or next to a sample) after rebinning
    narrow = _gauss(41, 20., 3.)
    spectra = [OneD.Spectrum(amp=wide, centerIndex=20), OneD.Spectrum(amp=narrow, centerIndex=20)]
    out = A.rebin_OneD(spectra, verbose=False)
    amp = np.asarray(out[0].amp)
    assert 0 < np.isnan(amp).sum() <= 2 and np.isfinite(amp[len(amp) // 2])
    z = [OneD.Spectrum(amp=wide, frequencies=np.linspace(100e9, 101e9, 41), velOrFreq='freq',
                       z=1., fEmLine=201e9, name='a'),
         OneD.Spectrum(amp=narrow, frequencies=np.linspace(100e9, 101e9, 41), velOrFreq='freq',
                       z=0.5, fEmLine=151e9, name='b')]
    out = A.regridFromZ1D(z, regridMethod='scaleToMin', verbose=False)
    amp = np.asarray(out[0].amp)
    assert len(amp) != 41
    assert 0 < np.isnan(amp).sum() <= 6 and np.isfinite(amp).sum() >= len(amp) - 6


def test_rebin_cubes_keeps_blanks_local(workdir):
    coords, names, truth = line_cubes(n=2, noise=0.)
    truth[0]['cube'].arr[..., 4] = np.nan             # blank channel, far from the line
    truth[0]['cube'].write(names[0])
    li.clear_cache()
    out = A.rebin_CubesSpectra(coords, imagenames=names, widths=[4., 2.], verbose=False)
    first = li.readCube(out[0])                       # 50 -> 25 channels: samples 0, 2, 4...
    assert first.shape[3] == 25
    assert np.isnan(first[..., 2]).all()
    assert np.isfinite(np.delete(first, 2, axis=3)).all()


def test_rebin_cubes_region_spectra_keep_blank_channels(workdir, monkeypatch):
    coords, names, truth = line_cubes(n=2, noise=0.)
    truth[0]['cube'].arr[..., 4] = np.nan
    truth[0]['cube'].write(names[0])
    li.clear_cache()
    seen = []
    fit = fitTools.GaussFit

    def spy(fctToFit=None, **kw):
        seen.append(np.array(fctToFit, dtype=float))
        return fit(fctToFit=fctToFit, **kw)
    monkeypatch.setattr(fitTools, 'GaussFit', spy)
    A.rebin_CubesSpectra(coords, imagenames=names, regionSize=4, verbose=False)
    assert np.isnan(seen[0][4]) and np.isfinite(np.delete(seen[0], 4)).all()   # not 0
    assert np.isfinite(seen[1]).all()


def test_rebin_to_rest_upsampling_with_blanks():
    stamp = np.ones((4, 4, 1, 20))
    stamp[1, 2, 0, 7] = np.nan
    c = LineStacker.Coord(0., 0., spec=1.)
    new, changed = A.RebinToRest(c, stamp, c.z, zToScale=0.5)
    assert changed and new.shape == (4, 4, 1, 27)    # round(20 * 2 / 1.5)
    # new channel j samples old channel 0.75 j: the last one (19.5) is beyond the data
    assert np.isnan(new[..., -1]).all()
    flat = new[..., :-1].reshape(16, -1)
    assert np.isfinite(np.delete(flat, 1 * 4 + 2, axis=0)).all()
    # pixel (1, 2): the new channels between old channels 6 and 8 (6.75, 7.5)
    assert list(np.flatnonzero(np.isnan(flat[1 * 4 + 2]))) == [9, 10]


# -- SNR grade -----------------------------------------------------------------------

def test_maximize_snr_ignores_blank_noise_channels():
    rng = np.random.default_rng(1)
    s = _gauss(80, 40., 3., 5.) + rng.normal(0., 0.1, 80)
    clean = A.maximizeSNR(s)
    s[[2, 70, 75]] = np.nan
    grade = A.maximizeSNR(s)
    assert np.isfinite(grade) and grade > 10
    assert abs(grade / clean - 1.) < 0.2
    flat = np.full(80, np.nan)
    flat[38:43] = 1.
    assert A.maximizeSNR(flat) == 0


# -- Spectrum axes ---------------------------------------------------------------------

@pytest.mark.parametrize('kw', [dict(velocities=[0.]),
                                dict(frequencies=[1e9], velOrFreq='freq')])
def test_spectrum_short_axis_is_a_value_error(kw):
    with pytest.raises(ValueError):
        OneD.Spectrum(amp=np.ones(5), centerIndex=2, **kw)


def test_spectrum_axis_checks_come_before_bins():
    with pytest.raises(ValueError, match='finite'):
        OneD.Spectrum(amp=np.ones(3), velocities=[np.nan, 1., 2.], centerIndex=1)
    sp = OneD.Spectrum(amp=np.ones(3), frequencies=[3., 2., 1.], velOrFreq='freq', centerIndex=1)
    assert sp.freqBin == -1.


# -- distortion signatures -------------------------------------------------------------

def test_distortion_signatures_with_different_keys_differ():
    a = {'sip': None, 'cpdis1': None}
    assert li._same_distortion(a, dict(a))
    assert not li._same_distortion(a, {'sip': None})
    assert not li._same_distortion(a, {'sip': None, 'cpdis2': None})

"""Continuum stacking with LineStacker.image, and mixed FITS/CASA inputs."""
import warnings

import numpy as np
import pytest

from conftest import Cube, FEM, casa_style, needs_casa

import LineStacker
import LineStacker.image as li
import LineStacker.analysisTools as A
from LineStacker import _io


def sources(n=5, seed=0, nchan=None, noise=0.01, writer='fits', line_chan=None, nan_chans=()):
    """n images/cubes with a point source each; returns coord files, names,
    cubes, true continuum fluxes (sum over channels)."""
    rng = np.random.default_rng(seed)
    cf, names, cubes, flux = [], [], [], []
    for i in range(n):
        cube = Cube(nx=50, ny=50, nchan=nchan)
        x, y = float(rng.integers(20, 30)), float(rng.integers(20, 30))
        f = rng.uniform(1, 2)
        if nchan is None or nchan == 1:
            cube.add_line(x, y, 0, f, 1.)
        else:
            # flat continuum plus a line in a few channels
            X, Y = np.meshgrid(np.arange(50), np.arange(50), indexing='ij')
            spat = np.exp(-((X - x)**2 + (Y - y)**2) / (2 * 1.5**2))
            cube.arr += f * spat[:, :, None]
            if line_chan is not None:
                cube.add_line(x, y, line_chan, 3 * f, 2.)
        cube.add_noise(noise, rng)
        for ch in nan_chans:
            cube.arr[:, :, ch] = np.nan
        name = 'src%d.%s' % (i, 'image' if writer == 'casa' else 'fits')
        cube.write(name, casa=(writer == 'casa'))
        ra, dec = casa_style(cube.sky(x, y))
        f_true = np.nansum(cube.arr[int(x), int(y), :])
        with open('src%d_coord.txt' % i, 'w') as fh:
            if nchan and nchan > 1:
                z = FEM / (232e9 - (line_chan or 20) * 7.6e7) - 1
                fh.write('%s, %s, %r\n' % (ra, dec, z))
            else:
                fh.write('%s, %s\n' % (ra, dec))
        cf.append('src%d_coord.txt' % i)
        names.append(name)
        cubes.append(cube)
        flux.append(f_true)
    return cf, names, cubes, np.array(flux)


def no_linestacker_warning(record):
    return not [w for w in record if 'continuum stacking' in str(w.message)]


def test_2d_images_are_stacked_as_continuum_without_warning(workdir):
    cf, names, _, flux = sources(nchan=None)
    coords = LineStacker.readCoords(cf, lineON=False)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        res = li.stack(coords, imagenames=names, stampsize=16, outfile='cont.fits')
    assert no_linestacker_warning(rec)
    assert res[0] == pytest.approx(np.mean(flux), abs=0.01)
    assert res[2].shape == (16, 16, 1, 1) and len(res[1]) == 1
    from astropy.io import fits
    h = fits.getheader('cont.fits')
    assert h['NAXIS'] == 2 and 'CTYPE3' not in h and h['NSTACK'] == 5
    assert li.getFlux('cont.fits') == pytest.approx(res[0], rel=1e-6)
    # without lineON=False the coordinate files have only 2 columns: plain
    # readCoords refuses them, so lineON=False is the documented way


def test_single_channel_cubes_are_continuum_without_warning(workdir):
    cf, names, _, flux = sources(nchan=1)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        res = li.stack(LineStacker.readCoords(cf, lineON=False), imagenames=names, stampsize=16, outfile='c1.fits')
    assert no_linestacker_warning(rec)
    assert res[0] == pytest.approx(np.mean(flux), abs=0.01)
    from astropy.io import fits
    h = fits.getheader('c1.fits')
    assert h['CTYPE3'] == 'FREQ' and h['NAXIS3'] == 1


def test_spectral_cubes_collapsed_with_warning(workdir):
    cf, names, cubes, flux = sources(nchan=30, line_chan=15, noise=0.)
    coords = LineStacker.readCoords(cf)                   # numeric z: ignored in continuum mode
    with pytest.warns(UserWarning, match='collapsing 5 cube'):
        res = li.stack(coords, imagenames=names, stampsize=12, continuum=True, outfile='coll.fits')
    # reference: stack of the channel-summed images (default collapse)
    pix = LineStacker.getPixelCoords(coords, names)
    stamps = []
    for c, cube in zip(pix, cubes):
        img = cube.arr.astype('f4').sum(axis=2)
        bx, by = int(np.floor(c.ra - 6 + 0.5)), int(np.floor(c.dec - 6 + 0.5))
        stamps.append(img[bx:bx + 12, by:by + 12])
    np.testing.assert_allclose(res[2][:, :, 0, 0], np.mean(stamps, axis=0), rtol=1e-5)
    # only the collapsed version of each cube is kept in memory
    keys = [k for k in _io._data_cache if k[0].endswith('src0.fits')]
    assert [k[-1] for k in keys] == ['sum']
    from astropy.io import fits
    h = fits.getheader('coll.fits')
    assert h['NAXIS3'] == 1 and abs(h['CDELT3']) == pytest.approx(30 * 7.6e7)
    assert h['COLLAPSE'] == 'sum'
    assert h['CRVAL3'] == pytest.approx(232e9 - 14.5 * 7.6e7)


def test_coordinates_without_spectral_value_imply_continuum(workdir):
    cf, names, _, flux = sources(nchan=20, noise=0.)
    with open(cf[0]) as fh:
        assert len(fh.read().split(',')) == 3
    for f in cf:                                           # drop the redshift column
        with open(f) as fh:
            ra, dec, _ = fh.read().split(',')
        with open(f, 'w') as fh:
            fh.write('%s,%s\n' % (ra, dec))
    coords = LineStacker.readCoords(cf, lineON=False)
    with pytest.warns(UserWarning, match='coordinates have no spectral value'):
        res = li.stack(coords, imagenames=names, stampsize=12, outfile=None)
    assert res[0] == pytest.approx(np.mean(flux), rel=1e-5)


@pytest.mark.parametrize('how', ['mean', 'median', 'sum'])
def test_collapse_ignores_blank_channels(workdir, how):
    cf, names, cubes, _ = sources(n=3, nchan=12, noise=0.05, nan_chans=(0, 5, 6))
    coords = LineStacker.readCoords(cf)
    with pytest.warns(UserWarning):
        res = li.stack(coords, imagenames=names, stampsize=10, continuum=True, collapse=how, outfile=None)
    pix = LineStacker.getPixelCoords(coords, names)
    func = {'mean': np.nanmean, 'median': np.nanmedian, 'sum': np.nansum}[how]
    ref = []
    for c, cube in zip(pix, cubes):
        img = func(cube.arr.astype('f4'), axis=2)
        bx, by = int(np.floor(c.ra - 5 + 0.5)), int(np.floor(c.dec - 5 + 0.5))
        ref.append(img[bx:bx + 10, by:by + 10])
    np.testing.assert_allclose(res[2][:, :, 0, 0], np.mean(ref, axis=0), rtol=1e-5, atol=1e-6)


def test_line_mode_refuses_images_without_spectral_axis(workdir):
    cf, names, _, _ = sources(n=2, nchan=20)
    Cube(nx=50, ny=50, nchan=None).write('flat.fits')
    coords = LineStacker.readCoords(cf)
    coords[1].image, coords[1].imagenumber = 'flat.fits', 1
    with pytest.raises(ValueError, match='continuum=True'):
        li.stack(coords, imagenames=[names[0], 'flat.fits'], fEm=FEM, stampsize=8, outfile=None)


def test_continuum_weighting_and_median(workdir):
    cf, names, _, flux = sources(nchan=None, noise=0.05, seed=3)
    coords = LineStacker.readCoords(cf, lineON=False)
    plain = li.stack(coords, imagenames=names, stampsize=16, outfile=None)
    for w in ('sigma2', 'sigma'):
        r = li.stack(coords, imagenames=names, stampsize=16, weighting=w, maskradius=4, outfile=None)
        assert r[0] == pytest.approx(plain[0], abs=0.05)
    med = li.stack(coords, imagenames=names, stampsize=16, method='median', outfile=None)
    assert med[0] == pytest.approx(np.median(flux), abs=0.06)

    # primary-beam weights: callable(dx, dy, nu) -> response
    def pb(dx, dy, nu):
        return np.exp(-(dx**2 + dy**2) / (2 * (2e-5)**2))
    r = li.stack(coords, imagenames=names, stampsize=16, weighting='pb', primarybeam=pb, outfile=None)
    prep = li._prepare(coords, names, stampsize=16)
    w = []
    for c, inf in zip(prep.coords, prep.infos):
        scale = np.deg2rad(inf.pixscale_deg)
        crpix = inf.celestial.wcs.crpix - 1
        w.append(pb((c.ra - crpix[0]) * scale, (c.dec - crpix[1]) * scale, None)**2)
    ref = np.tensordot(w, prep.data[:, :, :, 0, 0], axes=(0, 0)) / np.sum(w)
    np.testing.assert_allclose(r[2][:, :, 0, 0], ref, rtol=1e-6)
    with pytest.raises(ValueError, match='primarybeam'):
        li.stack(coords, imagenames=names, stampsize=16, weighting='pb', outfile=None)


def test_continuum_bootstrap_and_estimator(workdir):
    cf, names, _, flux = sources(nchan=None, noise=0.1, seed=4)
    coords = LineStacker.readCoords(cf, lineON=False)
    boot = A.bootStraping_Cube(coords, imagenames=names, stampsize=8, nRandom=50, save='all', seed=1)
    assert len(boot) == 51 and boot[1].shape == (8, 8, 1, 1)
    amp = A.bootStraping_Cube(coords, imagenames=names, stampsize=8, nRandom=50, save='amp', seed=1)
    assert amp[1] == pytest.approx(boot[1].sum(), rel=1e-9)
    with pytest.raises(ValueError, match='needs line stacking'):
        A.bootStraping_Cube(coords, imagenames=names, stampsize=8, nRandom=5, save='ampAndWidth')
    est = A.stack_estimator(coords, imagenames=names, nRandom=30, seed=2, n_jobs=2)
    assert len(est[2]) == 30 and est[2][0].shape == (1,)
    assert abs(est[1]) < 0.2


def test_continuum_estimator_on_collapsed_cubes(workdir):
    cf, names, _, _ = sources(n=3, nchan=10, noise=0.2, seed=6)
    coords = LineStacker.readCoords(cf)
    with pytest.warns(UserWarning):
        est = A.stack_estimator(coords, imagenames=names, nRandom=20, continuum=True, seed=1)
    assert est[2][0].shape == (1,) and abs(est[1]) < 0.2


@needs_casa
def test_mixed_fits_and_casa_inputs(workdir):
    cf, names_f, _, _ = sources(n=4, nchan=30, line_chan=15, noise=0.05, seed=7)
    coords = LineStacker.readCoords(cf)
    a = li.stack(coords, imagenames=names_f, fEm=FEM, stampsize=8, N_chans=9, outfile=None)
    from casatools import image
    mixed = list(names_f)
    for k in (1, 3):                              # same data as CASA images
        out = names_f[k].replace('.fits', '.image')
        ia = image()
        ia.fromfits(outfile=out, infile=names_f[k], overwrite=True)
        ia.done()
        mixed[k] = out
    _io.clear_cache()
    b = li.stack(coords, imagenames=mixed, fEm=FEM, stampsize=8, N_chans=9, outfile=None)
    np.testing.assert_array_equal(a[2], b[2])


@needs_casa
def test_casa_continuum_images(workdir):
    cf, names, _, flux = sources(nchan=None, writer='casa', seed=8)
    res = li.stack(LineStacker.readCoords(cf, lineON=False), imagenames=names, stampsize=16,
                   outfile='cont_stack.image')
    assert res[0] == pytest.approx(np.mean(flux), abs=0.01)
    from casatools import image
    ia = image()
    ia.open('cont_stack.image')
    shape = ia.shape()
    ia.done()
    assert list(shape) == [16, 16]
    assert li.getFlux('cont_stack.image') == pytest.approx(res[0], rel=1e-5)


def test_verbose_reports_written_files(workdir, capsys, monkeypatch):
    import os
    cf, names, _, _ = sources(n=3, nchan=None, seed=9)
    coords = LineStacker.readCoords(cf, lineON=False)
    assert 'read 3 positions from 3 coordinate files' in capsys.readouterr().out
    li.stack(coords, imagenames=names, stampsize=8, outfile='v.fits', returnWeights=True)
    out = capsys.readouterr().out
    assert os.path.abspath('v.fits') in out and 'continuum stack of 3 sources' in out
    assert os.path.abspath('weightFile.csv') in out
    li.stack(coords, imagenames=names, stampsize=8, outfile='q.fits', verbose=False)
    A.bootStraping_Cube(coords, imagenames=names, stampsize=8, nRandom=20, seed=1, verbose=False)
    assert capsys.readouterr().out == ''
    monkeypatch.setattr(LineStacker, 'verbose', False)
    LineStacker.readCoords(cf, lineON=False)
    li.stack(coords, imagenames=names, stampsize=8, outfile='q.fits')
    A.stack_estimator(coords, imagenames=names, nRandom=5, seed=1)
    assert capsys.readouterr().out == ''
    li.stack(coords, imagenames=names, stampsize=8, outfile='q.fits', verbose=True)   # per-call wins
    assert os.path.abspath('q.fits') in capsys.readouterr().out

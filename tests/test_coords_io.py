import numpy as np
import pytest

from conftest import Cube, FEM, casa_style

import LineStacker
import LineStacker.image as li
import LineStacker.analysisTools as A
import LineStacker.OneD_Stacker as OneD


def test_readcoords_formats(workdir):
    cube = Cube()
    sky = cube.sky(20, 21)
    rad, decd = float(sky.ra.deg), float(sky.dec.deg)
    ra_c, dec_c = casa_style(sky)
    lines = {
        'casa_coord.txt': '%s, %s, 6.5\n' % (ra_c, dec_c),
        'colon_coord.txt': '%s, %s, 6.5\n' % (ra_c, dec_c.replace('.', ':', 2)),
        'deg_coord.txt': '%r, %r, 6.5, 2.0\n' % (rad, decd),
        'hms_coord.txt': '%s, %s, 6.5\n' % (sky.ra.to_string('hourangle', sep='hms'), sky.dec.to_string(sep='dms', alwayssign=True)),
        'space_coord.txt': '# comment line\n\n%r   %r   6.5\n' % (rad, decd),
        'tab_coord.txt': '%r\t%r\t6.5\n' % (rad, decd),
    }
    for name, text in lines.items():
        with open(name, 'w') as fh:
            fh.write(text)
        c = LineStacker.readCoords([name])
        assert len(c) == 1, name
        assert c[0].spec == 6.5 and c[0].z == 6.5
        sep = LineStacker._to_skycoord(c).separation(sky).arcsec[0]
        assert sep < 0.01, (name, sep)
    assert LineStacker.readCoords(['deg_coord.txt'])[0].weight == 2.0


def test_readcoords_pixel_and_continuum(workdir):
    with open('pix_coord.txt', 'w') as fh:
        fh.write('10.5, 12, 3\n11, 13, 4, 0.5\n')
    c = LineStacker.readCoords(['pix_coord.txt'], type='pixel')
    assert [(x.x, x.y, x.spec, x.weight, x.type) for x in c] == [(10.5, 12., 3., 1., 'pixel'), (11., 13., 4., .5, 'pixel')]
    with open('cont_coord.txt', 'w') as fh:
        fh.write('53.1, -27.8\n53.2, -27.7, 3.0\n')
    c = LineStacker.readCoords(['cont_coord.txt'], lineON=False)
    assert c[0].spec == 'continuum' and c[1].weight == 3.0


def test_readcoords_bad_file_message(workdir):
    with open('bad_coord.txt', 'w') as fh:
        fh.write('not, a, number\n')
    with pytest.raises(Exception) as err:
        LineStacker.readCoords(['bad_coord.txt'])
    assert 'bad_coord.txt' in str(err.value)


def test_image_matching_by_file_order(workdir):
    """Coordinates whose image name is not in imagenames are matched by
    coordinate-file order (the documented examples rely on this)."""
    names, cf = [], []
    for i in range(3):
        cube = Cube(nchan=30)
        cube.add_line(20, 20, 15, 1. + i, 2.)
        cube.write('img%d.fits' % i)
        names.append('img%d.fits' % i)
        ra, dec = casa_style(cube.sky(20, 20))
        with open('obj%d_coord.txt' % i, 'w') as fh:
            fh.write('%s, %s, %r\n' % (ra, dec, FEM / (232e9 - 15 * 7.6e7) - 1))
        cf.append('obj%d_coord.txt' % i)
    coords = LineStacker.readCoords(cf)        # default images: obj?.fits, not on disk
    res = li.stack(coords, imagenames=names, fEm=FEM, stampsize=8, N_chans=9, outfile=None)
    assert res[0] == pytest.approx(2., rel=1e-6)
    # mixing physical and pixel coordinates
    mixed = [coords[0], LineStacker.Coord(20., 20., spec=coords[1].spec, image='img1.fits', type='pixel'), coords[2]]
    res2 = li.stack(mixed, imagenames=names, fEm=FEM, stampsize=8, N_chans=9, outfile=None)
    np.testing.assert_allclose(res2[2], res[2], rtol=1e-6)


def test_extract_many_matches_single_at_edges(workdir):
    cube = Cube(nx=30, ny=30, nchan=12)
    cube.arr = np.random.default_rng(0).normal(size=cube.arr.shape)
    cube.write('e.fits')
    from LineStacker import _io
    arr, info = _io.read_cube('e.fits')
    xs = np.array([5., 10.3, 20.7]); ys = np.array([6., 14.5, 22.2])
    edges = 0
    for centre, step in ((2., 1.), (10., -1.), (5.4, 0.7), (1., 1.3)):
        many = li._extract_many(arr, info, xs, ys, centre, step, 5, 15)
        for r in range(3):
            one, inside = li._extract(arr, info, xs[r], ys[r], centre, step, 5, 15)
            np.testing.assert_array_equal(many[r], one)
        edges += not inside.all()
    assert edges == 3


def test_failed_fits_become_nan(monkeypatch):
    import LineStacker.tools.fit as fitTools
    ims = [OneD.Spectrum(amp=np.exp(-(np.arange(30) - 15.)**2 / 8), centerIndex=15, velocities=np.arange(30.))
           for _ in range(4)]
    real = fitTools.GaussFit
    calls = {'n': 0}

    def flaky(*a, **k):
        calls['n'] += 1
        if calls['n'] % 2 == 0:
            raise RuntimeError('Optimal parameters not found')
        return real(*a, **k)
    monkeypatch.setattr(fitTools, 'GaussFit', flaky)
    with pytest.warns(UserWarning, match='did not converge'):
        res = A.bootstraping_OneD(ims, nRandom=5, save='ampAndWidth', seed=0)
    assert len(res) == 6 and sum(np.isnan(r[0]) for r in res) == 3


def test_round_half_up():
    assert [LineStacker.round_half_up(v) for v in (0.5, 1.5, 2.5, -0.5, -2.5, 2.4999)] == [1, 2, 3, -1, -3, 2]


def test_astropy_numpy_mismatch_is_detected():
    msg = LineStacker._astropy_numpy_problem('2.4.0', '7.1.1', '/venv/site-packages/astropy')
    assert 'astropy 7.1.1' in msg and 'NumPy 2.4.0' in msg and "astropy>=7.2" in msg
    assert LineStacker._astropy_numpy_problem('2.5.4', '7.2.0', 'x') is None
    assert LineStacker._astropy_numpy_problem('2.3.5', '6.1.7', 'x') is None
    assert LineStacker._astropy_numpy_problem('1.26.4', '5.3', 'x') is None
    assert LineStacker._astropy_numpy_problem('2.4.0rc1', '7.1.0', 'x') is not None
    assert LineStacker._astropy_numpy_problem() is None      # this environment is fine

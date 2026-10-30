"""Catalogue (readCatalog) and array (makeCoords) inputs."""
import os

import numpy as np
import pytest

from conftest import Cube, FEM, casa_style

import LineStacker
import LineStacker.image as li
import LineStacker.analysisTools as A


def sample(n=6, seed=0, folder='.'):
    """n cubes with one line each; returns coord files, cube names, rows
    (ra_str, dec_str, z, ra_deg, dec_deg)."""
    rng = np.random.default_rng(seed)
    os.makedirs(folder, exist_ok=True)
    cf, names, rows = [], [], []
    for i in range(n):
        cube = Cube(nchan=40)
        ch = int(rng.integers(12, 28))
        z = FEM / (232e9 - ch * 7.6e7) - 1
        x, y = float(rng.integers(14, 26)), float(rng.integers(14, 26))
        cube.add_line(x, y, ch, rng.uniform(1, 3), 2.5)
        cube.add_noise(0.05, rng)
        name = os.path.join(folder, 'cube%02d.fits' % i)
        cube.write(name)
        sky = cube.sky(x, y)
        ra, dec = casa_style(sky)
        with open(os.path.join(folder, 'cube%02d_coord.txt' % i), 'w') as fh:
            fh.write('%s, %s, %r\n' % (ra, dec, z))
        cf.append(os.path.join(folder, 'cube%02d_coord.txt' % i))
        names.append(name)
        rows.append((ra, dec, z, float(sky.ra.deg), float(sky.dec.deg)))
    return cf, names, rows


KW = dict(fEm=FEM, stampsize=8, N_chans=11, outfile=None)


def reference(cf, names):
    return li.stack(LineStacker.readCoords(cf), imagenames=names, **KW)[2]


def test_csv_catalogue_matches_coordinate_files(workdir):
    cf, names, rows = sample()
    ref = reference(cf, names)
    with open('catalog.csv', 'w') as fh:
        fh.write('ra,dec,z,image\n')
        for (ra, dec, z, _, _), name in zip(rows, names):
            fh.write('%s,%s,%r,%s\n' % (ra, dec, z, name))
    coords = LineStacker.readCatalog('catalog.csv')
    assert coords[0].image == os.path.abspath(names[0])
    np.testing.assert_array_equal(li.stack(coords, **KW)[2], ref)                    # no imagenames
    np.testing.assert_array_equal(li.stack(coords, imagenames=names, **KW)[2], ref)  # explicit, same order
    np.testing.assert_array_equal(li.stack(coords, imagenames=names[::-1], **KW)[2], ref)  # any order


def test_catalogue_formats_and_aliases(workdir):
    from astropy.table import Table
    import astropy.units as u
    cf, names, rows = sample(n=4, seed=1)
    ref = reference(cf, names)
    ra = [r[3] for r in rows]
    dec = [r[4] for r in rows]
    z = [r[2] for r in rows]
    t = Table({'RAJ2000': ra * u.deg, 'DEJ2000': dec * u.deg, 'redshift': z, 'cube': names})
    t.write('cat.fits', overwrite=True)
    t.write('cat.ecsv', overwrite=True)
    t.write('cat.vot', format='votable', overwrite=True)
    for path in ('cat.fits', 'cat.ecsv', 'cat.vot'):
        np.testing.assert_allclose(li.stack(LineStacker.readCatalog(path), **KW)[2], ref, rtol=1e-6, atol=1e-7)
    # units are honoured: the same positions in arcmin
    t2 = Table({'ra': (np.array(ra) * 60) * u.arcmin, 'dec': (np.array(dec) * 60) * u.arcmin, 'z': z, 'image': names})
    np.testing.assert_allclose(li.stack(LineStacker.readCatalog(t2), **KW)[2], ref, rtol=1e-6, atol=1e-7)
    # space separated text with a commented header, explicit column names
    with open('cat.txt', 'w') as fh:
        fh.write('# alpha_deg delta_deg zz file\n')
        for r, name in zip(rows, names):
            fh.write('%r %r %r %s\n' % (r[3], r[4], r[2], name))
    c = LineStacker.readCatalog('cat.txt', ra='alpha_deg', dec='delta_deg', z='zz', image='file')
    np.testing.assert_allclose(li.stack(c, **KW)[2], ref, rtol=1e-6, atol=1e-7)
    # radians
    trad = Table({'ra': np.deg2rad(ra), 'dec': np.deg2rad(dec), 'z': z, 'image': names})
    np.testing.assert_allclose(li.stack(LineStacker.readCatalog(trad, unit='rad'), **KW)[2], ref, rtol=1e-6, atol=1e-7)
    # dict of columns
    d = {'ra': ra, 'dec': dec, 'z': z, 'image': names}
    np.testing.assert_allclose(li.stack(LineStacker.readCatalog(d), **KW)[2], ref, rtol=1e-6, atol=1e-7)


def test_pandas_dataframe(workdir):
    pd = pytest.importorskip('pandas')
    cf, names, rows = sample(n=3, seed=2)
    df = pd.DataFrame({'ra': [r[3] for r in rows], 'dec': [r[4] for r in rows], 'z': [r[2] for r in rows],
                       'image': names})
    np.testing.assert_allclose(li.stack(LineStacker.readCatalog(df), **KW)[2], reference(cf, names),
                               rtol=1e-6, atol=1e-7)


def test_weights_and_relative_paths(workdir):
    cf, names, rows = sample(n=4, seed=3, folder='cubes')
    with open('cubes/catalog.csv', 'w') as fh:            # names relative to the catalogue
        fh.write('ra,dec,z,weight,image\n')
        for k, (r, name) in enumerate(zip(rows, names)):
            fh.write('%s,%s,%r,%r,%s\n' % (r[0], r[1], r[2], k + 1., os.path.basename(name)))
    coords = LineStacker.readCatalog('cubes/catalog.csv')
    assert [c.weight for c in coords] == [1., 2., 3., 4.]
    assert all(os.path.exists(c.image) for c in coords)
    res = li.stack(coords, **KW)
    plain = LineStacker.readCoords(cf)
    for c, w in zip(plain, [1., 2., 3., 4.]):
        c.weight = w
    np.testing.assert_array_equal(res[2], li.stack(plain, imagenames=names, **KW)[2])
    with open('catalog2.csv', 'w') as fh:                 # names relative to the working directory
        fh.write('ra,dec,z,image\n')
        for r, name in zip(rows, names):
            fh.write('%s,%s,%r,%s\n' % (r[0], r[1], r[2], name))
    assert all(os.path.exists(c.image) for c in LineStacker.readCatalog('catalog2.csv'))
    c3 = LineStacker.readCatalog('catalog2.csv', imagedir='/somewhere')
    assert c3[0].image == os.path.join('/somewhere', names[0])


def test_make_coords_inputs(workdir):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    cf, names, rows = sample(n=5, seed=4)
    ref = reference(cf, names)
    ra = np.array([r[3] for r in rows])
    dec = np.array([r[4] for r in rows])
    z = np.array([r[2] for r in rows])
    for coords in (LineStacker.makeCoords(ra, dec, z=z, image=names),
                   LineStacker.makeCoords(list(ra), list(dec), z=list(z), image=names),
                   LineStacker.makeCoords(ra * u.deg, dec * u.deg, z=z, image=names),
                   LineStacker.makeCoords(SkyCoord(ra, dec, unit='deg'), z=z, image=names),
                   LineStacker.makeCoords([r[0] for r in rows], [r[1] for r in rows], z=z, image=names),
                   LineStacker.makeCoords(np.deg2rad(ra), np.deg2rad(dec), z=z, image=names, unit='rad')):
        np.testing.assert_allclose(li.stack(coords, **KW)[2], ref, rtol=1e-6, atol=1e-7)
    with pytest.raises(ValueError, match='z has 3 values'):
        LineStacker.makeCoords(ra, dec, z=z[:3], image=names)
    with pytest.raises(ValueError, match='dec is missing'):
        LineStacker.makeCoords(ra)


def test_single_cube_with_many_sources(workdir):
    """Deep-field case: one cube, many sources, no image column."""
    rng = np.random.default_rng(5)
    cube = Cube(nx=80, ny=80, nchan=60)
    pos, zs = [], []
    for k in range(6):
        x, y, ch = 15 + 10 * k, float(rng.integers(20, 60)), int(rng.integers(15, 45))
        cube.add_line(x, y, ch, 2., 2.)
        sky = cube.sky(x, y)
        pos.append((float(sky.ra.deg), float(sky.dec.deg)))
        zs.append(FEM / (232e9 - ch * 7.6e7) - 1)
    cube.write('field.fits')
    coords = LineStacker.makeCoords([p[0] for p in pos], [p[1] for p in pos], z=zs)
    res = li.stack(coords, imagenames=['field.fits'], **KW)
    assert res[0] == pytest.approx(2., rel=1e-5)
    same = LineStacker.makeCoords([p[0] for p in pos], [p[1] for p in pos], z=zs, image='field.fits')
    np.testing.assert_array_equal(li.stack(same, **KW)[2], res[2])
    # several cubes but no cube per source: refuse to guess
    Cube(nx=80, ny=80, nchan=60).write('other.fits')
    with pytest.raises(ValueError, match='not in imagenames'):
        li.stack(coords, imagenames=['field.fits', 'other.fits'], **KW)
    with pytest.raises(ValueError, match='does not name its cube'):
        li.stack(coords, **KW)


def test_catalogue_cube_must_be_listed(workdir):
    cf, names, rows = sample(n=3, seed=6)
    coords = LineStacker.makeCoords([r[3] for r in rows], [r[4] for r in rows], z=[r[2] for r in rows],
                                    image=names)
    # a cube missing from imagenames is an error, never matched by position in the list
    with pytest.raises(ValueError, match='not in imagenames'):
        li.stack(coords, imagenames=names[:2], **KW)


def test_continuum_catalogue_and_analysis_tools(workdir):
    cf, names, rows = sample(n=5, seed=7)
    with open('cont.csv', 'w') as fh:
        fh.write('ra,dec,image\n')
        for r, name in zip(rows, names):
            fh.write('%r,%r,%s\n' % (r[3], r[4], name))
    coords = LineStacker.readCatalog('cont.csv')
    assert all(c.spec == 'continuum' for c in coords)
    with pytest.warns(UserWarning, match='coordinates have no spectral value'):
        res = li.stack(coords, stampsize=8, outfile=None)
    assert res[2].shape == (8, 8, 1, 1)
    line = LineStacker.makeCoords([r[3] for r in rows], [r[4] for r in rows], z=[r[2] for r in rows], image=names)
    boot = A.bootStraping_Cube(line, fEm=FEM, stampsize=8, N_chans=11, nRandom=10, seed=1)
    assert len(boot) == 11
    est = A.stack_estimator(line, fEm=FEM, N_chans=11, nRandom=5, seed=1, lowerLimit=2e-5, upperLimit=3e-5)
    assert len(est[2]) == 5
    assert len(A.rebin_CubesSpectra(line, outdir='rebinned', widths=[2., 2.5, 3., 2., 2.])) == 5
    assert sorted(os.listdir('rebinned'))[0] == 'cube00_SpectralRebinned.fits'


def test_catalogue_errors(workdir):
    with open('bad.csv', 'w') as fh:
        fh.write('x,y,z\n1,2,3\n')
    with pytest.raises(ValueError, match="no column 'ra'"):
        LineStacker.readCatalog('bad.csv')
    with open('miss.csv', 'w') as fh:
        fh.write('ra,dec,z,image\n53.1,-27.8,6.0,a.fits\n53.2,-27.8,,b.fits\n')
    with pytest.raises(ValueError, match='missing values'):
        LineStacker.readCatalog('miss.csv')
    with pytest.raises(FileNotFoundError):
        LineStacker.readCatalog('nope.csv')
    with pytest.raises(ValueError, match="no column 'redshift_x'"):
        LineStacker.readCatalog('bad.csv', ra='x', dec='y', z='redshift_x')


def test_read_cube(workdir):
    cube = Cube(nx=20, ny=30, nchan=12)
    cube.arr = np.random.default_rng(0).normal(size=cube.arr.shape)
    cube.write('r.fits')
    arr = li.readCube('r.fits')
    assert arr.shape == (20, 30, 1, 12)
    np.testing.assert_allclose(arr[:, :, 0, :], cube.arr, rtol=1e-6)
    np.testing.assert_allclose(li.readCube('r.fits', collapse='sum')[:, :, 0, 0], cube.arr.sum(axis=2), rtol=1e-5, atol=1e-5)

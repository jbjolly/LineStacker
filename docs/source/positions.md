# Stacking positions

Every stacking function takes a list of {class}`LineStacker.Coord`, one per
source. There are three ways to make it; all give the same coordinates.

## From a catalogue

{func}`LineStacker.readCatalog` reads the whole sample from one table, one
row per source:

```text
ra,dec,z,weight,image
53.1475,-27.7826,6.6036,1,cubes/source_01.fits
53.1612,-27.7791,6.5899,1,cubes/source_02.fits
03:32:41.31,-27.46.12.4,6.6210,2,cubes/source_03.image
```

```python
coords = LineStacker.readCatalog('sources.csv')
```

- **Formats**: CSV (with a header line), ECSV, FITS tables, VOTable and
  anything else `astropy.table` reads (`format=` if it cannot be guessed),
  or an astropy `Table`, a pandas `DataFrame` or a dict of columns.
- **Columns**: `ra`, `dec`, `z` (optional), `weight` (optional) and `image`.
  Names are matched without regard to case, and common alternatives are
  found: `RAJ2000`/`DEJ2000`, `ra_deg`/`dec_deg`, `redshift`, `zspec`,
  `cube`, `file`... Other names can be given:
  `readCatalog(cat, ra='alpha', z='zspec', image='file')`.
- **Positions** are in degrees (or `unit='rad'`), in the unit of the column
  when the table has units (FITS, ECSV, VOTable), or sexagesimal strings with
  RA in hours (`03:32:41.31`, `03h32m41.31s`, or CASA's dotted declination
  `-027.46.12.4`). Each value is read on its own, so one column can mix
  forms. The frame is ICRS unless `frame='fk5'` (J2000) or another
  equatorial frame is given.
- **z** is the spectral value of each source: its redshift for the default
  `spectralMethod='z'`, or the observed frequency, wavelength or velocity
  for `'obsSpec'`, or the channel for `'channel'`. Without a `z` column the
  positions are for continuum stacking.
- **image** names the cube of each source. Relative paths are looked up next
  to the catalogue file, then from the working directory; `imagedir=` sets
  another folder.

The cube of a source is never guessed from its position: it is the one the
catalogue names. If the catalogue has no `image` column, the stacking
function must be given a single cube (`imagenames=['field.fits']`), which
then holds every source.

## From Python arrays

{func}`LineStacker.makeCoords` builds the same coordinates from lists, NumPy
arrays, astropy quantities or a `SkyCoord`, for example in a notebook:

```python
import numpy as np
from astropy.coordinates import SkyCoord

coords = LineStacker.makeCoords(ra_deg, dec_deg, z=redshifts,
                                image=cube_names)
coords = LineStacker.makeCoords(SkyCoord(ra, dec, unit='deg'), z=redshifts,
                                image='field.fits')
coords = LineStacker.makeCoords(x_pix, y_pix, z=channels, image='cube.fits',
                                type='pixel')
```

A single value (`image='field.fits'`, `weight=1`) is used for every source.
A `SkyCoord` in another frame, such as Galactic, is converted.

## From coordinate files

{func}`LineStacker.readCoords` reads the format of LineStacker 1.x: one text
file per cube, one line per source, with `ra, dec, z` and optionally a
weight (`ra, dec` and an optional weight with `lineON=False` for continuum):

```text
03:30:34.79, -027.45.34.566, 7.22611675851
53.1475, -27.7826, 6.6036, 2.0
```

```python
coords = LineStacker.readCoords(['cube_1_coord.txt', 'cube_2_coord.txt'])
LineStacker.image.stack(coords, imagenames=['cube_1.image', 'cube_2.image'],
                        fEm=fEm)
```

Positions can be in degrees, in sexagesimal with RA in hours (colons,
spaces, `h m s` letters or CASA's dots), or pixels with `type='pixel'`.
Commas, spaces and tabs are tried as delimiters; empty lines and lines
starting with `#` are skipped. The sources of each file belong to the
image of the same name (`cube_1_coord.txt` to `cube_1.fits`, or the name
given with `images=`) if that image is in `imagenames`, and otherwise to the
image at the same position in `imagenames` as the file in the list.

## Which cube a source belongs to

The stacking functions take the list of cubes as `imagenames`. When every
coordinate names its cube (catalogue `image` column, `makeCoords(image=...)`),
`imagenames` can be left out: the cubes are those the coordinates name.
A coordinate whose cube is not in `imagenames` is an error, so a source is
never stacked from the wrong cube.

`coords=None` stacks the centre of every cube in `imagenames` (spatially and
spectrally).

## Weights

The `weight` of a coordinate multiplies its contribution to mean stacks;
the automatic weightings of {func}`LineStacker.image.stack` (`'sigma2'`,
`'1/A'`...) multiply it. Weights must be finite and not negative; a weight
of 0 leaves the source out of mean stacks. Medians ignore weights.

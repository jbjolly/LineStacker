# Stacking cubes

{func}`LineStacker.image.stack` stacks data cubes, FITS files or CASA
images (both can be mixed), around the sources of a coordinate list.

```python
import LineStacker
import LineStacker.image

coords = LineStacker.readCatalog('sources.csv')
stacked = LineStacker.image.stack(
    coords,
    fEm=1897420620253.1646,   # rest frequency of the line (Hz)
    stampsize=16,             # pixels
    N_chans=31,               # channels of the stack
    outfile='stack.fits')
peak, spectrum, cube, missing = stacked
```

It returns the central pixel at the central channel, the spectrum of the
central pixel, the whole stacked cube as an array ordered (x, y, Stokes,
channel), and the number of sources missing from each channel (their cube
does not cover it).

## Centring the lines

`spectralMethod` says what the spectral value of each coordinate (its `z`)
is:

| `spectralMethod` | spectral value | also needed |
|---|---|---|
| `'z'` (default) | redshift of the line | `fEm` (frequency cubes) or `lEm` (wavelength cubes) |
| `'obsSpec'` | observed frequency, wavelength or velocity of the line | |
| `'channel'` | channel of the line (0-based) | |

`fEm`, `lEm` and observed values are in the spectral unit of the first cube
(Hz for CASA images, Angstrom for MUSE cubes, for example). The cubes must
all have the same kind of spectral axis (all frequency, or all wavelength,
or all radio velocity); other units are converted (a cube in GHz stacks with
cubes in Hz). Cubes whose spectral axis runs the other way are read
backwards, so all spectra share the orientation of the first source's
cube. Non-linear spectral axes (`WAVE-LOG`, `FREQ-TAB`...) are refused for
line stacks, whatever `spectralMethod`, since the stack could not be given a
correct spectral axis; regrid such cubes to a linear axis first. Cubes in
different spectral reference frames (`SPECSYS`: LSRK, TOPOCENT...), or with
`SPECSYS` in some cubes only, give a warning: LineStacker does not convert
frames, so regrid such cubes to one frame first (with CASA's `imregrid`, or
`tclean(outframe=...)`).

### Intensity units

The cubes' pixel units (`BUNIT`) are compared before stacking. Units that
differ only by a scale are converted to the unit of the first cube that
has one (a cube in mJy/beam stacks with cubes in Jy/beam); recognised but
physically different units (K and Jy/beam, Jy/pixel and Jy/beam) are an
error. Cubes without `BUNIT`, or with a unit astropy cannot read, give a
warning and are stacked without conversion.

`N_chans` channels are stacked around the line centre (default: the number
of channels of the first cube). By default each source contributes its
nearest channel to each stack channel. With `regridFromZ=True` (needs
`spectralMethod='z'`) every spectrum is resampled by linear interpolation
so that all sources share the same rest-frame channel width: the width of
the source at the highest redshift (`regridMethod='scaleToMax'`, default) or
the lowest (`'scaleToMin'`). This corrects the narrowing of lines of the
same rest-frame width at higher redshift. A stack channel interpolated from
a blank channel is blank for that source; one that falls exactly on a
channel takes its value.

## Weights

With `method='mean'` each stamp is weighted, channel by channel, by the
weight of its coordinate times the weighting asked for:

| `weighting` | weight multiplied by |
|---|---|
| `None` (default) | 1: only the coordinates' weights |
| `'sigma2'` | 1/σ², σ being the standard deviation of the source's stamp |
| `'sigma2F'` | 1/σ² computed channel by channel |
| `'sigma'` | 1/σ (as for `'sigma2'`) |
| `'1/A'` | 1/A, A being the line amplitude: the maximum of the central-pixel spectrum. Only for lines visible source by source |
| `'pb'` | the square of the primary-beam response at the source (see below) |

σ is computed from the stamp pixels that are not blank and lie outside
`maskradius` pixels of the source; set `maskradius` so that the source does
not count as noise. For `'sigma'` the default mask radius is twice the beam's
major axis.
Sources (or channels) without valid noise pixels, or without a positive line
amplitude for `'1/A'`, get weight 0, with a warning. Weights must be finite
and not negative.

`method='median'` takes, for each pixel and channel, the median over the
sources whose cube covers that channel. It ignores weights.

### Primary-beam weighting

For primary-beam corrected images, whose noise grows as 1/PB away from the
pointing centre, weighting by PB² is inverse-variance weighting.
`primarybeam` gives the primary beam of each cube:

```python
# primary-beam images (e.g. tclean's .pb), one per cube, read at each
# source's position and line channel
LineStacker.image.stack(coords, imagenames=cubes, fEm=fEm,
                        weighting='pb', primarybeam=pb_images)
LineStacker.image.stack(coords, imagenames=cubes, fEm=fEm, weighting='pb',
                        primarybeam={'a.image': 'a.pb', 'b.image': 'b.pb'})

# a Gaussian beam of 12 m dishes, FWHM 1.13 λ/D, centred on the reference
# pixel of each image
pb = LineStacker.image.GaussianPrimaryBeam(dish=12.)
LineStacker.image.stack(coords, imagenames=cubes, fEm=fEm,
                        weighting='pb', primarybeam=pb)
```

A list gives one primary beam per cube, in the order of `imagenames`; one
path is enough for a single cube. Any function `f(dx, dy, nu)` returning the
response at offsets `dx`, `dy` (radians, east and north on the sky, from the
position of the image reference pixel) and frequency `nu` (Hz) can be used
too.

For line stacks, a primary-beam cube is read at the channel of each
source's line, so it needs the same kind of spectral axis as the science
cube (a single-plane primary beam can always be used). For continuum
stacks, the response is the mean of the primary-beam image's non-blank
planes at the source.

## Continuum

`continuum=True` stacks images instead of lines. Images without a spectral
axis, or with a single channel, are always stacked as continuum, and so are
coordinates without a spectral value (a catalogue without `z`). Cubes with a
spectral axis are collapsed over their channels in memory, with a warning:
`collapse='sum'` (default), `'mean'` or `'median'`; blank channels are
ignored. Only the collapsed images are kept in memory. The result is the
same list as for lines, with a single channel; the FITS output is 2D when
the inputs had no spectral axis.

## Output

`outfile` names the stacked cube: a name ending in `.fits` (or `.fits.gz`)
is written as FITS, any other name (`stack.image`) as a CASA image, which
needs casatools. `outfile=None` writes nothing. With verbose on (default)
LineStacker prints the full path of the file it wrote.

The header of the stacked cube:

- Pixels without valid data in any source are NaN (blank).
- **Sky axes** are offsets from the stacked position, (0, 0) at the stamp
  centre, with the pixel size and orientation of the first cube (a rotated
  grid is written with a PC matrix).
- **Spectral axis**: for `spectralMethod='z'`, the rest frame: the line at
  `fEm` (or `lEm`), channels of the mean rest-frame width, `RESTFRQ` (or
  `RESTWAV`) set, and `SPECSYS='SOURCE'` if the cubes give a frame;
  otherwise an axis centred on the mean observed value of the lines.
- **Beam**: `BMAJ`, `BMIN` are the means of the input beams and `BPA` their
  mean as an axis (179° and 1° give 0°): an approximation of the stack's
  beam when the input beams differ.
- `BUNIT` of the first cube that has one (other cubes are converted to it,
  see above); for a continuum stack of collapsed cubes it is per channel
  ('sum' adds channels), as the `COLLAPSE` keyword says. A continuum stack
  of cubes with non-linear spectral axes has no spectral axis.
- `NSTACK` (number of sources) and a `HISTORY` line with the method.

LineStacker aligns sources on their centres but does not reproject: cubes
with a different pixel size, orientation, projection (SIN, TAN...),
celestial frame (ICRS, FK5...) or distortion terms from the first one give a
warning, and their stamps are cut on their own grid. Reproject such cubes to
a common grid first (CASA's `imregrid`) to use the stacked image spatially.
Even the central-pixel spectrum is only comparable between cubes of similar
beams and pixel sampling.

Other outputs: `returnWeights=True` writes `weightFile.csv` (ra, dec, z and
the weight of each source, per channel if they differ), `saveSubCubes=True`
saves the stamps as a NumPy file (`outSubCubes.npy`, with their pixel
positions in `outCoordsFile.txt`), and `plotIt=True` plots the stacked
spectrum and the number of sources per channel.

## Reading cubes

{func}`LineStacker.image.readCube` reads a FITS file or a CASA image as a
NumPy array ordered (x, y, Stokes, channel), x along RA and y along Dec, as
LineStacker uses it, for plots and checks:

```python
cube = LineStacker.image.readCube('stack.fits')
spectrum = cube[8, 8, 0, :]
```

FITS cubes can have their data in any HDU (`hdul='auto'` takes the first
one with an image, `hdul=1` for MUSE cubes for instance) and their axes in
any order. Cubes are read once and cached; a cube is read again when its
file changes (checking costs about 0.1 ms per CASA image each time a cube is
used; `LineStacker.cache_check = False` trusts the cache instead).
`LineStacker.image.clear_cache()` frees the memory.

Blank pixels (NaN in FITS, masked pixels of CASA images) stay NaN. Stacks
leave them out pixel by pixel: a stacked pixel is the mean (or median) of
the sources that have a valid value there, and NaN if none has, as are
stack channels that no source's cube covers. Pixels that are exactly 0 are
data like any other.

# LineStacker 2.0: Python 3 / CASA 6 port

Tested with Python 3.13 (requires >= 3.9), NumPy 2.5, SciPy 1.18, astropy 8.0 and CASA 6.7.6
(casatools from PyPI), on the example data in `LineStacker/Example`
and on the synthetic test suite in `tests/` (183 tests; the 12 that need casatools
are skipped without it).

## Installing

    pip install -e .[casa]        # with CASA 6 for CASA images
    pip install -e .              # FITS cubes only
    python -m pytest tests        # optional: run the tests

To use LineStacker from the CASA prompt (monolithic CASA tarball), install
its dependencies with CASA's own pip, so they match CASA's NumPy:

    PYTHONNOUSERSITE=1 ~/casa-6.7.6-14-py3.12.el9/bin/pip3 install 'astropy>=7.2'

Mixing CASA with packages from another environment (e.g. a venv on
PYTHONPATH) breaks when versions disagree: CASA 6.7.6 bundles NumPy >= 2.4,
which astropy < 7.2 cannot use (`module 'numpy' has no attribute 'in1d'`).
`import LineStacker` warns when it detects this combination.

CASA is only imported to read or write CASA images (FITS-only use never
starts CASA). When LineStacker starts it, CASA's log goes to one file,
`<temporary dir>/linestacker-casa-<user>.log` (appended to), and CASA's
messages go there instead of the terminal, so no `casa-<date>.log` piles up
in the working directory. `LineStacker.casa_logfile` changes this (a path,
`None` for no log, `'casa'` for CASA's default); a `logfile` set in
`~/.casa/config.py` is used, and inside a CASA session (or after importing
casatools/casatasks yourself) CASA's own settings apply. CASA 5 (Python 2)
is no longer supported.

## One stacking module: `LineStacker.image`

`line_image` and the old continuum module `image` are merged into a single
`LineStacker.image`, which does what `line_image` did, plus continuum
stacking with the same functions:

    LineStacker.image.stack(coords, imagenames=names, fEm=fEm)      # line stack
    LineStacker.image.stack(coords, imagenames=names, continuum=True)

- Images without a spectral axis (or with one channel) are stacked as
  continuum automatically, without warning.
- `continuum=True` on cubes with a spectral axis issues a warning and stacks
  a collapsed version of each cube, made in memory (`collapse='sum'`
  default, or `'mean'`, `'median'`; blanked channels are ignored). The
  output header records it in `COLLAPSE`. Only the
  collapsed arrays are kept, not the full cubes.
- Coordinates without a spectral value (`readCoords(..., lineON=False)`)
  also select continuum stacking (same warning if cubes get collapsed).
- Line stacking with an image that has no spectral axis raises an error
  pointing to `continuum=True`.
- Continuum stacks return the same list as line stacks (the "spectrum" has
  one channel); FITS output is 2D when the inputs had no spectral axis.
- From the old continuum module: weightings `'sigma'` (1/σ, mask of two
  beams by default) and `'pb'` (squared primary-beam response) are available
  to all stacks; `getFlux` is kept. `primarybeam=` now takes the
  primary-beam image of each cube (FITS or CASA, e.g. tclean's `.pb`; one
  path, a list in the order of `imagenames`, or a dict `{cube: pb}`), read at
  each source's position and line channel, or a function `f(dx, dy, nu)`
  such as `LineStacker.image.GaussianPrimaryBeam(dish=12.)` (FWHM 1.13 λ/D,
  as LineStacker.uv). The `LineStacker.pb` models are gone (see below).
  Its `noise` and `boot_straping` are covered by
  `analysisTools.stack_estimator` and `bootStraping_Cube`, which accept
  `continuum=True`. Dropped: `psfmode='star'`, `calculate_peak_fluxes` and
  `calculate_flux_weights` (`.flux` images).
- `import LineStacker.line_image` no longer works: use `LineStacker.image`
  (same `stack` arguments).

## What the stacked cube contains

- Stamps are cut on each cube's own pixel grid, centred on the source;
  cubes are not reprojected. The output keeps the first cube's pixel size
  and orientation (a rotated grid is written with a PC matrix).
- All cubes must have the same kind of spectral axis (FREQ, WAVE, AWAV,
  VRAD...); stacking frequency with wavelength cubes, or using a non-linear
  axis (`WAVE-LOG`, `FREQ-TAB`...) in a line stack (any `spectralMethod`), is
  an error. Units are converted to the first cube's (Hz and GHz cubes stack
  together); `fEm`, `lEm` and observed values are in that unit. Different
  `SPECSYS` frames, or `SPECSYS` missing in some cubes, give a warning:
  LineStacker does not convert frames.
- Pixel units (`BUNIT`) are compared: scaled units are converted to the
  first cube's (mJy/beam to Jy/beam), recognised but physically different
  units (K and Jy/beam, Jy/pixel and Jy/beam) are an error, and cubes
  without `BUNIT`, or with units astropy cannot read, give a warning and
  are stacked without conversion.
- Cubes on a different grid from the first one (pixel size, orientation,
  projection, celestial frame, distortions) give a warning; stamps are not
  reprojected.
- BMAJ, BMIN are the mean of the input beams and BPA their axial mean (179°
  and 1° give 0°), an approximation of the stack's beam.
- Weights must be finite and >= 0 (a negative or NaN weight is an error).
  Every weighting multiplies the coordinates' own weights.

## 1D stacking: `Image` is now `Spectrum`

The object holding one spectrum to stack is `LineStacker.OneD_Stacker.Spectrum`
(was `Image`), with the same arguments and attributes. The 1D functions take
`spectra` instead of `Images`/`images`: `Stack(spectra, ...)`,
`bootstraping_OneD(spectra, ...)`, `subsample_OneD(spectra, ...)`,
`rebin_OneD(spectra)`, `regridFromZ1D(spectra, ...)`. Positional calls are
unchanged; calls passing `Images=` by keyword need updating.

Velocities are in km/s. A spectrum with velocities and a redshift
(`Spectrum(spectrum=..., z=z, fEmLine=f0)`) has velocities relative to the
rest frequency `f0` (as a CASA spectrum with `restfreq=f0`); its line is at
v = c z/(1+z) in the radio convention (default, as CASA and LineStacker.uv),
c z with `velConvention='optical'`, or the relativistic velocity with
`'relativistic'`. `velocity_to_frequency` and `frequency_to_velocity` do the
conversions.

## Stacking positions: catalogue, arrays or coordinate files

- `LineStacker.readCatalog(catalog)`: one table for the whole sample, one row
  per source: `ra`, `dec`, `z` (optional; without it the positions are for
  continuum stacking), `weight` (optional) and `image`, the cube of each
  source. Reads CSV, ECSV, FITS tables, VOTable (any astropy table format),
  an astropy `Table`, a pandas `DataFrame` or a dict. Column names are
  matched case-insensitively, common alternatives are recognised (`RAJ2000`,
  `DEJ2000`, `redshift`, `cube`, ...), and any name can be given
  (`readCatalog(cat, ra='alpha', z='zspec', image='file')`). Positions in
  degrees (or `unit='rad'`), with units from the table, or sexagesimal (RA in
  hours, CASA's dotted declinations accepted); each value is read on its own,
  so a column mixing forms is read correctly. Relative cube names are looked
  up next to the catalogue, then from the working directory (`imagedir=`
  overrides).
- `LineStacker.makeCoords(ra, dec, z=..., weight=..., image=...)`: the same
  from lists, numpy arrays, astropy quantities or a `SkyCoord` (notebooks).
- The cube of a source is never guessed: it comes from the catalogue/`image=`,
  or, when no cube is named, from a single cube given as `imagenames`. A cube
  missing from `imagenames` is an error. `imagenames` can be omitted when every
  source names its cube.
- `readCoords` (one file per cube) works as before.
- `LineStacker.image.readCube(name)` reads a FITS file or CASA image as a numpy
  array (x, y, stokes, channel), for plots and checks.
- `rebin_CubesSpectra(..., outdir=)` writes the rebinned cubes to a folder.

## Examples notebook

`LineStacker/Example/LineStacker_examples.ipynb` runs every example (1D
stacking, bootstrapping and subsampling; cube stacking from coordinate files,
a catalogue and arrays; line fits, rebinning and bootstrapping; continuum
stacking and random positions) on the example data in about 30 s, writing to
`Example/notebook_outputs/`. It is saved with its outputs (figures).

## Messages: `verbose`

`LineStacker.verbose` (default `True`) makes LineStacker report what it does:
positions read, the full path of every file written (stacked cube, weights,
sub-cubes, rebinned or regridded cubes) and progress bars. Each function that
prints takes `verbose=True/False` to override it for one call. With verbose
off LineStacker prints nothing; problems (zero-noise stamps, collapsed cubes,
failed fits...) are Python warnings and show either way.

## What works now that did not

- **CASA images as input.** `line_image` had been rewritten to read FITS only,
  so the shipped examples (CASA `.image` cubes) failed with
  `IsADirectoryError`. All modules now share one reader
  (`LineStacker/_io.py`) that takes FITS files or CASA images (exported
  through casatools; both kinds can be mixed in one stack), in any axis order (RA/Dec swapped, Stokes before or after
  the spectral axis, no Stokes, 2D images), with `CDELT` or `CD` matrices and
  frequency or wavelength (`FREQ`, `WAVE`, `AWAV`) axes. Cubes are read once
  and cached (`LineStacker._io.clear_cache()` frees them).
- **CASA 6 tools.** The shim imported `casalog` from `casatools`, which CASA 6
  does not export, so it failed on every real CASA 6 installation; its CASA 5
  fallback imported itself. Both fixed and checked against casatools 6.7.6.
- **Output cubes.** `outfile='x.fits'` writes FITS, any other name (e.g.
  `stackResult.image`) a CASA image with axes RA, Dec, Stokes, Frequency as
  before. The header is rebuilt cleanly instead of averaging every keyword of
  every input header.
- **NumPy 2.** `np.NaN`, `np.float`, `np.cast`, `np.rank`, list indexing and
  int-array in-place float maths were removed or broke `OneD_Stacker.Stack`,
  `congrid`, `regridFromZ1D`, `rebin_OneD` and spline regridding. Comparing a
  weights array with `'1/A'` raised "truth value is ambiguous".
- **Modules that used the removed Stacker `CoordList` model**
  (`image` (continuum stacking), `bootStraping_Cube`, `stack_estimator`,
  `rebin_CubesSpectra`, cube `regridFromZ`, and `randomCoords` /
  `randomizeCoords`, which only printed 'TODO') are ported to `Coord(ra, dec, spec)`. `Coord` also answers
  to `x`, `y`, `z` for older scripts.
- **regridFromZ** in `line_image.stack` was silently ignored; it now resamples
  each source to a common rest-frame channel (linear interpolation, fractional
  line centres), for frequency and wavelength cubes.
- `readCoords`: CASA `-027.40.59.866` declinations, blank/comment lines,
  delimiter fallback (comma, space, tab), `type='pixel'` (was ignored).
  Coordinate files are matched to `imagenames` by name, or else by file order.
- **Gaussian fits** (`LineStacker.tools.fit`) ignore NaN channels and raise
  `FitError` (a `RuntimeError`) for too few points, constant data or a fit
  that does not converge; the bootstraps and grades turn those into NaN or 0.
  `DoubleGaussPeak` is centred on the middle of its axis (it assumed an axis
  starting at 0).
- **Input checks** with clear messages: empty spectrum or image lists,
  negative or non-integer `nRandom`, `minSize`/`maxSize` of `subsample_OneD`,
  `regridMethod`, rebinning widths, `randomCoords` margins, empty cube names
  and non-finite `z`/weights in catalogues, equatorial frames only
  (`makeCoords(SkyCoord)` converts e.g. Galactic positions). Coordinates in
  different frames (ICRS, FK5...) are each converted.
- **Cached cubes** are read again when their file changes (for CASA images,
  any file of the image directory); the stale copy is dropped.
- Smaller crashes: `KeyError: 'freqlen'` when a line was near the upper edge
  of a cube, `coord.z` in obsSpec mode, `returnWeights`, `plotIt`,
  `stampsize=1` in `stack_estimator`, non-square images in `noise_estimator`,
  float slice indices in `rebin_CubesSpectra`, `freqBin` never rescaled when
  rebinning 1D spectra, `LineStacker.tools.GaussFit` typo in the 1D bootstrap.

## Changes in results (please check against your expectations)

1. **Stack centre moves by one channel** for FITS input. The FITS rewrite used
   the 1-based `CRPIX3` where the CASA code used CASA's 0-based reference pixel.
   The new code agrees with CASA's own `topixel`.
2. Line centres are rounded half away from zero, as Python 2's `round()` did
   (Python 3 rounds halves to even).
3. Median stacks ignore channels that fall outside a source's cube instead of
   counting them as zeros. Mean stacks already gave those channels zero weight.
4. `bootStraping_Cube` 'amp'/'ampAndWidth'/'outflow' now use the spectrum
   summed over the stamp, as documented. The code summed axes (1,2,3), i.e. a
   spatial profile.
5. Weights are computed on copies: the user's `Coord` objects are not
   modified, and `weighting='sigma2'` no longer compounds when the same
   coordinates are stacked twice.
6. Cubes whose spectral axis runs the other way (positive vs negative CDELT3)
   are flipped to match the first cube instead of being stacked mirrored.
7. Channels whose total weight is zero give NaN (blank) instead of
   ZeroDivisionError; sigma2F channels with zero noise get weight 0 instead
   of infinity.
8. Gaussian fits that do not converge are stored as `[nan, nan]` (with a
   warning) instead of aborting a whole bootstrap.
9. `hdul` defaults to `'auto'` (first HDU with an image: 0 for ALMA/CASA FITS,
   1 for MUSE). An explicit `hdul=1` still works.
10. `noise_estimator(maskCenter=[x, y])` uses x along RA and y along Dec as
    documented; it no longer blocks on `plt.show()` (`plotIt=False` disables
    the plot).
11. `maximizeSNR` (the default grade of `subsample_OneD`) returns 0 when the
    Gaussian fit converges on a line centred outside the spectrum (it can lock
    onto a stack edge). Its noise was then computed from a negative slice
    index, i.e. a few channels at the end, giving grades of ~1000 instead of
    ~4 that dominated the summed grades: the subsampling example no longer
    singled out the bright spectra. With the fix it does again (with 10 000
    iterations already). The fit's sigma is also used as an absolute value.
12. Random draws are vectorised, so a given `np.random.seed` gives different
    (statistically equivalent) samples than LineStacker 1.x. All Monte-Carlo
    tools also take `seed=`.
13. `weighting='sigma2F'` multiplies the coordinates' weights by 1/σ² (it
    replaced them, so catalogue weights were lost), and σ of each channel
    uses only the valid pixels outside the mask (masked and blank pixels
    were counted as zeros, which lowered σ). Unchanged for unit weights,
    no mask and no blank pixels, as in the examples.
14. `weighting='1/A'` gives weight 0 (with a warning) to sources whose
    amplitude is not positive, instead of a negative or infinite weight.
15. Weighted 1D stacks: channels where a spectrum's flux is NaN no longer
    count in the sum of weights (they pulled the mean towards 0).
16. Spectra given as amplitudes only, with an odd number of channels, lost
    their last channel in 'full' stacks (Python 3 rounding of the default
    velocity axis).
17. `Spectrum(z=..., fEmLine=...)` with velocities is centred on the line
    (`velToFreq` returned `fObs v/(v+c)`, not a frequency).
18. `regridFromZ` (writing cubes): regridded cubes have a rest-frame spectral
    axis (CDELT3 was the observed width, so the stack's CDELT3 was wrong) and
    keep the sky coordinates of the original cube (CRVAL was set to 0, 0).
    Its result is stacked as returned, with `spectralMethod='z'` and `fEm`.
19. Beam position angles are averaged as axes (BPA changes only when the
    input angles straddle ±90°).
20. Cubes in compatible but different intensity units (mJy/beam, Jy/beam)
    are converted to the first cube's unit before stacking (they were
    averaged as raw numbers).
21. `regridFromZ`: the coordinates returned for regridded cubes are the
    sources' own (fractional) pixels in the stamps, not the integer stamp
    centre, so their sky positions are exact.
22. 1D stacks leave out infinite fluxes, as they do NaN (one infinite
    channel made the stack infinite).

23. **Blank pixels stay NaN** (they were set to 0, as in LineStacker 1.x).
    Stacks leave them out pixel by pixel; stacked pixels without valid data,
    and stack channels that no source covers, are NaN instead of 0; pixels
    that are exactly 0 now count as data in the noise-based weightings and
    `noise_estimator`; collapsed images keep NaN where every channel is
    blank. `readCube` returns NaN for blank pixels. A blank channel blanks
    only the stack channels interpolated from it (a stack channel that falls
    exactly on a valid channel keeps its value). Rebinning and regridding
    never use blank channels: with linear interpolation (`rebin_OneD`,
    `rebin_CubesSpectra`) a new channel interpolated from a blank one is
    NaN; the cubic spline of `regridFromZ1D` and `regridFromZ` runs through
    each stretch of valid channels separately, so only new channels between
    a valid and a blank channel are NaN. Blanks thus stay local instead of
    spreading over the whole spectrum. The `maximizeSNR` grade measures the
    noise on valid channels only. 1D stacks also give NaN, instead of 0,
    for channels without data (no spectrum covers them, or only blank
    fluxes or zero weights), as cube stacks do.
24. **Regridding to a common rest-frame channel width** (`regridFromZ1D`,
    `Stack(regridFromZ=True)`, `regridFromZ`) samples old channel
    j (1+z_ref)/(1+z) for new channel j, so every source gets exactly the
    reference channel width (the number of channels was rounded and the
    width followed it, up to half a channel over the spectrum). New channels
    beyond the last old channel are NaN; the spline used to give 0 there,
    and a velocity of 0 at the end of the axis. Spectral axes are
    interpolated linearly (the spline bent them near the ends), and the
    channel widths, line centres (rounded, they were truncated), per-channel
    weights and `spectrum` arrays follow the new channels. Regridded cubes
    label their spectral axis with the exact, fractional, channel of the line
    (it was rounded to a whole channel).
25. `rebin_OneD` moves the line centre (`centerIndex`) with the rebinned
    channels (it kept the old channel, so line-centre stacks of rebinned
    spectra were misaligned unless `fit=True`), and compares line widths in
    channels (fitted widths in km/s mixed spectra with different channel
    widths).
26. The spectral axis of a 1D `Spectrum` must be evenly spaced (within 5% of
    a channel, so rounded values from text files are fine): every channel is
    located with one channel width, which is now the mean spacing (it was
    the first spacing, which rounding of the values could make off by a few
    percent). `weights='sigma2'` with `fit=True` left the line in the noise
    estimate of spectra with decreasing axes. A line centre outside its
    spectrum is now an error when stacking (it made very long stacks).

On the example data (no blank pixels) all stacks, with every weighting,
line and continuum, are bit-identical to the first 2.0 port; items 13-26
change results only for data that hit those cases (the examples do not
regrid spectra by redshift).

## Speed and parallelism

The 1D stacker, the bootstraps and the random-position estimator were
vectorised; this matters more than extra cores. Times on the example data,
2-core sandbox (expect better scaling with more cores):

| task | LineStacker 1.x | 2.0 serial | 2.0, 2 cores |
|---|---|---|---|
| 1D median bootstrap, 50 spectra, 10 000 draws | 263 s | 1.6 s | 0.9 s (threads) |
| 1D mean bootstrap, 100 000 draws | ~45 min (est.) | 0.8 s | |
| 1D bootstrap with Gaussian fits, 20 000 draws | | 18.0 s | 10.9 s (processes) |
| 1D subsample (maximizeSNR), 10 000 draws | ~190 s (est.) | 13.7 s | 8.2 s (processes) |
| cube stack, 50 CASA cubes (first / cached) | | 2.4 s / 0.09 s | |
| continuum stack of the same cubes (collapsed in memory) | | 1.6 s / 0.16 s | |
| cube mean bootstrap, 10 000 draws | broken | 0.2 s | |
| cube median bootstrap, 1 000 draws | broken | 8.5 s | 4.3 s (threads) |
| cube bootstrap with fits, 20 000 draws | broken | 11.7 s | 6.9 s (processes) |
| stack_estimator, 2 000 random stacks | broken | 0.2 s | |

`bootstraping_OneD`, `bootStraping_Cube`, `subsample_OneD` and
`stack_estimator` take `n_jobs` (-1 = all cores), `backend`
(`'process'` or `'thread'`; a sensible one is chosen per task) and `seed`.
Random numbers are drawn before the work is split, so results are identical
for any `n_jobs`. Worker processes start with `spawn` (safe with CASA and on
macOS), which costs ~0.7 s, so processes pay off for runs longer than a few
seconds. On macOS/Windows, scripts using processes need an
`if __name__ == '__main__':` guard (see the examples).

## Removed

- `LineStacker.modsub` and `LineStacker.pb`: the uv-domain code of Stacker
  (it needed the C++ `libstacker`). uv-domain stacking is
  `LineStacker.uv` (2.1); primary-beam weighting of image stacks takes
  primary-beam images or a function (see above).
- The vendored `interval` package (unused).
- CASA 5 support (`taskinit`): CASA 5 runs Python 2.
- Helpers working on the module state of the last stack:
  `image.calculate_sigma2_weights`, `calculate_sigma2_weights_spectral`,
  `calculate_amp_weights`, `myVeryOwnWeightedAverage` and the module
  variables (`image.data`, `skymap`, `unusedFrequencies`...), which kept
  every cube of the last stack in memory; use `weighting=` and the values
  `stack` returns. Also `analysisTools.randomizeSample` and
  `Coord.setWeightToX`/`setZeroWeightLeft`/`setZeroWeightRight` (unused).
- The documentation inside the package (`LineStacker/docs`): it is built
  from `docs/source` into `docs/` (the github.io site).

## Review fixes (October 2026)

A code review of the 2.0 port led to: the weighting fixes (items 13-15
above), the spectral-unit and axis checks, grid warnings and rotated output
WCS, primary-beam offsets on rotated grids (the cross terms of the pixel
matrix were ignored), the velocity conversion, `regridFromZ`'s cube labels,
NaN-aware fits, input validation, cache invalidation of modified CASA
images, closing every CASA tool, and per-call isolation of the inputs of
parallel jobs (two jobs in different threads could see each other's data).
Each fix has a test in `tests/test_weights_wcs.py` or
`tests/test_validation.py`.

A second review led to: non-linear spectral axes refused in every line
mode (the stack's header could not describe them), intensity-unit checks
and conversion, projection/frame/distortion checks, fractional pixels from
`regridFromZ`, infinite fluxes and one-channel spectra in 1D, `FitError`
for every fit failure (`DoublePeakFit` on short inputs, default widths from
valid channels only), the same "at least 2 valid pixels" rule for `'sigma'`
and `'sigma2'` as for `'sigma2F'`, integer subsample sizes, one image path
and limit checks in both `randomizeCoords`, `congrid(minusone=True)` with a
target size of 1, `None`/NaN cube names in catalogues, primary-beam cubes on
another spectral axis (an error), `LineStacker.cache_check` to trust the
cube cache, a lock around the thread variables set for worker processes,
and subprocess tests of CASA's log settings. Tests:
`tests/test_review2.py`.

A third round: blank pixels stay NaN (see "Changes in results", 23), SIP
and lookup-table distortions are compared term by term, the Gaussian
primary beam uses spherical sky offsets, `randomCoords(ncoords)` and
`randomizeCoords(maxTries)` must be integers, `noise_estimator` and
`stack_estimator` need `nRandom >= 1`, default fit widths use the median
channel spacing, and spectral axes of spectra must be finite and steadily
increasing or decreasing.

A fourth round finished the NaN handling: exact channels next to a blank
one are no longer blanked (`NaN * 0` is NaN) when stamps are cut,
rebinning and regridding share one NaN-aware resampling (blank channels
used by an interpolation blank its result, and nothing else),
`rebin_CubesSpectra(regionSize=...)` keeps channels blank in the whole
region as NaN for its width fits, `maximizeSNR` ignores blank channels in
its noise (and gives 0 with fewer than 2 valid ones), spectra with a
one-element velocity or frequency axis raise `ValueError`, and distortion
signatures with different terms compare as different. `regridFromZ` with
blanks now also works when it upsamples (it raised an error). Tests:
`tests/test_review4.py`.

A fifth round and a release review: splines that never use blank channels
(the blank mask approximated the spline's reach, which its prefilter
extends over the whole spectrum), evenly spaced 1D axes, `sigma2`
weights of fitted spectra with decreasing axes, NaN for 1D channels without
data, and the resampling fixes of items 24-25 (axes, centres, widths and
weights of rebinned and regridded spectra; exact channel widths and line
labels of regridded cubes). Also: line centres outside a spectrum,
`center='zero_vel'` without velocities, Gaussian-fit centring of frequency
spectra (it used the velocities), AIPS `FELO` axes treated as non-linear,
`rebin_CubesSpectra(regionSize=...)` regions outside the cube (and of the
size asked for: an odd size gave one pixel less, the region is now placed as
the stacker places stamps), a warning and weight 0 for blank primary-beam
responses (beyond the primary-beam cutoff, also for single-plane primary-beam
images, which gave a silent 0), dimension checks in the NaN-aware resampling, NaN-safe
`maximizeAmp` and 1D bootstrap amplitudes, and the missing docstrings of
`congrid`, `maximizeAmp` and `maximizeOutflow`. Tests: `tests/test_review5.py`.

The version stays 2.0.0 until the first 2.0 release; if a 2.0.0 was already
distributed, release this as 2.0.1.

## Example scripts

`LineStacker/Example/*/*.py` are ported (the originals were Python 2 and used
`taskinit`); `cube/continuum_example.py` is new. They read cubes with
`LineStacker.image.readCube` (no casatools calls of their own). The cube scripts accept the data folder as an argument instead of
the file dialog: `python basic_example.py data`. They use
`str.removesuffix('_coord.txt')`; the old `strip('_coords.txt')` removed
characters, not a suffix, and only worked by luck with the example names.

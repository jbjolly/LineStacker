# Analysis tools

{mod}`LineStacker.analysisTools` tests what a stack shows and prepares
cubes or spectra for stacking.

## Bootstrapping

Bootstrapping stacks random draws, with replacement, of the stacked sources.
The spread of the results shows how much the stack depends on which sources
it contains; a few bright sources dominating it give a skewed or
multi-peaked distribution.

```python
import LineStacker.analysisTools as A

boot = A.bootstraping_OneD(spectra, nRandom=10000, method='median', save='amp')
boot = A.bootStraping_Cube(coords, fEm=fEm, stampsize=16, N_chans=31,
                           nRandom=1000, save='ampAndWidth')
```

Both return a list of `nRandom + 1` entries, the first one for the stack of
all sources. `save` decides what is kept from each draw: `'all'` (the whole
stack), `'amp'` (its maximum), `'ampAndWidth'` (amplitude and width of a
Gaussian fit; fits that do not converge are stored as NaN), and for cubes
`'outflow'` (the six parameters of a two-Gaussian fit). For cubes the line
is the spectrum summed over the stamp.

## Subsampling

{func}`~LineStacker.analysisTools.subsample_OneD` stacks random subsets of
the spectra and grades every stack; each spectrum collects the grades of the
stacks it was part of. Spectra that make stacks better get higher totals.

```python
grades = A.subsample_OneD(spectra, nRandom=10000, maxTest=A.maximizeSNR)
best = sorted(grades, key=grades.get, reverse=True)[:10]
```

The grade function is {func}`~LineStacker.analysisTools.maximizeSNR` by
default (fitted amplitude over the standard deviation of the valid channels
outside the line; 0 when the fit fails or fewer than 2 such channels are
left),
{func}`~LineStacker.analysisTools.maximizeAmp` (maximum),
{func}`~LineStacker.analysisTools.maximizeOutflow` (amplitude × width² of the
fainter component of a two-Gaussian fit), or any function of a stacked
spectrum.
`minSize` and `maxSize` bound the sizes of the subsets.

## Random positions and noise

{func}`~LineStacker.analysisTools.stack_estimator` repeats a cube stack at
random positions, between `lowerLimit` and `upperLimit` radians from each
source (5 and 10 beams by default), keeping each source's redshift. The
distribution of these stacks is what noise alone gives:

```python
std, mean, stacks = A.stack_estimator(coords, fEm=fEm, nRandom=2000, N_chans=31)
```

{func}`~LineStacker.analysisTools.noise_estimator` measures the noise of
images or cubes from random pixels, optionally outside a mask, for the
whole cube or channel by channel. Blank (NaN) pixels are never drawn; pixels
that are exactly 0 are data. {func}`LineStacker.randomCoords` draws
random pixel positions anywhere in a set of images, and
{func}`LineStacker.randomizeCoords` moves positions randomly by one to a few
beams.

## Line widths and redshifts

Lines of the same rest-frame width are narrower at higher redshift, and
lines of different widths blur a stack. LineStacker can bring them to a
common width before stacking:

- {func}`~LineStacker.analysisTools.regridFromZ1D` (spectra) and
  `LineStacker.image.stack(..., regridFromZ=True)` (cubes, in memory)
  resample every source to the same rest-frame channel width;
  {func}`~LineStacker.analysisTools.regridFromZ` writes the regridded cubes
  (stamps around each source) to disk, to stack as returned with
  `spectralMethod='z'`.
- {func}`~LineStacker.analysisTools.rebin_OneD` and
  {func}`~LineStacker.analysisTools.rebin_CubesSpectra` rebin every source
  so that its line spans as many channels as the narrowest line, measured by
  Gaussian fits (or given as `widths`). This needs lines visible source by
  source.

Rebinned and regridded spectra keep their line centres, spectral axes and
per-channel weights in step with the new channels. Blank (NaN) channels are
never used to make a value: `rebin_OneD` and `rebin_CubesSpectra`
interpolate linearly, and a new channel interpolated from a blank channel is
blank; `regridFromZ1D` and `regridFromZ` interpolate with a cubic spline through each
stretch of valid channels separately, so only new channels between a valid
and a blank channel are blank. New channels beyond the last channel of a
source (when its channels are made narrower) are blank too.

## Speed, parallelism and reproducibility

The Monte-Carlo tools (`bootstraping_OneD`, `bootStraping_Cube`,
`subsample_OneD`, `stack_estimator`) take:

- `n_jobs`: the number of CPU cores to use (`-1` for all);
- `backend`: `'process'` or `'thread'`. A suitable one is chosen per task:
  threads for NumPy-bound work (mean and median stacks), processes for
  Python-bound work (Gaussian fits);
- `seed`: a seed for reproducible draws. With `seed=None` the draws come
  from NumPy's global random state, so `np.random.seed()` also makes runs
  reproducible.

All random draws are made before the work is split between workers, so the
results do not depend on `n_jobs` or `backend`. Worker processes are
started with `spawn`; scripts that use them on macOS or Windows need their
code under `if __name__ == '__main__':`. Two process-backed jobs started at
the same time from different threads run one after the other (each uses
all its workers).

## Fitting

{mod}`LineStacker.tools.fit` has the Gaussian fits the tools use:
{func}`~LineStacker.tools.fit.GaussFit`,
{func}`~LineStacker.tools.fit.DoubleGaussFit` and
{func}`~LineStacker.tools.fit.DoublePeakFit`. They ignore NaN channels and
raise {class}`~LineStacker.tools.fit.FitError` (a `RuntimeError`) when a fit
cannot be done or does not converge.

# Stacking spectra

{mod}`LineStacker.OneD_Stacker` stacks spectra that were already extracted.
Each spectrum is a {class}`LineStacker.OneD_Stacker.Spectrum`, which knows
where its line is; {func}`LineStacker.OneD_Stacker.Stack` aligns the spectra
on their lines and averages them channel by channel.

```python
import numpy as np
import LineStacker.OneD_Stacker as OneD

spectra = []
for i in range(100):
    data = np.loadtxt('data/spectra_%d' % i)        # columns: velocity, flux
    centre = np.loadtxt('data/central_velocity_%d' % i)
    spectra.append(OneD.Spectrum(spectrum=data, centralVelocity=centre,
                                 name=str(i)))

stack, nspectra = OneD.Stack(spectra)
```

`Stack` returns the stacked spectrum and the number of spectra in each of
its channels. Channels without data (no spectrum covers them, or only
blank fluxes or spectra of weight 0) are NaN.

## Spectra

A spectrum is given either as a two-column array, spectral values then flux
(shape `[N, 2]` or `[2, N]`): `Spectrum(spectrum=data)`; or as fluxes alone,
with or without their spectral values: `Spectrum(amp=flux,
velocities=v)`. `velOrFreq='vel'` (default) or `'freq'` says whether the
spectral values are velocities or frequencies. Velocities are in km/s.
The spectral values must increase or decrease evenly (within 5% of a
channel, so values rounded in a text file are fine): channels are located
with one channel width, the mean spacing. Resample spectra on other grids
(log-spaced, with gaps) to an even grid first.

The line centre comes from, in this order:

- `z=` and `fEmLine=`, the redshift and the rest frequency of the line. For
  frequency spectra the line is at `fEmLine/(1+z)`. For velocity spectra the
  velocities are relative to the rest frequency (as a CASA spectrum made with
  `restfreq=fEmLine`), and the line is at v = c z/(1+z) in the radio
  convention (`velConvention='radio'`, default, as CASA uses), c z with
  `'optical'`, or the relativistic velocity with `'relativistic'`;
- `centralVelocity=` or `centralFrequency=`, the observed line centre;
- `centerIndex=`, its channel;
- `fit=True`: a Gaussian fit of the spectrum, which overrides the above.

`weights` is the weight of the spectrum in mean stacks: a number (default 1),
an array of one weight per channel, `'sigma2'` (1/σ² of the spectrum, or of
the channels outside the line if `fit=True`), or `'1/A'` (one over the line
amplitude: the maximum, or the fitted amplitude). Channels whose flux is not
finite (NaN, inf) are left out. A spectrum needs at least two channels, and
as many spectral values as fluxes. `name` labels the spectrum, for subsampling for instance.

## Stack

```python
stack, nspectra = OneD.Stack(spectra, method='median', center='fit',
                             chansStack=41)
```

- `method`: `'mean'` (weighted) or `'median'` (unweighted).
- `center`: `'lineCenterIndex'` (default) uses the centre each spectrum was
  given; `'fit'` fits a Gaussian to every spectrum; `'zero_vel'` uses the
  channel at velocity 0; `'center'` the middle channel; an integer, that
  channel for all spectra.
- `chansStack`: `'full'` (default) covers every channel of every spectrum,
  aligned on the lines; an integer stacks that many channels around the line.
  The line centre of every spectrum must be one of its channels (an error
  says which spectrum it is not).
- `regridFromZ=True` resamples the spectra so that all have the same
  rest-frame channel width (see {func}`LineStacker.analysisTools.regridFromZ1D`);
  every spectrum then needs its `z`. Channels and line centres are scaled by
  (1+z)/(1+z_ref); with `'scaleToMin'` the last channels of the regridded
  spectra can fall beyond their data and are then NaN.

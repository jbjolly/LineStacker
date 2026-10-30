# LineStacker

LineStacker stacks the spectral lines of many faint sources to measure their
average emission: the line of each source is centred on its own redshift,
then all sources are averaged channel by channel. It stacks

- **data cubes**, FITS files or CASA images, around each source
  (`LineStacker.image`), for lines or continuum;
- **1D spectra** already extracted (`LineStacker.OneD_Stacker`);

and comes with tools to test what a stack shows: bootstrapping, subsampling
to find the sources that dominate it, stacks at random positions, noise
estimates, and rebinning or regridding of lines of different widths and
redshifts.

```python
import LineStacker
import LineStacker.image

# one row per source: ra, dec, z and the cube it lies in
coords = LineStacker.readCatalog('sources.csv')
fEm = 1897420620253.1646          # rest frequency of the line (Hz)
result = LineStacker.image.stack(coords, fEm=fEm, stampsize=16, N_chans=31)
# LineStacker: mean line stack of 30 sources (16x16 pixels, 31 channels)
#              written to .../stackResult.fits
spectrum = result[1]          # stacked spectrum of the central pixel
```

LineStacker 2.0 runs on Python 3 and CASA 6 (CASA is only needed for CASA
images: FITS cubes need no CASA). See [what changed in 2.0](changes.md) if
you used LineStacker 1.x.

```{toctree}
:maxdepth: 1
:caption: User guide

installation
principle
positions
cubes
spectra
analysis
```

```{include} _generated/uv_toc.md
```

```{toctree}
:maxdepth: 1
:caption: Examples and reference

examples
_generated/examples_notebook
api
changes
```

## Citing LineStacker

If LineStacker is useful for your work, please cite
Jolly, Knudsen & Stanley (2020),
[MNRAS 499, 3992](https://doi.org/10.1093/mnras/staa2908).
LineStacker started from Stacker
([Lindroos et al. 2015, MNRAS 446, 3502](https://doi.org/10.1093/mnras/stu2344)).

The code is on [GitHub](https://github.com/jbjolly/LineStacker), under the
GNU General Public License.

## Indices

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`

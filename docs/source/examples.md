# Examples

The repository comes with example data and scripts in `LineStacker/Example`,
and a notebook that runs all of them:
[Examples notebook](_generated/examples_notebook.ipynb).

## Notebook

`LineStacker/Example/LineStacker_examples.ipynb` runs every example on the
example data in about half a minute and writes its files to
`Example/notebook_outputs/`. Start it from the `Example` folder:

```bash
cd LineStacker/Example
jupyter notebook LineStacker_examples.ipynb
```

It covers 1D stacking, bootstrapping and subsampling of spectra; cube
stacking from coordinate files, from a CSV catalogue and from Python arrays;
line fits, rebinning and bootstrapping of cubes; and continuum stacking with
stacks at random positions.

## Scripts

Spectra (`Example/oneD`, run from that folder):

- `basic-example.py` stacks 100 spectra (`data/spectra_<i>`, velocity and
  flux columns) on their central velocities (`data/central_velocity_<i>`).
- `extended-example.py` takes 50 spectra, 10 of which have lines ten times
  brighter: a median stack, a bootstrap that shows the stack depends on
  which sources are drawn, and a subsampling that finds the bright spectra.

Cubes (`Example/cube`): 30 CASA cubes of simulated line emitters in `data`
and 50 in `dataExtended`, each with its coordinate file. Give the scripts
the data folder:

```bash
cd LineStacker/Example/cube
python basic_example.py data              # line stack of the 30 cubes
python extended_example.py dataExtended   # spectra, fits, rebinning, bootstrap
python continuum_example.py data          # continuum stack and random positions
```

Without a folder, the scripts open a file dialog to choose the coordinate
files.

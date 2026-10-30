<p align="center">
  <img src="https://jbjolly.github.io/LineStacker/_static/LineStacker-logo.svg" alt="LineStacker" width="600">
</p>

Python module to stack spectra in the image domain.

Install directly using pip: `pip install LineStacker`, or from a clone:
`pip install -e .` (add `[casa]` to also install CASA's Python packages).

Requirements:
- Python >= 3.9, numpy, scipy, astropy (astropy >= 7.2 with NumPy >= 2.4,
  e.g. with CASA 6.7's bundled NumPy)
- Cubes can be FITS files or CASA images. CASA images need the modular CASA 6
  casatools (`pip install casatools`, tested with CASA 6.7.6); FITS cubes
  need no CASA at all. CASA's log goes to one file in the temporary directory
  (`LineStacker.casa_logfile` changes it), not a new `casa-*.log` per run.
- matplotlib for the plotting options.

Stacking positions come from a catalogue table, one row per source with the
cube of each source in an `image` column (`LineStacker.readCatalog`: CSV,
ECSV, FITS, VOTable, astropy Table or pandas DataFrame), from Python arrays
(`LineStacker.makeCoords`), or from one coordinate file per cube
(`LineStacker.readCoords`).

`LineStacker.image.stack` stacks lines in cubes, or continuum with
`continuum=True` (images without a spectral axis are always stacked as
continuum; cubes are then summed over channels in memory).
`LineStacker.OneD_Stacker.Stack` stacks 1D spectra given as
`LineStacker.OneD_Stacker.Spectrum` objects.

LineStacker prints what it does, including the full path of every file it
writes (e.g. the stacked cube). `LineStacker.verbose = False` silences it;
every function that prints also accepts `verbose=True/False`.

The Monte-Carlo tools (bootstrapping, subsampling, random-position stacks)
accept `n_jobs=-1` to use all CPU cores. When using worker processes from a
script on macOS or Windows, put the script's code under
`if __name__ == '__main__':`.

`LineStacker/Example/LineStacker_examples.ipynb` runs all the examples on the
example data (run it from the `Example` folder).

Tests: `python -m pytest tests` (synthetic data; CASA tests are skipped when
casatools is not installed).

Documentation: https://jbjolly.github.io/LineStacker/ (also `docs/index.html`
and `docs/LineStacker.pdf`), built from `docs/source` with `make -C docs html`
(see `docs/README.md`). Changes from LineStacker 1.x: `PYTHON3_MIGRATION.md`.

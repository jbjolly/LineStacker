# Installation

LineStacker needs Python 3.9 or newer, NumPy, SciPy and astropy.

```bash
pip install LineStacker                 # FITS cubes and 1D spectra
pip install "LineStacker[casa]"         # also CASA images (casatools, casatasks)
pip install "LineStacker[plot]"         # matplotlib, for the plotting options
```

To install the development version, clone the
[repository](https://github.com/jbjolly/LineStacker) and install it in
editable mode:

```bash
git clone https://github.com/jbjolly/LineStacker.git
cd LineStacker
pip install -e ".[casa,plot]"
```

## CASA images

FITS cubes are read with astropy and need no CASA. To read and write CASA
images (`.image` directories), LineStacker uses CASA 6's `casatools`, the
modular CASA from PyPI (tested with CASA 6.7.6):

```bash
pip install casatools
```

`casatools` needs CASA's runtime data (the measures tables). By default it
downloads them into `~/.casa/data` the first time it is imported; see
[External data](https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html)
in the CASA documentation to set another place, for example without internet
access:

```python
# ~/.casa/config.py
measurespath = "/path/to/casa/data"
measures_auto_update = False
data_auto_update = False
```

### From the CASA prompt

LineStacker also runs inside the monolithic CASA (the CASA tarball with its
own Python). Install its dependencies with CASA's own `pip`, so they match
the NumPy that CASA ships:

```bash
PYTHONNOUSERSITE=1 /path/to/casa-6.7.6/bin/pip3 install "astropy>=7.2"
PYTHONNOUSERSITE=1 /path/to/casa-6.7.6/bin/pip3 install /path/to/LineStacker
```

Mixing CASA with packages from another environment (a virtual environment on
`PYTHONPATH`, say) breaks when the versions disagree: CASA 6.7 ships NumPy
2.4, which astropy older than 7.2 cannot use (the error is
`module 'numpy' has no attribute 'in1d'`). `import LineStacker` warns when it
finds that combination.

### CASA's log file

CASA is only started when LineStacker reads or writes a CASA image. When
LineStacker starts it, CASA's log goes to a single file in the temporary
directory, `<tmp>/linestacker-casa-<user>.log`, and CASA's messages go there
rather than to the terminal: no new `casa-<date>.log` appears in the working
directory at every run. Set `LineStacker.casa_logfile` before the first CASA
image is read to change this:

```python
import LineStacker
LineStacker.casa_logfile = 'stacking.log'   # this file
LineStacker.casa_logfile = None             # no log file
LineStacker.casa_logfile = 'casa'           # CASA's default, casa-<date>.log here
```

A `logfile` set in `~/.casa/config.py` is used as it is. Inside a CASA
session, or after importing `casatools` or `casatasks` yourself, CASA's own
settings apply.

## Checking the installation

The tests run on synthetic cubes and spectra (the CASA tests are skipped
when casatools is not installed):

```bash
pip install pytest
python -m pytest tests
```

The examples in `LineStacker/Example` run on the example data of the
repository; see [Examples](examples.md).

## Messages

LineStacker reports what it does: positions read, the full path of every
file it writes, progress bars. `LineStacker.verbose = False` silences it;
every function that prints also takes `verbose=True` or `verbose=False` for
one call. Possible problems (cubes collapsed for continuum stacking,
sources with no valid noise, cubes on different pixel grids...) are Python
warnings and show either way.

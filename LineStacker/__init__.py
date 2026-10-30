"""
Stacking positions (:class:`Coord`, :func:`readCatalog`, :func:`makeCoords`,
:func:`readCoords`), their conversion to pixels, random positions, and the
package settings (:data:`verbose`, :data:`casa_logfile`). The stacking itself
is in :mod:`LineStacker.image` (cubes) and :mod:`LineStacker.OneD_Stacker`
(spectra); :mod:`LineStacker.analysisTools` tests stacks.
"""
# -*- coding: utf-8; -*-
# LineStacker, Python module for stacking data cubes.
# Copyright (C) 2019  Jean-Baptiste Jolly
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301, USA.

#LineStacker is based on Stacker by Lukas Lindroos,
#Most of the main module (this one) is hence directly taken from Stacker
"""
    Library to stack interferometric data-cubes, or lines from 1D data.
"""

import copy as _copy
import os
import re
import numpy as np


__author__ = 'Jean-Baptiste Jolly'
__copyright__ = 'Copyright 2019'
__license__ = 'GPL'
__version__ = '2.0.0'
__maintainer__ = 'Jean-Baptiste Jolly'
__email__ = 'jbjolly@mpe.mpg.de'


def _major_minor(version):
    nums = re.findall(r'\d+', str(version))[:2]
    return tuple(int(n) for n in nums) + (0,) * (2 - len(nums))


def _astropy_numpy_problem(numpy_version=None, astropy_version=None, astropy_path=None):
    """Message if the astropy that would be imported cannot work with this
    NumPy, else None. astropy < 7.2 uses np.in1d, which NumPy 2.4 removed;
    this happens e.g. when CASA's bundled NumPy is mixed with an older astropy
    from another environment. Checked without importing astropy."""
    try:
        if astropy_version is None:
            from importlib.metadata import version
            from importlib.util import find_spec
            astropy_version = version('astropy')
            spec = find_spec('astropy')
            astropy_path = os.path.dirname(spec.origin) if spec and spec.origin else '?'
        numpy_version = numpy_version or np.__version__
    except Exception:            # astropy not installed: reported when needed
        return None
    if _major_minor(numpy_version) < (2, 4) or _major_minor(astropy_version) >= (7, 2):
        return None
    return ('astropy {0} (from {1}) cannot be used with NumPy {2}: astropy < 7.2 needs '
            'np.in1d, which NumPy 2.4 removed, so reading coordinates and cubes will fail with '
            "\"module 'numpy' has no attribute 'in1d'\". Install astropy >= 7.2 into the Python "
            'that runs LineStacker; for the CASA prompt (monolithic CASA): '
            "PYTHONNOUSERSITE=1 <casa>/bin/pip3 install 'astropy>=7.2'.").format(
                astropy_version, astropy_path, numpy_version)


def _warn_astropy_numpy():
    msg = _astropy_numpy_problem()
    if msg:
        import warnings
        warnings.warn(msg, stacklevel=3)


_warn_astropy_numpy()


#: When True (default), LineStacker reports what it does, the files it writes
#: (with their full path), what it reads, and progress bars. Set
#: ``LineStacker.verbose = False`` to silence it; every function that prints
#: also takes ``verbose=True/False`` to decide for one call. Warnings about
#: possible problems are shown either way.
verbose = True

#: Where CASA writes its log when LineStacker starts CASA itself (only to
#: read or write CASA images; FITS cubes never start CASA):
#:
#: - ``'auto'`` (default): one file, ``<temporary dir>/linestacker-casa-<user>.log``,
#:   instead of a new ``casa-<date>.log`` in the working directory at every
#:   run; a ``logfile`` set in ``~/.casa/config.py`` is used instead.
#: - a path: CASA's log goes there.
#: - ``None``: no log file.
#: - ``'casa'``: CASA's own default (``casa-<date>.log`` in the working directory).
#:
#: Set it before the first CASA image is read. Inside a CASA session, or once
#: casatools/casatasks were imported, CASA's own settings apply.
casa_logfile = 'auto'

#: Cubes are read once and kept in memory. When True (default) LineStacker
#: checks, every time it uses a cached cube, that its file has not changed
#: (about 0.1 ms per CASA image, a few microseconds per FITS file) and reads
#: it again if it has. Set ``LineStacker.cache_check = False`` to trust the
#: cache until ``LineStacker.image.clear_cache()``.
cache_check = True


def _verbose(flag=None):
    """A function's verbose argument, defaulting to LineStacker.verbose."""
    return bool(verbose) if flag is None else bool(flag)


def _info(message, flag=None):
    if _verbose(flag):
        print('LineStacker: ' + message, flush=True)


# CASA writes sexagesimal coordinates with dots as deg/min/sec separators
# (e.g. '-027.45.34.566' for declination, sometimes also '03.30.34.79' for RA).
# Astropy's angle parser doesn't understand that form -- it expects colons
# ('-027:45:34.566'), letter separators ('-027d45m34.566s'), or spaces. This
# helper normalises CASA's dot form to colons; anything else passes through.
_CASA_SEXAGESIMAL_RE = re.compile(r'^\s*([+-]?\d+)\.(\d+)\.(\d+(?:\.\d+)?)\s*$')


def _normalize_sexagesimal(s):
    """Convert CASA dot-separated sexagesimal ('DD.MM.SS.fff') to colon form
    ('DD:MM:SS.fff'). Returns the input unchanged if it doesn't match the
    CASA pattern (so plain degrees, 'HH:MM:SS', '00h42m30s', and space-
    separated forms all pass through untouched)."""
    if not isinstance(s, str):
        return s
    m = _CASA_SEXAGESIMAL_RE.match(s)
    if m:
        return '{0}:{1}:{2}'.format(m.group(1), m.group(2), m.group(3))
    return s


def round_half_up(x):
    """Round like Python 2's round() did (halves away from zero). Python 3's
    round() rounds halves to even, which would move some stack centres by one
    channel compared with LineStacker < 2.0."""
    return int(np.floor(np.abs(x) + 0.5) * np.sign(x))


class Coord:
    """
        Describes a stacking position on an image.

        Class used internally to represent coordinates. May describe a
        physical coordinate or a pixel coordinate.

        Init: Creates a coordinate. A pixel coordinate should always specify
        to which image it belongs. Physical coordinates are astropy angles
        (as produced by **LineStacker.readCoords**, **readCatalog** and
        **makeCoords**) or plain numbers in radians, in the equatorial frame
        given by **frame** (default ICRS).

        Parameters
        ---------
        ra
            x coordinate of stacking target (RA, or x pixel for pixel coordinates).
            Also available as ``coord.x``.
        dec
            y coordinate of stacking target (Dec, or y pixel). Also ``coord.y``.
        spec
            spectral coordinate of stacking target (redshift, observed spectral
            value or channel depending on spectralMethod). Also ``coord.z``.
        weight
            Weight of the source for mean stacking.
            Default is 1.
        image
            Cube associated to the source, either its path or its index in the
            list of images given to the stacking function.
            Automatically set with **LineStacker.readCoords**.
        type
            Either 'physical' or 'pixel'
        imagenumber
            Index of the cube in the image list (set by readCoords to the index
            of the coordinate file the source was read from).
        frame
            Equatorial frame of ra, dec: 'icrs' (default), 'fk5' (J2000), 'fk4'...

    """

    def __init__(self, ra, dec, spec='continuum', weight=1, image=0, type='physical',
                 imagenumber=-1, frame='icrs', z=None, index=None):
        self.ra = ra
        self.dec = dec
        self.spec = spec if z is None else z
        self.weight = weight
        self.image = image
        self.type = type
        self.imagenumber = imagenumber
        self.frame = frame
        self.index = index

    # Stacker-era attribute names, still used by older scripts.
    @property
    def x(self):
        return self.ra

    @x.setter
    def x(self, value):
        self.ra = value

    @property
    def y(self):
        return self.dec

    @y.setter
    def y(self, value):
        self.dec = value

    @property
    def z(self):
        return self.spec

    @z.setter
    def z(self, value):
        self.spec = value

    def copy(self):
        """Independent copy (weights arrays are copied too)."""
        new = _copy.copy(self)
        if isinstance(self.weight, np.ndarray):
            new.weight = self.weight.copy()
        elif isinstance(self.weight, list):
            new.weight = list(self.weight)
        return new

    def __str__(self):
        return '{0}, {1}'.format(self.ra, self.dec)

    def __repr__(self):
        return 'Coord({0}, {1}, spec={2}, weight={3}, image={4!r}, type={5!r})'.format(
            self.ra, self.dec, self.spec, self.weight, self.image, self.type)



def readCoordsNamesGUI():
    """
        Opens a file dialog (Tk) to select coordinate files; returns their paths.
    """
    import tkinter
    from tkinter import filedialog
    root = tkinter.Tk()
    root.withdraw()
    filez = filedialog.askopenfilenames(initialdir=".", title="Select coord files",
                                        filetypes=(("txt files", "*.txt"), ("all files", "*.*")))
    root.destroy()
    return list(filez)


class _CoordFormatError(RuntimeError):
    """A coordinate file could not be parsed with the current delimiter."""


def readCoords(coordfiles, unit='deg', lineON=True, delimiter='default', frame='icrs', images='default', type='physical',
               verbose=None):
    """
        Reads coordinate files (one file per cube, the format of LineStacker
        1.x) and returns a list of LineStacker.Coord.

        Parameters
        ---------
        coordfiles:
            List of path to coordinate files. Files in csv format. x and y should be in J2000. If using to stack line, redshift should be in the third column. A weight may be added in the last column to weight positions for mean stacking. If set to **None** stacking position is defined as the center of each cube.
            Coordinates can either be in degrees (e.g. "10.625, 41.2"), or in hour angle + deg, as hms dms format (e.g. "'00h42m30s', '+41d12m00s'"), or as hour angle + deg without hms nor dms (e.g. "00 42 30, +41 12 00" or CASA's "00:42:30, +41.12.00"). Empty lines and lines starting with # are ignored.
        unit:
            Unit of input coordinates. Allows two values, **'deg'** or **'rad'**.
        lineON
            Either **True** or **False**, should be set to **True** if doing line stacking. Default is **True**.
        delimiter:
            Delimiter used to separate coordinates, 'default' tries comma, space and tab (in that order).
        images
            List of images to which each coordfile is assigned, if set to 'default', image names are same as coordfiles, with '.fits' instead of '_coord.txt'.
            The stacking functions also accept coordinates whose image is not in their imagenames list: the n-th coordinate file is then matched to the n-th image.
        type
            Coordinate type, can be either 'physical' (default), or 'pixel'
        verbose
            report how many positions were read (default: LineStacker.verbose)

    """
    if coordfiles is None:
        _info('coordfile is None, stacking on the center', verbose)
        return None
    if isinstance(coordfiles, str):
        coordfiles = [coordfiles]
    if delimiter=='default':
        delimiters=[',',' ','\t']
        delimiterNames=['coma', 'space', 'tab']
    else:
        delimiters=[delimiter]
        delimiterNames=[repr(delimiter)]
    if images=='default':
        images=[coordname.replace('_coord.txt', '.fits') for coordname in coordfiles]
    elif not isinstance(images, (list, tuple)):
        images=[images]
    if len(images) < len(coordfiles):
        raise ValueError('{0} coordinate files but only {1} images given'.format(len(coordfiles), len(images)))

    errors = []
    for (j, delim) in enumerate(delimiters):
        try:
            coords = _readCoords(coordfiles, unit=unit, lineON=lineON, delimiter=delim,
                                 frame=frame, images=images, type=type)
        except _CoordFormatError as err:
            errors.append('{0}: {1}'.format(delimiterNames[j], err))
            continue
        _info('read {0} position{1} from {2} coordinate file{3}'.format(
            len(coords), 's' if len(coords) != 1 else '', len(coordfiles), 's' if len(coordfiles) != 1 else ''),
            verbose)
        return coords
    raise Exception("Could not read the coordinate files. If set to default delimiters can be comma, space, and tab. "
                    "Otherwise feed own delimiter to this function.\n" + '\n'.join(errors))


def _parse_angle_pair(ra, dec, unit, frame):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    uni = {'deg': u.deg, 'rad': u.rad}.get(unit)
    if uni is None:
        raise ValueError("unit must be 'deg' or 'rad'")
    ra_str = _normalize_sexagesimal(ra)
    dec_str = _normalize_sexagesimal(dec)
    try:
        return SkyCoord(ra=float(ra_str)*uni, dec=float(dec_str)*uni, frame=frame)
    except ValueError:
        pass
    try:
        return SkyCoord(ra=ra_str, dec=dec_str, frame=frame)
    except (u.UnitsError, ValueError):
        pass
    try:
        return SkyCoord(ra=ra_str, dec=dec_str, frame=frame, unit=(u.hourangle, u.deg))
    except ValueError as err:
        raise _CoordFormatError('Angles must be in a format recognized by astropy.coordinates, '
                                'with ra in hourangle and dec in deg ({0})'.format(err))


def _readCoords(coordfiles, unit='deg', lineON=True, delimiter=',', frame='icrs', images='', type='physical'):
    import csv

    coords = []
    for (i, coordi) in enumerate(coordfiles):
        with open(coordi, 'r') as fh:
            rows = list(csv.reader(fh, delimiter=delimiter))
        for lineno, row in enumerate(rows, 1):
            row = [field.strip() for field in row]
            row = [field for field in row if field != '']
            if len(row) == 0 or row[0].startswith('#'):
                continue
            if len(row) < 2:
                raise _CoordFormatError('{0}, line {1}: expected at least 2 columns, found {2}'.format(coordi, lineno, len(row)))
            lineMode = lineON
            if not lineMode and len(row) == 4:
                lineMode = True
                import warnings
                warnings.warn('{0}, line {1}: lineON=False but 4 columns found; reading it as x, y, z, w '
                              '(line mode)'.format(coordi, lineno), stacklevel=3)
            try:
                if lineMode:
                    if len(row) == 4:
                        weight = float(row[3])
                    elif len(row) == 3:
                        weight = 1.
                    else:
                        raise _CoordFormatError('{0}, line {1}: structure should be x, y, z (, w)'.format(coordi, lineno))
                    spec = float(row[2])
                else:
                    if len(row) == 3:
                        weight = float(row[2])
                    elif len(row) == 2:
                        weight = 1.
                    else:
                        raise _CoordFormatError('{0}, line {1}: too many columns, structure should be x, y, w '
                                                '(or x, y, spec, w and LineON set to True)'.format(coordi, lineno))
                    spec = 'continuum'
            except ValueError as err:
                raise _CoordFormatError('{0}, line {1}: {2}'.format(coordi, lineno, err))

            if type == 'physical':
                coord = _parse_angle_pair(row[0], row[1], unit, frame)
                coords.append(Coord(coord.ra, coord.dec, spec=spec, weight=weight, image=images[i],
                                    imagenumber=i, type='physical', frame=frame))
            elif type == 'pixel':
                try:
                    x, y = float(row[0]), float(row[1])
                except ValueError as err:
                    raise _CoordFormatError('{0}, line {1}: {2}'.format(coordi, lineno, err))
                coords.append(Coord(x, y, spec=spec, weight=weight, image=images[i],
                                    imagenumber=i, type='pixel'))
            else:
                raise ValueError("type must be 'physical' or 'pixel'")
    return coords


_HAS_UNIT_LETTER = re.compile(r'[hdms°\'"]', re.IGNORECASE)


_SEXAGESIMAL = re.compile(r'^([+-]?)(\d+)(?::|\s+)(\d+)(?::|\s+)(\d+(?:\.\d*)?)$')


def _per_source(values, n, name):
    """One value per source: scalars (and strings) are repeated n times."""
    if values is None:
        return [None] * n
    if isinstance(values, (str, os.PathLike)) or np.ndim(values) == 0:
        return [values] * n
    seq = list(values)
    if len(seq) != n:
        raise ValueError('{0} has {1} values but there are {2} sources'.format(name, len(seq), n))
    return seq


def _equatorial(frame):
    """Name of an equatorial astropy frame (ra/dec), or ValueError."""
    from astropy.coordinates import frame_transform_graph
    name = frame if isinstance(frame, str) else getattr(frame, 'name', None)
    cls = frame_transform_graph.lookup_name(name) if name else None
    if cls is None:
        raise ValueError('unknown coordinate frame {0!r}'.format(frame))
    names = cls().representation_component_names
    if 'ra' not in names or 'dec' not in names:
        raise ValueError("frame {0!r} has no RA/Dec: LineStacker coordinates are equatorial ('icrs', 'fk5' "
                         "for J2000, 'fk4'...). Convert the positions, e.g. SkyCoord(...).icrs".format(name))
    return name


def _angles_to_sky(ra, dec, unit, frame):
    """SkyCoord (1D) from numbers in unit, astropy angles/quantities, or
    sexagesimal strings (RA in hours; CASA's dotted declinations accepted)."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    uni = {'deg': u.deg, 'rad': u.rad}.get(unit)
    if uni is None:
        raise ValueError("unit must be 'deg' or 'rad'")
    def as_degrees(values, sexagesimal_unit, axis):
        if hasattr(values, 'unit'):
            return u.Quantity(values).reshape(-1).to_value(u.deg)
        arr = np.atleast_1d(np.asarray(values, dtype=object))
        try:
            return (np.asarray(arr, dtype=float) * uni).to_value(u.deg)
        except (TypeError, ValueError):
            pass
        strings = [_normalize_sexagesimal(str(v).strip()) for v in arr]
        parts = [_SEXAGESIMAL.match(s) for s in strings]
        if all(parts):                    # common case: 'DD:MM:SS.s' or 'DD MM SS.s'
            sign = np.array([-1. if p.group(1) == '-' else 1. for p in parts])
            dms = np.array([[float(p.group(2)), float(p.group(3)), float(p.group(4))] for p in parts])
            if np.all(dms[:, 1] < 60) and np.all(dms[:, 2] < 60):
                value = sign * (dms[:, 0] + dms[:, 1] / 60. + dms[:, 2] / 3600.)
                return value * (15. if sexagesimal_unit == u.hourangle else 1.)
        out = np.empty(len(arr))
        for k, v in enumerate(arr):       # strings: each value read on its own
            s = _normalize_sexagesimal(str(v).strip())
            try:
                out[k] = (float(s) * uni).to_value(u.deg)
                continue
            except ValueError:
                pass
            try:
                angle = Angle(s) if _HAS_UNIT_LETTER.search(s) else Angle(s, unit=sexagesimal_unit)
            except Exception as err:
                raise ValueError('source {0}: cannot read {1} {2!r} ({3})'.format(k, axis, str(v), err))
            out[k] = angle.deg
        return out

    from astropy.coordinates import Angle
    ra_deg = as_degrees(ra, u.hourangle, 'RA')
    dec_deg = as_degrees(dec, u.deg, 'Dec')
    if len(ra_deg) != len(dec_deg):
        raise ValueError('ra has {0} values but dec has {1}'.format(len(ra_deg), len(dec_deg)))
    return SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame=frame)


def makeCoords(ra, dec=None, z=None, weight=None, image=None, unit='deg', frame='icrs',
               type='physical', verbose=None):
    """
        Builds stacking coordinates directly from Python values (lists, numpy
        arrays, astropy quantities or a SkyCoord), e.g. in a notebook.

        Parameters
        ---------
        ra, dec
            Positions, one value per source. Numbers are in **unit** (or pixels
            if type='pixel'); astropy Angles/Quantities carry their own unit;
            strings are sexagesimal with RA in hours ('03:31:04.15',
            '-27:40:59.87', or CASA's '-027.40.59.866'). A SkyCoord can be
            given as ra, with dec=None.
        z
            Spectral value of each source: the redshift for spectralMethod='z'
            (the default of LineStacker.image.stack), the observed
            frequency/wavelength for 'obsSpec', the channel for 'channel'.
            None for continuum stacking.
        weight
            Weight of each source for mean stacks (default 1).
        image
            Cube of each source (path to a FITS file or CASA image), or one
            path for all sources. None when all sources lie in the single cube
            given to the stacking function.
        unit
            'deg' (default) or 'rad', for plain numbers.
        frame
            Equatorial frame of the positions: 'icrs' (default), 'fk5' (J2000), 'fk4'...
            A SkyCoord in any frame (e.g. Galactic) is converted to it.
        type
            'physical' (sky positions, default) or 'pixel'.

        Returns
        -------
        list of LineStacker.Coord
            to give to LineStacker.image.stack (imagenames can then be omitted
            when every source names its cube).
    """
    from astropy.coordinates import SkyCoord
    if type not in ('physical', 'pixel'):
        raise ValueError("type must be 'physical' or 'pixel'")
    frame = _equatorial(frame)
    if dec is None:
        if not isinstance(ra, SkyCoord):
            raise ValueError('dec is missing (or give a SkyCoord as ra)')
        if type == 'pixel':
            raise ValueError("a SkyCoord holds sky positions: use type='physical'")
        sky = ra.transform_to(frame)
        sky = sky.reshape((1,)) if sky.isscalar else sky.reshape(-1)
    elif type == 'physical':
        sky = _angles_to_sky(ra, dec, unit, frame)
    else:
        xs = np.atleast_1d(np.asarray(ra, dtype=float))
        ys = np.atleast_1d(np.asarray(dec, dtype=float))
        if len(xs) != len(ys):
            raise ValueError('ra has {0} values but dec has {1}'.format(len(xs), len(ys)))
    if type == 'physical':
        xs, ys = sky.ra, sky.dec          # extract once: indexing a SkyCoord is slow
        n = len(sky)
    else:
        n = len(xs)
    specs = _per_source(z, n, 'z')
    weights = _per_source(1. if weight is None else weight, n, 'weight')
    images = _per_source(image, n, 'image')
    for k in range(n):
        if specs[k] is not None and not np.isfinite(float(specs[k])):
            raise ValueError('source {0}: z is {1}; every source needs a finite spectral value '
                             '(or z=None for continuum)'.format(k, specs[k]))
        if not np.isfinite(float(weights[k])) or float(weights[k]) < 0:
            raise ValueError('source {0}: weight is {1}; weights must be finite and >= 0'.format(k, weights[k]))
        if images[k] is not None and not os.fspath(images[k]).strip():
            raise ValueError('source {0}: its image (cube) name is empty'.format(k))
    coords = []
    for k in range(n):
        spec = 'continuum' if specs[k] is None else float(specs[k])
        img = os.fspath(images[k]) if images[k] is not None else None
        if type == 'physical':
            px, py = xs[k], ys[k]
        else:
            px, py = float(xs[k]), float(ys[k])
        # imagenumber -1: the cube is matched by name only, never by order
        coords.append(Coord(px, py, spec=spec, weight=float(weights[k]), image=img, imagenumber=-1,
                            type=type, frame=frame))
    _info('made {0} coordinate{1}'.format(n, 's' if n != 1 else ''), verbose)
    return coords


_CATALOG_ALIASES = {
    'ra': ('ra', 'ra_deg', 'raj2000', 'ra_j2000', 'alpha'),
    'dec': ('dec', 'dec_deg', 'dej2000', 'decj2000', 'dec_j2000', 'delta'),
    'z': ('z', 'redshift', 'zspec', 'z_spec', 'spec'),
    'weight': ('weight', 'weights'),
    'image': ('image', 'imagename', 'cube', 'file', 'filename', 'path'),
}


def _read_table(catalog, format):
    """astropy Table and the folder of the catalogue file (None if not a file)."""
    from astropy.table import Table
    if isinstance(catalog, Table):
        return catalog, None
    if hasattr(catalog, 'columns') and hasattr(catalog, 'to_dict') and not isinstance(catalog, (str, os.PathLike)):
        return Table.from_pandas(catalog), None          # pandas DataFrame
    if isinstance(catalog, dict):
        return Table(catalog), None
    path = os.fspath(catalog)
    if not os.path.isfile(path):
        raise FileNotFoundError('catalogue not found: {0}'.format(path))
    import warnings
    lower = path.lower()
    binary = lower.endswith(('.fits', '.fit', '.fts', '.fits.gz', '.vot', '.votable', '.xml', '.hdf5', '.h5',
                             '.parquet'))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')              # unit-parsing chatter of some writers
        if format is not None:
            table = Table.read(path, format=format)
        elif binary:
            table = Table.read(path)                 # its own error is the useful one
        else:
            try:
                table = Table.read(path)
            except Exception as first:
                try:
                    table = Table.read(path, format='ascii')
                except Exception as second:
                    raise ValueError('cannot read the catalogue {0}: {1} (as a text table: {2}). Give format=, '
                                     "e.g. 'ascii.csv' or 'fits'.".format(path, first, second)) from second
    return table, os.path.dirname(os.path.abspath(path))


def _find_column(table, wanted, role, default):
    """Column of table for role ('ra', 'dec', ...): the name given, matched
    case-insensitively, and, when the default name was kept, its aliases."""
    if wanted is None:
        return None
    lower = {name.lower(): name for name in table.colnames}
    names = [wanted]
    if wanted == default:
        names += list(_CATALOG_ALIASES[role])
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    if wanted != default or role in ('ra', 'dec'):
        raise ValueError('catalogue has no column {0!r} for {1} (columns: {2})'.format(
            wanted, role, ', '.join(table.colnames)))
    return None


def _column_values(table, name, role):
    col = table[name]
    mask = np.ma.getmaskarray(col)
    if mask.any():
        raise ValueError('catalogue column {0!r} ({1}) has missing values in rows {2}'.format(
            name, role, list(np.where(mask)[0][:10])))
    if col.unit is not None and role in ('ra', 'dec'):
        return col.quantity
    data = np.asarray(col)
    if data.dtype.kind in 'SU':
        data = np.char.strip(data.astype(str))
    if role == 'image':
        def missing(v):
            if v is None or (isinstance(v, (float, np.floating)) and not np.isfinite(v)):
                return True
            return str(v).strip().lower() in ('', 'nan', 'none', '--')
        empty = [k for k, v in enumerate(data) if missing(v)]
        if data.dtype.kind not in 'SUO':
            data = data.astype(str)
        if empty:
            raise ValueError('catalogue column {0!r} (image) is empty in rows {1}'.format(name, empty[:10]))
    elif role in ('z', 'weight'):
        try:
            values = np.asarray(data, dtype=float)
        except ValueError:
            raise ValueError('catalogue column {0!r} ({1}) is not numeric'.format(name, role))
        bad = np.where(~np.isfinite(values))[0]
        if bad.size:
            raise ValueError('catalogue column {0!r} ({1}) has non-finite values in rows {2}'.format(
                name, role, list(bad[:10])))
        if role == 'weight' and np.any(values < 0):
            raise ValueError('catalogue column {0!r} has negative weights in rows {1}'.format(
                name, list(np.where(values < 0)[0][:10])))
    return data


def _resolve_cube(name, base, imagedir):
    name = os.path.expanduser(str(name).strip())
    if os.path.isabs(name):
        return os.path.normpath(name)
    if imagedir is not None:
        return os.path.normpath(os.path.join(imagedir, name))
    if base is not None:
        candidate = os.path.normpath(os.path.join(base, name))
        if os.path.exists(candidate) or not os.path.exists(name):
            return candidate
    return os.path.normpath(name)


def readCatalog(catalog, ra='ra', dec='dec', z='z', weight='weight', image='image', unit='deg',
                frame='icrs', imagedir=None, format=None, verbose=None):
    """
        Reads the stacking positions of all sources from one catalogue table,
        one row per source. The cube of each source is given in the table
        (image column); it is never guessed.

        Parameters
        ---------
        catalog
            Path to a table file: CSV (with a header line), ECSV, FITS table,
            VOTable or any format astropy.table reads (see **format**). Also
            accepted: an astropy Table, a pandas DataFrame or a dict of columns.
        ra, dec
            Column names of the positions (default 'ra', 'dec'; also found:
            RAJ2000/DEJ2000, ra_deg/dec_deg...). Numbers are in **unit** unless
            the column has a unit (FITS, ECSV, VOTable); strings are
            sexagesimal with RA in hours (CASA's '-027.40.59.866' accepted).
        z
            Column of the spectral value (default 'z', also 'redshift'): the
            redshift, or the observed frequency/wavelength or channel for
            spectralMethod 'obsSpec'/'channel'. If the catalogue has no such
            column (or z=None) the positions are for continuum stacking.
        weight
            Column of the weights (optional, default 1).
        image
            Column naming the cube of each source (FITS file or CASA image;
            default 'image', also 'cube', 'file'). Relative paths are looked
            up next to the catalogue first, then from the working directory.
            Without this column, give the stacking function a single cube.
        unit
            'deg' (default) or 'rad', for unit-less numeric positions.
        imagedir
            Folder to which relative cube names are relative (overrides the
            rule above).
        format
            astropy table format, e.g. 'ascii.csv', 'ascii.commented_header'
            (default: guessed from the file).
        verbose
            report what was read (default: LineStacker.verbose)

        Returns
        -------
        list of LineStacker.Coord
            for LineStacker.image.stack and the analysis tools; imagenames can
            be omitted there when the catalogue has an image column.
    """
    table, base = _read_table(catalog, format)
    if len(table) == 0:
        raise ValueError('the catalogue is empty')
    cols = {}
    for role, wanted, default in (('ra', ra, 'ra'), ('dec', dec, 'dec'), ('z', z, 'z'),
                                  ('weight', weight, 'weight'), ('image', image, 'image')):
        cols[role] = _find_column(table, wanted, role, default)
    values = {role: _column_values(table, name, role) for role, name in cols.items() if name is not None}
    images = None
    if 'image' in values:
        images = [_resolve_cube(v, base, imagedir) for v in values['image']]
    coords = makeCoords(values['ra'], values['dec'], z=values.get('z'), weight=values.get('weight'),
                        image=images, unit=unit, frame=frame, verbose=False)
    source = os.path.basename(os.fspath(catalog)) if isinstance(catalog, (str, os.PathLike)) else 'the table'
    ncubes = len(set(images)) if images is not None else 0
    _info('read {0} position{1}{2} from {3}{4}'.format(
        len(coords), 's' if len(coords) != 1 else '',
        ' in {0} cube{1}'.format(ncubes, 's' if ncubes != 1 else '') if images is not None else '',
        source, '' if 'z' in values else ' (no {0!r} column: continuum positions)'.format(z or 'z')), verbose)
    return coords


def _cube_list(coords, imagenames):
    """The cubes to use: imagenames as given, or else the cubes named by the
    coordinates (readCatalog image column, makeCoords image=...)."""
    if imagenames is not None and not isinstance(imagenames, (str, os.PathLike)) and len(imagenames) > 0:
        return [os.fspath(n) for n in imagenames]
    if isinstance(imagenames, (str, os.PathLike)):
        return [os.fspath(imagenames)]
    if coords is None or len(coords) == 0:
        raise ValueError('imagenames is empty: give the list of cubes to stack')
    names = []
    for c in coords:
        if not isinstance(c.image, (str, os.PathLike)):
            raise ValueError('imagenames was not given and the coordinate {0!r} does not name its cube. '
                             'Give imagenames, or the cube of each source (image column of readCatalog, '
                             'image= of makeCoords).'.format(c))
        name = os.fspath(c.image)
        if name not in names:
            names.append(name)
    return names


def _image_index(coord, imagenames):
    """Index in imagenames of the cube a coordinate belongs to."""
    img = coord.image
    if isinstance(img, (int, np.integer)) and not isinstance(img, bool):
        if 0 <= img < len(imagenames):
            return int(img)
    elif isinstance(img, (str, os.PathLike)):
        img = os.fspath(img)
        if img in imagenames:
            return imagenames.index(img)
        norm = os.path.normpath(os.path.abspath(img))
        for k, name in enumerate(imagenames):
            if os.path.normpath(os.path.abspath(name)) == norm:
                return k
    elif img is None and len(imagenames) == 1:
        return 0                       # positions without a cube, single cube given
    num = coord.imagenumber
    if num is not None and 0 <= num < len(imagenames):
        return int(num)
    if num is None or num < 0:
        raise ValueError('The cube {0!r} of coordinate {1!r} is not in imagenames {2}'.format(
            img, coord, imagenames if len(imagenames) <= 5 else imagenames[:5] + ['...']))
    raise ValueError('Cannot tell which image coordinate {0!r} belongs to: its image {1!r} is not '
                     'in imagenames and its imagenumber ({2}) is out of range'.format(coord, img, num))


def _to_skycoord(coords):
    """SkyCoord array (in the frame of the first coordinate) from physical
    Coords given as astropy angles or radians; coordinates in other frames
    are converted."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    def val(a):
        return a if hasattr(a, 'unit') else a * u.rad
    frames = [_equatorial(c.frame or 'icrs') for c in coords]
    target = frames[0]
    ra = u.Quantity([val(c.ra) for c in coords])
    dec = u.Quantity([val(c.dec) for c in coords])
    if all(f == target for f in frames):
        return SkyCoord(ra=ra, dec=dec, frame=target)
    out_ra, out_dec = ra.to(u.deg).copy(), dec.to(u.deg).copy()
    for f in set(frames) - {target}:
        sel = np.array([g == f for g in frames])
        conv = SkyCoord(ra=ra[sel], dec=dec[sel], frame=f).transform_to(target)
        out_ra[sel], out_dec[sel] = conv.ra.to(u.deg), conv.dec.to(u.deg)
    return SkyCoord(ra=out_ra, dec=out_dec, frame=target)


def getPixelCoords(coords, imagenames=None, hdul='auto'):
    """
        Creates pixel coordinate list from a physical coordinate
        list and a list of images.

        Parameters
        ---------
        coords
            A list of LineStacker.Coord coordinates.
        imagenames
            A list of images' paths (FITS files or CASA images). Optional when
            every coordinate names its cube (readCatalog, makeCoords).
        hdul
            index of the FITS HDU holding the cube; 'auto' (default) uses the
            first HDU containing an image (HDU 0 for ALMA/CASA FITS, 1 for MUSE).

        Returns
        -------
        list of LineStacker.Coord
            new coordinates of type 'pixel'; the original coordinates are not
            modified. Each new coordinate's ``index`` is the position of the
            original in ``coords``.
    """
    from LineStacker import _io
    imagenames = _cube_list(coords, imagenames)
    pixcoords = [None] * len(coords)
    groups = {}
    for n, coord in enumerate(coords):
        k = _image_index(coord, imagenames)
        if coord.type == 'pixel':
            new = coord.copy()
            new.image, new.imagenumber, new.index = imagenames[k], k, n
            pixcoords[n] = new
        else:
            groups.setdefault(k, []).append(n)
    for k, members in groups.items():
        info = _io.read_info(imagenames[k], hdul)
        sky = _to_skycoord([coords[n] for n in members])
        xs, ys = info.world_to_pixel(sky)
        for n, x, y in zip(members, np.atleast_1d(xs), np.atleast_1d(ys)):
            old = coords[n]
            new = old.copy()
            new.ra, new.dec = float(x), float(y)
            new.type, new.image, new.imagenumber, new.index = 'pixel', imagenames[k], k, n
            pixcoords[n] = new
    return pixcoords


def randomCoords(imagenames, ncoords=10, hdul='auto', margin=0, rng=None):
    """
        Random positions anywhere in a set of images: pixel coordinates whose
        spec is a random channel of their cube, to be stacked with
        spectralMethod='channel'.

        Parameters
        ---------
        imagenames
            A list of images paths (or one path).
        ncoords
            Number of random coordinates. Default is 10.
        margin
            Minimal distance (pixels) to the image edges, e.g. half the stampsize.
        rng
            numpy Generator or seed, for reproducible draws.
    """
    from LineStacker import _io
    if isinstance(imagenames, (str, os.PathLike)):
        imagenames = [imagenames]
    imagenames = [os.fspath(n) for n in imagenames]
    if not imagenames:
        raise ValueError('randomCoords needs at least one image')
    if isinstance(ncoords, bool) or int(ncoords) != ncoords or ncoords < 0:
        raise ValueError('ncoords must be an integer >= 0 (got {0!r})'.format(ncoords))
    ncoords = int(ncoords)
    if margin < 0:
        raise ValueError('margin must be >= 0 (got {0})'.format(margin))
    rng = np.random.default_rng(rng)
    infos = [_io.read_info(name, hdul) for name in imagenames]
    for name, inf in zip(imagenames, infos):
        if 2 * margin > min(inf.nx, inf.ny) - 1:
            raise ValueError('margin {0} leaves no room in {1} ({2}x{3} pixels): it must be at most {4}'.format(
                margin, name, inf.nx, inf.ny, (min(inf.nx, inf.ny) - 1) / 2.))
    out = []
    for _ in range(ncoords):
        k = int(rng.integers(len(imagenames)))
        inf = infos[k]
        x = rng.uniform(margin, inf.nx - 1 - margin)
        y = rng.uniform(margin, inf.ny - 1 - margin)
        chan = int(rng.integers(inf.nchan))
        out.append(Coord(x, y, spec=chan, image=imagenames[k], imagenumber=k, type='pixel'))
    return out


def _image_list(imagenames):
    """imagenames as a list of paths (one path is a list of one)."""
    if imagenames is None:
        return None
    if isinstance(imagenames, (str, os.PathLike)):
        return [os.fspath(imagenames)]
    return [os.fspath(n) for n in imagenames]


def randomizeCoords(coords, beam, maxBeamRange=5, imagenames=None, hdul='auto', rng=None, maxTries=100):
    """
        New positions at a random distance between beam and maxBeamRange*beam
        from each of the original coordinates (random direction).

        Parameters
        ---------
        coords
            list of original coordinates (LineStacker.Coord instances)
        beam
            beam size is radians,
            new random coordinates will be at a minimum distance beam
            from the original coordinates.
        maxBeamRange
            maximum distance from original coordinates
            at which new coordinates can be located, in units of beams.
            Default is 5.
        imagenames
            images the coordinates refer to; when given, new positions are
            kept inside their image.
        rng
            numpy Generator or seed, for reproducible draws.
    """
    from LineStacker import _io
    imagenames = _image_list(imagenames)
    if not np.isfinite(beam) or beam <= 0:
        raise ValueError('beam must be > 0 (radians), got {0}'.format(beam))
    if not (np.isfinite(maxBeamRange) and maxBeamRange >= 1):
        raise ValueError('maxBeamRange must be finite and >= 1 (in beams), got {0}'.format(maxBeamRange))
    if isinstance(maxTries, bool) or int(maxTries) != maxTries or maxTries < 1:
        raise ValueError('maxTries must be an integer >= 1 (got {0!r})'.format(maxTries))
    maxTries = int(maxTries)
    rng = np.random.default_rng(rng)
    out = []
    for coord in coords:
        info = None
        if imagenames is not None:
            info = _io.read_info(imagenames[_image_index(coord, imagenames)], hdul)
        for _ in range(maxTries):
            sep = rng.uniform(beam, maxBeamRange * beam)
            pa = rng.uniform(0, 2 * np.pi)
            new = coord.copy()
            if coord.type == 'pixel':
                if info is None:
                    raise ValueError('imagenames is needed to randomize pixel coordinates')
                step = sep / np.deg2rad(info.pixscale_deg)
                new.ra = coord.ra + step * np.cos(pa)
                new.dec = coord.dec + step * np.sin(pa)
                px, py = new.ra, new.dec
            else:
                import astropy.units as u
                sky = _to_skycoord([coord])[0].directional_offset_by(pa * u.rad, sep * u.rad)
                new.ra, new.dec = sky.ra, sky.dec
                px, py = info.world_to_pixel(sky) if info is not None else (0, 0)
            if info is None or (0 <= px <= info.nx - 1 and 0 <= py <= info.ny - 1):
                out.append(new)
                break
        else:
            raise Exception('Cannot find random position, change search region')
    return out

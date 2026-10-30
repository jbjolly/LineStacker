'''
Tools to test stacks and prepare data for stacking: bootstrapping,
subsampling, stacks at random positions, noise estimates, rebinning and
regridding of lines of different widths and redshifts.

The Monte-Carlo tools (bootstraping_OneD, bootStraping_Cube, subsample_OneD,
stack_estimator) accept ``n_jobs`` (number of CPU cores, -1 for all),
``backend`` ('process' or 'thread', see LineStacker.tools.parallel),
``seed`` and ``verbose``. All random draws happen before the work is split, so results do not
depend on n_jobs. With seed=None the draws come from numpy's global random
state, so np.random.seed() still makes runs reproducible.
'''

import copy
import pickle
import warnings

import numpy as np

import LineStacker
import LineStacker.tools
from LineStacker.tools import parallel


# ---------------------------------------------------------------------------
# random draws
# ---------------------------------------------------------------------------

class _Draws(object):
    """Random numbers from numpy's global state (seed=None) or a seeded Generator."""

    def __init__(self, seed=None):
        self.rng = None if seed is None else np.random.default_rng(seed)

    def integers(self, low, high, size=None):
        if self.rng is None:
            return np.random.randint(low, high, size=size)
        return self.rng.integers(low, high, size=size)

    def uniform(self, low, high, size=None):
        if self.rng is None:
            return np.random.uniform(low, high, size=size)
        return self.rng.uniform(low, high, size=size)

    def random(self, size):
        if self.rng is None:
            return np.random.random_sample(size)
        return self.rng.random(size)


def _n_random(nRandom, minimum=0):
    """Validated number of Monte-Carlo iterations."""
    if isinstance(nRandom, bool) or int(nRandom) != nRandom or nRandom < minimum:
        raise ValueError('nRandom must be an integer >= {0} (got {1!r})'.format(minimum, nRandom))
    return int(nRandom)


def _check_spectra(spectra, what='spectra'):
    import LineStacker.OneD_Stacker as OneD
    spectra = list(spectra)
    if not spectra:
        raise ValueError('no {0} given'.format(what))
    if not all(isinstance(sp, OneD.Spectrum) for sp in spectra):
        raise Exception('spectra must be LineStacker.OneD_Stacker.Spectrum objects, please initialize your spectra accordingly')
    return spectra


def _z_reference(spectra, regridMethod):
    zs = [sp.z for sp in spectra]
    if any(z is None for z in zs):
        raise ValueError('regridding from redshifts needs the z of every spectrum (Spectrum(z=...))')
    if regridMethod == 'scaleToMin':
        return min(zs)
    if regridMethod == 'scaleToMax':
        return max(zs)
    raise ValueError("regridMethod must be 'scaleToMin' or 'scaleToMax' (got {0!r})".format(regridMethod))


def _chunks(n_items, per_item_cost, n_jobs, budget=8e6):
    """Slices of range(n_items): small enough for memory, enough for workers."""
    size = max(1, int(budget // max(1, per_item_cost)))
    n_workers = LineStacker.tools.resolve_n_jobs(n_jobs)
    if n_workers > 1:
        size = min(size, max(1, -(-n_items // (4 * n_workers))))
    return [slice(i, min(n_items, i + size)) for i in range(0, n_items, size)]


def _picklable(obj):
    try:
        pickle.dumps(obj)
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# 1D bootstrapping
# ---------------------------------------------------------------------------

def _summarise_1d(stack, save):
    import LineStacker.tools.fit as fitTools
    if save == 'all':
        return stack
    if save == 'amp':
        return np.nanmax(stack) if np.any(np.isfinite(stack)) else np.nan
    if save == 'ampAndWidth':
        try:
            fited = fitTools.GaussFit(fctToFit=stack, returnInfos=True)
        except RuntimeError:               # fit did not converge
            return [np.nan, np.nan]
        return [fited[1][0], fited[1][2]]
    raise ValueError("save must be 'all', 'amp' or 'ampAndWidth'")


def _warn_failed_fits(saved, save):
    if save not in ('ampAndWidth', 'outflow'):
        return
    bad = sum(1 for s in saved if np.ndim(s) and np.size(s) and
              (np.any(np.isnan(np.asarray(s, dtype=float))) or np.all(np.asarray(s, dtype=float) == -100)))
    if bad:
        warnings.warn('{0} of {1} Gaussian fits did not converge; they are stored as {2}'.format(
            bad, len(saved), '[nan, nan]' if save == 'ampAndWidth' else '[-100]*6'))


def _quiet_fits():
    """Context ignoring scipy's 'covariance could not be estimated': the
    Monte-Carlo tools only use the fitted parameters."""
    from scipy.optimize import OptimizeWarning
    ctx = warnings.catch_warnings()
    ctx.__enter__()
    warnings.simplefilter('ignore', OptimizeWarning)
    return ctx


def _boot1d_worker(sl):
    import LineStacker.OneD_Stacker as OneD
    al = parallel.get_shared('aligned')
    idx = parallel.get_shared('idx')[sl]
    method, save = parallel.get_shared('method'), parallel.get_shared('save')
    stacks, _ = OneD.combine_indices(al, idx, method)
    ctx = _quiet_fits()
    try:
        return [_summarise_1d(OneD.crop(al, stacks[r], idx[r]), save) for r in range(len(idx))]
    finally:
        ctx.__exit__(None, None, None)


def bootstraping_OneD(    spectra,
                        nRandom=1000,
                        chansStack='full',
                        method='mean',
                        center='lineCenterIndex',
                        save='all',
                        seed=None,
                        n_jobs=1,
                        backend=None,
                        verbose=None):

    """
        Bootstrapping of a 1D stack: stacks of random draws, with replacement, of the spectra. See LineStacker.OneD_Stacker.Stack for the stacking parameters.

        Parameters
        ---------

        spectra
            a list of LineStacker.OneD_Stacker.Spectrum objects
        nRandom
            number of bootstrap iterations
        chansStack
            number of channels to stack, either a fixed number or 'full' for entire spectra
        method
            stacking method, 'mean' or 'median'
        center
            Method to center spectra.
            See LineStacker.OneD_Stacker.Stack for further information on centering methods.
        save
            data to save at each bootstrap iterations\n
            possible values are **'all'**, **'amp'** and **'ampAndWidth'**\n
            **'all'** saves the full stack at each bootstrap iteration (can use a lot of memory)\n
            **'amp'** saves the amplitude of the line (maximum value of the stack) at each bootstrap iteration, fastest\n
            **'ampAndWidth'** fits the line with a gaussian and saves the corresponding amplitude and width at each bootstrap iteration (fits that do not converge are stored as [nan, nan]), can be cpu expensive.
        seed
            seed for reproducible resampling (None uses numpy's global random state)
        n_jobs
            number of CPU cores to use (-1 for all). Worth it mostly for 'ampAndWidth' and for medians.
        backend
            'process' or 'thread', see LineStacker.tools.parallel
        verbose
            print what is done and a progress bar (default: LineStacker.verbose)

        Returns
        -------
        list
            nRandom+1 entries (as chosen by save), the first one for the stack
            of all spectra.
    """
    import LineStacker.OneD_Stacker as OneD

    spectra = _check_spectra(spectra)
    nRandom = _n_random(nRandom)
    if save not in ('all', 'amp', 'ampAndWidth'):
        raise ValueError("save must be 'all', 'amp' or 'ampAndWidth'")
    baseStack = OneD.Stack(spectra, chansStack=chansStack, method=method, center=center)
    saved = [_summarise_1d(baseStack[0], save)]

    # Centres do not depend on which spectra are drawn: align everything once.
    al = OneD.Aligned(spectra, chansStack)
    n = len(spectra)
    idx = _Draws(seed).integers(0, n, size=(nRandom, n))

    LineStacker._info('bootstrapping {0} spectra, {1} iterations'.format(n, nRandom), verbose)
    tasks = _chunks(nRandom, n * al.L, n_jobs)
    if backend is None and save != 'ampAndWidth':
        backend = 'thread'          # numpy-bound work: no need to copy data around
    results = parallel.map_chunks(_boot1d_worker, tasks, n_jobs=n_jobs, backend=backend,
                                  shared={'aligned': al, 'idx': idx, 'method': method, 'save': save},
                                  progress=LineStacker._verbose(verbose))
    for r in results:
        saved.extend(r)
    _warn_failed_fits(saved, save)
    return saved


# ---------------------------------------------------------------------------
# cube bootstrapping
# ---------------------------------------------------------------------------

def _stamp_spectrum(pixels):
    """Spectrum summed over a stacked stamp (x, y, stokes, chan): blank
    pixels left out, NaN for channels without any valid pixel."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        out = np.nansum(pixels, axis=(0, 1, 2))
    out[~np.isfinite(pixels).any(axis=(0, 1, 2))] = np.nan
    return out


def _summarise_cube_spectrum(spectrum, save):
    import LineStacker.tools.fit as fitTools
    if save == 'amp':
        return np.nanmax(spectrum) if np.any(np.isfinite(spectrum)) else np.nan
    if save == 'ampAndWidth':
        try:
            fited = fitTools.GaussFit(fctToFit=spectrum, returnInfos=True)
        except RuntimeError:               # fit did not converge
            return [np.nan, np.nan]
        return [fited[1][0], fited[1][2]]
    if save == 'outflow':
        try:
            fited = fitTools.DoubleGaussFit(fctToFit=spectrum, returnInfos=True)
            return [fited[1]]
        except RuntimeError:
            return [-100 for i in range(6)]
    raise ValueError("save must be 'all', 'amp', 'ampAndWidth' or 'outflow'")


def _bootcube_worker(sl):
    """One chunk of cube bootstraps."""
    ctx = _quiet_fits()
    try:
        return _bootcube_chunk(sl)
    finally:
        ctx.__exit__(None, None, None)


def _bootcube_chunk(sl):
    import LineStacker.image as li
    idx = parallel.get_shared('idx')[sl]
    method, save = parallel.get_shared('method'), parallel.get_shared('save')
    out = []
    if parallel.get_shared('reduced'):
        # mean stack summed over the stamp: linear, so sum the sources first
        S = parallel.get_shared('spectra')            # (nsrc, nchan) weighted spectra
        W = parallel.get_shared('weights')            # (nsrc, nchan)
        for r in range(len(idx)):
            cnt = np.bincount(idx[r], minlength=W.shape[0]).astype(float)
            den = cnt @ W
            with np.errstate(invalid='ignore', divide='ignore'):
                spec = np.where(den > 0, (cnt @ S) / den, np.nan)
            out.append(_summarise_cube_spectrum(spec, save))
        return out
    data, W, inside = parallel.get_shared('data'), parallel.get_shared('weights'), parallel.get_shared('inside')
    for r in range(len(idx)):
        sel = idx[r]
        pixels = li._combine(data[sel], W[sel], inside[sel], method)
        if save == 'all':
            out.append(pixels)
        else:
            out.append(_summarise_cube_spectrum(_stamp_spectrum(pixels), save))
    return out


def bootStraping_Cube(    coords,
                        stampsize = 32,
                        imagenames=None,
                        method = 'mean',
                        weighting = None,
                        maxmaskradius=0,
                        fEm = 0,
                        spectralMethod='z',
                        N_chans=30,
                        nRandom=1000,
                        save='amp',
                        lEm=-1,
                        hdul='auto',
                        maskradius=None,
                        seed=None,
                        n_jobs=1,
                        backend=None,
                        continuum=False,
                        collapse='sum',
                        verbose=None):#'all', 'amp', 'ampAndWidth'

    """
        Bootstrapping of a cube stack: stacks of random draws, with replacement, of the sources. See LineStacker.image.stack for the stacking parameters.

        Parameters
        ---------
        coords
            A list of LineStacker.Coord of all target coordinates
        stampsize
            size of target image in pixels
        imagenames
            Name of images to extract cubes from
        method
            stacking method, 'mean' or 'median'
        weighting
            weighting method to use if stacking method is mean\n
            possible values are 'sigma2', 'sigma2F', '1/A' or None (weights of the coords)\n
            see LineStacker.image.stack for a complete description of weighting methods\n
        maxmaskradius
            radius of the mask used to blank the centre pixels in weight calculation (also accepted as maskradius)
        fEm
            rest emission frequency of the line (or lEm for wavelength cubes)
        N_chans
            number of channels of the resulting stack
        nRandom
            number of bootstrap iterations
        save
            data to save at each bootstrap iterations\n
            possible values are **'all'**, **'amp'**, **'ampAndWidth'** and **'outflow'**\n
            **'all'** saves the full stack at each bootstrap iteration (can use a lot of memory)\n
            **'amp'** saves the amplitude (maximum value of the stack) of the line, at each bootstrap iteration, fastest\n
            **'ampAndWidth'** fits the line with a gaussian and saves the corresponding amplitude and width at each bootstrap iteration (fits that do not converge are stored as [nan, nan]), can be cpu expensive\n
            **'outflow'** fits the line with two gaussian components and saves the stack parameters at each bootstrap iteration, can be cpu expensive\n
            for 'amp', 'ampAndWidth' and 'outflow' the line is obtained by summing all pixels inside the stack stamp.
        seed, n_jobs, backend
            see the module description
        continuum, collapse
            continuum stacking, see LineStacker.image.stack ('amp' is then the
            stamp-summed flux; 'ampAndWidth' and 'outflow' need a spectral axis)
        verbose
            print what is done and a progress bar (default: LineStacker.verbose)

        Returns
        -------
        list
            nRandom+1 entries (as chosen by save), the first one for the stack
            of all sources.
    """
    import LineStacker.image as li

    if maskradius is None:
        maskradius = maxmaskradius
    if fEm == 0:
        fEm = -1
    nRandom = _n_random(nRandom)
    if save not in ('all', 'amp', 'ampAndWidth', 'outflow'):
        raise ValueError("save must be 'all', 'amp', 'ampAndWidth' or 'outflow'")
    prep = li._prepare(coords, imagenames, stampsize=stampsize, N_chans=N_chans,
                       spectralMethod=spectralMethod, fEm=fEm, lEm=lEm, hdul=hdul,
                       method=method, weighting=weighting, maskradius=maskradius,
                       continuum=continuum, collapse=collapse)
    if prep.continuum and save in ('ampAndWidth', 'outflow'):
        raise ValueError("save='{0}' fits a line profile; it needs line stacking, not continuum".format(save))
    W = prep.weights()
    stacked_im = li._combine(prep.data, W, prep.valid, method)
    if save == 'all':
        saved = [stacked_im]
    else:
        saved = [_summarise_cube_spectrum(_stamp_spectrum(stacked_im), save)]

    n = len(prep.coords)
    idx = _Draws(seed).integers(0, n, size=(nRandom, n))
    # Mean stacks summed over the stamp are linear: the sources can be summed
    # first, unless blank pixels make the sum depend on which sources are drawn.
    usable = np.isfinite(prep.data) | (W[:, None, None, None, :] == 0)
    reduced = method == 'mean' and save != 'all' and bool(usable.all())
    if reduced:
        spectra = np.sum(np.where(np.isfinite(prep.data), prep.data, 0.), axis=(1, 2, 3)) * W   # (nsrc, nchan)
        shared = {'reduced': True, 'spectra': spectra, 'weights': W}
        cost = n
    else:
        shared = {'reduced': False, 'data': prep.data, 'weights': W, 'inside': prep.valid}
        cost = prep.data[0].size * n
    shared.update({'idx': idx, 'method': method, 'save': save})
    if backend is None and save in ('all', 'amp'):
        backend = 'thread'
    LineStacker._info('bootstrapping {0} sources, {1} iterations'.format(n, nRandom), verbose)
    results = parallel.map_chunks(_bootcube_worker, _chunks(nRandom, cost, n_jobs, budget=5e7),
                                  n_jobs=n_jobs, backend=backend, shared=shared,
                                  progress=LineStacker._verbose(verbose))
    for r in results:
        saved.extend(r)
    _warn_failed_fits(saved, save)
    return saved


# ---------------------------------------------------------------------------
# resampling helpers
# ---------------------------------------------------------------------------

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    """
        Arbitrary resampling of an array to new dimension sizes (same number
        of dimensions), with the sampling points of IDL's congrid.

        Parameters
        ----------
        a
            array to resample
        newdims
            new size of every dimension
        method
            **'neighbour'** (closest value), **'nearest'** or **'linear'**
            (n 1-D interpolations with scipy.interpolate.interp1d), or
            **'spline'** (scipy.ndimage.map_coordinates, cubic)
        centre
            True: interpolation points at the centres of the bins; False
            (default): at their front edges
        minusone
            False (default): an axis of i samples resized to x is sampled by
            factors of i/x; True: by (i-1)/(x-1), which avoids extrapolating
            one sample beyond the input

        Blank (NaN) values are not handled specially here; LineStacker's
        rebinning functions handle them (see rebin_OneD).

        From https://scipy-cookbook.readthedocs.io/items/Rebinning.html
        (Copyright (c) 2001, 2002 Enthought, Inc.; Copyright (c) 2003-2017
        SciPy Developers. All rights reserved.)
    """
    import scipy.interpolate
    import scipy.ndimage

    a = np.asarray(a)
    if a.dtype not in [np.float64, np.float32]:
        a = a.astype(float)

    m1 = int(bool(minusone))
    if m1 and np.any(np.asarray(newdims) == 1):
        raise ValueError('[congrid] minusone=True needs new dimensions > 1 (got {0})'.format(list(newdims)))
    ofs = int(bool(centre)) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        raise ValueError("[congrid] dimensions error. "
                         "This routine currently only support "
                         "rebinning to the same number of dimensions.")
    newdims = np.asarray( newdims, dtype=float )
    if np.any(newdims < 1) or np.any(newdims != np.round(newdims)):
        raise ValueError('[congrid] new dimensions must be integers >= 1 (got {0})'.format(list(newdims)))
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = np.indices(newdims.astype(int))[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = np.array( dimlist ).round().astype(int)
        cd = np.stack([np.clip(cd[i], 0, old[i] - 1) for i in range(ndims)])
        newa = a[tuple( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = np.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [np.arange(i, dtype = float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + list(range( ndims - 1 ))
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        nslices = [ slice(0,int(j)) for j in list(newdims) ]
        newcoords = np.mgrid[tuple(nslices)].astype(float)

        newcoords_dims = list(range(newcoords.ndim))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        raise ValueError("Congrid error: Unrecognized interpolation type. "
                         "Currently only 'neighbour', 'nearest', 'linear' "
                         "and 'spline' are supported.")


def _sample_positions(old, new, centre=False, minusone=False):
    """Positions (in old samples) where congrid samples an axis of length
    old resized to new."""
    m1 = int(bool(minusone))
    ofs = int(bool(centre)) * 0.5
    return (old - m1) / float(new - m1) * (np.arange(new) + ofs) - ofs


def _congrid_nan(a, newdims, method='linear', centre=False, minusone=False):
    """congrid for arrays with blank (NaN or infinite) values, which are never
    used to make a value.

    'linear': an output is NaN when its interpolation uses a blank sample
    (the two nearest samples along each resampled axis); the others are
    exactly congrid's. 'spline': a cubic spline (as congrid's) through each
    run of valid samples between blanks, so values depend on valid data only;
    outputs between a valid sample and a blank one, and outputs beyond the
    last sample (congrid extrapolates when it enlarges an axis), are NaN.
    'neighbour', 'nearest': congrid (a blank is copied where it is picked)."""
    import scipy.ndimage
    a = np.asarray(a, dtype=float)
    dims = np.asarray(newdims, dtype=float)
    if dims.shape != (a.ndim,) or np.any(dims < 1) or np.any(dims != np.round(dims)):
        raise ValueError('[congrid] new dimensions must be {0} integers >= 1 (got {1})'.format(a.ndim, list(newdims)))
    if minusone and np.any(dims == 1):
        raise ValueError('[congrid] minusone=True needs new dimensions > 1 (got {0})'.format(list(newdims)))
    if method not in ('neighbour', 'nearest', 'linear', 'spline'):
        raise ValueError("[congrid] method must be 'neighbour', 'nearest', 'linear' or 'spline'")
    if method == 'spline':
        out = a
        for i in range(a.ndim):
            n_old, n_new = out.shape[i], int(newdims[i])
            if n_new == n_old and not minusone and not centre:
                continue                          # identity: the spline returns the samples
            out = _spline_axis(out, i, _sample_positions(n_old, n_new, centre, minusone))
        return out
    blank = ~np.isfinite(a)
    if method in ('neighbour', 'nearest') or not blank.any():
        return congrid(a, newdims, method=method, centre=centre, minusone=minusone)
    values = congrid(np.where(blank, 0., a), newdims, method=method, centre=centre, minusone=minusone)
    axes = []
    for i in range(a.ndim):
        if int(newdims[i]) == a.shape[i] and not minusone and not centre:
            axes.append(np.arange(a.shape[i], dtype=float))   # identity: exact samples
        else:
            axes.append(_sample_positions(a.shape[i], int(newdims[i]), centre, minusone))
    coords = np.meshgrid(*axes, indexing='ij')
    reach = scipy.ndimage.map_coordinates(blank.astype(float), coords, order=1, mode='nearest')
    values = np.asarray(values, dtype=float)
    values[reach > 1e-9] = np.nan
    return values


def _spline_axis(a, axis, pos):
    """Cubic-spline resampling of a along one axis at positions pos (in
    samples), separately through every run of finite samples; NaN outside
    the runs and beyond the ends of the axis."""
    import scipy.ndimage
    moved = np.moveaxis(a, axis, -1)
    n = moved.shape[-1]
    lines = moved.reshape(-1, n)
    out = np.full((lines.shape[0], len(pos)), np.nan)
    tol = 1e-9
    for k in range(lines.shape[0]):
        y = lines[k]
        ok = np.isfinite(y)
        if ok.all():
            runs = [(0, n - 1)]
        else:
            edges = np.diff(np.concatenate([[0], ok.astype(int), [0]]))
            runs = list(zip(np.where(edges == 1)[0], np.where(edges == -1)[0] - 1))
        for first, last in runs:
            sel = (pos >= first - tol) & (pos <= last + tol)
            if sel.any():
                local = np.clip(pos[sel], first, last) - first
                out[k, sel] = scipy.ndimage.map_coordinates(y[first:last + 1], [local], order=3, mode='constant')
    return np.moveaxis(out.reshape(moved.shape[:-1] + (len(pos),)), -1, axis)


def _linear_at(y, pos):
    """Linear interpolation of the 1D array y at positions pos (in samples):
    NaN beyond the ends and where a blank (NaN) sample would be used."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    out = np.full(len(pos), np.nan)
    inside = (pos >= -1e-9) & (pos <= n - 1 + 1e-9)
    p = np.clip(pos, 0., n - 1.)
    j0 = np.minimum(np.floor(p + 1e-9).astype(int), n - 1)
    j1 = np.minimum(j0 + 1, n - 1)
    t = np.clip(p - j0, 0., 1.)
    with np.errstate(invalid='ignore'):
        val = np.where(t <= 1e-9, y[j0], np.where(t >= 1. - 1e-9, y[j1], y[j0] * (1. - t) + y[j1] * t))
    out[inside] = val[inside]
    return out


def _axis_at(axis, pos):
    """Spectral axis at positions pos (in channels): linear interpolation,
    and linear extrapolation beyond the last channel."""
    axis = np.asarray(axis, dtype=float)
    out = np.interp(pos, np.arange(len(axis)), axis)
    beyond = pos > len(axis) - 1
    out[beyond] = axis[-1] + (pos[beyond] - (len(axis) - 1)) * (axis[-1] - axis[-2])
    return out


def _resample_spectrum(spectrum, ratio, method):
    """Copy of a Spectrum whose channels are `ratio` times narrower: new
    channel j samples old channel j/ratio, int(n*ratio) channels. Fluxes
    are interpolated (method 'linear' or 'spline') from valid channels only
    (NaN next to a blank channel and beyond the last channel), per-channel
    weights linearly, spectral axes linearly (extrapolated at the end); the
    line centre, channel widths and spectrum array follow."""
    new = copy.deepcopy(spectrum)
    n = len(spectrum.amp)
    m = int(n * ratio)
    if m < 2:
        raise ValueError('spectrum {0!r}: resampling {1} channels by {2:.3g} leaves fewer than 2'.format(
            spectrum.name, n, ratio))
    pos = np.arange(m) / float(ratio)
    amp = np.asarray(spectrum.amp, dtype=float)
    if method == 'spline':
        new.amp = _spline_axis(amp, 0, pos)
    else:
        new.amp = _linear_at(amp, pos)
    for name in ('velocities', 'frequencies'):
        ax = getattr(spectrum, name, [])
        if len(ax) == n:
            setattr(new, name, _axis_at(ax, pos))
    w = spectrum.weights
    if not isinstance(w, str) and np.ndim(w) > 0 and len(w) > 1:
        w = np.asarray(w, dtype=float)
        new.weights = _linear_at(w, np.arange(max(2, int(len(w) * ratio))) / float(ratio))
    if len(new.velocities) == m:
        new.velBin = (new.velocities[-1] - new.velocities[0]) / (m - 1.)
    if len(new.frequencies) == m:
        new.freqBin = (new.frequencies[-1] - new.frequencies[0]) / (m - 1.)
    sp = np.asarray(spectrum.spectrum)
    if sp.size:
        axis = new.velocities if spectrum.velOrFreq == 'vel' else new.frequencies
        pair = np.vstack([axis, new.amp])
        new.spectrum = pair.T if (sp.ndim == 2 and sp.shape[1] == 2) else pair
    if spectrum.centerIndex is not None:
        new.centerIndex = LineStacker.round_half_up(spectrum.centerIndex * ratio)
    if new.fit:
        new.fitAgain()
    return new


def rebin_OneD(spectra, verbose=None):
    """
        Rebin a list of spectra so that all lines span as many channels as
        the narrowest one (widths in channels, from Gaussian fits).
        **The lines must be visible source by source.**
        Fluxes are linearly interpolated; a rebinned channel interpolated
        from a blank (NaN) channel is NaN. The spectral axes, per-channel
        weights and line centres (centerIndex) follow the new channels.

        Parameters
        ----------
        spectra
            A list of LineStacker.OneD_Stacker.Spectrum objects
        verbose
            print what is done (default: LineStacker.verbose)

        Returns
        -------
        list of LineStacker.OneD_Stacker.Spectrum
            the rebinned spectra (copies).
    """

    import LineStacker.tools.fit as fitTools
    spectra = _check_spectra(spectra)
    widths = []
    for spectrum in spectra:
        # line widths in channels (rebinning changes the number of channels)
        if spectrum.fit:
            axis = spectrum.velocities if spectrum.velOrFreq == 'vel' else spectrum.frequencies
            step = abs(axis[-1] - axis[0]) / (len(axis) - 1.)
            widths.append(abs(spectrum.gaussfit[1][2]) / step)
        else:
            widths.append(abs(fitTools.GaussFit(fctToFit=spectrum.amp, returnInfos=True)[1][2]))
    widths = np.array(widths, dtype=float)
    if not np.all(np.isfinite(widths)) or np.min(widths) <= 0:
        raise ValueError('fitted line widths must be > 0 (got {0})'.format(list(widths)))
    narrowest = int(np.argmin(widths))
    allNewSpectra = []
    for (i, spectrum) in enumerate(spectra):
        if i == narrowest:
            allNewSpectra.append(copy.deepcopy(spectrum))
        else:
            allNewSpectra.append(_resample_spectrum(spectrum, widths[narrowest] / widths[i], 'linear'))
    LineStacker._info('rebinned {0} spectra to the narrowest line width'.format(len(allNewSpectra)), verbose)
    return allNewSpectra


def _cube_header(info, nchan, crpix0, cdelt, crval=None):
    """FITS header for an (x, y, stokes, chan) array cut from a cube described
    by info, with a new spectral axis."""
    sub = info.wcs.sub([info.lng_axis + 1, info.lat_axis + 1])
    h = sub.to_header()
    for k in list(h.keys()):
        if k.startswith(('WCSAXES', 'MJD', 'DATE')):
            del h[k]
    h['CTYPE3'] = info.spec_ctype or 'FREQ'
    if info.spec_cunit:
        h['CUNIT3'] = info.spec_cunit
    h['CRPIX3'] = crpix0 + 1.
    h['CRVAL3'] = info.spec_crval if crval is None else crval
    h['CDELT3'] = cdelt
    for key in ('BMAJ', 'BMIN', 'BPA', 'BUNIT', 'RESTFRQ', 'RESTFREQ', 'SPECSYS', 'OBJECT'):
        if key in info.header:
            h[key] = info.header[key]
    return h


_IMAGE_SUFFIXES = ('.fits.gz', '.fit.gz', '.fts.gz', '.fits', '.fit', '.fts', '.image', '.im', '.img')


def _derived_name(imagename, suffix, outdir=None):
    """Name of a cube derived from imagename: suffix added before the
    extension ('cube.fits.gz' -> 'cube<suffix>.fits'); CASA images stay CASA
    images."""
    import os
    from LineStacker import _io
    base = os.fspath(imagename).rstrip('/')
    ext = '.image' if _io.is_casa_image(base) else '.fits'
    stem = base
    for old in _IMAGE_SUFFIXES:
        if base.lower().endswith(old):
            stem = base[:-len(old)]
            break
    name = stem + suffix + ext
    if outdir is not None:
        name = os.path.join(outdir, os.path.basename(name))
    return name


def rebin_CubesSpectra(    coords,
                        imagenames=None,
                        regionSize=False,
                        widths=False,
                        outputName='_SpectralRebinned',
                        hdul='auto',
                        outdir=None,
                        verbose=None):

    """
        Rebin a list of image-cubes so that all lines have the width of the
        narrowest one (widths from Gaussian fits, or given).
        **The lines must be visible source by source.**
        A rebinned channel interpolated from a blank (NaN) channel is NaN. One coordinate per cube
        is enough (the first one of each cube is used).

        Parameters
        ----------
        coords
            A list of LineStacker.Coord (one per image).
        imagenames
            Name of images to rebin (FITS or CASA images)
        regionSize
            size (in pixels) of the square region whose summed spectrum is fitted (placed around the source as
            the stamps of LineStacker.image.stack)\n
            If set to False spectra will be extracted solely from the coord pixel
        widths
            widths of the lines\n
            If set to 'False' the spectra will be fitted with a gaussian to extract the width
        outputName
            Suffix of the rebinned cube, to add to the original cube name.
        outdir
            Folder for the rebinned cubes (created if needed); default: next to the originals.
        verbose
            print where the rebinned cubes are written, and a progress bar (default: LineStacker.verbose)

        Returns
        -------
        list of str
            the names of the rebinned cubes (same format as the inputs).
    """
    import os
    from LineStacker import _io
    imagenames = LineStacker._cube_list(coords, imagenames)
    if widths is not False:
        widths = np.asarray(widths, dtype=float).ravel()
        if len(widths) != len(imagenames):
            raise ValueError('{0} widths given for {1} cubes'.format(len(widths), len(imagenames)))
        if not np.all(np.isfinite(widths) & (widths > 0)):
            raise ValueError('line widths must be finite and > 0 (got {0})'.format(list(widths)))
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
    LineStacker._info('rebinning {0} cubes to the narrowest line width'.format(len(imagenames)), verbose)
    pix = LineStacker.getPixelCoords(coords, imagenames, hdul=hdul)
    first = {}
    for c in pix:
        first.setdefault(c.imagenumber, c)
    if widths is False:
        import LineStacker.tools.fit as fitTools
        allWidths = []
        for i, image in enumerate(imagenames):
            if i not in first:
                raise ValueError('no coordinate for image {0}'.format(image))
            cube, info = _io.read_cube(image, hdul)
            c = first[i]
            x, y = LineStacker.round_half_up(c.ra), LineStacker.round_half_up(c.dec)
            if regionSize:
                # regionSize x regionSize pixels, placed as the stacker places stamps
                size = int(regionSize)
                bx = int(np.floor(c.ra - size // 2 + 0.5))
                by = int(np.floor(c.dec - size // 2 + 0.5))
                if bx < 0 or by < 0 or bx + size > info.nx or by + size > info.ny:
                    raise ValueError('the {0}-pixel region around source {1} ({2}, {3}) is not inside {4}'.format(
                        regionSize, i, x, y, image))
                region = cube[bx:bx + size, by:by + size, 0, :]
                spectra = np.nansum(region, axis=(0, 1))
                spectra[~np.isfinite(region).any(axis=(0, 1))] = np.nan  # all-blank channels
            else:
                spectra = cube[x, y, 0, :]
            allWidths.append(abs(fitTools.GaussFit(fctToFit=spectra, returnInfos=True)[1][2]))
    else:
        allWidths = list(widths)
    if not np.all(np.isfinite(allWidths)) or np.min(allWidths) <= 0:
        raise ValueError('fitted line widths must be > 0 (got {0}); give widths='.format(allWidths))
    minWidth = np.min(allWidths)
    widthsRatios = minWidth / np.array(allWidths, dtype=float)
    rebinnedImageNames = []
    for (i, image) in enumerate(imagenames):
        LineStacker.tools.ProgressBar(i, len(imagenames), verbose)
        cube, info = _io.read_cube(image, hdul)
        nchan = cube.shape[3]
        newN = max(1, int(nchan * widthsRatios[i]))
        if newN != nchan:
            newArray = _congrid_nan(cube[:, :, 0, :], [cube.shape[0], cube.shape[1], newN])
            newArray = newArray.reshape(cube.shape[0], cube.shape[1], 1, newN)
        else:
            newArray = cube[:, :, :1, :]
        # congrid samples old channel j*nchan/newN for new channel j
        r = nchan / float(newN)
        header = _cube_header(info, newN, info.spec_crpix0 / r, info.spec_cdelt * r)
        outname = _derived_name(image, outputName, outdir=outdir)
        _io.write_cube(outname, newArray, header)
        rebinnedImageNames.append(outname)
    LineStacker._info('{0} rebinned cubes written: {1}'.format(
        len(rebinnedImageNames), _describe_paths(rebinnedImageNames)), verbose)
    return rebinnedImageNames


def _describe_paths(paths):
    """'<dir>/<first> ... <last>' summary of written files."""
    import os
    full = [os.path.abspath(p) for p in paths]
    if len(full) == 1:
        return full[0]
    folders = set(os.path.dirname(p) for p in full)
    if len(folders) == 1:
        return 'in {0} ({1} ... {2})'.format(folders.pop(), os.path.basename(full[0]), os.path.basename(full[-1]))
    return '{0} ... {1}'.format(full[0], full[-1])


def maximizeAmp(spectra):
    """Grade of a stacked spectrum: its maximum (blank channels ignored; 0 if
    every channel is blank)."""
    spectra = np.asarray(spectra, dtype=float)
    return float(np.nanmax(spectra)) if np.any(np.isfinite(spectra)) else 0.


def maximizeSNR(spectra):
    """Grade of a stacked spectrum: amplitude of a Gaussian fit over the
    standard deviation of the valid (finite) channels outside the line
    (2.35 sigma on each side). 0 if the fit fails, finds no line inside the
    spectrum, or leaves fewer than 2 valid channels outside it."""
    import LineStacker.tools.fit as fitTools
    spectra = np.asarray(spectra, dtype=float)
    try:
        fitos=fitTools.GaussFit(fctToFit=spectra, returnInfos=True)
    except RuntimeError:
        return 0
    amp, center, width=fitos[1]
    width = abs(width)                 # the fit can return a negative sigma
    # A Gaussian centred outside the spectrum (the fit locked onto an edge) is
    # not a line; its limits would also index the spectrum from the end.
    if not np.isfinite(center) or not (0 <= center <= len(spectra) - 1):
        return 0
    toSTD=[]
    leftLimSpectra=int(center-width*2.35)
    rightLimSpectra=int(center+width*2.35)
    if rightLimSpectra<len(spectra)-1:
        toSTD.append(spectra[rightLimSpectra:][:])
    if leftLimSpectra>0:
        toSTD.append(spectra[:leftLimSpectra][:])
    if len(toSTD)==2:
        toSTD=np.concatenate(toSTD)
    elif len(toSTD)==1:
        toSTD=toSTD[0]
    else:
        return 0
    noise = np.asarray(toSTD, dtype=float).ravel()
    noise = noise[np.isfinite(noise)]  # blank channels are not noise samples
    if noise.size < 2:
        return 0
    STDed = np.std(noise)
    if not np.isfinite(amp) or not np.isfinite(STDed) or amp == 0 or STDed <= 0:
        return 0
    return float(amp) / STDed


def maximizeOutflow(spectra):
    """Grade of a stacked spectrum for a second (e.g. outflow) component:
    amplitude x width**2 of the fainter component of a two-Gaussian fit
    (blank channels ignored; 0 if the fit fails)."""
    import LineStacker.tools.fit as fitTools
    try:
        doubleFitos=fitTools.DoubleGaussFit(fctToFit=spectra, returnInfos=True)[1]
    except RuntimeError:
        return 0
    if doubleFitos[0]>doubleFitos[1]:
        area=doubleFitos[1]*doubleFitos[5]**2
    else:
        area=doubleFitos[0]*doubleFitos[4]**2
    return area


def _subsample_worker(sl):
    ctx = _quiet_fits()
    try:
        return _subsample_chunk(sl)
    finally:
        ctx.__exit__(None, None, None)


def _subsample_chunk(sl):
    import LineStacker.OneD_Stacker as OneD
    members = parallel.get_shared('members')[sl]
    sizes = parallel.get_shared('sizes')[sl]
    maxTest = parallel.get_shared('maxTest')
    spectra = parallel.get_shared('spectra')
    al = parallel.get_shared('aligned')
    kwargs = parallel.get_shared('kwargs')
    out = np.zeros(len(sizes))
    if al is None:                       # regridFromZ: stack each subset from scratch
        for r in range(len(sizes)):
            sub = [spectra[k] for k in members[r, :sizes[r]]]
            out[r] = maxTest(OneD.Stack(sub, verbose=False, **kwargs)[0])
        return out
    method = kwargs.get('method', 'mean')
    stacks, _ = OneD.combine_indices(al, members, method)   # padded with the empty row
    for r in range(len(sizes)):
        out[r] = maxTest(OneD.crop(al, stacks[r], members[r, :sizes[r]]))
    return out


def subsample_OneD(    spectra,
                    nRandom=10000,
                    maxTest=maximizeSNR,
                    seed=None,
                    n_jobs=1,
                    backend=None,
                    minSize=3,
                    maxSize=None,
                    verbose=None,
                    **kwargs):
    """
        Stacks random subsets of the spectra and grades every stack with maxTest; each spectrum collects the grades of the stacks it was part of.

        Parameters
        ----------

        spectra
            A list of LineStacker.OneD_Stacker.Spectrum objects
        nRandom
            Number of iterations of the Monte-Carlo process
        maxTest
            function test to grade the sources build your own, or use existing: **maximizeAmp**, **maximizeSNR**, **maximizeOutflow**.
            Default is to maximize the signal to noise ratio. With backend='process' a custom function must be defined at module level (not a lambda) so it can be sent to the workers.
        minSize, maxSize
            range of the random subsample sizes (default 3 to all spectra)
        seed, n_jobs, backend
            see the module description. Gaussian-fit based grades (maximizeSNR, maximizeOutflow) benefit from n_jobs>1.
        verbose
            print what is done and a progress bar (default: LineStacker.verbose)
        **kwargs
            passed to LineStacker.OneD_Stacker.Stack (method, center, chansStack...)

        Returns
        -------
        dict
            {spectrum name: summed grade} (spectra without a name are named by their index).
    """
    import LineStacker.OneD_Stacker as OneD
    spectra = _check_spectra(spectra)
    nRandom = _n_random(nRandom)
    n = len(spectra)
    for name, value in (('minSize', minSize), ('maxSize', maxSize)):
        if value is not None and (isinstance(value, bool) or int(value) != value):
            raise ValueError('{0} must be an integer (got {1!r})'.format(name, value))
    hi = n if maxSize is None else int(maxSize)
    if not (1 <= int(minSize) <= hi <= n):
        raise ValueError('subsample sizes must satisfy 1 <= minSize <= maxSize <= number of spectra '
                         '(got minSize={0}, maxSize={1}, {2} spectra)'.format(minSize, maxSize, n))
    LineStacker._info('subsampling {0} spectra, {1} iterations'.format(len(spectra), nRandom), verbose)
    spectraGrades={}
    for (i,spectrum) in enumerate(spectra):
        if spectrum.name=='':
            spectrum.name=str(i)
        spectraGrades[str(spectrum.name)]=0
    draws = _Draws(seed)
    sizes = draws.integers(int(minSize), hi + 1, size=nRandom)
    perm = np.argsort(draws.random((nRandom, n)), axis=1)
    members = np.full((nRandom, int(sizes.max()) if nRandom else 0), n)   # n = empty row
    for r in range(nRandom):
        members[r, :sizes[r]] = perm[r, :sizes[r]]

    if kwargs.get('regridFromZ', False):
        al = None
    else:
        OneD._set_centers(spectra, kwargs.get('center', 'lineCenterIndex'))
        al = OneD.Aligned(spectra, kwargs.get('chansStack', 'full'))
        # empty extra row used to pad subsets of different sizes
        al.X = np.vstack([al.X, np.full((1, al.L), np.nan)])
        al.X0 = np.vstack([al.X0, np.zeros((1, al.L))])
        al.W = np.vstack([al.W, np.zeros((1, al.L))])
        al.M = np.vstack([al.M, np.zeros((1, al.L), bool)])
        al.lo, al.hi = np.append(al.lo, 10**9), np.append(al.hi, -10**9)

    if backend in (None, 'process') and LineStacker.tools.resolve_n_jobs(n_jobs) > 1 and not _picklable(maxTest):
        warnings.warn('maxTest cannot be sent to worker processes (lambda or local function); using threads')
        backend = 'thread'
    shared = {'members': members, 'sizes': sizes, 'maxTest': maxTest, 'spectra': spectra,
              'aligned': al, 'kwargs': kwargs}
    tasks = _chunks(nRandom, n * (al.L if al is not None else 100), n_jobs, budget=2e6)
    results = parallel.map_chunks(_subsample_worker, tasks, n_jobs=n_jobs, backend=backend, shared=shared,
                                  progress=LineStacker._verbose(verbose))
    grades = np.concatenate(results) if results else np.zeros(0)
    for r in range(nRandom):
        for k in members[r, :sizes[r]]:
            spectraGrades[str(spectra[k].name)] += grades[r]
    return spectraGrades


def noise_estimator(    imagenames,
                        nRandom=10000,
                        maskradius=0,
                        maskCenter=None,
                        continuum=True,
                        chanToNoise='random',
                        hdul='auto',
                        plotIt=True,
                        seed=None,
                        verbose=None):
    """
        Estimates the noise of images or cubes from randomly drawn pixels
        (blank, i.e. NaN or infinite, pixels are never drawn; pixels that are
        exactly zero are data).

        Parameters
        ---------
        imagenames
            Path to image (FITS or CASA image), or list of paths
        nRandom
            Number of iterations
        maskradius
            Radius of the mask (in pixels)
        maskCenter
            Center location of the mask (in pixels, [x, y]; default [0, 0]). NB: CENTER OF IMAGE IS DEFINED AS [0,0]
        continuum
            If set to True points will be taken at random on the any of the spectral channels. Otherwise all channels are probed independently.
        chanToNoise
            If set to a specific (int) value then noise will be extracted from the specified channel. NB: only works if continuum is set to True.
        plotIt
            show the (masked) first image
        seed
            seed for reproducible draws (None uses numpy's global random state)
        verbose
            print what is done and a progress bar (default: LineStacker.verbose)

        Returns
        -------
        float or numpy array
            the standard deviation of the drawn pixels (continuum=True), or
            one per channel.
    """
    from LineStacker import _io
    import os
    maskCenter = [0, 0] if maskCenter is None else list(maskCenter)
    nRandom = _n_random(nRandom, minimum=1)
    draws = _Draws(seed)
    if isinstance(imagenames, (str, os.PathLike)) or not isinstance(imagenames, (list, tuple)):
        imagenames=[imagenames]
    imagenames = [os.fspath(n) for n in imagenames]
    if not imagenames:
        raise ValueError('no image given')

    if not continuum and len(imagenames) > 1:
        nchans = [_io.n_channels(name, hdul) for name in imagenames]
        if len(set(nchans)) > 1:
            raise ValueError('continuum=False estimates the noise channel by channel, so all images need '
                             'the same number of channels (found {0})'.format(sorted(set(nchans))))

    toSTD={}
    for (imNumber, imagename) in enumerate(imagenames):
        cube, info = _io.read_cube(imagename, hdul)
        allData = np.array(cube[:, :, :1, :], dtype=float)

        LineStacker.tools.ProgressBar(imNumber, len(imagenames), verbose)
        if maskradius!=0:
            X = np.arange(0, allData.shape[0])-int(allData.shape[0]/2.)
            Y = np.arange(0, allData.shape[1])-int(allData.shape[1]/2.)
            X,Y = np.meshgrid(X,Y, indexing='ij')
            keep = ((X-maskCenter[0])**2+(Y-maskCenter[1])**2) > maskradius**2
            allData[~keep] = np.nan
            if not np.isfinite(allData).any():
                raise Exception('no pixel left outside the mask in {0}: mask radius bigger than the image?'.format(imagename))
        if imagename==imagenames[0] and plotIt:
            import matplotlib.pyplot as plt
            plt.figure()
            plt.title(str(imagename))
            plt.imshow(np.nansum(allData, axis=(2,3)).T, origin='lower')
            plt.show(block=False)
        if continuum:
            if chanToNoise=='random':
                values = allData[:, :, 0, :]
            else:
                values = allData[:, :, 0, int(chanToNoise)]
            values = values[np.isfinite(values)]
            if values.size == 0:
                raise Exception('no valid (non-blank) pixel to draw from in {0}'.format(imagename))
            toSTD[imagename] = values[draws.integers(0, values.size, size=nRandom)]
        else:
            nchan = allData.shape[3]
            toSTD[imagename] = np.zeros((nRandom, nchan))
            for f in range(nchan):
                values = allData[:, :, 0, f]
                values = values[np.isfinite(values)]
                if values.size == 0:
                    raise Exception('channel {0} of {1} is blank'.format(f, imagename))
                toSTD[imagename][:, f] = values[draws.integers(0, values.size, size=nRandom)]

    if continuum:
        return np.std(np.concatenate([toSTD[u] for u in toSTD]))
    return np.std(np.concatenate([toSTD[u] for u in toSTD], axis=0), axis=0)


def _estimator_worker(sl):
    import LineStacker.image as li
    pos = parallel.get_shared('positions')[sl]          # (k, nsrc, 2)
    prep = parallel.get_shared('prep')
    W = parallel.get_shared('weights')
    method, save = parallel.get_shared('method'), parallel.get_shared('save')
    cubes = parallel.get_shared('cubes') or li._cubes(prep, parallel.get_shared('hdul'))
    s = prep.stampsize
    centre = int(s / 2)
    k, nchan = len(pos), prep.valid.shape[1]
    inside = prep.valid           # which channels exist does not depend on position

    def stamps_of(i):
        c = prep.coords[i]
        return li._extract_many(cubes[c.imagenumber], prep.infos[i], pos[:, i, 0], pos[:, i, 1],
                                prep.centres[i], prep.steps[i], s, nchan) * prep.scales[i]
    if method == 'mean':
        num = np.zeros((k, s, s, 1, nchan))
        den = np.zeros((k, s, s, 1, nchan))
        for i in range(len(prep.coords)):
            st = stamps_of(i)
            ok = np.isfinite(st)
            num += np.where(ok, st, 0.) * W[i]
            den += ok * W[i]
        with np.errstate(invalid='ignore', divide='ignore'):
            stacks = np.where(den > 0, num / den, np.nan)
    else:
        allst = np.stack([stamps_of(i) for i in range(len(prep.coords))], axis=1)
        allst = np.where(inside[None, :, None, None, None, :], allst, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            stacks = np.nanmedian(allst, axis=1)
    if save == 'central_pixel':
        return list(stacks[:, centre, centre, 0, :])
    return list(stacks)


def stack_estimator(    coords,
                        nRandom=100,
                        imagenames=None,
                        stampsize=1,
                        method='mean',
                        N_chans=30,
                        lowerLimit='default',
                        upperLimit='default',
                        save='central_pixel',
                        fEm=0,
                        spectralMethod='z',
                        lEm=-1,
                        hdul='auto',
                        seed=None,
                        n_jobs=1,
                        backend='thread',
                        verbose=None,
                        **kwargs):
    """
        Performs stacks at random positions a set number of times. Allows to probe the relevance of stack through stacking random positions as a Monte Carlo process.
        The random positions keep each source's redshift; positions whose stamp
        would fall outside the image are drawn again.

        Parameters
        ----------
        coords
            A list of LineStacker.Coord of all target coordinates.
        nRandom
            Number of iterations
        imagenames
            Name of imagenames to stack
        stampsize
            Size of the stamp to stack (because only the central pixel is extracted anyway, this should be kept to default to enhance efficiency)
        method
            Method for stacking, see LineStacker.image.stack
        N_chans
            Number of channels in the stack
        lowerLimit
            Lower spatial limit (distance from stacking position, radians) to randomize new stacking position.
            Default is 5 beams
        upperLimit
            Upper spatial limit (distance from stacking position, radians) to randomize new stacking position.
            Default is 10 beams
        save
            either 'central_pixel' or 'all'. If 'central_pixel' only the spectrum of central pixel is saved.
        spectralMethod
            Method to select the central frequency of the stack. See LineStacker.image.stack.
        seed, n_jobs, backend
            see the module description ('thread' is the default backend here: the cubes stay shared)
        verbose
            print what is done and a progress bar (default: LineStacker.verbose)
        **kwargs
            weighting, maskradius, regridFromZ... passed to LineStacker.image.stack's preparation

        Returns
        -------
        list
            [standard deviation, mean, all random stacks].
    """
    import LineStacker.image as li
    if fEm == 0:
        fEm = -1
    nRandom = _n_random(nRandom, minimum=1)
    if save not in ('central_pixel', 'all'):
        raise ValueError("save must be 'central_pixel' or 'all'")
    prep = li._prepare(coords, imagenames, stampsize=stampsize, N_chans=N_chans, spectralMethod=spectralMethod,
                       fEm=fEm, lEm=lEm, hdul=hdul, method=method, **kwargs)
    W = prep.weights()
    if lowerLimit == 'default' or upperLimit == 'default':
        bmaj = prep.infos[0].beam_major_deg
        if bmaj is None:
            raise ValueError('no beam (BMAJ) in {0}: give lowerLimit and upperLimit in radians'.format(prep.infos[0].path))
        beam = np.deg2rad(bmaj)
        lowerLimit = beam * 5. if lowerLimit == 'default' else lowerLimit
        upperLimit = beam * 10. if upperLimit == 'default' else upperLimit

    draws = _Draws(seed)
    n = len(prep.coords)
    half = int(stampsize / 2.)
    positions = np.zeros((nRandom, n, 2))
    for i, c in enumerate(prep.coords):
        inf = prep.infos[i]
        scale = np.deg2rad(inf.pixscale_deg)
        todo = np.arange(nRandom)
        for attempt in range(1000):
            dr = draws.uniform(lowerLimit, upperLimit, size=todo.size) / scale
            phi = draws.uniform(0, 2 * np.pi, size=todo.size)
            x = c.ra + dr * np.cos(phi)
            y = c.dec + dr * np.sin(phi)
            bx = np.floor(x - half + 0.5)
            by = np.floor(y - half + 0.5)
            ok = (bx >= 0) & (by >= 0) & (bx + stampsize <= inf.nx) & (by + stampsize <= inf.ny)
            positions[todo[ok], i, 0] = x[ok]
            positions[todo[ok], i, 1] = y[ok]
            todo = todo[~ok]
            if todo.size == 0:
                break
        else:
            raise Exception('no random position {0:.3g}-{1:.3g} rad ({2:.1f}-{3:.1f} pixels) away from source {4} '
                            'keeps a {5}-pixel stamp inside {6} ({7}x{8} pixels). Use smaller lowerLimit/upperLimit '
                            '(radians) or larger images.'.format(lowerLimit, upperLimit, lowerLimit / scale,
                                                                 upperLimit / scale, i, stampsize, inf.path,
                                                                 inf.nx, inf.ny))

    shared = {'positions': positions, 'prep': prep, 'weights': W, 'method': method, 'save': save, 'hdul': hdul,
              # threads share memory: hand them the cubes (processes read them)
              'cubes': li._cubes(prep, hdul) if (backend == 'thread' or
                                                 LineStacker.tools.resolve_n_jobs(n_jobs) == 1) else None}
    cost = prep.data.size
    LineStacker._info('stacking {0} sets of random positions around {1} sources'.format(nRandom, n), verbose)
    results = parallel.map_chunks(_estimator_worker, _chunks(nRandom, cost, n_jobs, budget=2e7),
                                  n_jobs=n_jobs, backend=backend, shared=shared,
                                  progress=LineStacker._verbose(verbose))
    dist = [d for r in results for d in r]
    # channels without data are NaN in every random stack: left out
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return [np.nanstd(dist), np.nanmean(dist), dist]


def randomizeCoords(    coords,
                        beam=0,
                        lowerLimit='default',
                        upperLimit='default',
                        imagenames=None,
                        hdul='auto',
                        seed=None):
    """
        New positions at a random distance (uniform between lowerLimit and upperLimit) and in a random direction from each coordinate.

        Parameters
        ---------
        coords
            A list of LineStacker.Coord coordinates
        lowerLimit
            Lower limit of the distance (radians) to the original position.\n
            Default is 5 beams
        upperLimit
            Upper limit of the distance (radians) to the original position.\n
            Default is 10 beams
        beam
            beam size (radians), needed if lowerLimit or upperLimit are set to 'default'
        imagenames
            needed for pixel coordinates (to convert distances to pixels)
    """
    if (lowerLimit == 'default' or upperLimit == 'default') and not beam > 0:
        raise ValueError('give beam (radians, > 0) for the default limits, or lowerLimit and upperLimit')
    if lowerLimit=='default':
        lowerLimit=beam*5.
    if upperLimit=='default':
        upperLimit=beam*10.
    if not (0 <= lowerLimit <= upperLimit) or not np.isfinite(upperLimit):
        raise ValueError('the limits must satisfy 0 <= lowerLimit <= upperLimit (radians), got {0}, {1}'.format(
            lowerLimit, upperLimit))
    imagenames = LineStacker._image_list(imagenames)
    rng = np.random.default_rng(seed)
    out = []
    for coord in coords:
        dr = rng.uniform(lowerLimit, upperLimit)
        pa = rng.uniform(0, 2 * np.pi)
        new = coord.copy()
        if coord.type == 'pixel':
            from LineStacker import _io
            if imagenames is None:
                raise ValueError('imagenames is needed to randomize pixel coordinates')
            info = _io.read_info(imagenames[LineStacker._image_index(coord, imagenames)], hdul)
            step = dr / np.deg2rad(info.pixscale_deg)
            new.ra, new.dec = coord.ra + step * np.cos(pa), coord.dec + step * np.sin(pa)
        else:
            import astropy.units as u
            sky = LineStacker._to_skycoord([coord])[0].directional_offset_by(pa * u.rad, dr * u.rad)
            new.ra, new.dec = sky.ra, sky.dec
        out.append(new)
    return out


def RebinToRest(    coord,
                    stamp,
                    z,
                    zToScale=0):
    """Resample the channels of an (x, y, stokes, chan) stamp of a source at
    redshift coord.z so that its channels are (1+z)/(1+zToScale) times
    narrower: new channel j samples old channel j*(1+zToScale)/(1+z), with a
    cubic spline through the valid channels (blank channels are never used;
    new channels next to a blank one, or beyond the last channel, are NaN).
    Returns (new stamp, changed)."""
    zRatio = (1. + coord.z) / (1. + zToScale)
    if zRatio == 1:
        return stamp, False
    newN = max(2, LineStacker.round_half_up(stamp.shape[-1] * zRatio))
    newStamp = _spline_axis(np.asarray(stamp[:, :, 0, :], dtype=float), 2, np.arange(newN) / zRatio)
    return newStamp[:, :, None, :], True


def regridFromZ(        coords,
                        imagenames=None,
                        stampsize=32,
                        N_chans=0,
                        fEm=0,
                        writeImage=True,
                        regridMethod='scaleToMin',
                        hdul='auto',
                        outdir='./regridedImages',
                        verbose=None):

    """
        Cube regridding function correcting the observed linewidth change due to redshift.
        Spectral dimensions are rebinned to a new_size=old_size*zRatio where zRatio=(1+z)/(1+zRef)
        where z is the observed line redshift and zRef is either the smaller (**scaleToMin**) or bigger (**scaleToMax**) redshift among the studied lines.
        Spectra are resampled with a cubic spline through the valid channels
        (separately between blank channels, which are never used): new
        channels between a valid and a blank channel, or beyond the last
        channel, are NaN.

        LineStacker.image.stack(..., regridFromZ=True) corrects the line widths in memory instead, without writing
        cubes (with linear interpolation).

        Parameters
        ----------
        coords
            A list of LineStacker.Coord of all target coordinates.
        imagenames
            Name of images to extract flux from.
        stampsize
            size of target image in pixels
        N_chans
            unused, kept for compatibility
        fEm
            rest emission frequency of the line (frequency cubes only)
        writeImage
            boolean, if set to **True** regridded images will be written to disk (in **outdir**)
        regridMethod
            Method of regridding, can be set either to **'scaleToMin'** or **'scaleToMax'**.
            In the first case all spectra are regridded to match the smallest redshift (over-gridding), all spectra are regridded to match the highest redshift in the other case (under-gridding).
        outdir
            folder for the regridded cubes (default ./regridedImages; emptied first)
        verbose
            print where the regridded cubes are written, and a progress bar (default: LineStacker.verbose)

        Returns
        -------
        tuple or list
            If writeImage is **True**: (imagenames, coords), the regridded cubes and their coordinates. Stack them as they are,
            with LineStacker.image.stack(coords, imagenames=imagenames, spectralMethod='z', fEm=fEm): a regridded cube's
            spectral axis is labelled in the rest frame (its coordinate has z = 0, channels of rest-frame width), and sources
            already at the reference redshift keep their original cube and redshift.
            Regridded cubes are stamps of stampsize pixels around each source, with the sky coordinates of the original cube.
            If writeImage is **False**: the list of regridded stamps.
    """
    import os
    import shutil
    from LineStacker import _io

    imagenames = LineStacker._cube_list(coords, imagenames)
    pix = LineStacker.getPixelCoords(coords, imagenames, hdul=hdul)
    LineStacker._info('regridding {0} sources using their redshifts'.format(len(pix)), verbose)

    zs = [float(cc.spec) for cc in pix]
    if regridMethod=='scaleToMin':
        zToScale=min(zs)
    elif regridMethod=='scaleToMax':
        zToScale=max(zs)
    else:
        raise ValueError("regridMethod must be 'scaleToMin' or 'scaleToMax'")
    if writeImage and fEm in (0, -1, None):
        raise ValueError('regridFromZ needs fEm (rest frequency of the line) to label the regridded cubes')

    if writeImage:
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)
        os.makedirs(outdir)
        regridedImagenames=[]
        regridedCoords=[]
    else:
        returnStamp=([0 for i in range(len(pix))])

    s = stampsize
    for (coordNumber, coord) in enumerate(pix):
        LineStacker.tools.ProgressBar(coordNumber, len(pix), verbose)
        imagename = imagenames[coord.imagenumber]
        cube, info = _io.read_cube(imagename, hdul)
        if info.spectral_kind != 'frequency':
            raise ValueError('regridFromZ writes frequency cubes only; use LineStacker.image.stack(regridFromZ=True) '
                             'for wavelength cubes')
        bx = int(np.floor(coord.ra - int(s / 2.) + 0.5))
        by = int(np.floor(coord.dec - int(s / 2.) + 0.5))
        if bx < 0 or by < 0 or bx + s > info.nx or by + s > info.ny:
            raise Exception('source {0} is too close to the edge of {1} for stampsize {2}'.format(coordNumber, imagename, s))
        tempStamp = np.array(cube[bx:bx + s, by:by + s, :1, :], dtype=float)

        if not writeImage:
            returnStamp[coordNumber]=RebinToRest(coord, tempStamp, coord.z, zToScale=zToScale)[0]
            continue
        if coord.z == zToScale:
            regridedImagenames.append(imagename)
            tempCoord = coord.copy()
            tempCoord.image, tempCoord.imagenumber = imagename, coordNumber
            regridedCoords.append(tempCoord)
            continue
        newStamp, changed = RebinToRest(coord, tempStamp, coord.z, zToScale=zToScale)
        if not changed:
            regridedImagenames.append(imagename)
            tempCoord = coord.copy()
            tempCoord.image, tempCoord.imagenumber = imagename, coordNumber
            regridedCoords.append(tempCoord)
            continue
        # new channel j samples old channel j/zRatio (RebinToRest)
        zRatio = (1. + coord.z) / (1. + zToScale)
        obsStep = info.spec_cdelt / zRatio                   # observed width of a new channel
        lineChannel = info.channel_of(fEm / (1. + coord.z)) * zRatio   # the line, in new channels
        # The new spectral axis is labelled in the rest frame (z = 0): the line
        # at fEm, channels of rest-frame width obsStep * (1 + z).
        header = _cube_header(info, newStamp.shape[3], lineChannel, obsStep * (1. + coord.z), crval=fEm)
        # the stamp keeps the sky coordinates of the original cube
        header['CRPIX1'] = header['CRPIX1'] - bx
        header['CRPIX2'] = header['CRPIX2'] - by
        outImName = _derived_name(imagename, '_source-'+str(coordNumber)+'_ZRegrid', outdir=outdir)
        _io.write_cube(outImName, newStamp, header)
        regridedImagenames.append(outImName)
        # the source's own (fractional) pixel in the stamp
        regridedCoords.append(LineStacker.Coord(coord.ra - bx, coord.dec - by, spec=0, image=outImName,
                                                imagenumber=coordNumber, weight=coord.weight, type='pixel'))
    if writeImage:
        written = [n for n in regridedImagenames if n.startswith(outdir)]
        LineStacker._info('{0} regridded cubes written to {1}{2}'.format(
            len(written), os.path.abspath(outdir),
            '' if len(written) == len(regridedImagenames) else
            ' ({0} sources at the reference redshift keep their original cube)'.format(len(regridedImagenames) - len(written))),
            verbose)
        return regridedImagenames, regridedCoords
    return returnStamp


def regridFromZ1D(spectra, regridMethod='scaleToMin', verbose=None):

    """
        1D regridding function correcting the observed linewidth change due to redshift.
        Spectra are rebinned to a new_size=old_size*zRatio where zRatio=(1+z)/(1+zRef)
        where z is the observed line redshift and zRef is either the smaller (**scaleToMin**) or bigger (**scaleToMax**) redshift among the studied lines.
        Fluxes are resampled with a cubic spline through the valid channels
        (separately between blank channels, which are never used): new
        channels between a valid and a blank channel, or beyond the last
        channel, are NaN. Spectral axes, per-channel weights and line
        centres (centerIndex) follow the new channels.

        Parameters
        ----------
        spectra
            List of spectra (LineStacker.OneD_Stacker.Spectrum objects).
        regridMethod
            Method of regridding, can be set either to **'scaleToMin'** or **'scaleToMax'**.
            In the first case all spectra are regridded to match the smallest redshift (over-gridding), all spectra are regridded to match the highest redshift in the other case (under-gridding).
        verbose
            print what is done and a progress bar (default: LineStacker.verbose)

        Returns
        -------
        list of LineStacker.OneD_Stacker.Spectrum
            the regridded spectra (copies).
    """

    spectra = _check_spectra(spectra)
    zToScale = _z_reference(spectra, regridMethod)

    allNewSpectra=[]

    LineStacker._info('regridding {0} spectra using their redshifts'.format(len(spectra)), verbose)
    for (i, spectrum) in enumerate(spectra):
        LineStacker.tools.ProgressBar(i, len(spectra), verbose)
        if spectrum.z == zToScale:
            allNewSpectra.append(copy.deepcopy(spectrum))
        else:
            zRatio = (1. + spectrum.z) / (1. + zToScale)
            allNewSpectra.append(_resample_spectrum(spectrum, zRatio, 'spline'))
    return allNewSpectra

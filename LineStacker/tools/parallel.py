"""
Small helper to spread Monte-Carlo iterations (bootstraps, subsamples, random
stacks) over several CPU cores.

Every function in LineStacker that accepts ``n_jobs`` draws all its random
numbers up front, in the calling process, and only distributes the
deterministic work. Results are therefore identical whatever ``n_jobs`` or
``backend`` is, and reproducible with ``seed``.

``backend='process'`` (default) runs workers in separate processes and is the
right choice for Python-heavy work such as Gaussian fits. ``backend='thread'``
avoids starting processes and copying data; it helps when the work is mostly
large numpy operations (e.g. median stacks of big cubes).

When using processes from a script on macOS or Windows, protect the script's
entry point::

    if __name__ == '__main__':
        main()

otherwise each worker re-runs the script while starting.
"""

import contextvars
import threading

import numpy as np

from LineStacker.tools import resolve_n_jobs

_DEFAULT_BACKEND = 'process'
_BLAS_THREAD_VARS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                     'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')

# Large read-only inputs of the running map_chunks call. A context variable,
# so that map_chunks calls running at the same time in different threads
# each see their own; worker processes set it once (_init_worker), so the
# inputs are not pickled with every task.
_shared = contextvars.ContextVar('linestacker_shared', default=None)
# Process jobs set BLAS thread variables while their workers start: one at a
# time, so concurrent jobs do not restore each other's values.
_env_lock = threading.Lock()


def set_default_backend(backend):
    """Set the backend used when a function is called with backend=None
    ('process' or 'thread')."""
    global _DEFAULT_BACKEND
    if backend not in ('process', 'thread'):
        raise ValueError("backend must be 'process' or 'thread'")
    _DEFAULT_BACKEND = backend


def _init_worker(shared):
    _shared.set(dict(shared))


def get_shared(name):
    """Access, from a worker function, an object passed through ``shared``."""
    current = _shared.get()
    if current is None:
        raise RuntimeError('get_shared() is only available inside a function run by map_chunks')
    return current[name]


def chunk_slices(n, nchunks):
    """Split range(n) into at most nchunks contiguous slices."""
    nchunks = max(1, min(n, nchunks))
    edges = np.linspace(0, n, nchunks + 1).astype(int)
    return [slice(edges[i], edges[i + 1]) for i in range(nchunks) if edges[i + 1] > edges[i]]


def map_chunks(func, tasks, n_jobs=1, backend=None, shared=None, progress=True):
    """
        Apply ``func`` to every element of ``tasks`` and return the results in
        order. ``func`` must be a module-level function (so it can be sent to
        worker processes). Large inputs common to all tasks go in ``shared``
        (a dict); workers read them with :func:`get_shared`.

        Parameters
        ---------
        func
            function taking one task and returning its result
        tasks
            list of tasks (usually slices of iterations)
        n_jobs
            number of workers; 1 runs serially, -1 uses all cores
        backend
            'process' or 'thread' (None uses the default, see set_default_backend)
        shared
            dict of read-only objects made available to func via get_shared
        progress
            show LineStacker's progress bar
    """
    from LineStacker.tools import ProgressBar
    tasks = list(tasks)
    n_jobs = resolve_n_jobs(n_jobs)
    shared = shared or {}
    backend = backend or _DEFAULT_BACKEND
    if n_jobs == 1 or len(tasks) <= 1:
        token = _shared.set(dict(shared))
        try:
            out = []
            for i, t in enumerate(tasks):
                out.append(func(t))
                if progress:
                    ProgressBar(i, len(tasks))
            return out
        finally:
            _shared.reset(token)

    if backend == 'process':
        with _env_lock:
            return _run_pool(func, tasks, n_jobs, backend, shared, progress)
    return _run_pool(func, tasks, n_jobs, backend, shared, progress)


def _run_pool(func, tasks, n_jobs, backend, shared, progress):
    from LineStacker.tools import ProgressBar
    import concurrent.futures as cf
    import os
    saved_env = {}
    if backend == 'process':
        # Each worker would otherwise start one BLAS thread per core: share
        # the cores out instead (workers inherit the environment when spawned).
        per_worker = str(max(1, (os.cpu_count() or 1) // n_jobs))
        for var in _BLAS_THREAD_VARS:
            saved_env[var] = os.environ.get(var)
            os.environ[var] = per_worker
    token = None
    if backend == 'thread':
        # each task runs in a copy of this context, which holds `shared`
        token = _shared.set(dict(shared))
        executor = cf.ThreadPoolExecutor(max_workers=n_jobs)
    elif backend == 'process':
        import multiprocessing as mp
        # 'spawn' everywhere: forking a process that holds CASA's C++ tools
        # (and macOS system libraries) is unsafe.
        executor = cf.ProcessPoolExecutor(max_workers=n_jobs, mp_context=mp.get_context('spawn'),
                                          initializer=_init_worker, initargs=(shared,))
    else:
        raise ValueError("backend must be 'process' or 'thread'")
    try:
        with executor:
            if backend == 'thread':
                futures = [executor.submit(contextvars.copy_context().run, func, t) for t in tasks]
            else:
                futures = [executor.submit(func, t) for t in tasks]
            done = 0
            for _ in cf.as_completed(futures):
                if progress:
                    ProgressBar(done, len(futures))
                done += 1
            return [f.result() for f in futures]
    finally:
        if token is not None:
            _shared.reset(token)
        for var, value in saved_env.items():
            if value is None:
                os.environ.pop(var, None)
            else:
                os.environ[var] = value

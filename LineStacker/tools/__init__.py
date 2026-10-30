def _bar(fraction, final=False):
    g = 10 if final else int(fraction * 10)
    done = ("\033[0;36;46m" + " " + '\x1b[0m') * g
    todo = '' if final else ' ' * (10 - g)
    s = '|{0}|'.format(done + todo)
    s += "\033[1;33m" + ' LOADING ' + '\x1b[0m'
    if final:
        s += '\033[1;31m' + ' 100%' + '\x1b[0m' + '\n\n'
    else:
        s += "\033[0;32m" + str(int(fraction * 100)) + '% ' + '\x1b[0m'
    return '\r' + s


def ProgressBar(n, total, verbose=None):
    """
        Draws a progress bar (at most about 100 redraws).

        Parameters
        ----------
        n
            current step
        total
            total number of steps
        verbose
            draw the bar (default: LineStacker.verbose)
    """
    import LineStacker
    if total <= 0 or not LineStacker._verbose(verbose):
        return
    if n == total - 1:
        print(_bar(1., final=True), end="", flush=True)
        return
    # Redraw at most ~100 times, whatever the number of steps.
    step = max(1, total // 100)
    if n % step == 0:
        print(_bar(float(n) / total), end="", flush=True)


def resolve_n_jobs(n_jobs):
    """Number of workers for an n_jobs argument (None/1 = serial, -1 = all cores)."""
    import os
    if n_jobs is None:
        return 1
    n = int(n_jobs)
    ncpu = os.cpu_count() or 1
    if n < 0:
        n = max(1, ncpu + 1 + n)
    return max(1, n)

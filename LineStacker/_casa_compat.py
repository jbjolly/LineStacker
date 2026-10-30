"""
CASA 6 tools, used only to read and write CASA images.

LineStacker imports casatools the first time it opens a CASA image (FITS
cubes never need CASA). :func:`new_tool` returns a fresh, private instance of
a tool, so concurrent reads (threads) never share one.

CASA's log: when LineStacker is the one starting CASA (casatools not yet
imported, i.e. not inside a CASA session), CASA's log goes to one file in the
temporary directory instead of a new ``casa-<date>.log`` in the working
directory each time, and CASA's messages go to that file rather than the
terminal; see ``LineStacker.casa_logfile``. A ``logfile`` set in
``~/.casa/config.py`` is respected, and so is everything when CASA was
started before (CASA prompt, or ``import casatools`` / ``casatasks`` first).
"""

import os
import re
import sys
import tempfile

_DEFAULT_LOG = re.compile(r'^casa-\d{8}-\d{6}\.log$')


def default_logfile():
    """CASA log file used when LineStacker starts CASA: one file per user in
    the temporary directory (appended to)."""
    try:
        import getpass
        user = getpass.getuser()
    except Exception:
        user = str(os.getuid()) if hasattr(os, 'getuid') else 'user'
    user = re.sub(r'[^A-Za-z0-9_.-]', '_', user)
    return os.path.join(tempfile.gettempdir(), 'linestacker-casa-{0}.log'.format(user))


def _configure_log():
    """Point CASA's log where LineStacker.casa_logfile says, before casatools
    is imported. Returns True if CASA's messages should then go to that log
    only (not the terminal); False when CASA is already running or
    casa_logfile is 'casa' (CASA's settings are kept)."""
    if 'casatools' in sys.modules:
        return False
    import LineStacker
    setting = getattr(LineStacker, 'casa_logfile', 'auto')
    if setting is None:
        setting = os.devnull
    if setting == 'casa':
        return False
    try:
        from casaconfig import config
    except ImportError:
        return False
    if setting == 'auto':
        current = getattr(config, 'logfile', None)
        casa_default = (isinstance(current, str) and _DEFAULT_LOG.match(os.path.basename(current))
                        and os.path.dirname(os.path.abspath(current)) == os.getcwd())
        if getattr(config, 'nologfile', False):
            config.logfile = os.devnull
        elif current is None or casa_default:
            config.logfile = default_logfile()
        # else: logfile set in ~/.casa/config.py, keep it
        return True
    config.logfile = os.fspath(setting)
    return True


def _casatools():
    quiet = _configure_log()
    try:
        import casatools
    except ImportError as exc:
        raise ImportError(
            "Reading or writing CASA images needs CASA 6's casatools "
            "(`pip install casatools`, see https://casadocs.readthedocs.io/), "
            "or convert the cubes to FITS: FITS cubes do not need CASA.") from exc
    if quiet:
        # casatools' own log sink (writing to config.logfile) becomes the
        # global one: tool messages go to the log file, not to the terminal,
        # as in a CASA session.
        sink = getattr(casatools, 'logger', None)
        if sink is not None:
            try:
                sink.setglobal(True)
            except Exception:
                pass
    return casatools


def new_tool(name):
    """A fresh, private instance of a CASA tool ('image', 'regionmanager',
    'quanta', ...). Close it with .done() when finished."""
    return getattr(_casatools(), name)()


__all__ = ['new_tool', 'default_logfile']

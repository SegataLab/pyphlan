import bz2
import sys
import gzip


def openr(fn, mode="r"):
    if fn is None:
        return sys.stdin

    if fn.endswith(".bz2"):
        if sys.version_info[0] < 3:
            return bz2.BZ2File(fn)
        else:
            return bz2.open(fn, 'rt')
    elif fn.endswith(".gz"):
        if sys.version_info[0] < 3:
            return None  # need to check if gzip is different in Python2
        else:
            return gzip.open(fn, 'rt')
    else:
        return open(fn, mode)


def openw(fn):
    if fn is None:
        return sys.stdout

    if fn.endswith(".bz2"):
        if sys.version_info[0] < 3:
            return bz2.BZ2File(fn, "w")
        else:
            return bz2.open(fn, 'wt')
    elif fn.endswith(".gz"):
        if sys.version_info[0] < 3:
            return None  # need to check if gzip is different in Python2
        else:
            return gzip.open(fn, 'wt')
    else:
        return open(fn, "w")


def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

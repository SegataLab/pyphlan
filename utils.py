import bz2

def openr( fn ):
    return bz2.BZ2File(fn) if fn.endswith(".bz2") else open(fn)

def openw( fn ):
    return bz2.BZ2File(fn,"w") if fn.endswith(".bz2") else open(fn,"w")


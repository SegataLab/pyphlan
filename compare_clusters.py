#!/usr/bin/env python

import sys
import collections
import utils
import numpy as np

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert core gene files to core gene summaries\n')

    p.add_argument( '-1', default=None, type=str,
            help=   "Cluster file 1")
    p.add_argument('-2', default=None, type=str,
            help=   "Cluster file 2")

    return vars( p.parse_args() )

def stats( cl ):
    lens = [len(c) for c in cl]
    avg, stdev = np.mean(lens), np.std(lens)
    singletons = len([1 for l in lens if l == 1])
    return avg, stdev, len(lens), singletons, lens

if __name__ == "__main__":
    args = read_params( sys.argv )

    with utils.openr(args['1']) as inp:
        cl1 = [set(l.strip().split('\t')) for l in inp]
    with utils.openr(args['2']) as inp:
        cl2 = [set(l.strip().split('\t')) for l in inp]

    print stats(cl1)
    print stats(cl2)




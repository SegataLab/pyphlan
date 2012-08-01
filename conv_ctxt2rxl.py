#!/usr/bin/env python

import sys
import collections

try:
    import argparse as ap
    import bz2 
    import random
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert Usearch ".uc" files in tab-delimited'
            ' files with the seed as first field followed by the other IDs\n')

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input uc file [stdin if not present]")
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output txt file compresse if fiven with bz2 extension\n"
                    "[stdout if not present]")
    p.add_argument('--subsample', metavar="Subsampling rate",
            default=1.0, type=float )
    p.add_argument('-n', metavar="Minimum number of matching taxa",
            default=0, type=int )

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    openr = bz2.BZ2File if args['ctxt'] and args['ctxt'].endswith(".bz2") else open
    valin = []
    with (openr(args['ctxt']) if args['ctxt'] else sys.stdin) as inp:
        for l in inp:
            tset = set([int(a) for a in l.strip().split('\t')])
            if len(tset) < args['n']:
                continue
            valin.append(tset)
        #valin = [set([int(a) for a in l.strip().split('\t')]) for l in inp]
    all_t = set()
    for v in valin:
        all_t |= v

    res = {}
    for t in all_t:
        #if len(t) < args['n']:
        #    continue
        res[t] = [int(t in v) for v in valin]

    openw = bz2.BZ2File if args['txt'].endswith(".bz2") else open
    with (openw(args['txt'],"w") if args['txt'] else sys.stdout) as out:
        n = len(res.values()[0])
        n_s = int(float(n)*args['subsample'])
        out.write( str(len(res))+" "+str(n_s)+"\n" )
        indok = set(random.sample( list(range(n)), n_s))

        for k,v in res.items():
            out.write( str(k)[1:]+" "*(9-len(str(k)[1:]))+"".join([str(s) for i,s in enumerate(v) if i in indok]) +"\n" )
#!/usr/bin/env python

import sys
import collections

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description="Merge two core txt files"
            "where the second are subclusters of the first" )

    p.add_argument( 'ctxt1', nargs='?', default=None, type=str,
            help=   "the first ctxt")
    p.add_argument('ctxt2', nargs='?', default=None, type=str,
            help=   "the second ctxt")
    p.add_argument('out_ctxt', nargs='?', default=None, type=str,
            help=   "the output txt file compressef if given with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    openr = bz2.BZ2File if args['ctxt1'].endswith(".bz2") else open
    ctxt1 = dict([v[0],v[1:] if len(v) > 1 else []) 
                for v in (l.strip().split('\t') for l in openr(args['ctxt1']))])
    openr = bz2.BZ2File if args['ctxt2'].endswith(".bz2") else open
    ctxt2 = dict([v[0],v[1:] if len(v) > 1 else []) 
                for v in (l.strip().split('\t') for l in openr(args['ctxt2']))])

    mtxt = dict()
    for k,v in ctxt1.items():
        nv = v
        for kk in [k]+v:
            if kk in ctxt2:
                nv += ctxt2[kk]
        mtxt[k] = nv

    openw = bz2.BZ2File if args['txt'].endswith(".bz2") else open
    with (openw(args['txt'],"w") if args['txt'] else sys.stdout) as out:
        for k,v in sorted(mtxt.items(),key=lambda x:-len(x[1])):
            out.write( "\t".join([k]+list(v)) +"\n" )

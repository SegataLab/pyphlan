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

    p.add_argument( 'ctxt1', default=None, type=str,
            help=   "the first ctxt")
    p.add_argument('--out_ctxt', default=None, type=str,
            help=   "the output txt file compressef if given with bz2 extension\n"
                    "[stdout if not present]")
    p.add_argument('ctxt2', nargs='*', default=None, type=str,
            help=   "the second ctxt")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    print args['ctxt1']
    print args['ctxt2']

    openr = bz2.BZ2File if args['ctxt1'].endswith(".bz2") else open
    mtxt = dict([([v[0],(v[1:] if len(v) > 1 else [])]) for v in (l.strip().split('\t') for l in openr(args['ctxt1']))])
    
    for f in args['ctxt2']:
        openr = bz2.BZ2File if f.endswith(".bz2") else open
        ctxt = dict([([v[0],(v[1:] if len(v) > 1 else [])]) for v in (l.strip().split('\t') for l in openr(f))])

        for k,v in mtxt.items():
            nv = v
            for kk in [k]+v:
                if kk in ctxt:
                    nv += ctxt[kk]
            mtxt[k] = nv

    openw = bz2.BZ2File if args['out_ctxt'].endswith(".bz2") else open
    with (openw(args['out_ctxt'],"w") if args['out_ctxt'] else sys.stdout) as out:
        for k,v in sorted(mtxt.items(),key=lambda x:-len(x[1])):
            out.write( "\t".join([k]+list(v)) +"\n" )

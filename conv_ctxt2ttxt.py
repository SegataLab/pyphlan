#!/usr/bin/env python

import sys
import collections
import utils
try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert core gene txt file'
            ' substituting gene IDs with genomes IDs\n')

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input ctxt file [stdin if not present]")
    p.add_argument('--g2t', metavar="Mapping file from genes to taxa",
            default=None, type=str )
    p.add_argument('--t2g', metavar="Mapping file from taxa to genes",
            default=None, type=str )
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output gtxt file, compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    #openr = bz2.BZ2File if args['ctxt'] and args['ctxt'].endswith(".bz2") else open
    #with (openr(args['ctxt']) if args['ctxt'] else sys.stdin) as inp:
    with utils.openr( args['ctxt'] ) as inp:
        valin = [[int(a) for a in l.strip().split('\t')] for l in inp]
    if not args['g2t'] and not args['t2g']:
        sys.stdout.write("Error one of --t2g and --g2t must be provided\n")
        sys.eit(0)
    g2t = {}


    if args['g2t']:
        with open( args['g2t'] ) as inp:
            g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
    elif args['t2g']:
        with open( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[int(g)] = int(ll[0])
    
    #valout = [[g2t[k]]+[g2t[vv] for vv in v] for k,v in valin]
    for j,v in enumerate(valin):
        #for i,vv in enumerate(v):
        #    vn = set(g2t[vv] for vv in v)
        valin[j] = set(g2t[vv] for vv in v)

    #openw = bz2.BZ2File if args['txt'] and args['txt'].endswith(".bz2") else open
    with utils.openw(args['txt']) as out:
        for v in sorted(valin,key=lambda x:-len(x)):
            out.write( "\t".join([str(s) for s in v]) +"\n" )

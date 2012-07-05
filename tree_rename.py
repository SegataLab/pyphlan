#!/usr/bin/env python

import sys

try:
    import argparse as ap
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

try:
    import phylophlan as ppa
except ImportError:
    sys.stderr.write( "phylophlan.py not found\n" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Extract a subtree from a tree')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree [stdin if not present]")
    p.add_argument('outtree', nargs='?', default=None, type=str,
            help=   "the out file [stdout if not present]")
    p.add_argument('-f', default=None, type=str,
            help=   "file containing the list of clades for lca and ltcs " ) 
    p.add_argument('-n', default=None, type=str,
            help=   "name of the new root " )

    st = ['root_name','lca','ltcs']
    p.add_argument('-s', choices=st, default='root_name', type=str,
            help=  "set the strategy the select the root of the new tree")

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    ctree = ppa.PpaTree( args['intree'] )
    ctree.rename( strategy = args['s'], n = args['n'], fn = args['f'] )
    ctree.export( args['outtree'] )


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
    p = ap.ArgumentParser(description='Extract sub tree from atree')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree [stdin if not present]")
    p.add_argument('outtree', nargs='?', default=None, type=str,
            help=   "the output PhyloXML tree [stdout if not present]")
    p.add_argument('-f', metavar="File containing clade names",
            default=None, type=str )
    p.add_argument( '-n', metavar="N leaves for longest_internal_edge_n, or N for n_anc", 
                    default=1, type=int)

    reroot_st = ['name','lca','ltcs','longest_edge','longest_internal_edge','longest_internal_edge_n',"n_anc"]
    #reroot_st = ['lca','ltcs','midpoint','longest_edge','longest_internal_edge','longest_internal_edge_n']
    p.add_argument('-s', choices=reroot_st, default='longest_internal_edge', type=str,
            help=  "select rerooting strategy")

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    tree = ppa.PpaTree( args['intree'] )
    tree.subtree( strategy = args['s'], fn = args['f'], n = args['n'] )
    tree.export( args['outtree'] )

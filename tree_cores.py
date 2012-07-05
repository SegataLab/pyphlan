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
    p = ap.ArgumentParser(
            description='Screen clusters of taxa for core genes')

    p.add_argument( 'intree', nargs='?', default=None, type=str,
            help=   "the input tree [stdin if not present]")
    p.add_argument('outfile', nargs='?', default=None, type=str,
            help=   "the output core file [stdout if not present]")
    p.add_argument('-f', metavar="File containing sets of taxa",
            default=None, type=str )

    return vars( p.parse_args() )


if __name__ == "__main__":
    args = read_params( sys.argv )
    tree = ppa.PpaTree( args['intree'] )
    cores = tree.find_cores(args['f'])

    with open( args['outfile'], "w" ) as outf:
        for k,v in sorted(cores.items(),key=lambda x:x[0]):
            for vv in v:
                outf.write( "\t".join( [str(s) for s in [k]+list(vv)]) +"\n" )

    #ctree.export_cores( args['outfile']  )
    """
    ctree.reroot( strategy = args['s'], tf = args['f'] )
    ctree.export( args['outtree'] )
    """

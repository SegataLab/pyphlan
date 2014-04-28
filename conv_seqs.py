#!/usr/bin/env python

import sys
import collections
import utils
try:
    import argparse as ap
    import bz2 
    import random
    from Bio import SeqIO 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='')

    p.add_argument( '-i', nargs='?', default=None, type=str,
            help=   "the input file [stdin if not present]")
    p.add_argument( '--if', nargs='?', default='fastq', type=str,
            help=   "the input format [default fastq]")
    p.add_argument('--of', nargs='?', default='fasta', type=str,
            help=   "the output format [default fastqq]")

    """
    p.add_argument('--subsample', metavar="Subsampling rate",
            default=1.0, type=float )
    p.add_argument('-n', metavar="Minimum number of matching taxa",
            default=0, type=int )
    p.add_argument('-p', metavar="Prefix for taxon names",
            default="", type=str )
    """
    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )


    with utils.openw( None ) as outf:
        for seq in SeqIO.parse( utils.openr( args['i'], mode = 'rb' ), args['if']):
            SeqIO.write(seq, outf, args['of']) 

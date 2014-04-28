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

    p.add_argument('-q', metavar="Subsampling rate",
            default=10, type=int )
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
        for seq in SeqIO.parse( utils.openr( None ),'fastq' ):
            minp, maxp = None, None
            for i,v in enumerate(seq.letter_annotations["phred_quality"]):
                if minp is None and v >= args['q']:
                    minp = i
            for i,v in enumerate(seq.letter_annotations["phred_quality"][::-1]):
                if maxp is None and v >= args['q']:
                    maxp = len(seq.letter_annotations["phred_quality"]) - i
            SeqIO.write(seq[minp:maxp], outf, 'fastq') 

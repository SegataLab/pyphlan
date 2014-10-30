#!/usr/bin/env python 
from __future__ import with_statement

import sys
import argparse
import os
import textwrap
from collections import namedtuple as nt
import random as rnd
rnd.seed(1982)
import utils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np

def read_params(args):
    parser = argparse.ArgumentParser(description='Display the length of each fasta entry')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output txt file [stdout if not present]")
    arg( '-q', action='store_true', 
         help="set this for fastq")
    arg( '-t','--total', action='store_true', help="Print only the sum of the length of all sequences\n")
    arg( '-s','--stat', action='store_true', help="Print only the statistics about the length of sequences\n")
    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    ltot = 0
    lvec= []
    with utils.openw( par['out_f'] ) as outf:
        for r in SeqIO.parse( utils.openr(par['inp_f']), "fastq" if par['q'] else "fasta"):
            l = len(r.seq)
            if par['total']:
                ltot += l
            elif par['stat']:
                lvec.append(l);
            else:
                outf.write( "\t".join( [r.id,str(l)] ) + "\n" ) 
        if par['total']:
            outf.write( str(ltot) + "\n" )
        if par['stat']:
            outf.write("number_of_reads\t" + str(len(lvec))  + "\n")
            outf.write("minimum_read_length\t" + str(np.min(lvec))  + "\n")
            outf.write("median_read_length\t" + str(np.median(lvec))  + "\n")
            outf.write("mean_read_length\t" + str(np.mean(lvec)) + "\n")
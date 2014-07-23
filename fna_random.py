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


bin2nuc = { "00":"A", "01":"T", "10":"C", "11":"T"}

def read_params(args):
    parser = argparse.ArgumentParser(description='Display the length of each fasta entry')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")
    arg( '-t','--total', action='store_true', help="Print only the sum of the length of all sequences\n") 
    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    ltot = 0
    with utils.openw( par['out_f'] ) as outf:
        dna = ""
        for r in utils.openr(par['inp_f']):
            for l in r.strip():
                sn = "{0:b}".format(ord(l))
                if len(sn) % 2:
                    sn += "0"
                for i,s in enumerate( sn ):
                    if not i % 2:
                        dna += bin2nuc[sn[i:i+2]]
        outf.write( dna + "\n" ) 

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
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def read_params(args):
    parser = argparse.ArgumentParser(description='Split/Select/Randomize/Subsample a multi fasta file')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")

    parser.add_argument('-s', metavar='Subsequene to look for', required = True, type = str )

    return vars(parser.parse_args())

if __name__ == '__main__':
    par = read_params(sys.argv)

    ss = par['s'].lower()
    ssr = Seq(par['s']).reverse_complement().lower()
    f = os.path.basename(par['inp_f']).split(".")[0]
    with utils.openw( par['out_f'] ) as outf:
        for r in SeqIO.parse( utils.openr(par['inp_f']), "fasta"):
            rl = r.seq.lower()
            if ss in rl or ssr in rl:
                outf.write( f + "\t" + str(r.id) + "\n" )
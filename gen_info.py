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
    arg( '-n', default=None, type=str ) 
    return vars(parser.parse_args())


def N50( lens ):
    ll = sorted(lens,reverse=True)
    val =sum(lens)*0.5
    suml, i = 0, 0
    n50 = 0
    while suml < val:
        suml += ll[i]
        n50 = ll[i]
        i += 1
    return n50

if __name__ == '__main__':
    par = read_params(sys.argv)

    ltot = 0
    lvec= []


    with utils.openw( par['out_f'] ) as outf:
        lens = []
        cn, gn, tn, an = 0, 0, 0, 0
        for r in SeqIO.parse( utils.openr(par['inp_f']), "fastq" if par['q'] else "fasta"):
            l = len(r.seq)
            lens.append(l)
            an += r.seq.count('a') + r.seq.count('A')
            tn += r.seq.count('t') + r.seq.count('T')
            gn += r.seq.count('g') + r.seq.count('G')
            cn += r.seq.count('c') + r.seq.count('C')
            
        
        info = ["total_length","n_contigs","N50","median_length","mean_length","CG_perc"]
        if par['n']:
            info = ['name'] + info

        outf.write( "#"+"\t".join(info) +"\n")

        res = {}
        res['name'] = par['n']
        res['total_length'] = sum( lens )
        res['n_contigs'] = len( lens )
        res['N50'] = N50( lens )
        res['median_length'] = np.median( lens )
        res['mean_length'] = round(np.mean( lens ), 0)
        res['CG_perc'] = round( 100.0 * float(cn+gn) / float(an+tn+gn+cn), 2)

        outf.write( "\t".join([str(res[v]) for v in info]) +"\n")


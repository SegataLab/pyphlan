#!/usr/bin/env python
from Bio import SeqIO
import argparse as ap
import sys

def read_params(args):
    p = ap.ArgumentParser(description = 'fastax_len_filter.py Parameters\n')
    p.add_argument('--min_len', required = True, default = None, type = int)
    p.add_argument('--max_len', required = False, default = None, type = int)
    p.add_argument('-t', required = False, default = 'fastq', type = str)
    return vars(p.parse_args()) 

if __name__ == '__main__':
    args = read_params(sys.argv)
    min_len = args['min_len']
    max_len = args['max_len']
    with sys.stdout as outf:
        for r in SeqIO.parse(sys.stdin, args['t']):
            lenr = len(r)
            if ( max_len is not None ) and ( lenr > max_len ):
                continue
            if lenr >= min_len:
                SeqIO.write(r, outf, args['t'])

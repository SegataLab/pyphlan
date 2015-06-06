#!/usr/bin/env python
from Bio import SeqIO
import argparse as ap
import sys

def read_params():
    p = ap.ArgumentParser(description = 'fastq2fasta.py Parameters\n')
    p.add_argument('--ifn', required = False, default = None, type = str)
    p.add_argument('--ofn', required = False, default = None, type = str)
    return vars(p.parse_args())

    
if __name__ == '__main__':
    args = read_params()
    if args['ifn'] == None:
        ifile = sys.stdin
    else:
        ifile = open(args['ifn'], 'r')
    if args['ofn'] == None:
        ofile = sys.stdout
    else:
        ofile = open(args['ofn'], 'w')
    
    for r in SeqIO.parse(ifile, "fastq"):
        SeqIO.write(r, ofile, "fasta")

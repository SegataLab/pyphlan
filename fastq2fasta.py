#!/usr/bin/env python
from Bio import SeqIO
import argparse as ap
import sys

def read_params(args):
	p = ap.ArgumentParser(description = 'fastq2fasta.py Parameters\n')
	p.add_argument('--ifn', required = True, default = None, type = str)
	p.add_argument('--ofn', required = True, default = None, type = str)
	return vars(p.parse_args())
	
if __name__ == '__main__':
	args = read_params(sys.argv)
	with open(args['ofn'], 'w') as ofile:
		for r in SeqIO.parse(args['ifn'], "fastq"):
			SeqIO.write(r, ofile, "fasta")

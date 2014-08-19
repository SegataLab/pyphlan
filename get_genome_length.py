#!/usr/bin/env python

import sys
import Bio
from Bio import SeqIO	
import argparse as ap


def get_genome_length(ifn):
	length = 0
	for contig in SeqIO.parse(ifn, 'fasta'):
		length += len(contig.seq)
	return length

def read_params(args):
	p = ap.ArgumentParser(description = 'get_genome_length.py Parameters\n')
	p.add_argument('--ifn', required = True, default = None, type = str)
	return vars(p.parse_args())
	
if __name__ == '__main__':
	args = read_params(sys.argv)
	print get_genome_length(args['ifn'])

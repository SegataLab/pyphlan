#!/usr/bin/env python


import sys
import argparse as ap
from Bio import Phylo
import os
import pyphlan as ppa


def read_params():
    p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter,
                          description='Append a label in front of the current node labels')

    p.add_argument('intree', nargs='?', default=None, type=str,
                   help="the input tree, stdin if not present")
    p.add_argument('-f', default='newick', type=str,
                   choices=['newick', 'nexus', 'nexml', 'phyloxml', 'cdao'],
                   help="inp/out -put format")
    p.add_argument('-l', required=True, type=str,
                   help="the set of leaves (file or comma-separated values) to use for computing the partial branch length")

    args = vars(p.parse_args())

    if os.path.isfile(args['l']):
        args['l'] = [r.strip() for r in open(args['l']) if r.strip()]
    else:  # the input list is given as comma-separated values
        args['l'] = args['l'].split(',')

    if not args['intree']:
        args['intree'] = sys.stdin

    return args


if __name__ == "__main__":
    args = read_params()
    tree = Phylo.read(args['intree'], args['f'])

    print(ppa.partial_branch_length(tree.clade, args['l']))

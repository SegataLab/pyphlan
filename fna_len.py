#!/usr/bin/env python


from __future__ import with_statement
import sys
import argparse
import utils
from Bio import SeqIO
import numpy as np
import os


def read_params():
    p = argparse.ArgumentParser(description='Display the length of each fasta entry')
    p.add_argument('inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
        help="the input fna file [stdin if not present]")
    p.add_argument('out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
        help="the output txt file [stdout if not present]")
    p.add_argument('-q', action='store_true',
        help="set this for fastq")
    p.add_argument('-t', '--total', action='store_true', default=False,
        help="Print only the sum of the length of all sequences\n")
    p.add_argument('-s', '--stat', action='store_true', default=False,
        help="Print only the statistics about the length of sequences\n")

    return vars(p.parse_args())


if __name__ == '__main__':
    par = read_params()
    lvec = []
    samplename = None

    if par['inp_f']:
        samplename = par['inp_f'][par['inp_f'].rfind('/')+1:]
        samplename = samplename[:samplename.rfind('.')]
    else:
        samplename = 'stdin'
        samplename += '_fastq' if par['q'] else '_fasta'

    with utils.openw(par['out_f']) as outf:
        for r in SeqIO.parse(utils.openr(par['inp_f']), "fastq" if par['q'] else "fasta"):
            l = len(r.seq)

            if par['stat'] or par['total']:
                lvec.append(l)
            else:
                outf.write("\t".join([r.id, str(l)]) + "\n")

        if par['total']:
            outf.write(str(np.sum(lvec)) + "\n")

        if par['stat']:
            outf.write('\t'.join(["#samplename", "n_of_bases", "n_of_reads", "min_read_len", "median_read_len", "mean_read_len", "max_read_len"]) + '\n')
            outf.write('\t'.join([str(a) for a in [samplename, np.sum(lvec), len(lvec), np.min(lvec), np.median(lvec), np.mean(lvec), np.max(lvec)]]) + '\n')

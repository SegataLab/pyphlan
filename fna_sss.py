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

def read_params(args):
    parser = argparse.ArgumentParser(description='Split/Select/Randomize/Subsample a multi fasta file')
    arg = parser.add_argument
    arg( 'inp_f', metavar='INPUT_FILE', nargs='?', default=None, type=str,
         help="the input fna file [stdin if not present]")
    arg( 'out_f', metavar='OUTPUT_FILE', nargs='?', default=None, type=str,
         help="the output fna file [stdout if not present]")

    parser.add_argument('--select', action='store_true', help="Select fna entries\n")
    
    parser.add_argument('--subsample', metavar='Random subsample probability', 
        default=None, type = float )
    parser.add_argument('--randomize', help='Randomize output order. Not very memory efficient right now!', action='store_true' )
    parser.add_argument('--split', metavar='nsplits', default=1, type = int )

    parser.add_argument('--ids', metavar='s', default="", type=str, 
        help="the list of entries to select (separated by ::: if given as string otherwise as \\n if a file is given)")

    return vars(parser.parse_args())

class read:
    def __init__( self, n = None, seq = None ):
        self.n = n
        self.seq = "" if seq is None else seq

    def __str__( self ):
        return "\n".join( [">"+self.n] + textwrap.wrap(self.seq, 120) ) +"\n" 

    def __repr__( self ):
        return str(self)

class reader(object):
    def __init__( self, fn ):
        self.ret = False
        #openr = bz2.BZ2File if bool(fn) and fn.endswith(".bz2") else open
        self.inp = utils.openr(fn) if bool(fn) else sys.stdin
        self.cr = read( )

    def __iter__(self):
        return self

    def next( self ):
        while True:
            try:
                l = self.inp.next().strip()
                if l[0] == '>':
                    ret = None
                    if self.cr.n and len(self.cr.seq) > 0:
                        ret = read(self.cr.n,self.cr.seq)
                    self.cr.n = l[1:]# hash(l[1:])
                    self.cr.seq = ""
                    if ret:
                        return ret
                else:
                    self.cr.seq += l
            except StopIteration:
                if self.cr.n and len(self.cr.seq) > 0:
                    ret = read(self.cr.n,self.cr.seq)
                    self.cr.n = None
                    return ret
                raise StopIteration 

def sss( par ):
    subsample = bool(par['subsample']) 
    select = bool(par['select'])
    randomize = bool(par['randomize'])
    if bool(par['out_f']):
        n = par['split']
        #openw = bz2.BZ2File if par['out_f'].endswith(".bz2") else open
        if n == 1:
            out_stream = [utils.openw( par['out_f'], "w")]
        else:
            out_stream = [utils.openw( par['out_f']+str(r).zfill(len(str(n)))+".fna"+(".bz2" if par['out_f'].endswith(".bz2") else ""), "w") for r in range(n)]
    else:
        out_stream = [sys.stdout] # larger buffer?

    if select:
        if os.path.exists(par['ids']):
            #openr = bz2.BZ2File if par['ids'].endswith(".bz2") else open 
            es = [s.strip().split('\t')[0] for s in utils.openr(par['ids'])]
        else:
            es = [(s.split("$")[1] if s.count("$") else s) for s in  par['ids'].split(":::")]
        es = set(es)

    all_reads = []
    nstreams = len( out_stream )

    p = par['subsample']
    reads = reader( par['inp_f'] )
    cind = 0 
    for r in reads:
        if select and r.n not in es:
            continue
        if subsample and rnd.random() > p:
            continue
        if randomize:
            all_reads.append( r )
            continue
        out_stream[cind].write(  str(r)  )
        cind = (cind + 1) % nstreams
    
    if randomize:
        rnd.shuffle(all_reads)
        step = len(all_reads) / nstreams 
        for i,r in enumerate(all_reads):
            out_stream[cind].write( str(r) )
            if not i % step:
                cind = (cind + 1) % nstreams

    for o in out_stream:
        o.close()

if __name__ == '__main__':
    params = read_params(sys.argv)
    sss( params )


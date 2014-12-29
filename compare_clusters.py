#!/usr/bin/env python

import sys
import collections
import utils
import numpy as np
from Bio import SeqIO 
import subprocess
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert core gene files to core gene summaries\n')

    p.add_argument( '-1', default=None, type=str,
            help=   "Cluster file 1")
    p.add_argument('-2', default=None, type=str,
            help=   "Cluster file 2")
    p.add_argument('-i', default=None, type=str,
            help=   "Input sequences")

    return vars( p.parse_args() )

def get_clusters( clf ):
    with utils.openr(clf) as inp:
        return dict([(i,{'cl':set(l.strip().split('\t'))}) for i,l in enumerate(inp)])
    return None

def add_sequences( cl, ffn ):
    for i,c in cl.items():
        seqs = [ffn[idf] for idf in c['cl']]
        cl[i]['ffn'] = seqs
        print sorted([len(s) for s in seqs])
        sys.stdin.readline()
    return cl

def edit_distance( x, y ):

    start, stop = None, None
    for i,(a,b) in enumerate(zip(x,y)):
        if a != '-' and b != '-':
            start = i
            break
    for i,(a,b) in enumerate(zip(x,y)[::-1]):
        if a != '-' and b != '-':
            stop = i
            break

    eql = len([1 for a,b in zip(x[start:-stop],y[start:-stop]) if a == b])

    return float(eql)/(len(x)-(stop+start))

def add_alignments( cl ):
    for i,c in cl.items():
        cline = MuscleCommandline(clwstrict=True)
        child = subprocess.Popen( str(cline),
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  universal_newlines=True,
                                  shell=(sys.platform!="win32")) 
        SeqIO.write(c['ffn'], child.stdin, "fasta")
        child.stdin.close()
        align = AlignIO.read(child.stdout, "clustal")
        cl[i]['aln'] = align

        for a in align:
            for b in align:
                print edit_distance(str(a.seq),str(b.seq))
                print a.id
                print str(a.seq)
                print b.id
                print str(b.seq)
        sys.exit()

    return cl
        


    
def stats( cl ):
    lens = [len(c['cl']) for c in cl.values()]
    avg, stdev = np.mean(lens), np.std(lens)
    singletons = len([1 for l in lens if l == 1])
    print avg, stdev, len(lens), singletons
    for l in sorted(set(lens),reverse=True):
        print l, lens.count(l)

if __name__ == "__main__":
    args = read_params( sys.argv )

    cl1 = get_clusters( args['1'] )
    cl2 = get_clusters( args['2'] )

    stats(cl1)
    stats(cl2)

    ffn = SeqIO.to_dict(SeqIO.parse( args['i'], "fasta"))

    cl1 = add_sequences( cl1, ffn ) 
    cl2 = add_sequences( cl2, ffn )

    cl2 = add_alignments( cl2 )

    """
    equal = []
    for c1 in cl1:
        for c2 in cl2:
            if len(c1) == len(c2) and len(c1) == len(c1 & c2):
                equal.append(c1)
    print len(equal)

    for c1 in cl1:
        for c2 in cl2:
            if not (len(c1) == len(c2) and len(c1) == len(c1 & c2)):
                if c1 <= c2:
                    print c1, " contained in ", c2
    """

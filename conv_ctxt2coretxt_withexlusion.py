#!/usr/bin/env python

import sys
import collections
import utils
try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Convert core gene txt file'
            ' substituting gene IDs with genomes IDs\n')

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input ctxt file [stdin if not present]")
    p.add_argument('--g2t', metavar="Mapping file from genes to taxa",
            default=None, type=str )
    p.add_argument('--t2g', metavar="Mapping file from taxa to genes",
            default=None, type=str )
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output gtxt file, compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")
    #p.add_argument('-n', default=1, type=int )
    p.add_argument('--ncores', default=1, type=int )
    p.add_argument('-c', default=1, type=int )
    p.add_argument('--out_taxa', default=0, type=int )
    p.add_argument('--max_tot_repeats', default=0, type=int )
    p.add_argument('--removed_taxa', required = True, type=str )
    p.add_argument('--sk', action='store_true' )

    return vars( p.parse_args() )


def k_subsets_i(n, k):
    '''
    Yield each subset of size k from the set of intergers 0 .. n - 1
    n -- an integer > 0
    k -- an integer > 0
    '''
    # Validate args
    if n < 0:
        raise ValueError('n must be > 0, got n=%d' % n)
    if k < 0:
        raise ValueError('k must be > 0, got k=%d' % k)
    # check base cases
    if k == 0 or n < k:
        yield set()
    elif n == k:
        yield set(range(n))
    else:
        # Use recursive formula based on binomial coeffecients:
        # choose(n, k) = choose(n - 1, k - 1) + choose(n - 1, k)
        for s in k_subsets_i(n - 1, k - 1):
            s.add(n - 1)
            yield s
        for s in k_subsets_i(n - 1, k):
            yield s

def k_subsets(s, k):
    '''
    Yield all subsets of size k from set (or list) s
    s -- a set or list (any iterable will suffice)
    k -- an integer > 0
    '''
    s = list(s)
    n = len(s)
    for k_set in k_subsets_i(n, k):
        yield set([s[i] for i in k_set])

def count_ncores( pang, to_excl, tot ):
    ntot = tot - len(to_excl)
    n = 0
    for k,v in pang.items():
        if len(v) < ntot:
            continue
        if len( v - to_excl ) >= ntot:
            n += 1
    return n

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )
    
    gint = str if args['sk'] else int

    if not args['g2t'] and not args['t2g']:
        sys.stdout.write("Error one of --t2g and --g2t must be provided\n")
        sys.exit(0)
    g2t = {}


    if args['g2t']:
        with utils.openr( args['g2t'] ) as inp:
            #g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
            for l in inp:
                f,t = l.strip().split('\t')
                g2t[gint(f)] = gint(t)
    elif args['t2g']:
        with utils.openr( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[gint(g)] = gint(ll[0])

    pangenome = {}
    all_genomes = set()
    with utils.openr( args['ctxt'] ) as inp:
        for l in inp:
            line = l.strip().split('\t')
            genomes = set([g2t[al] for al in line])
            pangenome[line[0]] = genomes
            all_genomes |= genomes
    #print all_genomes 
    #print pangenome
    nn = count_ncores( pangenome, set(), len(all_genomes) )
    run = 1
    torem = set()
    while nn < args['ncores']:
        nn = 0
        for g in k_subsets( all_genomes, run ):
            newn = count_ncores( pangenome, g, len(all_genomes) )
            if newn > nn:
                nn = newn
                torem |= g 
        run += 1
    #print nn

    with utils.openw(args['removed_taxa']) as out:
        for a in list(torem):
            out.write( a + "\n" )

    ntaxa = len(all_genomes) - len(torem)

    with utils.openw(args['txt']) as out:
        with utils.openr( args['ctxt'] ) as inp:
            for l in inp:
                valin = [gint(a) for a in l.strip().split('\t') if g2t[a] not in torem]
                if len(valin) < ntaxa:
                    continue
                ko = False
                genes2taxa = collections.defaultdict(list)
                for vv in valin:
                    genes2taxa[g2t[vv]].append( vv )
                    if len( genes2taxa[g2t[vv]] ) > args['c']:
                        ko = True
                        break
                    #totlen += len(genes2taxa[g2t[vv]])
                    #genes2taxa[g2t[vv]] = genes2taxa[g2t[vv]][:1]

                if not ko:
                    totlen = 0
                    for k1, v1 in genes2taxa.items():
                        totlen += len( v1 )
                        genes2taxa[k1] = v1[:1]
                    if totlen - ntaxa > args['max_tot_repeats']:
                        ko = True

                if len(genes2taxa) < ntaxa:
                    continue

                if ko:
                    continue
                valin = list(valin) 
   
                if args['out_taxa']:
                    #out.write( "\t".join([str(g2t[s]) for s in valin]) +"\n" ) 
                    out.write( "\t".join([str(s) for s in genes2taxa.keys()]) +"\n" ) 
                else:    
                    out.write( "\t".join([str(s[0]) for s in genes2taxa.values()]) +"\n" )




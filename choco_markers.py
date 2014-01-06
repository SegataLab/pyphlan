#!/usr/bin/env python

import sys
import collections
import utils
#import cPickle as pickle
import pickle
import pyphlan as ppa  
from Bio import SeqIO

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Profile ChocoPhlAn genes\n')

    p.add_argument( '--sam', required = True, default=None, type=str )
    p.add_argument( '--centroids', required = True, default=None, type=str )
    #p.add_argument( '--centroids_cores', action='store_true' )
    p.add_argument( '--centroids_ffn', required = True, default=None, type=str )
    p.add_argument( '--lens', required = True, default=None, type=str )
    p.add_argument( '--g2c', required = True, default=None, type=str )
    p.add_argument( '--taxonomy', required = True, default=None, type=str )
    p.add_argument( '--out_markers', required = True, default=None, type=str )
    p.add_argument( '--out_ml', required = True, default=None, type=str )
    p.add_argument( '--out', required = True, default=None, type=str )
    p.add_argument( '--out_mpa_pkl', required = True, default=None, type=str )
    p.add_argument( '--out_m2c', required = True, default=None, type=str )
    p.add_argument( '--out_summary', required = True, default=None, type=str )
    p.add_argument( '--min_n_markers', default=10, type=int )
    p.add_argument( '--min_n_markers_strains', default=20, type=int )
    p.add_argument( '--top_n_markers', default=200, type=int )
    p.add_argument( '--score_th', default=100, type=int )
    p.add_argument( '--include_strains', default=False, action = 'store_true' )
    return vars( p.parse_args() )
    
if __name__ == "__main__":
    args = read_params( sys.argv )
    
    tree = ppa.PpaTree( args['taxonomy'],  lev_sep = "|" )

    clades2terms = ppa.clades2terms( tree.tree )
    
    clades2taxa = dict([(clade.full_name,{'taxa':set([taxon.name[3:] for taxon in taxa]),'genes':set()}) for clade,taxa in clades2terms.items()])
    #ttab = 1 if args['centroids_cores'] else 2
    #gtab = 0 if args['centroids_cores'] else 1
    ttab = 2
    gtab = 1
    genes2taxa, all_genes, all_reads = {}, set(), set()
    with open(args['centroids']) as inp:
        for line in (l.split('\t') for l in inp):
            taxon = line[ttab]
            gene = line[gtab]
            genes2taxa[line[gtab]] = taxon
            all_genes.add( gene )
            clades2taxa[taxon]['genes'].add(gene)


    contigs2genomes = {}
    with open( args['g2c'] ) as inp:
        lines = (list(i.strip().split('\t')) for i in inp)
        for line in lines: 
            genome = line[0]
            for l in line[1:]:
                contigs2genomes[l] = genome
    
    reads2lens = {}
    #genes2nreads = collections.defaultdict( int )
    with open( args['lens'] ) as inp:
        lines = (i.strip().split('\t') for i in inp)
        for r,l in lines:
            #reads2lens[r] = int(l)
            gene = "_".join(r.split("_")[:-1])
            if gene in all_genes:
                #genes2nreads[gene] += 1
                all_reads.add( r )
                reads2lens[r] = int(l)


    #genes2nreads = collections.defaultdict( int )
    #for r in reads2lens:
    #    gene = "_".join(r.split("_")[:-1])
    #    if gene in all_genes:
    #        genes2nreads[gene] += 1
    #        all_reads.add( r )


    reads2ghits = {}
    inp = bz2.BZ2File( args['sam'] )
    lines = (list(l.split('\t')) for l in inp)
    lines_strip = ([l[0],l[1],int(l[2])] for l in lines if l[0] in all_reads)
    for fr,to,val in lines_strip:
        to_genome = contigs2genomes[to]
        valf = float(val)/float(reads2lens[fr])
        if fr not in reads2ghits:
            reads2ghits[fr] = {to_genome: valf}
        elif to_genome not in reads2ghits[fr]:
            reads2ghits[fr][to_genome] = valf
        elif to_genome in reads2ghits[fr] and valf > reads2ghits[fr][to_genome]:
            reads2ghits[fr][to_genome] = valf
    inp.close()

    genes2genomes = {}
    for r,h in reads2ghits.items():
        gene = "_".join(r.split("_")[:-1])
        if gene not in genes2genomes:
            genes2genomes[gene] = {}
        
        for hh,vv in h.items():
            if hh not in genes2genomes[gene]:
                genes2genomes[gene][hh] = set()
            genes2genomes[gene][hh].add( vv )

    ffn = SeqIO.to_dict(SeqIO.parse(args['centroids_ffn'], "fasta"))
    ffn_len = dict([(k,len(v)) for k,v in ffn.items()])

    res = {}
    for gene,taxon in genes2taxa.items():
        if gene not in genes2genomes:
            continue
        int_genomes = set(clades2taxa[taxon]['taxa'])
        genomes_hit = set(genes2genomes[gene])
        int_genomes_hit = int_genomes & genomes_hit
        ext_genomes_hit = genomes_hit - int_genomes
       
        n_int_genomes, n_genomes_hit, n_int_genomes_hit, n_ext_genomes_hit = len(int_genomes), len(genomes_hit), len(int_genomes_hit), len(ext_genomes_hit)

        len_pen = 0 
        if ffn_len[gene] < 250:
            len_penlen_pen = 6
        if ffn_len[gene] < 500:
            len_penlen_pen = 4
        if ffn_len[gene] < 750:
            len_penlen_pen = 2
        
        score = n_ext_genomes_hit + len_pen
        score += 5.0*float( n_int_genomes-n_int_genomes_hit ) / float(n_int_genomes) 

        if taxon not in res:
            res[taxon] = {}

        if score > args['score_th']:
            continue

        if "t__" in taxon and score > 0:
            continue

        res[taxon][gene] = { 'score' : score, 'ext' : ext_genomes_hit }


    selected_markers = set()
    for taxon, markers in res.items():
        for i, (marker, marker_v) in enumerate(sorted( markers.items(), key = lambda x: x[1]['score'] )):
            if i >= args['top_n_markers']:
                del res[taxon][marker]
            #res[taxon][marker]['seq'] = ffn[marker]

    
    markers_ffn = []
    to_pkl = {'markers':{}}
    with open( args['out_summary'], "w" ) as outf:
        with open( args['out'], "w" ) as out:
            for taxon, markers in res.items():
                if not args['include_strains']  and "t__" in taxon:
                    continue
                if args['include_strains']  and "t__" in taxon:
                    if len(markers) < args['min_n_markers_strains']:
                        continue
                if "t__"  not in taxon and len(markers) < args['min_n_markers']:
                    continue
                
                outf.write( "\t".join( [str(taxon), str(len(markers))] ) +"\n" )

                for marker, marker_v in markers.items():
                    out.write( "\t".join([  marker,
                                            #taxon,
                                            genes2taxa[marker].split("|")[-1],
                                            str( ffn_len[marker]),
                                            str(marker_v['score']),
                                            str(len(marker_v['ext'])),
                                            str(",".join(marker_v['ext'])),
                                            ])
                                            +"\n" )
                    to_pkl['markers'][marker] = {   'taxon':taxon,
                                                    'clade':genes2taxa[marker].split("|")[-1],
                                                    'len': ffn_len[marker],
                                                    'score': marker_v['score'],
                                                    'ext': marker_v['ext'] }

                    res[taxon][marker]['seq'] = ffn[marker]
                    markers_ffn.append( ffn[marker] )
                    selected_markers.add( ffn[marker].id  )
    
    SeqIO.write( markers_ffn, args['out_markers'], "fasta")

    to_pkl['taxonomy'] = [l.strip() for l in open(args['taxonomy'])]

    #with open(args['out_mpa_pkl'], 'wb') as out:
    #   bz2.compress(pickle.dump(to_pkl, out, pickle.HIGHEST_PROTOCOL))
    out = bz2.BZ2File(args['out_mpa_pkl'],"wb")
    pickle.dump(to_pkl, out, pickle.HIGHEST_PROTOCOL)
    out.close()




    #with open( args['out_ml'], "w" ) as outf:
    with open( args['out_ml'], "w" ) as out_ml:
        with open( args['out_m2c'], "w" ) as out_m2c:
            for k,v in ffn_len.items():
                if not args['include_strains'] and "t__" in genes2taxa[k]:
                    continue
                if k not in selected_markers:
                    continue
                #outf.write( "\t".join([k,str(v)]) + "\n" )
                out_ml.write( "\t".join([k,str( ffn_len[k]  )]) + "\n" ) 
                out_m2c.write( "\t".join([k, genes2taxa[k].split("|")[-1]  ]) + "\n" ) 
                #out_m2c.write( "\t".join([k,str( ffn_len[k]  )]) + "\n" ) 




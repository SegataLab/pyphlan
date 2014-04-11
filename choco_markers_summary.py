#!/usr/bin/env python

import sys
import collections
import utils
import cPickle as pickle
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

    p.add_argument( '--inp', required = True, default=None, type=str )
    p.add_argument( '--taxonomy', required = True, default=None, type=str )
    return vars( p.parse_args() )
    
if __name__ == "__main__":
    args = read_params( sys.argv )
    
    tree = ppa.PpaTree( args['taxonomy'],  lev_sep = "|" )

    with open( args['inp'] ) as inp:
        markers = dict([(line[0],int(line[1])) for line in (l.strip().split('\t') for l in inp)])


    def traverse( clade, father = None ):
        if len(clade.clades) == 1:
            if clade.full_name in markers:
                print clade.full_name, markers[clade.full_name]
                traverse(clade.clades[0],markers[clade.full_name])
            else:
                if father is None:
                    print clade.full_name, 0
                    traverse(clade.clades[0])
                else:
                    print clade.full_name, father
                    traverse(clade.clades[0],father)
        else:
            if clade.full_name in markers:
                print clade.full_name, markers[clade.full_name] 
            else:
                if father is None:
                    print clade.full_name, 0
                else:
                    print clade.full_name, father
            for c in clade.clades:
                traverse(c)

    traverse(tree.tree.root)





    """   


    clades2terms = ppa.clades2terms( tree.tree )
    
    clades2taxa = dict([(clade.full_name,{'taxa':set([taxon.name[3:] for taxon in taxa]),'genes':set()}) for clade,taxa in clades2terms.items()])
    
    genes2taxa, all_genes, all_reads = {}, set(), set()
    with open(args['centroids']) as inp:
        for line in (l.split('\t') for l in inp):
            taxon = line[2]
            gene = line[1]
            genes2taxa[line[1]] = taxon
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
        
        score = n_ext_genompes_hit + len_pen
        score += 5.0*float( n_int_genomes-n_int_genomes_hit ) / float(n_int_genomes) 

        if taxon not in res:
            res[taxon] = {}

        if score > 100:
            continue

        res[taxon][gene] = { 'score' : score, 'ext' : ext_genomes_hit }


    selected_markers = set()
    for taxon, markers in res.items():
        for i, (marker, marker_v) in enumerate(sorted( markers.items(), key = lambda x: x[1]['score'] )):
            if i >= args['top_n_markers']:
                del res[taxon][marker]
            #res[taxon][marker]['seq'] = ffn[marker]

    
    markers_ffn = []

    with open( args['out_summary'], "w" ) as outf:
        for taxon, markers in res.items():
            if "t__" in taxon:
                continue
            if len(markers) < args['min_n_markers']:
                continue
            outf.write( "\t".join( [str(taxon), str(len(markers))] ) +"\n" )
            for marker, marker_v in markers.items():
                res[taxon][marker]['seq'] = ffn[marker]
                markers_ffn.append( ffn[marker] )
                selected_markers.add( ffn[marker].id  )
    
    SeqIO.write( markers_ffn, args['out_markers'], "fasta")

    #with open( args['out_ml'], "w" ) as outf:
    with open( args['out_ml'], "w" ) as out_ml:
        with open( args['out_m2c'], "w" ) as out_m2c:
            for k,v in ffn_len.items():
                if "t__" in genes2taxa[k]:
                    continue
                if k not in selected_markers:
                    continue
                #outf.write( "\t".join([k,str(v)]) + "\n" )
                out_ml.write( "\t".join([k,str( ffn_len[k]  )]) + "\n" ) 
                out_m2c.write( "\t".join([k, genes2taxa[k].split("|")[-1]  ]) + "\n" ) 
                #out_m2c.write( "\t".join([k,str( ffn_len[k]  )]) + "\n" ) 


    """

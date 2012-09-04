#!/usr/bin/env python

import sys
import collections
import utils
import pandas
import StringIO

try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Generate a taxonomic txt hierarchy'
            ' with genome assignments from IMG taxa table\n')

    
    p.add_argument( 'img', nargs='?', default=None, type=str,
            help=   "the input uc file [stdin if not present]")
    p.add_argument('txt', nargs='?', default=None, type=str,
            help=   "the output txt file compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")
    p.add_argument('-d', metavar="Domain",
            default='Mic', choices=['Mic','Vir','Euk'] )

    return vars( p.parse_args() )

qm = "?"
def get_qm( s ):
    if s is None:
        return qm
    if not s:
        return qm
    if type(s) != str:
        return qm
    if s in ['Unclassified','unclassified']:
        return qm
    if s in ['sp','sp.','Sp','Sp.','spp','spp.','Spp','Spp.']:
        return qm
    return s.replace("Candidatus ","").replace(" ","_").replace(".","_").replace(",","_")


if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )

    tax_lev = "dpcofgs"
    tax_lev_exp = ['Domain','Phylum','Class','Order','Family','Genus','Species']

    with utils.openr(args['img'],"rU") as inp:
        table = pandas.read_table( inp, sep='\t', index_col=0)
 
    if args['d'] == 'Mic': 
        table = table[table['Domain'].isin(['Bacteria','Archaea'])]
        table = table[table['Gene Count'] > 250.0]
        table = table[table['Genome Size'] > 50000.0]
    elif args['d'] == 'Vir':
        table = table[table['Domain'].isin(['Viruses'])]
        table = table[table['Gene Count'] > 0.0]
    elif args['d'] == 'Euk':
        table = table[table['Domain'].isin(['Eukaryota'])]
        able = table[table['Gene Count'] > 1000.0]
        table = table[table['Genome Size'] > 500000.0]

    table = table.reindex(columns=tax_lev_exp+['Genome Name'])

    with utils.openw(args['txt']) as out:
        for i,t in table.iterrows():
            out.write( "\t".join( [ #str(-int(i)),
                                  ".".join( ["__".join([taxl, get_qm(t[taxle])]) 
                                      for taxl,taxle in  zip(tax_lev,tax_lev_exp)] + ["t__"+str(-int(i) if i.is_integer() else "Nan")]
                                        ),
                                  #str(t['Genome Name'])
                                  ] ) + "\n")
    
    

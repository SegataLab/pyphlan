#!/usr/bin/env python
import sys
import collections
import utils
import cPickle as pickle
import pyphlan as ppa  
import argparse as ap

def read_params( args ):
    p = ap.ArgumentParser(description='Reduce Sam for ChocoPhlAn\n')
    p.add_argument( '--sam', required = False, default=None, type=str )
    p.add_argument( '--out', required = False, default=None, type=str )
    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    with utils.openr( args['sam'] ) as inp:

        lines = (l.split('\t') for l in inp)
        out = utils.openw( args['out'] )
        for l in lines:
            out.write(l[0]+"\t"+l[2]+"\t"+l[11][5:] + "\n" )
            #out.write(l[0]+"\t"+l[2]+"\t"+[ll for ll in l[3:] if "AS:i:" in ll][0][5:] + "\n" )
        out.close()



#!/usr/bin/env python

import sys
import collections
import cPickle as pickle
import utils
import bz2
try:
    import argparse as ap
    import bz2 
except ImportError:
    sys.stderr.write( "argparse not found" )
    sys.exit(-1)

def read_params( args ):
    p = ap.ArgumentParser(description='Add markerness'
            ' to core gene file\n')

    p.add_argument( 'ctxt', nargs='?', default=None, type=str,
            help=   "the input ctxt file [stdin if not present]")
    p.add_argument('--b6o', metavar="The outfmt6 file for the cores",
            default=None, type=str )
    p.add_argument('-n', metavar="Total number of target sets (total targets in the b6o file if unspecified)",
            default=None, required = True, type=int )
    p.add_argument('--g2t', metavar="Mapping file from genes to taxa",
            default=None, type=str )
    p.add_argument('--t2g', metavar="Mapping file from taxa to genes",
            default=None, type=str )
    p.add_argument('--pkl', metavar="The compressed pickled output",
            default=None, type=str )
    p.add_argument('mtxt', nargs='?', default=None, type=str,
            help=   "the output mtxt file, compressed if fiven with bz2 extension\n"
                    "[stdout if not present]")

    return vars( p.parse_args() )

if __name__ == "__main__":
    args = read_params( sys.argv )
    uc2cl = collections.defaultdict( set )
   
    if not args['g2t'] and not args['t2g']:
        sys.stdout.write("Error one of --t2g and --g2t must be provided\n")
        sys.exit(0)
    g2t = {}
    if args['g2t']:
        with open( args['g2t'] ) as inp:
            g2t = dict(([a for a in l.strip().split('\t')] for l in inp))
            #g2t = dict(([int(a) for a in l.strip().split('\t')] for l in inp))
    elif args['t2g']:
        with open( args['t2g'] ) as inp:
            for ll in (l.strip().split('\t') for l in inp):
                for g in ll[1:]:
                    g2t[g] = ll[0]
    
    all_fr = set()
    with utils.openr( args['ctxt'] ) as inp:
        valin = (l.strip().split('\t') for l in inp)

        g2c = collections.defaultdict( set )
        
        if args['b6o']:
            inp_mat = ((a,b) for a,b in (l.rstrip('\n').split("\t")[:2] for l in utils.openr(args['b6o'])))
    
            #all_targets = set()
            for fr,to in inp_mat:
                all_fr.add( fr )
                #all_targets.add( to )
                if fr != to:
                    g2c[fr].add( to )

        n = args['n'] # if args['n'] else len(all_targets)
        n = float(n)
   
        out_dict = collections.defaultdict( dict ) 
        all_fr_taxa = set([g2t[tt] for tt in all_fr])

        with utils.openw(args['mtxt']) as out:
            last,lastv = "",[]
            outbuf = []
            gt = None

            outtmp = [] 

            for v in valin:
                gt = v[0]
                if last == gt:
                    outtmp.append( v )
                    continue
                    
                if len(outtmp) == 1:
                    #print outtmp
                    outbuf.append( outtmp[0] )
                elif len(outtmp) > 1:
                    #print outtmp
                    outtmp2 = [o for o in outtmp if "_sp_" not in o[1] and "_bacterium_" not in o[1] and "_unclass" not in o[1]]
                    if len(outtmp2) == 1:
                        outbuf.append( outtmp2[0] )
                    elif len(outtmp2) > 1:
                        outtmp3 = sorted(outtmp2,key=lambda x:-int(x[2]))
                        cur = int(outtmp3[0][2])
                        other = sum([int(vv[2]) for vv in outtmp3[1:]])
                        #print cur, other, outtmp3[0]
                        if float(cur)/float(other) > 2.5 and float(other) < 10:
                            outbuf.append( outtmp3[0] )
                
                last = gt
                outtmp = [v]

            if last and last != gt:
                outbuf.append( lastv )
            """
            for v in valin:
                gt = v[0]
                if last == gt:
                    lastv = ""
                    continue
                if lastv:
                    outbuf.append( lastv )
                last = gt
                lastv = v
            if last and last != gt:
                outbuf.append( lastv )
            """
            for v in outbuf:
                fr = v[0]
                frt = g2t[fr]
                targett = list(set([g2t[s] for s in g2c[fr]])-all_fr_taxa) if args['b6o'] else []
                nu = len(g2c[fr]) - len(targett) if args['b6o'] else 0
                targets = ":".join([str(s) for s in targett]) if targett else "-"
                uniqueness = round(float(len(targett)) / n,3)
                out.write( "\t".join([str(g2t[fr])]+v+[str(nu),str(len(targett)),str(uniqueness),targets]) +"\n" )
                
                toadd = {   'gid' : v[0], 'tid' : g2t[fr], 'n_cores' : v[2], 'coreness' : float(v[4]), 
                            'uniqueness' : float(uniqueness),
                            'targets' : targett }
                out_dict[v[1]][v[0]] = toadd

            if args['pkl']:
                outf = open( args['pkl'], 'wb' ) 
                # pickle.loads( bz2.decompress( open("/tmp/a.pkl").read() ) )
                outf.write( bz2.compress(pickle.dumps(out_dict,pickle.HIGHEST_PROTOCOL)) )
                outf.close()


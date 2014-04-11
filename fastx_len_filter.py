from Bio import SeqIO
import sys


with sys.stdout as outf:
    for r in SeqIO.parse( sys.stdin, "fastq"):
        if len(r) > 100:
            SeqIO.write(r, outf, "fasta")



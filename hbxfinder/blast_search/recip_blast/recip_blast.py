#!/usr/bin/env python3
from subprocess import call as unix
import argparse

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 recip_blast.py


parse = argparse.ArgumentParser()

parse.add_argument("--infile",type=str, help="input file to run reciprocal BLAST on",required=True)
parse.add_argument("--outfile",type=str, help="output file name from reciprocal BLAST run",required=True)

args = parse.parse_args()


#Run reciprocal blast
print('Running reciprocal BLASTx search...')
unix('blastx -query ' + args.infile + ' -db ../../../hbx_data/homeobox -evalue 1e-5 -num_threads 4 -seg yes -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out ' + args.outfile, shell=True)


#Parse reciprocal BLASTx output
print('\n')
print('Parsing reciprocal BLAST output...')
sp_assem_list = {}
with open(args.outfile) as f:
    for line in f:
        sp = line.split('.')[0]
        sp_assem = line.split('|')[0]
        if sp_assem not in sp_assem_list:
            sp_assem_list[sp_assem] = sp


for assemb, species in sp_assem_list.items():
    print(species)
    with open(args.outfile) as f, open(assemb + '_hbx.fasta', 'w') as outF:
        for line in f:
            line = line.strip()
            sp = line.split('.')[0]
            sp_assem = line.split('|')[0]
            query = line.split('\t')[0]
            query_hit = query.split('|')[1:]
            query_hit = '|'.join(query_hit)
            subject = line.split('\t')[1:]
            perc_ident = float(subject[2])
            subject = '\t'.join(subject)
            if assemb == sp_assem:
                outF.write(query_hit + '\t' + subject + '\n')

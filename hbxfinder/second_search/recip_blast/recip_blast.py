#!/usr/bin/env python3
from subprocess import call as unix
import argparse
import os

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 recip_blast.py --gene <gene name>

parse = argparse.ArgumentParser()

parse.add_argument("--gene",type=str, help="name of homeobox gene to obtain results for",required=True)
parse.add_argument("--group",type=str, help="group of animals to use as seed search i.e. vertebrate or invertebrate",required=True)

args = parse.parse_args()

#Run reciprocal blast
print('Running reciprocal BLASTx search...')
unix('blastx -query genome_' + args.gene + '_recipBlast.fasta -db ../../../hbx_data/' + args.group + '_data/family_data/' + args.gene + ' -evalue 1e-5 -num_threads 5 -seg yes -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out recipBlast_' + args.gene + '.fa', shell=True)


#Parse reciprocal BLASTx output
print('\n')
print('Parsing reciprocal BLAST output...')
sp_list = []
with open('recipBlast_' + args.gene + '.fa') as f:
    for line in f:
        sp = line.split('|')[0]
        if sp not in sp_list:
            sp_list.append(sp)

os.makedirs('species_hbx', exist_ok=True)
for species in sp_list:
    print(species)
    with open('recipBlast_' + args.gene + '.fa') as f, open('species_hbx/' + species + '_' + args.gene + '.fasta', 'w') as outF:
        for line in f:
            line = line.strip()
            sp = line.split('|')[0]
            query = line.split('\t')[0]
            query_hit = query.split('|')[1:]
            query_hit = '|'.join(query_hit)
            subject = line.split('\t')[1:]
            perc_ident = float(subject[2])
            subject = '\t'.join(subject)
            if species == sp:
                outF.write(query_hit + '\t' + subject + '\n')

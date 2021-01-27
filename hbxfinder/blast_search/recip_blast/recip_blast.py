#!/usr/bin/env python3
from subprocess import call as unix

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 recip_blast.py


#Run reciprocal blast
print('Running reciprocal BLASTx search...')
unix('blastx -query best_assemb_recipBlast.fasta -db ../../../hbx_data/homeobox -evalue 1e-5 -num_threads 10 -seg yes -max_target_seqs 1 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out recipBlast_homeodomains.fa', shell=True)


#Parse reciprocal BLASTx output
print('\n')
print('Parsing reciprocal BLAST output...')
sp_assem_list = {}
with open('recipBlast_homeodomains.fa') as f:
    for line in f:
        sp = line.split('.')[0]
        sp_assem = line.split('|')[0]
        if sp_assem not in sp_assem_list:
            sp_assem_list[sp_assem] = sp


for assemb, species in sp_assem_list.items():
    print(species)
    with open('recipBlast_homeodomains.fa') as f, open(assemb + '_hbx.fasta', 'w') as outF:
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

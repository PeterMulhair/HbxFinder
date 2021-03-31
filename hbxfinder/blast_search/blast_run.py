#!/usr/bin/env python3
import glob
from joblib import Parallel, delayed
from subprocess import call as unix
from Bio import SeqIO
import argparse
import os

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 blast_run.py --path <path to genomes to blast against> 

parse = argparse.ArgumentParser()

parse.add_argument("--path",type=str, help="path to genomes in fasta format",required=True)

args = parse.parse_args()

finished_genome = []
for blastout in glob.glob("*blastoutput.fa"):
        genome = blastout.split('.blast')[0]
        finished_genome.append(genome)

db_list = []
sp_name_list = []
for fasta in glob.glob(args.path + "*.fasta"):
        sp_fasta = fasta.split('/')[-1].split('.fa')[0]
        if sp_fasta not in finished_genome:
                sp_name_list.append(sp_fasta)
                db_list.append(fasta)

#Do tBLASTn search using genome fasta file
os.makedirs('genome_blastdb', exist_ok=True)
def hbx_blast(fasta):
        sp_file = fasta.split('/')[-1]
        sp_name = sp_file.split('.fast')[0]
        print('Running BLAST;', sp_name)
        #Make blast database from genome fasta files
        unix('makeblastdb -dbtype nucl -in ' + fasta + ' -out genome_blastdb/' + sp_name, shell=True)
        #Search for hbx genes in unannotated genomes
        unix('tblastn -query ../../raw/hbx_data/homeobox.fasta -db genome_blastdb/' + sp_name + ' -evalue 1 -seg yes -max_target_seqs 5000 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out ' + sp_name + '.blastoutput.fa', shell=True)

Parallel(n_jobs=5)(delayed(hbx_blast)(sp_assem) for sp_assem in db_list)


#Parse tBLASTn output to get sequences for reciprocal BLASTx search
print('\n')
print('Parsing BLAST output...')
#os.mkdir('recip_blast')
assemb_sp = {}
for blast_out in sp_name_list:
        blast_output = blast_out + '.blastoutput.fa'
        sp = blast_out.split('.')[0]
        sp_assem = blast_out.split('.blastout')[0]
        assemb_sp[sp_assem] = sp
        

outF = open('recip_blast/best_assemb_recipBlast.fa', 'w')
for assemb, sp in assemb_sp.items():
        contig_nuc_assem = {}
                
        sp_assem = assemb + '.fasta'
        with open(args.path + sp_assem) as f:
                for record in SeqIO.parse(f, 'fasta'):
                        contig = record.id
                        seq = str(record.seq)
                        contig_nuc_assem[contig] = seq
                                
                                
        print(assemb)
        sp_blast = assemb + '.blastoutput.fa'
        with open(sp_blast) as f:
                for line in f:
                        lines = line.split('\t')
                        geneID = lines[0]
                        hbx_gene = geneID.split('|')[1]
                        scaff = lines[1]
                        scaff_nuc_seq = contig_nuc_assem[scaff]
                        
                        perc_ident = lines[3]
                        qstart = lines[5]
                        qend = lines[6]
                        
                        blast_hit = []
                        if int(lines[8]) < int(lines[9]):
                                sstart = int(lines[8]) - 10
                                send = int(lines[9]) + 10
                                for i in range(sstart, send + 1):
                                        try:
                                                nuc = scaff_nuc_seq[i]
                                                blast_hit.append(nuc)
                                                orientation = 'R'
                                        except:
                                                continue
                        else:
                                sstart = int(lines[9]) - 10
                                send = int(lines[8]) + 10
                                for i in range(sstart, send + 1):
                                        try:
                                                nuc = scaff_nuc_seq[i]
                                                blast_hit.append(nuc)
                                                orientation = 'L'
                                        except:
                                                continue
                        blast_hit = ''.join(blast_hit)
                        
                        outF.write('>' + assemb + '|' + hbx_gene +  '|' + scaff + '|' + str(sstart + 10) + '|' + str(send - 10) + '|' + orientation + '\n' + blast_hit + '\n')
                        
outF.close()

#!/usr/bin/env python3
import sys
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
parse.add_argument("--group",type=str, help="group of animals to use as seed search i.e. vertebrate or invertebrate",required=True)

args = parse.parse_args()

finished_genome = []
for blastout in glob.glob("*blastoutput.tsv"):
        genome = blastout.split('.')[0]
        finished_genome.append(genome)

db_list = []
sp_name_list = []
for fasta in glob.glob(args.path + "*"):
        if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
                sp_fasta = fasta.split('/')[-1].split('.')[0]
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
        unix('tblastn -query ../../data_hbx/' + args.group + '_data/homeobox.fasta -db genome_blastdb/' + sp_name + ' -evalue 1 -seg yes -max_target_seqs 5000 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out ' + sp_name + '.blastoutput.tsv', shell=True)

if len(db_list) >= 1:        
        Parallel(n_jobs=45)(delayed(hbx_blast)(sp_assem) for sp_assem in db_list)
else:
        print('No new genomes to annotate')
        sys.exit()

        
#Parse tBLASTn output to get sequences for reciprocal BLASTx search
print('\n')
print('Parsing BLAST output...')
assemb_sp = {}
for blast_out in db_list:
        blast_output = blast_out + '.blastoutput.tsv'
        sp = blast_out.split('/')[-1].split('.f')[0]
        sp_assem = blast_out.split('.blastout')[0]
        assemb_sp[sp_assem] = sp
        

outF = open('recip_blast/blast_parsed_output.fasta', 'w')
for assemb, sp in assemb_sp.items():
        contig_nuc_assem = {}
                
        with open(assemb) as f:
                for record in SeqIO.parse(f, 'fasta'):
                        contig = record.id
                        seq = str(record.seq)
                        contig_nuc_assem[contig] = seq
                                
                                
        print(sp)
        sp_blast = assemb.split('/')[-1] + '.blastoutput.tsv'
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
                        
                        outF.write('>' + sp + '|' + hbx_gene +  '|' + scaff + '|' + str(sstart + 10) + '|' + str(send - 10) + '|' + orientation + '\n' + blast_hit + '\n')
                        
outF.close()

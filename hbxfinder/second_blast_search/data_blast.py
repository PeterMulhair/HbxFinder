#!/usr/bin/env python3
import argparse
import glob
from collections import defaultdict
from Bio import SeqIO
from joblib import Parallel, delayed
from subprocess import call as unix
import os

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 data_blast.py --gene <Hbx gene class> --path <path to genomes to blast against>

parse = argparse.ArgumentParser()

parse.add_argument("--gene",type=str, help="name of homeobox gene to obtain results for",required=True)
parse.add_argument("--path",type=str, help="path to genomes in fasta format",required=True)

args = parse.parse_args()

#Parse results from first BLAST search to run second sequence search
os.mkdir('genome_hbx')
os.mkdir('genome_hbx/' + args.gene)
def genome_data(cluster):
    print('Parsing first BLAST run...')
    with open('../blast_search/recip_blast/hbx_contigs/' + cluster + '_cluster.txt') as f:
        lines = f.read()
        sp_cluster = lines.split('>')
        for sp in sp_cluster:
            if not sp:
                continue
            else:
                contig_gene_region = defaultdict(list)
                sp_info = sp.split('\n')
                sp_name = sp_info[0]
                for line in sp_info[1:]:
                    if not line:
                        continue
                    else:
                        gene_info = line.split('\t')
                        contig = gene_info[0]
                        gene_pos1 = gene_info[1].split(', ')[0].strip('(')
                        gene_pos2 = gene_info[1].split(', ')[1].strip(')')
                        contig_gene_region[contig].append(int(gene_pos1))
                        contig_gene_region[contig].append(int(gene_pos2))

                sp_genome = (args.path + sp_name + '.fasta')
                with open(sp_genome) as fg, open('genome_hbx/' + args.gene + '/' + sp_name + '_blastRegion.fasta', 'w') as outF:
                    for record in SeqIO.parse(fg, 'fasta'):
                        contigID = record.id
                        seq = str(record.seq)
                        
                        if contigID in contig_gene_region:
                            #print(sp_name, contigID)
                            blast_hit = []
                            reg = contig_gene_region[contigID]
                            reg_start = min(reg) - 1000
                            reg_end = max(reg) + 1000
                            cluster_range = range(reg_start, reg_end)
                            for i in cluster_range:
                                try:
                                    nuc = seq[i]
                                    blast_hit.append(nuc)
                                except:
                                    continue
                            blast_hit = ''.join(blast_hit)
                            
                            outF.write('>' + contigID + '\n' + blast_hit + '\n')

genome_data(args.gene)



sp_blast_list = []
for fasta in glob.glob("genome_hbx/" + args.gene + "/*.fasta"):
    sp_blast_list.append(fasta)

os.mkdir('genome_blastdb')
os.mkdir('genome_blastdb/' + args.gene)
os.mkdir('hbx_blastoutput')
def second_blast(fasta):
    #Make blast database from fasta files
    sp_fasta = fasta.split('/')[-1]
    blastdb = sp_fasta.split('.')[0]
    unix('makeblastdb -dbtype nucl -in ' + fasta + ' -out genome_blastdb/' + args.gene + '/' + blastdb, shell=True)
    #Search for hbx genes; output format 6
    print('\n')
    print('Running BLAST;',blastdb)
    unix('tblastn -query ../../raw/hbx_data/family_data/' + args.gene + '.fasta -db genome_blastdb/' + args.gene + '/' + blastdb + ' -evalue 1e-1 -seg yes -max_target_seqs 5000 -outfmt "6 qseqid sseqid evalue pident bitscore qstart qend qlen sstart send slen" -out hbx_blastoutput/' + blastdb + '_' + args.gene + '.blastoutput.fa', shell=True)

Parallel(n_jobs=20)(delayed(second_blast)(sp) for sp in sp_blast_list)



def seq_blast(fasta):
    #Search for hbx genes; output format 5
    sp_fasta = fasta.split('/')[-1]
    blastdb = sp_fasta.split('.')[0]
    unix('tblastn -query ../../raw/hbx_data/family_data/' + args.gene + '.fasta -db genome_blastdb/' + args.gene + '/' + blastdb + ' -evalue 1 -seg yes -max_target_seqs 5000 -outfmt "5" -out blast_seq_out/' + blastdb + '_' + args.gene +  '.blastoutput.xml', shell=True)

    
Parallel(n_jobs=20)(delayed(seq_blast)(sp) for sp in sp_blast_list)





#Parse tBLASTn output to get sequences for reciprocal BLASTx search
print('\n')
print('Parsing new BLAST output...')
assemb_sp = {}
for blast_out in glob.glob('hbx_blastoutput/*' + args.gene + '.blastoutput.fa'):
    blast_file = blast_out.split('/')[-1]
    sp = blast_file.split('.')[0]
    sp_assem = blast_file.split('.blastout')[0]
    sp = sp.split('_' + args.gene)[0]
    sp_assem = sp_assem.split('_' + args.gene)[0]
    assemb_sp[sp_assem] = sp

outF = open('recip_blast/genome_' + args.gene + '_recipBlast.fasta', 'a+')
for assemb, sp in assemb_sp.items():
    contig_nuc_assem = {}
    sp_assem = assemb + '.fasta'
    #print(assemb, sp)
    for blastRegion_file in glob.glob('genome_hbx/' + args.gene + '/' + sp + '*.fasta'):
        with open(blastRegion_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                contig = record.id
                seq = str(record.seq)
                contig_nuc_assem[contig] = seq
            
            
    print(assemb)
    #print(contig_nuc_assem.keys())
    sp_blast = 'hbx_blastoutput/' + assemb + '_' + args.gene + '.blastoutput.fa'
    with open(sp_blast) as f:
        for line in f:
            lines = line.split('\t')
            geneID = lines[0]
            hbx_gene = geneID.split('|')[1]
            scaff = lines[1]
            #print(scaff)
            scaff_nuc_seq = contig_nuc_assem[scaff]
            
            perc_ident = lines[3]
            qstart = lines[5]
            qend = lines[6]
            sstart = int(lines[8])
            send = int(lines[9])
            subject_range = (sstart, send)
            
            blast_hit = []
            if sstart < send:
                for i in range(sstart, send):
                    nuc = scaff_nuc_seq[i]
                    blast_hit.append(nuc)
                    orientation = 'R'
            else:
                for i in range(send, sstart):
                    nuc = scaff_nuc_seq[i]
                    blast_hit.append(nuc)
                    orientation = 'L'
            blast_hit = ''.join(blast_hit)
            
            outF.write('>' + assemb + '|' + hbx_gene +  '|' + scaff + '|' + str(sstart) + '|' + str(send) + '|' + orientation + '\n' + blast_hit + '\n')
            
outF.close()

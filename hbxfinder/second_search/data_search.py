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

finished_genome = []
for blastparse in glob.glob('genome_hbx/' + args.gene + '/*fasta'):
    genome = blastparse.split('/')[-1].split('_blastRegion')[0]
    finished_genome.append(genome)
    
os.makedirs('genome_hbx', exist_ok=True)
os.makedirs('genome_hbx/' + args.gene, exist_ok=True)
gene_orientation = {}
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
                        gene_hit_name = gene_info[2]
                        gene_pos1 = gene_info[1].split(', ')[0].strip('(')
                        gene_pos2 = gene_info[1].split(', ')[1].strip(')')

                        if int(gene_pos1) < int(gene_pos2):
                            gene_len = int(gene_pos1), int(gene_pos2)
                            gene_orientation[gene_len] = 'R'
                        else:
                            gene_len = int(gene_pos2), int(gene_pos1)
                            gene_orientation[gene_len] = 'L'
                        contig_gene_region[contig].append(gene_len)
                        
                sp_genome = (args.path + sp_name + '.fasta')
                with open(sp_genome) as fg, open('genome_hbx/' + args.gene + '/' + sp_name + '_blastRegion.fasta', 'w') as outF:
                    for record in SeqIO.parse(fg, 'fasta'):
                        contigID = record.id
                        seq = str(record.seq)
                        
                        if contigID in contig_gene_region:
                            blast_hit = []
                            reg = contig_gene_region[contigID]
                            
                            for gene in reg:
                                blast_hit = []
                                gene_start = gene[0] - 1000
                                gene_end = gene[1] + 1000
                                gene_range = range(gene_start, gene_end + 1)
                                for i in gene_range:
                                    try:
                                        nuc = seq[i]
                                        blast_hit.append(nuc)
                                    except:
                                        continue
                                blast_hit = ''.join(blast_hit)
                                outF.write('>' + contigID + '|' + str(gene_start) + '|' + str(gene_end) + '\n' + blast_hit + '\n')


genome_data(args.gene)

sp_blast_list = []
for fasta in glob.glob("genome_hbx/" + args.gene + "/*.fasta"):
    sp_fasta = fasta.split('/')[-1].split('_blastRegion')[0]
    if sp_fasta not in finished_genome:       
        sp_blast_list.append(fasta)


mmseq_outfile = []
os.makedirs('hbx_mmseqoutput', exist_ok=True)
def second_blast(fasta):
    sp_fasta = fasta.split('/')[-1]
    sp_assem = sp_fasta.split('.')[0]
    #Search for hbx genes with mmseq easy-search
    print('\n')
    print('Running MMseqs;',sp_assem)
    unix('mmseqs easy-search ../../raw/hbx_data/family_data/' + args.gene + '.fasta ' + fasta + ' hbx_mmseqoutput/' + sp_assem + '_' + args.gene + '.m8 tmp --spaced-kmer-pattern 1101111 -k 6 -a -e 1 --num-iterations 2', shell=True)
    mmseq_outfile.append(sp_assem + '_' + args.gene + '.m8')
        
Parallel(n_jobs=1)(delayed(second_blast)(sp) for sp in sp_blast_list)



#Parse MMseq search output to get sequences for reciprocal BLASTx search
print('\n')
print('Parsing MMseqs output...')

outF = open('recip_blast/genome_' + args.gene + '_recipBlast.fasta', 'a+')
for mmseq_out in mmseq_outfile:
    #mmseq_file = mmseq_out.split('/')[-1]
    sp = mmseq_out.split('_' + args.gene)[0]

    contig_nuc_assem = {}
    for blastRegion_file in glob.glob('genome_hbx/' + args.gene + '/' + sp + '*.fasta'):
        with open(blastRegion_file) as f:
            for record in SeqIO.parse(f, 'fasta'):
                contig = record.id
                seq = str(record.seq)
                contig_nuc_assem[contig] = seq
            
    print(sp)
    with open('hbx_mmseqoutput/' + mmseq_out) as f:
        for line in f:
            lines = line.split('\t')
            geneID = lines[0]
            hbx_gene = geneID.split('|')[1]
            scaff = lines[1]
            scaff_nuc_seq = contig_nuc_assem[scaff]
            contig = scaff.split('|')[0]
            hbx_start = int(scaff.split('|')[1]) + 1000
            hbx_end = int(scaff.split('|')[2]) - 1000
            hbx_len = hbx_start, hbx_end
            hbx_orientation = gene_orientation[hbx_len]
                        
            perc_ident = float(lines[2])*100
            qstart = lines[6]
            qend = lines[7]
            if int(lines[8]) < int(lines[9]):
                sstart = int(lines[8])
                send = int(lines[9])
            else:
                sstart = int(lines[9])
                send = int(lines[8])

            query_start_diff = int(qstart) - 1
            query_start_diff = query_start_diff*3
            query_end_diff = 60 - int(qend)
            query_end_diff = query_end_diff*3

            sub_cor_start = sstart - query_start_diff
            sub_cor_end = send + query_end_diff

            subject_range = range(sub_cor_start - 1, sub_cor_end + 1)

            blast_hit = []
            for i in subject_range:
                try:
                    nuc = scaff_nuc_seq[i]
                    blast_hit.append(nuc)
                except:
                    continue
            blast_hit = ''.join(blast_hit)

            outF.write('>' + sp + '|' + hbx_gene +  '|' + contig + '|' + str(hbx_start) + '|' + str(hbx_end) + '|' + '\n' + blast_hit + '\n')
            
outF.close()

import os
import glob
import json
import shutil
import argparse
from Bio import Align, SeqIO
from collections import defaultdict
from subprocess import call as unix

# Author: Peter Mulhair
# Date: 09/02/2020
# Usage python3 get_seq_data.py --gene <homeobox gene to search for>

aligner = Align.PairwiseAligner()
parse = argparse.ArgumentParser()

parse.add_argument("--gene",type=str, help="name of homeobox gene to parse",required=True)
parse.add_argument("--group",type=str, help="group of animals to use as seed search i.e. vertebrate or invertebrate",required=True)

args = parse.parse_args()

intron_list = ['Pb', 'Ro', 'Abd-B', 'lab']
#intron_list = ['Msx', 'NK1', 'NK4', 'NK6', 'Emx']

##Check if blast is installed localy
if shutil.which('sixpack') is None:
    sys.exit('sixpack not installed locally')
else:
    continue

##Open dictionary of Homeodomain classes and their genes
with open('../../../hbx_data/' + args.group + '_data/hbx_naming.json') as f:
    hbx_naming_dict = json.load(f)

hbx_dict = hbx_naming_dict[args.gene]

##Create dictionary of hbx gene IDs and AA sequences
hbx_gene_seq_dict = {}
with open('../../../hbx_data/' + args.group + '_data/family_data/' + args.gene + '.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        header = record.description
        gene = header.split('|')[1].strip()
        geneName = hbx_dict[gene]
        seq = str(record.seq)
        hbx_gene_seq_dict[geneName] = seq
    
#Create dictionary of the nucleotide sequences for the hbx genes
genome_hbx_seqs = defaultdict(list)
with open('../recip_blast/genome_' + args.gene + '_recipBlast.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        header = record.description
        sp = header.split('|')[0]
        gene = header.split('|')[1]
        geneName = hbx_dict[gene]
        contig = header.split('|')[2]
        gene_start = header.split('|')[3]
        gene_end = header.split('|')[4]
        seq = str(record.seq)
        header = sp + '|' + geneName + '|' + contig + '|' + gene_start + '|' + gene_end + '|'
        genome_hbx_seqs[header].append(seq)
        
##Parse cluster file and output nucleotide sequences
count=0
gene_nuc_data = {}
os.mkdir('temp_seqs/' + args.gene)
with open('../recip_blast/hbx_clusters/' + args.gene + '_cluster.txt') as f, open(args.gene + '_HD_nuc.fasta', 'w') as outF:
    lines = f.read()
    sp_homeobox = lines.split('>')
    for cluster in sp_homeobox:
        if not cluster:
            continue
        else:
            genes = cluster.split('\n')
            spName = genes[0].split('_i')[0]
            for gene in genes[1:]:
                if not gene:
                    continue
                else:
                    gene_info = gene.split('\t')
                    contig = gene_info[0]
                    gene_start = gene_info[1].split(',')[0].strip('(')
                    gene_end = gene_info[1].split(', ')[1].strip(')')
                    gene_end = int(gene_end) + 1
                    gene_end = str(gene_end)
                    geneName = gene_info[2]
                    
                    #Get the nucleotide sequence for the gene from the recip blast fasta file
                    gene_header_search = spName + '|' + geneName + '|' + contig + '|' + gene_start + '|' + gene_end + '|'
                    #print(gene_header_search)
                    gene_nuc_seq = genome_hbx_seqs[gene_header_search]
                    #print(gene_nuc_seq)
                    gene_header = spName + '_' + geneName + '_' + contig + '_' + gene_start + '_' + gene_end
                    try:
                        outF.write('>' + gene_header + '\n' + gene_nuc_seq[0] + '\n')
                    except:
                        #print(gene_header)
                        count+=1
                        continue
                    if '_' in contig:
                        contig = contig.split('_')[0]
                    
                    if geneName in intron_list:
                        if len(set(gene_nuc_seq)) > 1:
                            with open('temp_seqs/' + args.gene + '/' + spName + '_' + geneName + '_' + contig + '_' + gene_start + '_' + gene_end + '.1.fasta', 'w') as outF1:
                                outF1.write('>' + gene_header + '\n' + gene_nuc_seq[0] + '\n')
                            with open('temp_seqs/' + args.gene + '/' + spName + '_' + geneName + '_' + contig + '_' + gene_start + '_' + gene_end + '.2.fasta', 'w') as outF2:
                                outF2.write('>' + gene_header + '\n' + gene_nuc_seq[1] + '\n')
                        else:
                            with open('temp_seqs/' + args.gene + '/' + spName + '_' + geneName + '_' + contig + '_' + gene_start + '_' + gene_end + '.fasta', 'w') as outF1:
                                outF1.write('>' + gene_header + '\n' + gene_nuc_seq[0] + '\n')
                    else:
                        with open('temp_seqs/' + args.gene + '/' + spName + '_' + geneName + '_' + contig + '_' + gene_start + '_' + gene_end + '.fasta', 'w') as outF3:
                            outF3.write('>' + gene_header + '\n' + gene_nuc_seq[0] + '\n')
#print(count, 'genes with naming err')
##Translate nucleotide hbx sequences
os.chdir('temp_seqs/' + args.gene)
for fasta in glob.glob("*.fasta"):
    unix('sixpack -sequence ' + fasta + ' -outfile ' + fasta + '.sixpack -outseq ' + fasta + '.sixpack.fa', shell=True)

outF = open(args.gene + '_HD_AA.fasta', 'w')
sorted_fasta = []
for fa in glob.glob("*fa"):
    sorted_fasta.append(fa)
sorted_fasta = sorted(sorted_fasta)

for fa in sorted_fasta:
    fa_info = fa.split('.fasta')[0]
    if len(fa_info.split('_')) == 6:
        spName = '_'.join(fa_info.split('_')[:1])
        geneID = fa_info.split('_')[2]
        contig = fa_info.split('_')[3]
        pos = '_'.join(fa_info.split('_')[4:5])
    else:
        #print(fa)
        spName = '_'.join(fa_info.split('_')[:1])
        geneID = fa_info.split('_')[3]
        contig = fa_info.split('_')[4]
        pos = '_'.join(fa_info.split('_')[5:6])
        
    AA_seq = hbx_gene_seq_dict[geneID]
    sp_gene_seq_idents = {}
    with open(fa) as f:
        for record in SeqIO.parse(f, 'fasta'):
            seq = str(record.seq)
            alignments = aligner.align(AA_seq, seq)
            
            alignment = alignments[0]
            align_score = alignment.score
            
            sp_gene_seq_idents[seq] = align_score

    best_sp_gene_seq = max(sp_gene_seq_idents.values())
    for gene_seq, score in sp_gene_seq_idents.items():
        if score == best_sp_gene_seq:
            if gene_seq[-1] == 'X':
                outF.write('>' + fa_info + '\n' + gene_seq[:-1] + '\n')
            else:
                outF.write('>' + fa_info + '\n' + gene_seq[:-1] + '\n')
                
unix('mv ' + args.gene + '_HD_AA.fasta ../../', shell=True)
os.chdir('../../')



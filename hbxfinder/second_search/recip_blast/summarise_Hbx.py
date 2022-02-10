#!/usr/bin/env python3
import argparse
from collections import defaultdict, OrderedDict
import json
import os
import glob

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 summarise_Hbx.py --taxa <files_to_consider> --gene <Hbx gene class>


parse = argparse.ArgumentParser()

parse.add_argument("--taxa",type=str, help="name of species HOX results to parse",required=True)
parse.add_argument("--gene",type=str, help="name of homeobox gene to obtain results for",required=True)
parse.add_argument("--path",type=str, help="path to genome fasta files",required=True)

args = parse.parse_args()

GCA2sp = {}
for fasta in glob.glob(args.path + '*fasta'):
    genome = fasta.split('/')[-1]
    GCA = genome.split('.')[0]
    sp = genome.split('.')[1].split('_')[1]
    GCA2sp[GCA] = sp
    
sp_file = args.taxa.split('/')[-1]
sp_name = sp_file.split('_' + args.gene)[0]
species = GCA2sp[sp_name]
print('Parsing', species, args.gene, 'gene content...')

##Open dictionary of Homeodomain classes and their genes
with open('../../../raw/hbx_data/hbx_naming.json') as f:
    hbx_naming_dict = json.load(f)


##Parse blastoutput to get gene positions
seq_ranges = []
seq_range_perc_ident = defaultdict(list)
seq_range_sub_gene = {}
contig_seq_ranges = defaultdict(list)
contig_seq_range_sub_gene = defaultdict(dict)
contig_seq_range_perc_ident = defaultdict(dict)
with open(args.taxa) as f:
    for line in f:
        lines = line.split('\t')
        query = lines[0]
        perc_ident = lines[3]
        subject = lines[1]
        subject_gene = subject.split('|')[1]

        query_info = query.split('|')
        contig = query_info[1]
        
        if int(query_info[2]) < int(query_info[3]):
            seq_start = int(query_info[2])
            seq_end = int(query_info[3])
        else:
            seq_start = int(query_info[3])
            seq_end = int(query_info[2])
            
        seq_range = range(seq_start, seq_end)
        seq_ranges.append(seq_range)
        contig_seq_ranges[contig].append(seq_range)
        seq_range_sub_gene[seq_range] = subject_gene
        seq_range_perc_ident[seq_range] = perc_ident
        
        contig_seq_range_sub_gene[contig] = seq_range_sub_gene
        contig_seq_range_perc_ident[contig] = seq_range_perc_ident


##Remove overlapping hits for each hbx in the class by retaining largest hit - write gene locations to file
hbx_gene_dict = hbx_naming_dict[args.gene]

os.makedirs('hbx_clusters', exist_ok=True)
outF = open('hbx_clusters/' + args.gene + '_cluster.txt', 'a+')
sp_assem = args.taxa.split('/')[-1].split('_' + args.gene)[0]
species_name = GCA2sp[sp_assem]
#print(sp_assem)
outF.write('>' + sp_assem + '_' + species_name + '\n')
for contig, seq_range_list in contig_seq_ranges.items():

    seq_ranges_dict = defaultdict(list)
    for seq_range in set(seq_range_list):
        seq_ranges_dict[len(seq_range)].append(seq_range)
        

    ordered_range_list = []
    seq_ranges_dict = OrderedDict(sorted(seq_ranges_dict.items()))
    for k, v in reversed(seq_ranges_dict.items()):
        for range_seqs in v:
            ordered_range_list.append(range_seqs)

    subset_gene_list = []
    for i in range(0, len(ordered_range_list)):
        overlap = False
        range_i = ordered_range_list[i]
        
        for j in range(0, len(ordered_range_list)):
            range_j = ordered_range_list[j]
            if range_i != range_j:
                set_range_i = set(range_i)
                
                if set_range_i.intersection(range_j):
                    overlap = True
                    range_longest = max(range_i, range_j, key=len)
                    if range_longest not in subset_gene_list:
                        subset_overlap = False
                        for ranges in subset_gene_list:
                            if set(range_longest).intersection(ranges):
                                subset_overlap = True
                        if subset_overlap == False:
                            subset_gene_list.append(range_longest)
                        
        if overlap == False:
            subset_gene_list.append(range_i)

    for ranges in sorted(subset_gene_list, key=lambda r: r.start):
        range_genes = contig_seq_range_sub_gene[contig]
        range_perc = contig_seq_range_perc_ident[contig]
        range_gene = range_genes[ranges]

        if range_genes[ranges] in hbx_gene_dict:
            geneID = range_genes[ranges]
            hbx_geneID = hbx_gene_dict[geneID]
            #print(contig, ranges, len(ranges), hbx_geneID, range_perc[ranges])
            gene_start = ranges[0]
            gene_end = ranges[-1]
            gene_positions = gene_start, gene_end
            outF.write(contig + '\t' + str(gene_positions) + '\t' + hbx_geneID + '\t' + str(range_perc[ranges]) + '\n') 
        
        
outF.close()

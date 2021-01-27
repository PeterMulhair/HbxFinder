#!/usr/bin/env python3
from Bio import SearchIO
from collections import defaultdict
import argparse
import json

# Author: Peter Mulhair
# Date: 15/01/2020
# Usage python3 parse_blast_xml.py --taxa <xml blastoutput file> --gene <homeobox gene to search for>

parse = argparse.ArgumentParser()

parse.add_argument("--taxa",type=str, help="name of species blast file in xml format",required=True)
parse.add_argument("--gene",type=str, help="name of homeobox gene to parse",required=True)

args = parse.parse_args()

species_name = args.taxa.split('_' + args.gene)[0]
print(species_name)


##Open dictionary of Homeodomain classes and their genes
with open('../../../hbx_data/hbx_naming.json') as f:
    hbx_naming_dict = json.load(f)


hbx_dict = hbx_naming_dict[args.gene]
Homeobox_list = hbx_dict.values()

intron_hbx_genes = ['Pb', 'Ro']

sp_gene_ranges_dict = {}
for homeo_gene in Homeobox_list:
    
    with open('../recip_blast/hbx_clusters/' + args.gene + '_cluster.txt') as f:
        lines = f.read()
        sp_homeobox = lines.split('>')
        for cluster in sp_homeobox: 
            if not cluster:
                continue
            else:
                genes = cluster.split('\n')
                spName = genes[0]
                for gene in genes[1:]:
                    if not gene:
                        continue
                    else:
                        gene_info = gene.split('\t')
                        contig = gene_info[0]
                        gene_start = gene_info[1].split(',')[0].strip('(')
                        gene_end = gene_info[1].split(',')[1].strip(')')
                        geneName = gene_info[2]

                        if spName == species_name:
                            if geneName == homeo_gene:
                                sp_contig = contig
                                sp_gene_start = gene_start.strip()
                                sp_gene_end = gene_end.strip()
                                sp_geneName = geneName
                                sp_geneName_range = sp_geneName + '_' + sp_gene_start + '_' + sp_gene_end

                                sp_gene_ranges_dict[sp_geneName_range] = sp_contig
                                #print(contig, gene_start, gene_end, geneName)
                                
                        else:
                            break

contig_blast_hit = defaultdict(dict)
for gene_name_ranges, contig in sp_gene_ranges_dict.items():
    blast_hit_dict = defaultdict(list)
    homeo_geneID = gene_name_ranges.split('_')[0]
    sp_gene_start = gene_name_ranges.split('_')[1]
    sp_gene_end = gene_name_ranges.split('_')[2]
    
    results = args.taxa

    blast_Result = SearchIO.parse(results, 'blast-xml')
    for res in blast_Result:
        query_ID = res.id
        gene_ID = query_ID.split('|')[1]
        if gene_ID in hbx_dict:
            homeo_ID = hbx_dict[gene_ID]

            if homeo_ID == homeo_geneID:
                sp = query_ID.split('|')[0]
                for hsp in res:              
                    hit_contig = hsp.id

                    if hit_contig == contig:
                        for hits in hsp:
                            hit_start = hits.hit_start
                            hit_start_range = range(int(hit_start) - 5000, int(hit_start) + 5000)
                            hit_end = hits.hit_end
                            hit_end_range =range(int(hit_end) - 5000, int(hit_end) + 5000)

                            hit_seq = hits.hit.seq
                            
                            if (int(sp_gene_start) in hit_start_range) and (int(sp_gene_end) in hit_end_range):
                                blast_hit_dict[gene_name_ranges].append(hit_seq)

                        contig_blast_hit[contig].update(blast_hit_dict)
                                


with open(args.gene + '_AA_seq.fasta', 'a+') as outF:
    for contig, blast_hit in contig_blast_hit.items():
        for genName_range, AA_seq in blast_hit.items():
            hbx_name = genName_range.split('_')[0]
            header = species_name + '|' + hbx_name + '|' + '_'.join(genName_range.split('_')[1:]) + '|' + contig

            outF.write('>' + header + '\n')

            if len(set(AA_seq)) > 1:
                seq_list = []
                for seqs in set(AA_seq):
                    if len(seqs) == 60:
                        seq_list.append(seqs)

                if len(seq_list) < 1:
                    if hbx_name in intron_hbx_genes:
                        for AA_seq_exons in set(AA_seq):
                            outF.write(str(AA_seq_exons) + '\n')
                    else:
                        seq_max = max(set(AA_seq), key=len)
                        outF.write(str(seq_max) + '\n')

                else:
                    for AA_seqs in seq_list:
                        outF.write(str(AA_seqs) + '\n')


            else:
                for seq in set(AA_seq):
                    outF.write(str(seq) + '\n')


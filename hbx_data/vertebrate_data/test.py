import json
import glob
from collections import defaultdict

class_gene_names = defaultdict(dict)
for fas in glob.glob('family_data/*fasta'):
    assign_gene_name = {}
    hbx_class = fas.split('/')[-1].split('.')[0]
    with open(fas) as f:
        for line in f:
            if line.startswith('>'):
                line = line.strip('\n')
                lines = line.split('|')
                gene_name = lines[1]
                assign_name = lines[2]
                if assign_name == 'unassigned':
                    assign_gene_name[gene_name] = gene_name
                else:
                    assign_gene_name[gene_name] = assign_name

    class_gene_names[hbx_class] = assign_gene_name

#print(class_gene_names['Other'])
with open('hbx_naming.json', 'w') as outf:
    json.dump(class_gene_names, outf)

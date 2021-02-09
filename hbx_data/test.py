import json
from Bio import SeqIO

#hbx_dict = {'HOX': {'abd-A': 'abd-A', 'Abd-A': 'abd-A', 'Abd-B': 'Abd-B', 'Antp': 'Antp', 'Dfd': 'Dfd', 'eve': 'eve', 'Eve': 'eve', 'ftz': 'ftz', 'Ftz': 'Ftz', 'lab': 'lab', 'Lab': 'lab', 'pb': 'Pb', 'Pb': 'Pb', 'ptl': 'Antp', 'ro': 'Ro', 'Ro': 'Ro', 'Scr': 'Scr', 'Ubx': 'Ubx', 'Utx': 'Ubx', 'zen': 'zen', 'Zen': 'zen', 'zen1': 'zen', 'zen2': 'zen', 'ShxA': 'ShxA', 'ShxB': 'ShxB', 'ShxC': 'ShxC', 'ShxD': 'ShxD'}, 'NK': {'C15': 'Tlx', 'Dr2': 'Msx', 'Dr': 'Msx', 'Abox': 'Abox', 'Bap': 'NK3', 'Bari': 'Bari', 'Ems': 'Emx', 'Hgtx': 'NK6', 'lbe': 'Lbx', 'Lbx': 'Lbx', 'NK7': 'NK7', 'Prrx': 'Prrx', 'slou': 'NK1', 'Slou': 'NK1', 'Tin': 'NK4', 'Vnd': 'NK2', 'Hmx': 'Hmx'},'PRD': {'al': 'Arx', 'Pph13': 'Arx', 'Al': 'Arx', 'CG34340': 'Drgx', 'hbn': 'Hbn', 'otp': 'Otp', 'oc': 'Otx', 'Oc1': 'Otx', 'Oc2': 'Otx', 'sv': 'Pax2/5/8', 'Sv': 'Pax2/5/8', 'gsb': 'Pax3/7', 'gsb-n': 'Pax3/7', 'prd': 'Pax3/7', 'Gsb': 'Pax3/7', 'Gsbn': 'Pax3/7', 'Prd': 'Pax3/7', 'ey': 'Pax4/6', 'eyg': 'Pax4/6', 'toe': 'Pax4/6', 'toy': 'Pax4/6', 'Ey': 'Pax4/6', 'Eyg': 'Pax4/6', 'Toy': 'Pax4/6', 'PHDP': 'Phox', 'Ptx1': 'Pitx', 'CG32532': 'Prop', 'CG9876': 'Prop', 'Rx': 'Rax', 'repo': 'Repo', 'CG34367': 'Shox', 'OdsH': 'Uncx', 'unc-4': 'Uncx', 'Unc4': 'Uncx', 'tup': 'Vsx', 'Vsx1': 'Vsx', 'Vsx2': 'Vsx'}, 'LIM': {'Tup': 'Isl', 'Lim1': 'Lhx1/5', 'ap': 'Lhx2/9', 'Ap1': 'Lhx2/9', 'Ap2': 'Lhx2/9', 'Lim3': 'Lhx3/4', 'Awh': 'Lhx6/8', 'Awh1': 'Lhx6/8', 'Awh2': 'Lhx6/8', 'CG32105': 'Lmx', 'CG4328': 'Lmx', 'Lmxa': 'Lmx', 'Lmxb': 'Lmx'}, 'POU': {'nub': 'Pou2', 'pdm2': 'Pou2', 'Nub': 'Pou2', 'vvl': 'Pou3', 'Vvl': 'Pou3', 'acj6': 'Pou4', 'Acj6': 'Pou4', 'Acj6-l': 'Pou4', 'pdm3': 'Pou6', 'Pdm3a': 'Pou6', 'Pdm3b': 'Pou6', 'Pdm3': 'Pou6'}, 'CUT': {'dve': 'Cmp', 'Dve': 'Cmp', 'ct': 'Cux', 'Ct': 'Cux', 'onecut': 'Onecut'}, 'TALE': {'ara': 'Irx', 'caup': 'Irx', 'mirr': 'Irx', 'Ara': 'Irx', 'Mirr': 'Irx', 'hth': 'Meis', 'Hth': 'Meis', 'CG11617': 'Mkx', 'exd': 'Pbx', 'Exd': 'Pbx', 'achi': 'Tgif', 'vis': 'Tgif', 'Tgif1': 'Tgif', 'Tgif2': 'Tgif'}, 'SINE': {'so': 'Six1/2', 'So': 'Six1/2', 'Optix': 'Six3/6', 'Six4': 'Six4/5'}, 'ZF': {'Zfh1': 'Zeb', 'zfh1': 'Zfhx', 'zfh2': 'Zfhx', 'Zfh2': 'Zfhx'}, 'PROS': {'pros': 'Prox', 'Pros': 'Prox'}, 'CERS': {'Lag1': 'Cers', 'Lag1a': 'Cers', 'Lag1b': 'Cers'}}

#with open('hbx_naming.json', 'w') as fp:
#    json.dump(hbx_dict, fp)


with open('hbx_naming.json') as f:
    hbx_dict = json.load(f)

homeobox_dict = {}
with open('homeobox.fasta') as f:
    for record in SeqIO.parse(f, 'fasta'):
        ID = record.description
        seq = str(record.seq)
        homeobox_dict[ID] = seq


for k, v in hbx_dict.items():
    print(k)
    with open('family_data/' + k + '.fasta', 'w') as outF:
        for gene in set(v.values()):
            for genes, seqs in homeobox_dict.items():
                if '|' + gene + '|' in genes:
                    outF.write('>' + genes + '\n' + seqs + '\n')
    

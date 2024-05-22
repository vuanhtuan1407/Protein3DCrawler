import json
import math

import numpy as np
import pandas as pd
from Bio import PDB, SeqIO
from tqdm import tqdm

from crawler import SP6_PATH, load_cache


def extract_fasta(fasta_file=SP6_PATH):
    prots = []
    records = SeqIO.parse(fasta_file, 'fasta')
    for record in records:
        prot_id, organism, label, partition = str(record.id).split('|')
        prots.append({
            'prot_id': prot_id,
            'organism': organism,
            "label": label,
            "partition": partition
        })
    return prots


def extract_pdb(prot_id, use_lib=True):
    if use_lib:
        prot = []
        parser = PDB.PDBParser()
        structure = parser.get_structure(prot_id, f'./data/pdb/AF-{prot_id}-F1-model_v4.pdb')
        for model in structure:
            for chain in model:
                for residue in chain:
                    ca = residue['CA']
                    prot.append({
                        "residue": residue.resname,
                        "coord": ca.coord.tolist()
                    })
                    if len(prot) == 70:
                        break
        return prot
    else:
        with open(f'./data/pdb/AF-{prot_id}-F1-model_v4.pdb', 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith('ATOM'):
                    line = line.rstrip().split()[:70]

        # Pending...


def build_3d_dataset():
    prots = extract_fasta()
    done_ids, _ = load_cache()
    for prot in tqdm(prots):
        prot_id = prot['prot_id']
        if prot_id in done_ids:
            prot['prot_3d'] = extract_pdb(prot_id=prot_id)
        else:
            prot['prot_3d'] = []
    with open('./data/3d_dataset.json', 'w') as f:
        json.dump(prots, f)


def is_closed_enough(residue1, residue2):
    coord1 = residue1['coord']
    coord2 = residue2['coord']

    dis = math.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (coord1[2] - coord2[2]) ** 2)
    return dis < 9.0


def extract_3d_dataset():
    PROT3D_PATH = './data/3d_dataset.json'
    AA_INFO_PATH = './data/smiles_string_aa.csv'
    df = pd.read_csv(AA_INFO_PATH)
    prots = []
    with open(PROT3D_PATH, 'r') as f:
        data = json.load(f)
        for prot in data:
            from_list, to_list = [], []
            adj_matrix = np.zeros((20, 20), dtype=np.int8)
            if len(prot['prot_3d']) != 0:
                prot_3d = prot['prot_3d']
                for i, r1 in enumerate(prot_3d[:-1]):
                    for j, r2 in enumerate(prot_3d[i + 1:]):
                        ret = is_closed_enough(r1, r2)
                        if ret:
                            idx1 = df.index[df['abbreviations'] == r1['residue']].tolist()[0]
                            idx2 = df.index[df['abbreviations'] == r2['residue']].tolist()[0]
                            from_list.append(idx1)
                            to_list.append(idx2)
                            from_list.append(idx2)
                            to_list.append(idx1)
                            adj_matrix[idx1, idx2] = 1
                            adj_matrix[idx2, idx1] = 1

                prots.append({
                    "prot_id": prot['prot_id'],
                    "kingdom": prot['organism'],
                    "label": prot['label'],
                    "from_list": from_list,
                    "to_list": to_list,
                    "adj_matrix": adj_matrix.tolist()
                })

    with open(PROT3D_PATH, 'w') as f:
        json.dump(prots, f)


if __name__ == '__main__':
    # build_3d_dataset()

    extract_3d_dataset()

    # P82290
    # A0A0A1I6E7
    # prot = extract_pdb('P82290')
    # with open('./test.json', 'w') as f:
    #     json.dump(prot, f)

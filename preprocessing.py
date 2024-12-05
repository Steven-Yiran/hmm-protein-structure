from Bio import SeqIO
from Bio.PDB import PDBParser, DSSP
import numpy as np

# Amino acid to integer mapping
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_to_int = {aa: idx for idx, aa in enumerate(amino_acids)}

# Function to extract sequences and secondary structures from PDB files
def parse_pdb(pdb_id, pdb_path):
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, f"{pdb_path}/{pdb_id}.pdb")
    model = structure[0]
    dssp = DSSP(model, f"{pdb_path}/{pdb_id}.pdb")

    seq = ''
    sec_struc = ''
    for key in dssp.keys():
        aa = dssp[key][1]
        ss = dssp[key][2]
        if aa in amino_acids:
            seq += aa
            if ss == 'H':
                sec_struc += 'H'
            elif ss == 'E':
                sec_struc += 'E'
            elif ss in ['T', 'S', 'G']:
                sec_struc += 'T'
            else:
                sec_struc += 'C'  # Coil or other

    return seq, sec_struc

# Function to encode sequences and secondary structures
def encode_sequences(sequences, secondary_structures):
    ss_to_states = {
        'H': [1, 2, 3, 4],
        'E': [5, 6, 7],
        'T': [8, 9, 10],
        'C': [11, 12, 13]
    }

    encoded_sequences = []
    encoded_state_sequences = []
    sequence_lengths = []

    for seq, ss in zip(sequences, secondary_structures):
        obs_seq = [aa_to_int[aa] for aa in seq]
        state_seq = []
        for s in ss:
            state_seq.extend(ss_to_states[s])
        encoded_sequences.append(obs_seq)
        encoded_state_sequences.append(state_seq)
        sequence_lengths.append(len(obs_seq))

    return encoded_sequences, encoded_state_sequences, sequence_lengths
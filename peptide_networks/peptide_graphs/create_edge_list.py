import pandas as pd
import os
from Levenshtein import distance
import argparse
from tqdm import tqdm
import re
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from quantiprot.metrics.aaindex import get_aa2volume

aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

def levenshtein_distance(seq_1, seq_2):
    return distance(seq_1, seq_2)

def blosum_distance(seq1, seq2):
    return aligner.align(seq1, seq2).score
       
def biophysical_distance(seq1,seq2):
    """ Returns the weighted euclidean distance of biophysical parameters of peptides """
    pp_seq1 = ProteinAnalysis(seq1)
    pp_seq2 = ProteinAnalysis(seq2)

    # peptide parameters
    ip_seq1 = pp_seq1.isoelectric_point()
    ip_seq2 = pp_seq2.isoelectric_point()
    helix_fraction_seq1 = pp_seq1.secondary_structure_fraction()[0]
    helix_fraction_seq2 = pp_seq2.secondary_structure_fraction()[0]
    hp_seq1 = pp_seq1.gravy()
    hp_seq2 = pp_seq2.gravy()
    mw_seq1 = pp_seq1.molecular_weight()
    mw_seq2 = pp_seq2.molecular_weight()
    vol_seq1 =  sum(get_aa2volume(pp_seq1).mapping.values())
    vol_seq2 = sum(get_aa2volume(pp_seq2).mapping.values())
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)

    # coefficients
    k_len = 1
    k_hp = 1
    k_ip = 1
    k_mw = 1
    k_helix = 1
    k_vol = 1 
    
    # weighted normalized euclidean_distance
    def normalized_euclidean(v1, v2):
        return sum(((p-q)/(p+q))**2 if p != 0 and q != 0 else 0 for p, q in zip(v1, v2)) ** .5
    
    vec1 = [ip_seq1, helix_fraction_seq1, hp_seq1, mw_seq1, vol_seq1, len_seq1]
    vec2 = [ip_seq2, helix_fraction_seq2, hp_seq2, mw_seq2, vol_seq2, len_seq2]
    k_vec = [k_ip, k_helix, k_hp, k_mw, k_vol, k_len]
    w_vec1 = np.multiply(vec1, k_vec)
    w_vec2 = np.multiply(vec2, k_vec)
    result =  normalized_euclidean(w_vec1, w_vec2)
    return result

def amino_acid_distance(seq1,seq2):
    
    pp_seq1 = ProteinAnalysis(seq1)
    pp_seq2 = ProteinAnalysis(seq2)
    aa_seq1 = list(pp_seq1.get_amino_acids_percent().values())
    aa_seq2 = list(pp_seq2.get_amino_acids_percent().values())
    def euclidean(v1, v2):
        return sum((p-q)**2 if p != 0 and q != 0 else 0 for p, q in zip(v1, v2)) ** .5
    return euclidean(aa_seq1, aa_seq2) 

def create_adjacency_matrix(df, matrix):
    """
    Creates an adjanceny matrix. 
    This can probably be made quicker.
    """
    similarity_matrix = []
    column_and_index = []
    Peptides = {'aa':df['Peptide'].values, 'Area': df['Area'].values, 'Accession':df['Accession'].values}
    for i, seq1 in tqdm(zip(range(len(Peptides['aa'])), Peptides['aa']), total=len(Peptides['aa'])):
        similarity_for_seq1 = []
        area_from = Peptides['Area'][i]
        accession_from = Peptides['Accession'][i]
        column_and_index.append(seq1)
        for k in range(i+1):
            seq2 = Peptides['aa'][k]
            area_to = Peptides['Area'][k]
            accession_to = Peptides['Accession'][k]
            if(matrix=='levenshtein'):
                distance = levenshtein_distance(seq1,seq2)
            elif(matrix=='biophysical'):
                distance = biophysical_distance(seq1,seq2)
            elif(matrix=='blosum'):
                distance = blosum_distance(seq1,seq2)
            elif(matrix=='aa'):
                distance = amino_acid_distance(seq1,seq2)
            similarity_for_seq1.append((distance,area_from,area_to,accession_from, accession_to))    
        similarity_matrix.append(similarity_for_seq1)

    similarity_df = pd.DataFrame(similarity_matrix)
    similarity_df.index=column_and_index
    similarity_df.columns=column_and_index
    return similarity_df

def adjacency_matrix_to_edge_list(adjacency_matrix):
    return adjacency_matrix.stack().reset_index()


def main(args):
    filepath = args.filepath
    matrix = args.matrix
    sample_threshold = args.threshold
    df = pd.read_excel(filepath, engine='openpyxl')
    df = df.sample(frac = sample_threshold, random_state = 42)
    df.dropna(subset=['Accession'], inplace=True) # remove unidentifiable peptides
    df['Peptide'] = df['Peptide'].apply(lambda x: re.sub("[^a-zA-Z]+", "", x))
    df['Accession'] = df['Accession'].apply(lambda x: str(x).split(':')[0])
    df['Accession'] = df['Accession'].apply(lambda x: str(x).split('|')[2]) # most prominent protein
    
    area_col = [col for col in df.columns if col.startswith('Area')]
    df.rename(columns={area_col[0]:'Area'}, inplace=True)
    
    df = df[['Peptide','Area','Accession']]
    df.replace(0, np.nan, inplace=True)
    df.dropna(inplace=True)
    df['Area'] = df['Area'].apply(lambda x: np.log10(x))

    adjacency_matrix = create_adjacency_matrix(df, matrix)
    edge_list = adjacency_matrix_to_edge_list(adjacency_matrix)
    
    edge_list.columns = ['from', 'to','distance,area_from,area_to,accession_from,accession_to']
    edge_list[['distance', 'area_from', 'area_to', 'accession_from','accession_to']] = pd.DataFrame(edge_list['distance,area_from,area_to,accession_from,accession_to'].tolist(), index=edge_list.index)
    edge_list.drop(columns=['distance,area_from,area_to,accession_from,accession_to'], inplace=True)
    file = re.sub('\D', '', filepath)
    print(edge_list)
    edge_list_file_name = f'data/edge_lists/{matrix}/{file}_edge_list.gz'
    print(f'Writing: {edge_list_file_name}')
    edge_list.to_csv(edge_list_file_name, index=False, compression='gzip')    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    parser.add_argument("matrix", type=str, choices=['blosum', 'biophysical', 'levenshtein', 'aa'], default='levenshtein', help="distance matrix")
    parser.add_argument("threshold", type=float, default=1, help="subset fraction")
    args = parser.parse_args()
    main(args)

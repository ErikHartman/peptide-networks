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
from Bio.SeqUtils.ProtParam import ProtParamData
from quantiprot.metrics.aaindex import get_aa2volume, get_aa2hydropathy



def levenshtein_distance(seq_1, seq_2):
    return distance(seq_1, seq_2)

def blosum_distance(seq1, seq2):
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    return aligner.align(seq1, seq2).score
       
def custom_distance(seq1,seq2):
    """ Returns the weighted euclidian distance of biophysical parameters of peptides """
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
    
    # weighted euclidian_distance
    vec1 = [ip_seq1, helix_fraction_seq1, hp_seq1, mw_seq1, vol_seq1, len_seq1]
    vec2 = [ip_seq2, helix_fraction_seq2, hp_seq2, mw_seq2, vol_seq2, len_seq2]
    k_vec = [k_ip, k_helix, k_hp, k_mw, k_vol, k_len]
    euclidian_distance = np.linalg.norm(vec1-vec2)
    result = np.multiply(euclidian_distance,k_vec)

    return result

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
        area = Peptides['Area'][i]
        accession = Peptides['Accession'][i]
        column_and_index.append(seq1)
        for k in range(i+1):
            seq2 = Peptides['aa'][k]
            if(matrix=='levenshtein'):
                distance = levenshtein_distance(seq1,seq2)
            elif(matrix=='biophysical'):
                distance = custom_distance(seq1,seq2)
            elif(matrix=='blosum'):
                distance = blosum_distance(seq1,seq2)
            similarity_for_seq1.append((distance,area,accession)) # accession and area based on "from"        
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
    df = pd.read_excel(filepath, engine='openpyxl')
    
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
    
    edge_list.columns = ['from', 'to','distance,area,accession']
    edge_list[['distance', 'area', 'accession']] = pd.DataFrame(edge_list['distance,area,accession'].tolist(), index=edge_list.index)
    edge_list.drop(columns=['distance,area,accession'], inplace=True)
    file = re.sub('\D', '', filepath)
    edge_list_file_name = f'data/edge_lists/{matrix}/{file}_edge_list.gz'
    print(f'Writing: {edge_list_file_name}')
    edge_list.to_csv(edge_list_file_name, index=False, compression='gzip')    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create edge list from file based on distance matrix')
    parser.add_argument("filepath", type=str, help="Path to file")
    parser.add_argument("matrix", type=str, choices=['blosum', 'biophysical', 'levenshtein'], default='levenshtein', help="distance matrix")
    args = parser.parse_args()
    main(args)

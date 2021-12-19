import pandas as pd
import os
from Levenshtein import distance
import argparse
from tqdm import tqdm
import re
import numpy as np

"""
Creates edge lists for wanted matrices.
"""
parser = argparse.ArgumentParser()
parser.add_argument("filepath", type=str, help="Path to file")
args = parser.parse_args()
filepath = args.filepath


def levenshtein_distance(str_1, str_2):
    return distance(str_1, str_2)

def create_adjacency_matrix(df):
    """
    Creates adjacency matrix using the Levenshtein algorithm. 
    Only computes the upper triangle to reduce complexity.
    """
    similarity_matrix = []
    column_and_index = []
    Peptides = {'aa':df['Peptide'].values, 'Area': df['Area'].values, 'Accession':df['Accession'].values}
    for i, seq1 in tqdm(zip(range(len(Peptides['aa'])), Peptides['aa']), total=len(Peptides['aa'])):
        similarity_for_seq1 = []
        column_and_index.append(seq1)
        for k in range(i+1):
            similarity_for_seq1.append((levenshtein_distance(seq1,Peptides['aa'][k]), Peptides['Area'][i], Peptides['Accession'][i]))
        similarity_matrix.append(similarity_for_seq1)
    
    similarity_df = pd.DataFrame(similarity_matrix)
    similarity_df.index=column_and_index
    similarity_df.columns=column_and_index
    return similarity_df


def main():
    df = pd.read_excel(filepath, engine='openpyxl')
    df.dropna(subset=['Accession'], inplace=True)
    df['Peptide'] = df['Peptide'].apply(lambda x: re.sub("[^a-zA-Z]+", "", x))
    df['Accession'] = df['Accession'].apply(lambda x: str(x).split(':')[0])
    df['Accession'] = df['Accession'].apply(lambda x: str(x).split('|')[2])
    
    area_col = [col for col in df.columns if col.startswith('Area')]
    df.rename(columns={area_col[0]:'Area'}, inplace=True)
    df = df[['Peptide','Area','Accession']]
    df.replace(0, np.nan, inplace=True)
    df.dropna(inplace=True)
    df['Area'] = df['Area'].apply(lambda x: np.log10(x))
    
    similarity_df = create_adjacency_matrix(df)
    similarity_df = similarity_df.stack().reset_index()
    similarity_df.columns = ['from', 'to','distance,Area,Accession']
    similarity_df[['distance', 'Area', 'Accession']] = pd.DataFrame(similarity_df['distance,Area,Accession'].tolist(), index=similarity_df.index)
    similarity_df.drop(columns=['distance,Area,Accession'], inplace=True)
    file = filepath.split('/')[-1]
    edge_file_name = f'./edge_lists/{file}_edge_list.csv'
    print(f'Writing: {edge_file_name}')
    similarity_df.to_csv(edge_file_name)    

if __name__ == "__main__":
    main()

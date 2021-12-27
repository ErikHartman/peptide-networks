import pandas as pd
from peptide_networks.preprocessing import hash_util

def main():
    df = pd.read_csv('data/temp/peptide_sample_31.txt', sep='\t')
    clean_df= pd.read_csv('data/clean_data/peptide_sample_inf/peptide_sample_31.csv')
    df.columns=['size_hash', 'UniProt ID','Start','End']
    df['size'] = df['size_hash'].apply(lambda x: int(x.split('_')[0]))
    df['hash'] = df['size_hash'].apply(lambda x: int(x.split('_')[1]))
    merged_df = df.merge(clean_df, on=['size','hash'], how='inner')
    merged_df['Intensity s31'] = merged_df['Area']
    merged_df = merged_df[['Peptide','UniProt ID','Start','End','Intensity s31']]
    merged_df.to_csv('data/temp/peptigram_input_31.csv', index=False)
    

if __name__ == '__main__':
    main()
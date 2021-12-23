import pandas as pd
import glob
import re
import os
from hash_util import HashString


def main():
    PATH = '../../data/raw_data/'
    CLEAN_PATH = '../../data/clean_data/'
    columns_to_keep = ['Peptide', '-10lgP', 'ppm', 'Accession','Length']
    for filepath in glob.glob(PATH + '*/*'):
        folder = filepath.split('\\')[-2]
        filenm = filepath.split('\\')[-1]
        if not os.path.exists(CLEAN_PATH + folder):
            os.mkdir(CLEAN_PATH + folder)
        
        df = pd.read_excel(filepath)
        area_column_name = list(filter(lambda v: re.match('Area*', v), df.columns))
        df['Peptide'] = df['Peptide'].apply(lambda x: re.sub(r'[()+-0123456789]',"",x))
        df = df[df.columns.intersection(columns_to_keep + area_column_name)]
        df.rename(columns = {area_column_name[0] :'Area', 'Length':'size'}, inplace=True)
        
        df['hash'] = df['Peptide'].apply(lambda x: HashString(x).getHash())
        
        df = df.set_index(['size', 'hash'])
        df.sort_values(['size', 'hash'], inplace=True)
        
        df.to_csv(CLEAN_PATH + folder + '/' + filenm[:-4] + 'csv')
        print(f'{filenm} done')
        
        
        


if __name__ == '__main__':
    main()
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 09:15:09 2021

@author: matti
"""

import pandas as pd
import glob
import re
from pseAAC import GetPseudoAAC

if __name__ == '__main__':
    AAP = ['Hydrophobicity','hydrophilicity', 'pK1', 'pK2', 'pI']
    for file in glob.glob('../data/*inf/*.xlsx'):
        print('Processing: ', file)
        df = pd.read_excel(file)
        print(df.shape)
        df.Peptide = df['Peptide'].apply(lambda x: re.sub(r'[()+-0123456789]',"",x))
        vals = df.Peptide.map(lambda x: GetPseudoAAC(x,5, AAP=AAP)).tolist()
        df = pd.concat([df,pd.DataFrame(vals)],axis=1)
        print(df.shape, df.isnull().sum())
        df.to_csv(file[:-5] + "pseAAC.csv")



# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:12:35 2021

@author: matti
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob
from sklearn.preprocessing import normalize, MinMaxScaler
from sklearn.decomposition import PCA

def norm(series):
    m = series.min()
    ma = series.max()
    return series.apply(lambda x: (x-m)/(ma-m))

if __name__ == '__main__':
    dfs = []
    i = 0
    for path in glob.glob('../data/*inf/*pseAAC.csv'):
        tdf = pd.read_csv(path)
        tdf['type'] = np.ones(tdf.shape[0])*('ninf' in path)
        tdf['nr'] = np.ones(tdf.shape[0])*i 
        i += 1
        dfs.append(tdf)
    df = pd.concat(dfs, axis=0)
    datacolumns = df.select_dtypes(include='number')
    datacolumns = datacolumns.loc[:,datacolumns.isnull().sum() == 0].columns
    
    mean_data = df[datacolumns].groupby(['type','nr']).mean()
    mean_data = mean_data.apply(norm, axis=0)
    mean_data = mean_data[[x for x in mean_data.columns if 'AAC' in x]]
    """
    for column in mean_data:
        sns.histplot(x=column, data=mean_data, hue='type')
        plt.title(column)
        plt.show()
    """
    pca = PCA()
    pca_data = pca.fit_transform(mean_data)
    print(pca.explained_variance_ratio_)
    mean_data['0'] = pca_data[0]
    mean_data['1'] = pca_data[1]
    mean_data['2'] = pca_data[2]
    sns.scatterplot('0', '1', data=mean_data, hue='type')
    plt.show()
    sns.scatterplot('0', '2', data=mean_data, hue='type')
    plt.show()
    sns.scatterplot('1', '2', data=mean_data, hue='type')
    plt.show()
    
    
    
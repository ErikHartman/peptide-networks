# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 09:15:09 2021

@author: matti
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as stats

inf = pd.read_excel('../data/peptide_sample_inf/peptide_sample_31.xlsx')
inf['lens'] = inf.Peptide.apply(len)
con = pd.read_excel('../data/peptide_sample_ninf/peptide_sample_13.xlsx')
con['lens'] = con.Peptide.apply(len)

data = [inf, con]
numerical = [x == np.int64 or x == np.float64 for x in inf.dtypes]
for i, col in enumerate(inf.columns):
    if numerical[i]:
        try:
            plt.hist(inf[col], bins=20, alpha=0.3, density=True)
            plt.hist(con[col], bins=20, alpha=0.3, density=True)
            plt.title(col)
            plt.show()
            print(col, stats.ttest_ind(inf[col], con[col]))
        except KeyError:
            plt.show()
            pass
    


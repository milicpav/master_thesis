# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 23:13:30 2018

@author: pmilicka
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def boxplots(col_name, dfs, dataset_names, title = ""):
    plt.figure(figsize = (8, 6))
    n_dfs = len(dfs)
    n_problems = len(dfs[0])
    vals = np.array([dfs[0][col_name].values])
    for df_idx in range(1, n_dfs):
        vals = np.append(vals, [dfs[df_idx][col_name].values], axis = 0)
    vals = vals.transpose()
    plt.boxplot(vals,
                sym = "b+",
                patch_artist = True,
                labels = dataset_names,
                meanline = True,
                showmeans= True)
    plt.grid(axis = "y")
    plt.xlabel("Datasets")
    plt.ylabel(col_name)
    plt.title(title)
    plt.xticks = dataset_names
    
df02 = pd.read_csv("datasets/RG10/rg10_02/x.csv", header = 0, sep = ";")
df04 = pd.read_csv("datasets/RG10/rg10_04/x.csv", header = 0, sep = ";")
df06 = pd.read_csv("datasets/RG10/rg10_06/x.csv", header = 0, sep = ";")
df08 = pd.read_csv("datasets/RG10/rg10_08/x.csv", header = 0, sep = ";")

for col_name in list(df04):
    if "int" in str(df04[col_name].dtype)\
        or "float" in str(df04[col_name].dtype):
        boxplots(col_name, 
                 [df08, df06, df04, df02], 
                 dataset_names = ["rg10_08", "rg10_06", "rg10_04", "rg10_02"],
                 title = col_name
                 )

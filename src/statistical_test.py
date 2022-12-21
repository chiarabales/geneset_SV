import pandas as pd
import numpy as np
import scipy.stats as ss
import os


import SOUG as sg




def create_fisher_exact_test(relevant_genes, gen, which_ordering):
    p_values = []
    
    all_genes_in_pathway = set({})
    for i in range(gen.shape[0]):
        all_genes_in_pathway = all_genes_in_pathway.union(set(gen.iloc[i][1:]))
    all_genes_in_pathway = all_genes_in_pathway.difference({np.nan})

    name_significant = 'significant_' + relevant_genes + '.genes'
    name_all = 'all_' + relevant_genes + '.genes'
    
    _data_input = os.path.join('..', 'data/__PHENOTYPES')
    path_significant = os.path.join(_data_input, name_significant)
    path_all = os.path.join(_data_input, name_all)
    
    genes_all = pd.read_csv(path_all, header = None)
    
    genes_significant = pd.read_csv(path_significant, header = None)
    genes_all = np.asarray(genes_all[0]).tolist()
    genes_significant = np.asarray(genes_significant[0]).tolist()

    
    genes_significant = set(genes_significant)
    genes_all = set(genes_all)
    U = all_genes_in_pathway.union(genes_all)
    for j in which_ordering[:]:
        genes_in_pathway = set(gen.iloc[which_ordering[j]][1:]).difference({np.nan})
        
        matrix = [[0,0], [0,0]]
        matrix[0][0] = len(genes_in_pathway.intersection(genes_significant))
        matrix[0][1] = len(genes_significant.difference(genes_in_pathway))
        matrix[1][0] = len(genes_in_pathway.difference(genes_significant))
        matrix[1][1] = len(U.difference(genes_in_pathway).difference(genes_significant))
        
        p_values.append(ss.fisher_exact(matrix)[1])
    return p_values
import pandas as pd
import numpy as np
import scipy.stats as ss
import os

import SOUG as sg
import evaluation_measures as em
import data_preprocessing as dp

# ______________________________________________________________________________________________________________
#
# PUNISHED ORDERINGS
# ______________________________________________________________________________________________________________

def order_punished(df, type_punish = 'normal', max_ranking = None):
    
    ordering = []
    
    # data pre-processing
    # ___________________________________________________________________________
    
    data_features = df.copy()
    mydata = df.copy()
    mydata = mydata.drop(columns = 'pathway')
    mydata = np.asarray(mydata)
    mydata = mydata.astype(int)
    
    # first ranked element
    # ___________________________________________________________________________

    SV_features = sg.calculate_svs(mydata)
    data_features['SV'] = SV_features
    which = data_features['SV'].argmax()
    ordering.append(which)

    # initialization of punishment 
    # ___________________________________________________________________________

    dim = np.shape(mydata)[0]
    punish = np.zeros((dim, dim))
    y = mydata[which]
    for j in range(dim):
        x = mydata[j]
        if j not in ordering:
            punish[:, which][j] = em.jaccard_distance(x,y)
    max_svs = []
    
    # recomputation of punishment
    # ___________________________________________________________________________

    mydata[which] = np.zeros(np.shape(mydata)[1])
    
    if max_ranking == None:
        max_ranking = np.shape(mydata)[0]
        
    for j in range(1, max_ranking):
        SV_features = sg.calculate_svs(mydata) 
        SV_sum = np.sum(SV_features)
        
        p_max = np.max(np.sum(punish, axis = 1))
        SV_max = np.max(SV_features)
        
        max_svs.append(SV_max)
        
        for o in ordering:
            SV_features[o] = -100  
        if type_punish == 'normal':
            rescale = 1
        elif type_punish == 'rescaled':
            rescale = SV_max/p_max
        data_features['SV'] = SV_features - (np.sum(punish, axis = 1)*rescale)
        which = data_features['SV'].argmax()
        ordering.append(which)
        y = mydata[which]
        for j in range(dim):
            x = mydata[j]
            if j not in ordering:
                punish[:, which][j] = em.jaccard_distance(x,y)
        mydata[which] = np.zeros(np.shape(mydata)[1])
        
    return ordering

def order_punished_artificial_pathway(df, type_punish = 'normal', max_ranking = None):
    
    ordering = []
    
    # data pre-processing
    # ___________________________________________________________________________
    
    data_features = df.copy()
    mydata = df.copy()
    mydata = mydata.drop(columns = 'pathway')
    mydata = np.asarray(mydata)
    mydata = mydata.astype(int)

    # first ranked element
    # ___________________________________________________________________________
    
    SV_features = sg.calculate_svs(mydata)
    data_features['SV'] = SV_features
    
    which = data_features['SV'].argmax()
    ordering.append(which)

    # initialization of punishment 
    # ___________________________________________________________________________
    
    dim = np.shape(mydata)[0]
    punish = np.zeros(dim)
    
    genes_selected = data_features.loc[ordering]
    genes_selected.loc[len(ordering), 'pathway'] = 'artificial'
    for g in genes_selected.columns:
        if not g in ['pathway', 'SV']:
            if genes_selected[g].sum() > 0:
                genes_selected.loc[len(genes_selected)-1, g] = 1
            else:
                genes_selected.loc[len(genes_selected)-1, g] = 0
    genes_for_punish = genes_selected[genes_selected.columns.difference(['pathway', 'SV'])]
    y = np.asarray(genes_for_punish.loc[len(genes_for_punish)-1])
    
    for j in range(dim):
        x = mydata[j]
        if j not in ordering:
            punish[j] = em.jaccard_distance(x,y)            
    max_svs = []
    
    # recomputation of punishment
    # ___________________________________________________________________________

    mydata[which] = np.zeros(np.shape(mydata)[1])
        
    if max_ranking == None:
        max_ranking = np.shape(mydata)[0]
    
    for j in range(1, max_ranking):
    
        SV_features = sg.calculate_svs(mydata) 
        SV_sum = np.sum(SV_features)
        
        SV_max = np.max(SV_features)
        p_max = np.max(punish)
        
        if type_punish == 'normal':
            rescale = 1
        elif type_punish == 'rescaled':
            rescale = SV_max/p_max
        
        for o in ordering:
            SV_features[o] = -100   
        data_features['SV'] = SV_features - punish*rescale
        which = data_features['SV'].argmax()
        ordering.append(which)
        
        genes_selected = data_features.loc[ordering]
        genes_selected.loc[len(ordering), 'pathway'] = 'artificial'
        for g in genes_selected.columns:
            if not g in ['pathway', 'SV']:
                if genes_selected[g].sum() > 0:
                    genes_selected.loc[len(genes_selected)-1, g] = 1
                else:
                    genes_selected.loc[len(genes_selected)-1, g] = 0
        genes_for_punish = genes_selected[genes_selected.columns.difference(['pathway', 'SV'])]
        y = np.asarray(genes_for_punish.loc[len(genes_for_punish)-1])
        
        for j in range(dim):
            x = mydata[j]
            if j not in ordering:
                punish[j] = em.jaccard_distance(x,y)
        mydata[which] = np.zeros(np.shape(mydata)[1])
       
    return ordering


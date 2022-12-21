
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import os
import pickle

import pandas as pd
import random
from statsmodels.stats.multitest import multipletests

import genomics_old as g
import statistical_test as st


import click


@click.command()
@click.option('--_alpha', type=float)


def main(_alpha):
    phenotypes = [
        'HEIGHT',
        'BLOOD_PLATELET_COUNT',
        'STANDING_HEIGHT',
        'BLOOD_RED_COUNT',
        'HEEL_TSCORE',
        'BLOOD_WHITE_COUNT',
        'BLOOD_EOSINOPHIL_COUNT',
        'SITTING_HEIGHT',
        'TRUNK_FAT_FREE_MASS',
        'TRUNK_PREDICTED_MASS',
        'WHOLE_BODY_FAT_FREE_MASS',
        'WHOLE_BODY_WATER_MASS',
        'SYSTOLIC_BLOOD_PRESSURE',
        'BASAL_METABOLIC_RATE',
        'IMPENDANCE_OF_WHOLE_BODY',
        'BMI',
        'COMPARATIVE_HEIGHT_SIZE_AT_AGE_10',
        'ARM_PREDICTED_MASS_R',
        'ARM_FAT_FREE_MASS_R',
        'LEG_FAT_FREE_MASS_R',
        'LEG_PREDICTED_MASS_R',
        'LUNG_FEV1FVC_RATIO',
        'IMPENDANCE_OF_LEG_R',
        'IMPENDANCE_OF_ARM_R',
        'LUNG_FVC',
        'WEIGHT',
        'WHRATIO',
        'HAIR_PIGMENT',
        'HIP_CIRCUMFERENCE',
        'TRUNK_FAT_MASS',
        'WHOLE_BODY_FAT_MASS',
        'ARM_FAT_MASS_R',
        'LEG_FAT_MASS_R',
        'TRUNK_FAT_PERCENTAGE',
        'BODY_FAT_PERCENTAGE',
        'ARM_FAT_PERCENTAGE_R',
        'LEG_FAT_PERCENTAGE_R',
        'CARDIOVASCULAR_DISEASE'
    ]
    
    
    alpha = 0.1
    

    print('sto facendo qualcos')
    
    fig, ax = plt.subplots(len(phenotypes), 4)
    fig.set_size_inches((16,3.5*len(phenotypes)))

    name_ =  '_' + str(round(((1-alpha)*100)))
    graph_name = 'BONFERRONI_pvalues' + name_ + '.pdf'

    _data_output = os.path.join('..', 'data/_FISHEREXACTTEST' )
    output_where_graph = os.path.join(_data_output, graph_name)

    _genesets = ['KEGG', 'CGN', 'CM', 'TFT_LEGACY']

    f_BONFERRONI = []

    for geneset in range(len(_genesets)):

        _data_output = os.path.join('..', 'data/' + _genesets[geneset])
        _data_input = os.path.join('..', 'data/' + _genesets[geneset])
        _data_gsea = os.path.join('..', 'data/GSEA/' + _genesets[geneset])

        output_where = os.path.join(_data_output, _genesets[geneset] + '_SV_ordering.pickle')
        SV_ordering = pickle.load(open(output_where, "rb" ) )

        output_where = os.path.join(_data_output, _genesets[geneset] + '_ordering_normal.pickle')
        ordering_normal = pickle.load(open(output_where, "rb" ) )
        output_where = os.path.join(_data_output, _genesets[geneset] + '_ordering_normal_rescaled.pickle')
        ordering_normal_rescaled = pickle.load(open(output_where, "rb" ) )
        output_where = os.path.join(_data_output, _genesets[geneset] + '_ordering_artificial.pickle')
        ordering_artificial = pickle.load(open(output_where, "rb" ) )
        output_where = os.path.join(_data_output, _genesets[geneset] + '_ordering_artificial_rescaled.pickle')
        ordering_artificial_rescaled = pickle.load(open(output_where, "rb" ) )

        path = os.path.join(_data_input, _genesets[geneset] + '.csv')
        if _genesets[geneset] == 'KEGG':
            gen = pd.read_csv(path, sep = ';', header = None, low_memory = False)
        else:
            gen = pd.read_csv(path, sep = ',', header = None, low_memory = False)
        gen = gen.drop(columns = 1)

        ordering_random = list(range(gen.shape[0]))
        random.shuffle(ordering_random)
        for phenotype in range(len(phenotypes)):
        # DO THINGS FOR FISHER EXACT TEST
            p_random = np.asarray(st.create_fisher_exact_test(phenotypes[phenotype], gen, ordering_random))
            p_SV_ordering = np.asarray(st.create_fisher_exact_test(phenotypes[phenotype], gen, SV_ordering))
            p_ordering_normal = np.asarray(st.create_fisher_exact_test(phenotypes[phenotype], gen, ordering_normal))
            p_ordering_normal_rescaled = np.asarray(st.create_fisher_exact_test(phenotypes[phenotype], gen, ordering_normal_rescaled))
            p_ordering_artificial = np.asarray(st.create_fisher_exact_test(phenotypes[phenotype], gen, ordering_artificial))
            p_ordering_artificial_rescaled = np.asarray(st.create_fisher_exact_test(phenotypes[phenotype], gen, ordering_artificial_rescaled))
            total_genes = int(len(p_random))
            print('-------------------------------------------\n', phenotypes[phenotype])


            max_accepted = int(len(p_random)/4)
            count_random = []
            count_SV_ordering = []
            count_ordering_normal = []
            count_ordering_normal_rescaled = []
            count_ordering_artificial = []
            count_ordering_artificial_rescaled = []
            for max_accepted in range(int(len(p_random))):
                bf_p = alpha/(max_accepted+1)
                count_random.append(len(p_random[np.where(p_random[:max_accepted] < bf_p)]))
                count_SV_ordering.append(len(p_SV_ordering[np.where(p_SV_ordering[:max_accepted] < bf_p)]))
                count_ordering_normal.append(len(p_ordering_normal[np.where(p_ordering_normal[:max_accepted] < bf_p)]))
                count_ordering_normal_rescaled.append(len(p_ordering_normal_rescaled[np.where(p_ordering_normal_rescaled[:max_accepted] < bf_p)]))
                count_ordering_artificial.append((len(p_ordering_artificial[np.where(p_ordering_artificial[:max_accepted] < bf_p)])))
                count_ordering_artificial_rescaled.append(len(p_ordering_artificial_rescaled[np.where(p_ordering_artificial_rescaled[:max_accepted] < bf_p)]))
            selected_pathways = len(p_random[np.where(p_random < alpha/len(p_random))])

            # CREATE NOISE
            fdrA = np.random.normal(0, alpha, len(count_random))
            fdrB = np.random.normal(0, alpha, len(count_random))
            fdrC = np.random.normal(0, alpha, len(count_random))
            fdrD = np.random.normal(0, alpha, len(count_random))
            fdrE = np.random.normal(0, alpha, len(count_random))
            fdrF = np.random.normal(0, alpha, len(count_random))

            # PLOT FOR COUNT REGARDING FISHER_EXACT TEST
            ax[phenotype, geneset].plot(np.arange(len(count_random)), count_SV_ordering + fdrA[:], '.', label = 'SV') 
            ax[phenotype, geneset].plot(np.arange(len(count_random)), count_ordering_normal + fdrB[:], '.', label = 'PO') 
            ax[phenotype, geneset].plot(np.arange(len(count_random)), count_ordering_normal_rescaled + fdrC[:], '.', label = 'POR') 
            ax[phenotype, geneset].plot(np.arange(len(count_random)), count_ordering_artificial + fdrD[:], '.', label = 'AO') 
            ax[phenotype, geneset].plot(np.arange(len(count_random)), count_ordering_artificial_rescaled + fdrE[:], '.', label = 'AOR') 
            ax[phenotype, geneset].hlines(selected_pathways, 0, max_accepted, colors=None, linestyles='--', label='all', linewidth=7.0)
            if phenotype == 0:
                ax[phenotype, geneset].set_title('%s' % _genesets[geneset], fontsize = 16)
            if geneset == 0:
                ax[phenotype, geneset].set_ylabel('%s' % phenotypes[phenotype], fontsize = 14)
            if phenotype != len(phenotypes)-1:
                ax[phenotype, geneset].set_xticks([])

            ax[phenotype, geneset].tick_params(axis="x", labelsize=14)
            ax[phenotype, geneset].tick_params(axis="y", labelsize=14)
            y_lim_max = max(21, selected_pathways+5)

            ax[phenotype, geneset].set_ylim([-1,y_lim_max])

            if phenotype == len(phenotypes)-1:
                ax[phenotype, geneset].set_xlabel('# pathways selected', fontsize=14)
                if geneset == len(_genesets)-1:
                    ax[phenotype, geneset].legend(fontsize = 14)
    fig.savefig(output_where_graph, bbox_inches = 'tight')
    
if __name__ == '__main__':
    main()
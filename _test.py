import numpy as np 
import pandas as pd

import SOUG as sg
import genomics as g
import data_preprocessing as dp

import random
import click
import os
import pickle

@click.command()
@click.option('--geneset', type=str)
@click.option('--ordering', type=str)

def main(geneset, ordering):
 
    _data_output = os.path.join('..', 'data/' + geneset)
    _data_input = os.path.join('..', 'data/' + geneset)
    _plot_path = os.path.join('..', 'plot_genesets')
    path = os.path.join(_data_input, geneset + '.csv')
    if geneset == 'KEGG':
        gen = pd.read_csv(path, sep = ';', header = None)
    else:
        gen = pd.read_csv(path, sep = ',', header = None)
    gen = gen.drop(columns = 1)
    data = dp.transform_gen(gen, geneset)
    max_ranking = np.shape(data)[0]

    if ordering == 'SV':
        mydata = data.copy()
        mydata = mydata.drop(columns = 'set')
        mydata = np.asarray(mydata)
        SV_pathways = sg.calculate_svs(mydata)
        data_new = data.copy()
        data_new['SV_original'] = SV_pathways
        data_new = data_new.sort_values(by = 'SV_original', ascending=False)
        order = []
        for row in data_new[['set', 'SV_original']].index:
            order.append(row)

    elif ordering == 'PO':
        order = g.PO(data, 'N', max_ranking = max_ranking)

    elif ordering == 'POR':
        order = g.PO(data, 'R', max_ranking = max_ranking)

    elif ordering == 'AO':
        order = g.AO(data, 'N', max_ranking = max_ranking)

    elif ordering == 'AOR':
        order = g.AO(data,  'R', max_ranking = max_ranking)
            
    return order

if __name__ == '__main__':
    main()

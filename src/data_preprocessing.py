import numpy as np
import os
import pandas as pd

def create_columns(gen):
    columns = []
    for j in np.arange(gen.shape[0]):
        if j % 50 == 0:
            print(j, len(columns))
        for i in np.arange(2, gen.shape[1]+1):
            if str(gen.loc[j][i]) != 'nan':
                if (gen.loc[j][i] not in columns):
                    columns.append(gen.loc[j][i])
    return columns


def load_data(geneset):
    
    # load data and create the columns file
    
    _data_input = os.path.join('..', 'data')

    path = os.path.join(_data_input, geneset + '.csv')
    gen = pd.read_csv(path, sep = ';', header = None)
    gen = gen.drop(columns = 1)
    
    columns = []
    for j in np.arange(gen.shape[0]):
        if j % 50 == 0:
            print(j, len(columns))
        for i in np.arange(2, gen.shape[1]+1):
            if str(gen.loc[j][i]) != 'nan':
                if (gen.loc[j][i] not in columns):
                    columns.append(gen.loc[j][i])
    
    
    columns_output = 'cols_' + geneset + '.npy'
    
    path_columns = os.path.join(_data_input, columns_output)
    np.save(path_columns, columns)
    
    return gen


def transform_gen(gen_set):
    
#     # loading data 
#     _data_input = os.path.join('..', 'data')
#     file_input = gen_set + '.csv'
#     path = os.path.join(_data_input, file_input)
#     gen = pd.read_csv(path, sep = ';', header = None)
#     gen = gen.drop(columns = 1)
    
    gen = load_data(gen_set)
    
    # loading columns
    _data_cols = os.path.join('..', 'data')
    name_columns = 'cols_' + gen_set + '.npy'
    cols_where = os.path.join(_data_cols, name_columns)
    cols = np.load(cols_where)
    
    
#     if gen_set == 'KEGG':
#         cols_where = os.path.join(_data_cols, 'cols_KEGG.npy')
#         cols = np.load(cols_where)
        
#     elif gen_set == 'CGN':
#         cols_where = os.path.join(_data_cols, 'cols_CGN.npy')
#         cols = np.load(cols_where)

#     elif gen_set == 'CM':
#         cols_where = os.path.join(_data_cols, 'cols_CM.npy')
#         cols = np.load(cols_where)

#     elif gen_set == 'TFT_LEGACY':
#         cols_where = os.path.join(_data_cols, 'cols_TFT_LEGACY.npy')
#         cols = np.load(cols_where)
          
    columns = cols.tolist()
    
    
    ar = np.zeros(shape = (gen.shape[0],len(columns)))
    data = pd.DataFrame(ar, columns=columns)
    
    for j in range(gen.shape[0]):

        s = str(gen[0][j])
        data['set'] = 0
        data['set'] = gen[0]
        for i in range(2, gen.shape[1]):
            s = str(gen[i][j])
            if s != 'nan':
                data[s][j] = 1    
    return data



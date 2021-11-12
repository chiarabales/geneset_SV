"""
Created on Feb 11 20:43:02 2021

@author: Chiara
"""

import numpy as np
import random


random.seed(42)

def set_from_matrix(A):
    sets = []
    for i in range(np.shape(A)[1]):
        sets.append([])
        for j in range(np.shape(A)[0]):
            if A[j, i] != 0:
                sets[i].append(j)
    sets[:] = (value for value in sets if value != [])
    return sets

def create_soug(A):  
    game_dict = {}
    sets = set_from_matrix(A)
    for i in sets:
        key = str(i)
        if key in game_dict.keys():
            game_dict[key] = 1 + game_dict[key]
        else:
            game_dict[key] = 1
    return game_dict

def calculate_svs(A):
    players = np.shape(A)[0]
    svs = np.zeros(players)
    game_dict = create_soug(A)
    T_sets = game_dict.keys()
    l = 0
    for T in T_sets:
        l += game_dict[T]
    if l == 0:
        l = 1
    for feature in range(players):
        for T in T_sets:
            T_set = eval(T)
            if feature in T_set:
                svs[feature] += game_dict[T] / len(T_set)
    return svs/l

def jaccard_distance(x, y):
    j_min = 0
    j_max = 0
    for i in range(len(x)):
        j_min += np.minimum(x[i], y[i])
        j_max += np.maximum(x[i], y[i])
    return j_min/j_max

import numpy as np

# ______________________________________________________________________________________________________________

# DIFFERENT FUNCTIONS TO DEFINE THE DISTANCE BETWEEN PATHWAYS 
# ______________________________________________________________________________________________________________

def jaccard_distance(x, y):
    
    # ___________________________________________________________________________
    #
    # pairwise jaccard distance among two sets
    # ___________________________________________________________________________
    
    j_min = 0
    j_max = 0
    for i in range(len(x)):
        j_min += np.minimum(x[i], y[i])
        j_max += np.maximum(x[i], y[i])
    return j_min/j_max


def jaccard_rate(data, order, which_order):
    
    # ___________________________________________________________________________
    #
    # averaged pairwise jaccard distance within family of sets
    # ___________________________________________________________________________
    
    df = data.loc[which_order[:order]]
    num_genes = np.shape(data)[1]
    jr = 0
    for i in which_order[:order]:
        for j in which_order[:order]:
            x = np.asarray(df.loc[i][1:num_genes])
            y = np.asarray(df.loc[j][1:num_genes])
            if i != j:
                jr = jr + jaccard_distance(x,y)
            else:
                jr = jr
    return jr/(order*(order-1))



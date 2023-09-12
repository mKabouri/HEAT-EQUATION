import numpy as np
from math import *
from random import randint



# s: symetrique, d: définie, p: positive, c: creuse.
# On utilise une propriété sur les matrices symétriques positives:
#"Une matrice M est symétrique définie positive ssi les valeurs propres de M sont positives"
def sdpc_matrix_generator(dimension, nb_termes):
    A = np.zeros((dimension, dimension))
    for i in range(nb_termes//2):
        m = randint(0, dimension-1)
        n = randint(0, dimension-1)
        value = randint(1, 10)
        if m == n:
            if 0 < m:
                m = m - 1
            elif m == 0:
                m = m + 1
        A[m][n] = value
        A[n][m] = value
    for i in range(dimension):
        A[i][i] = randint(1, 10)
    L = np.linalg.eigvals(A)	
    for i in range(len(L)):
        if(L[i] < 0):
            return sdpc_matrix_generator(dimension, nb_termes)
    return A


if __name__=='__main__':
    print(sdpc_matrix_generator(4, 4))
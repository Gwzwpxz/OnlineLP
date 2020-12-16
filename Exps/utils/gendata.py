import numpy as np
from scipy import sparse

# Data generator
def gendata(m, n, typeid, tightness=0.25, bord=1.0, sps=1.0):
    
    if typeid == 0:
        if sps < 1.0:
            A = np.round(sparse.rand(m, n, sps, format="csc"), 2)
            b = np.sum(A, axis=1).reshape(m, 1) * (0.25 / n**(1 - bord))
            c = np.sum(A, axis=0).reshape(n, 1) / m + np.random.rand(n, 1) * 5
        else:
            A = np.random.randint(1, 1000, (m, n)) / 100
            b = np.sum(A, axis=1).reshape(m, 1) * (0.25 / n**(1 - bord))
            c = np.sum(A, axis=0).reshape(n, 1) / m + np.random.rand(n, 1) * 5
    else:   
        A = np.random.randn(m, n) + 1
        b = (np.random.rand(m, 1) + 1) / 3 * (n**bord) 
        c = np.sum(A, axis=0).reshape(n, 1) - np.random.rand(n, 1) * m
    
    return {"A": A, "b": b, "c": c}    
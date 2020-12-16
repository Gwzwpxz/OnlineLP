# Import dependencies
import numpy as np
import scipy.sparse
from scipy.io import savemat, loadmat
from gurobipy import *


# Online Algorithm
def fastLP(A, b, c, K, Method):
    m = A.shape[0]
    n = A.shape[1]
    
    # It is worth considerinvg whether it is better to exclude K here
    # stepsize = 1 / np.sqrt(n * K)
    
    # Initialize dual solution
    if Method == "M":
        y = np.ones((m, 1)) / np.exp(1)
    else:
        y = np.zeros((m, 1))
        
    # Initialize resource
    d = b / n
    
    # Initialize primal solution
    x = np.zeros((n, 1))
    
    # Start dual descent
    for i in range(K):
        
        p = np.random.permutation(n)
        
        for j in p:
            
            stepsize = 1 / np.sqrt(n * (i + 1))
            
            if type(A) == scipy.sparse.csc.csc_matrix:
                aa = A[:, j].todense().reshape(m, 1)
            else:
                aa = A[:, j].reshape(m, 1)
            xk = (c[j] > np.dot(aa.T, y))
            
            if Method  == "M":
                y = np.multiply(y, np.exp(- stepsize * (d - aa * xk)))
            else:
                y = y - stepsize * (d - aa * xk)
                y = np.maximum(y, 0.0)
            
            x[j] += xk[0][0]
            
    obj = np.dot(c.T, x / K)
            
    return {"x": x / K, "y": y, "obj": obj}


# Rounding
def rounding(A, b, c, x):
    
    m = b.size
    n = c.size
    p = np.random.permutation(n)
    rdx = np.zeros((n, 1))
    
    for i in p:
        aa = A[:, i].reshape(m, 1)
        isround = (np.random.rand() <= x[i])
        if isround and (np.min(b - aa) >= 0):
            rdx[i] = 1
            b = b - aa
    
    obj = np.dot(c.T, rdx)
    
    return {"rdx": rdx, "obj": obj}

# Gurobi LP
def GRBLP(A, b, c):
    
    model = Model()
    x = model.addMVar(c.size, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS)
    constr = model.addMConstr(A, x, GRB.LESS_EQUAL, b.squeeze())
    model.setMObjective(Q=None, c=c.squeeze(), constant=0.0, sense=GRB.MAXIMIZE)
    model.update()
    model.optimize()
    optdual = model.getAttr(GRB.Attr.Pi, model.getConstrs())
    optx = model.getAttr(GRB.Attr.X, model.getVars())
    time = model.getAttr(GRB.Attr.Runtime)
    obj = model.getAttr(GRB.Attr.ObjVal)
    
    return {"x": optx, "y": optdual, "time": time, "model": model, "obj": obj}

# Gurobi MIP with initialization
def GRBMIP(A, b, c, initX=None):
    
    model = Model()
    x = model.addMVar(c.size, vtype=GRB.BINARY)
    
    # Set initial solution
    if initX is not None:
        for i in range(c.size):
            x[i].setAttr(GRB.Attr.Start, initX[i])
    
    constr = model.addMConstr(A, x, GRB.LESS_EQUAL, b.squeeze())
    model.setMObjective(Q=None, c=c.squeeze(), constant=0.0, sense=GRB.MAXIMIZE)
    model.setParam(GRB.Param.MIPGap, 0.01)
    model.update()
    model.optimize()
    optx = model.getAttr(GRB.Attr.X, model.getVars())
    time = model.getAttr(GRB.Attr.Runtime)
    obj = model.getAttr(GRB.Attr.ObjVal)
    
    return {"x": optx, "time": time, "model": model, "obj": obj}


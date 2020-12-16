# This is a script testing relation between regret and K
import os
import sys
sys.path.append("..\\utils\\")
from matplotlib.pyplot import plot
from matplotlib.pyplot import savefig
from algs import *
from gendata import *

approx = np.zeros(20)
regret = np.zeros(20)

objextr = lambda x: x["obj"]
timeextr = lambda x: x["time"]
dualextr = lambda x: x["y"]

m = 5
n = 100
typeid = 0
tightness = 0.5
bord = 1

data = gendata(m, n, typeid, tightness, bord)
A = data["A"]
b = data["b"]
c = data["c"]

lpres = GRBLP(A, b, c)

counter = 0

for i in np.ceil(np.logspace(1, 3, 20)).astype(int):
    approxres = fastLP(A, b, c, i + 1, "M")
    approx[counter] = objextr(approxres) / objextr(lpres)
    regret[counter] = - objextr(approxres) + objextr(lpres)
    print("Approximation ratio: {0}".format(approx[counter]))
    counter += 1
    

plot(regret * np.sqrt(np.logspace(1, 3, 20)))
savefig("regret_K_type_{0}_m_{1}_n_{2}_tight_{3}_bord_{4}.png".format(typeid, m, n, tightness * 100, bord))
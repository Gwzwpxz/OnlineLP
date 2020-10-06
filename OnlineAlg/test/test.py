import os
import numpy as np
from scipy import sparse
from scipy.io import savemat
import gurobipy as gp
from gurobipy import GRB

def GenData(m, n, bord, tightness):
    
    A = np.random.randint(1, 1000, (m, n)) / 1000

    # A = np.multiply(A, (sparse.rand(m, n, 0.0375) > 0).todense())
    c = np.zeros(n)

    sumA2 = np.sum(A, axis=0) / m
    sumA2 = np.asarray(sumA2).squeeze()
    for i in range(n):
        c[i] = sumA2[i] + np.random.rand() * 0.5

    c = c / np.max(c) * (np.random.rand() + 4) / 5
    
    A = sparse.csc_matrix(A)
    b = np.zeros(m)
    if (bord == "1_2"):
        sumA = tightness * np.sum(A, axis=1) / (n ** 0.5)
        sumA.flatten()
        for i in range(m):
            b[i] = sumA[i]
    elif (bord == "1_3"):
        sumA = tightness * np.sum(A, axis=1) / (n ** (2/3))
        sumA.flatten()
        for i in range(m):
            b[i] = sumA[i]

    return {"A": A, "b": b, "c": c}


common_header_gurobi = """
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gurobi_c.h"

int main(int argc, const char * argv[]) {

    int error = 0;
    double objval = 0.0;
    double soltime = 0.0;
    
    GRBenv *env = NULL;
    GRBmodel *model = NULL;
    
"""

common_err_gurobi = """
    if (error != 0) {
        goto QUIT;
    }
    
"""

common_body_gurobi = """
    
    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_LOGTOCONSOLE, 0);
    
    if (error != 0) {
        goto QUIT;
    }

    error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MIPGAP, 0.01);

    if (error != 0) {
        goto QUIT;
    }

	error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_TIMELIMIT, 2000.0);

    if (error != 0) {
        goto QUIT;
    }

    
    error = GRBoptimize(model);

    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);

    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBgetdblattr(model, GRB_DBL_ATTR_RUNTIME, &soltime);

    if (error != 0) {
        goto QUIT;
    }

    printf("\\nObj:%f \\n", objval);
    printf("%f ms\\n", soltime * 1000);
    
QUIT:
    if (error) {
        printf("ERROR: %s\\n", GRBgeterrormsg(env));
        exit(1);
    }

    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}
"""

common_body_gurobi_lp = """
    
    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_LOGTOCONSOLE, 0);
    
    if (error != 0) {
        goto QUIT;
    }

    error = GRBsetdblparam(GRBgetenv(model), GRB_DBL_PAR_MIPGAP, 0.01);

    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBoptimize(model);

    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);

    if (error != 0) {
        goto QUIT;
    }
    
    error = GRBgetdblattr(model, GRB_DBL_ATTR_RUNTIME, &soltime);

    if (error != 0) {
        goto QUIT;
    }

    printf("\\nObj:%f \\n", objval);
    printf("%f ms\\n", soltime * 1000);

    
QUIT:
    if (error) {
        printf("ERROR: %s\\n", GRBgeterrormsg(env));
        exit(1);
    }

    GRBfreemodel(model);
    GRBfreeenv(env);
    return 0;
}
"""

common_header_ollp = """
#include "ollpinterface.h"
#include <stdio.h>
#include <stdlib.h>
"""

common_body_ollp = """
int main(int argc, const char * argv[]) {
    
    int retcode = OLLP_RETCODE_OK;
    double *sol = (double *) malloc(nCol * sizeof(double));  
    double *binsol = (double *) malloc(nCol * sizeof(double));  
    
    ollp_prob *prob = NULL;
    double ObjVal = 0.0;
    double ObjValbin = 0.0;
    double soltimeA = 0.0;
    double soltimeB = 0.0;
    double Xmax = 1.0;
    
    retcode = OLLP_CreateProb(&prob);
    
    if (retcode != OLLP_RETCODE_OK) {
      goto exit_cleanup;
    }
    
    retcode = OLLP_LoadLP(prob, nCol, nInEqRow, boostingParam, 0, Xmax, colObj, 
                            inEqMatBeg, inEqMatIdx, inEqMatElem, inEqRHS);
                            
    if (retcode != OLLP_RETCODE_OK) {
      goto exit_cleanup;
    }
    
    retcode = OLLP_Solve(prob);
    
    if (retcode != OLLP_RETCODE_OK) {
      goto exit_cleanup;
    }
    
    OLLP_GetTime(prob, &soltimeA);
    OLLP_GetObjValbin(prob, &ObjValbin);

    OLLP_DeleteProb(&prob);

    retcode = OLLP_CreateProb(&prob);
    if (retcode != OLLP_RETCODE_OK) {
      goto exit_cleanup;
    }
    
    retcode = OLLP_LoadLP(prob, nCol, nInEqRow, boostingParam, 1, Xmax, colObj, 
                            inEqMatBeg, inEqMatIdx, inEqMatElem, inEqRHS);
                            
    if (retcode != OLLP_RETCODE_OK) {
      goto exit_cleanup;
    }
    
    retcode = OLLP_Solve(prob);
    
    if (retcode != OLLP_RETCODE_OK) {
      goto exit_cleanup;
    }
    
    OLLP_GetTime(prob, &soltimeB);
    OLLP_GetObjValbin(prob, &ObjVal);

    OLLP_GetObjVal(prob, &ObjVal);
    
    
    /* 
    printf("\\nRelaxed Solution: \\n");
    for (int i = 0; i < nCol; ++i) {
        printf("%f, ", sol[i]);
    } 
    */
    //printf("\\nRelaxed LP Objective value:  %f \\n", ObjVal);
    
    /*
    printf("\\nMLP Solution: \\n");
    for (int i = 0; i < nCol; ++i) {
        printf("%f, ", binsol[i]);
    }
    */
    
    printf("\\nMIPObj:%f \\n", ObjValbin);
    printf("Obj:%f \\n", ObjVal);
    printf("%f ms \\n", (soltimeA + soltimeB) * 500);
    
    return 0;
    
exit_cleanup:
    if (retcode != OLLP_RETCODE_OK) {
      printf("Errcode %d: Failed to solve problem.\\n", retcode);
    }

    OLLP_DeleteProb(&prob);
}

"""

commom_header_test = \
"""#! /bin/bash
for i in {1..100}
do
    gcc -O2 -DNDEBUG -I. -o gurobi_test_lp_$i.out gurobi_test_lp_$i.c ./libgurobi90.so

done

for i in {1..100}
do
    gcc -O2 -DNDEBUG -I. -o ollpdata_$i.out ollpdata_$i.c -Wl,-rpath . -L. -lollp -lm

done

for i in {1..100}
do
    gcc -O2 -DNDEBUG -I. -o gurobi_test_mlp_$i.out gurobi_test_mlp_$i.c ./libgurobi90.so

done


for i in {1..100}
do
"""

numTestCase = 100

"""

REMEMBER TO ADJUST XMAX !!!

"""
all_test_size_m_n = [(5, 100), (8, 1000), (16, 2000), (32, 4000), 
						(64, 10000), (128, 100000), (1024, 1000000)]

Xmax = 1

"""

REMEMBER TO ADJUST XMAX !!!

"""

mrange = [5, 8, 16, 32]
nrange = [100, 1000, 2000, 4000]
boosting_range = [50, 1000]
# sps_range = [False]
bord_range = ["1_2", "1_3"]
fbdlist = [(4096, 1000000, 500, "1_2")]
tightness_range = [0.25, 0.50, 0.75]

for all_m_n_combinations in all_test_size_m_n:
    (test_size_m, test_size_n) = all_m_n_combinations
    for boosting_param in boosting_range:
        for tightness in tightness_range:
            for bord in bord_range:

                if (test_size_m, test_size_n, boosting_param, bord) in fbdlist:
                    break

                for p in range(numTestCase):
                    data = GenData(test_size_m, test_size_n, bord, tightness)
                    A = data["A"]
                    b = data["b"]
                    c = data["c"]

                    A_s = sparse.csc_matrix(A)                    
                    
                    with open('./ollpdata_{0}.c'.format(p + 1), 'w') as f:
                        c_str = ", ".join([str(c[i]) for i in np.arange(c.size)])
                        b_str  = ", ".join([str(b[i]) for i in np.arange(b.size)])
                        Ap_str = ", ".join([str(A_s.indptr[i]) for i in np.arange(A_s.indptr.size)])
                        Ai_str = ", ".join([str(A_s.indices[i]) for i in np.arange(A_s.indices.size)])
                        Ax_str = ", ".join([str(A_s.data[i]) for i in np.arange(A_s.data.size)])
                        print('Generating OLLP testing script ollpdata_{}.c...'.format(p + 1))
                       
                        f.write(common_header_ollp)
                        f.write("int    nCol = {0:d};\n".format(test_size_n))
                        f.write("int    nInEqRow = {0:d};\n".format(test_size_m))
                        f.write("double colObj[{0:d}] = {{{1:s}}};\n".format(c.size, c_str))
                        f.write("int    inEqMatBeg[{0:d}] = {{{1:s}}};\n".format(A_s.indptr.size, Ap_str))
                        f.write("int    inEqMatIdx[{0:d}] = {{{1:s}}};\n".format(A_s.indices.size, Ai_str))
                        f.write("double inEqMatElem[{0:d}] = {{{1:s}}};\n".format(A_s.data.size, Ax_str))
                        f.write("double inEqRHS[{0:d}] = {{{1:s}}};\n".format(b.size, b_str))
                        f.write("int    boostingParam = {0:d}; \n".format(boosting_param))
                        f.write(common_body_ollp)
                        
                    print('Generating LP model file ollp_model_lp_{}.lp...'.format(p + 1))    
                    model = gp.Model()
                    x = model.addVars(test_size_n, lb=0, ub=Xmax, obj=0.0, name="x", vtype=GRB.CONTINUOUS)
                    model.addMConstrs(A_s, None, GRB.LESS_EQUAL, b, "INEQCONSTR")
                    model.setObjective(gp.quicksum(c[i] * x[i] for i in range(test_size_n)) , GRB.MAXIMIZE)
                    model.update()
                    model.write("ollp_model_lp_{0:d}.lp".format(p + 1))
                    
                    print('Generating MLP model file ollp_model_mlp_{}.lp...'.format(p + 1))    
                    model = gp.Model()
                    x = model.addVars(test_size_n, lb=0, ub=Xmax, obj=0.0, name="x", vtype=GRB.INTEGER)
                    model.addMConstrs(A_s, None, GRB.LESS_EQUAL, b, "INEQCONSTR")
                    model.setObjective(gp.quicksum(c[i] * x[i] for i in range(test_size_n)) , GRB.MAXIMIZE)
                    model.update()
                    model.write("ollp_model_mlp_{0:d}.lp".format(p + 1))
                    
                    with open('./gurobi_test_lp_{0}.c'.format(p + 1), 'w') as g:
                        print('Generating Gurobi testing file gurobi_test_lp_{}.c...'.format(p + 1))
                        g.write(common_header_gurobi)
                        g.write('    error = GRBloadenv(&env, "gurobi_test_lp_{}.log");\n'.format(p + 1))
                        g.write(common_err_gurobi)
                        g.write('    error = GRBreadmodel(env, "ollp_model_lp_{}.lp", &model);\n'.format(p + 1))
                        g.write(common_body_gurobi_lp)
                    
                    with open('./gurobi_test_mlp_{0}.c'.format(p + 1), 'w') as h:
                        print('Generating Gurobi testing file gurobi_test_mlp_{}.c...'.format(p + 1))
                        h.write(common_header_gurobi)
                        h.write('    error = GRBloadenv(&env, "gurobi_test_mlp_{}.log");\n'.format(p + 1))
                        h.write(common_err_gurobi)
                        h.write('    error = GRBreadmodel(env, "ollp_model_mlp_{}.lp", &model);\n'.format(p + 1))
                        h.write(common_body_gurobi)

                with open('./Log_LP_m_{}_n_{}_tightness_{}_boo_{}_bord_{}.txt'.format(test_size_m, test_size_n, tightness, boosting_param, bord), 'w') as x:
                        x.write("Test Result: \n")              

                with open('./Log_MLP_m_{}_n_{}_tightness_{}_boo_{}_bord_{}.txt'.format(test_size_m, test_size_n, tightness, boosting_param, bord), 'w') as x:
                        x.write("Test Result: \n")              

                with open('./Log_OLLP_m_{}_n_{}_tightness_{}_boo_{}_bord_{}.txt'.format(test_size_m, test_size_n, tightness, boosting_param, bord), 'w') as x:
                        x.write("Test Result: \n")              

                with open('run.sh', 'w') as y:
                    y.write(commom_header_test)
                    y.write("\n   ./gurobi_test_lp_$i.out >> Log_LP_m_{}_n_{}_tightness_{}_boo_{}_bord_{}.txt\n".format(test_size_m, test_size_n, tightness, boosting_param, bord))
                    y.write("done ")
                    y.write("\nfor i in {1..100}")
                    y.write("\ndo\n ")
                    y.write("\n   ./ollpdata_$i.out >> Log_OLLP_m_{}_n_{}_tightness_{}_boo_{}_bord_{}.txt\n".format(test_size_m, test_size_n, tightness, boosting_param, bord))
                    y.write("done\n")
                    y.write("for i in {1..100}")
                    y.write("\ndo\n ")
                    y.write("\n   ./gurobi_test_mlp_$i.out >> Log_MLP_m_{}_n_{}_tightness_{}_boo_{}_bord_{}.txt\n".format(test_size_m, test_size_n, tightness, boosting_param, bord))
                    y.write("done \n")
                    y.write("./cleanup.sh")

                os.system("chmod +x *.sh")
                os.system("./run.sh")


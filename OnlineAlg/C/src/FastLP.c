//
//  FastLP.c
//  OLLP
//
//  Created by 高文智 on 2020/7/18.
//  Copyright © 2020 高文智. All rights reserved.
//

#include <string.h>
#include <time.h>
#include "FastLP.h"
#include "Sort.h"
#include "cs.h"
/*  Online Linear Programming Computation Routine */

// b: n^{1/2} ~ n^{1/3}
// -1 <= A_ij <= 1 (positive)
// c: arbitrary    (positive)
// m * n  m > 1e04; n > 1e05~1e06
// time measure + objective value

OLLP_int OLLP_FASTLP
(
    OLLP_int *inEqMatBeg,     // Column pointer in CSC format of inequality constraint matrix
    OLLP_int *inEqMatIdx,     // Row indices in CSC format of inequality constraint matrix
    double *inEqMatElem,      // Data in CSC format of inequality constraint matrix
    OLLP_int nInEqRow,        // Number of rows of inequality matrix
    double *inEqMatRHS,       // Right-hand-side vector of inequality constraints
    OLLP_int nCol,            // Dimension of optimization variable
    double *Colobj,           // Linear coefficient of objective function
    OLLP_int boostingParam,   // Boosting parameter of online LP algorithm
    double *sol,              // Solution to OLLP
    double *y,                // Price vector
    double *binsol,           // Solution to Online mixed integer program
    OLLP_int CheckInnerFeas,  // Choose whether to check inner feasibility
    double Xmax               // An optional parameter for solving general integer program
)
{
    OLLP_int retcode = OLLP_RETCODE_OK;
        
    /* Initialize size information */
    OLLP_int m = nInEqRow;
    OLLP_int n = nCol;
    
    int nUpdate = 0;
    int *lastUpdate = NULL;
    lastUpdate = (int *) malloc(m * sizeof(int));
    for (int i = 0; i < m; ++i) {
        lastUpdate[i] = 0;
    }
    
    
    // Get inequality RHS vector b and copy it to resource
    double *b = inEqMatRHS;
    double *resource = NULL;
    resource = (double *) malloc(m * sizeof(double));
    
    // Get linear coefficient c in objective function
    double *c = Colobj;
    
    // Define and initialize vector d: d = b / n
    double *d = NULL;
    d = (double *) malloc(m * sizeof(double));
    
    for (int i = 0; i < m; ++i) {
        d[i] = b[i] / n;
        resource[i] = b[i];
    }
    
    // Initialize stepsize
    double step = 1 / sqrt(boostingParam * n);
    
    // Initialize solution and price vector with 0
    memset(sol, 0, sizeof(double) * n);
    
    if (binsol) {
        memset(binsol, 0, sizeof(double) * n);
    }
    
    memset(y, 0, sizeof(double) * m);
    
    // Initialize remaining inventory: br = K * b
    double *br = NULL;
    br = (double *) malloc(m * sizeof(double));
    
    for (int i = 0; i < m; ++i) {
        br[i] = b[i] * boostingParam;
    }
    
    /* Start computation routine */
    /* Start the outer loop */
    
    // Declare random permutation
    csi *p = NULL;
    
    // Initialize auxiliary variables
    double xk = 0;
    double minbr = 0;
    
    for (int k = 0; k < boostingParam; ++k) {
        
        // Generate random permutation
        srand((unsigned int) time(NULL) + (unsigned int) rand());
        p = cs_randperm(n, rand() +(unsigned int) time(NULL));
        
        /* Start inner loop */
        for (int i = 0; i < n; ++i) {
            
            // j here stands for ii
            OLLP_int j = p[i];
            xk = 0.0;
            
            for (int q = inEqMatBeg[j]; q < inEqMatBeg[j + 1]; ++q) {
                int s = inEqMatIdx[q];
                double a = inEqMatElem[q];
                /* Update y[s] till this point */
                y[s] = y[s] - step * d[s] * (nUpdate - lastUpdate[s]);
                y[s] = OLLP_MAX(y[s], 0.0);
                lastUpdate[s] = nUpdate;
            
                /* Compute xk */
                xk = xk + a * y[s];
            }
            
            if (Xmax == 1.0) {
                xk = (c[j] >= xk);
            } else {
                xk = OLLP_MIN(Xmax, floor(c[j] / xk));
            }
            
        
            // Update the dual solution: y=max(0,y+step*(xk*aa-d));
            nUpdate++;
            if (xk) {
                for (int q = inEqMatBeg[j]; q < inEqMatBeg[j + 1]; ++q) {
                    int s = inEqMatIdx[q];
                    double a = inEqMatElem[q];
                    y[s] = y[s] + step * (a - d[s]);
                    y[s] = OLLP_MAX(y[s], 0.0);
                    lastUpdate[s] = nUpdate;
                }
            }
            
            // Update the remaining inventory and primal solution
            if (xk) {
                minbr = 0;
                for (int q = inEqMatBeg[j]; q < inEqMatBeg[j + 1]; ++q) {
                    minbr = OLLP_MIN(minbr, br[inEqMatIdx[q]] - inEqMatElem[inEqMatIdx[q]]);
                }
            } else {
                minbr = - 1;
            }
            
            if ((!CheckInnerFeas) || (minbr >= 0)) {
                for (int q = inEqMatBeg[j]; q < inEqMatBeg[j + 1]; ++q) {
                    br[inEqMatIdx[q]] = br[inEqMatIdx[q]] - inEqMatElem[q];
                }
                sol[j] = sol[j] + xk;
            }
        }
    }
    
    for (int i = 0; i < n; ++i) {
        sol[i] = sol[i] / boostingParam;
    }
    
    if (binsol) {
        if (Xmax == 1.0) {
            /* Perform rounding after permutation to get an approximate solution */
            
            OLLP_FREE(p);
            p = cs_randperm(n, rand() +(unsigned int) time(NULL));
            
            for (int i = 0; i < n; ++i) {
                int SetTrue = 0;
                
                /* Generate Bernoulli random variable */
                 srand((unsigned int) time(NULL) + (unsigned int) rand());
                 double randnum = (double) rand() / RAND_MAX;
                 if (randnum <= sol[p[i]]) {
                    SetTrue = 1;
                }
                
                /*
                if (sol[p[i]] >= 0.5) {
                    SetTrue = 1;
                }
                */
                
                if (SetTrue) {
                    for (int s = inEqMatBeg[p[i]]; s < inEqMatBeg[p[i] + 1]; ++s) {
                        if (resource[inEqMatIdx[s]] < inEqMatElem[s]) {
                            SetTrue = 0;
                        }
                    }
                    
                    if (SetTrue) {
                        for (int s = inEqMatBeg[p[i]]; s < inEqMatBeg[p[i] + 1]; ++s) {
                            resource[inEqMatIdx[s]] -= inEqMatElem[s];
                        }
                        binsol[p[i]] = 1.0;
                    }
                }
            }
        } else {
            
            /* Perform sorting after permutation to get an approximate solution */
            
            /* Sort solution in descending order */
            OLLP_int *idx = NULL;
            idx = (OLLP_int *) malloc(sizeof(OLLP_int) * n);
            
            double *sortsol = (double *) malloc(sizeof(double) * n);
            
            for (int i = 0; i < n; ++i) {
                sortsol[i] = sol[i];
                idx[i] = i;
            }
            
            // Sort the solution array
            QSort(idx, sortsol, 0, n);
            
            for (int i = n - 1; i >= 0; --i) {
                double SetTrue = 0;

                double rdx = 0;
                double intpart = floor(sol[idx[i]]);
                double digitpart = sol[idx[i]] - intpart;
                
                srand((unsigned int) time(NULL) + (unsigned int) rand());
                double randnum = (double) rand() / RAND_MAX;

                if (randnum <= digitpart) {
                    SetTrue = 1;
                }

                rdx = intpart + SetTrue;

                /* Remember to memset to clear non-zero values in aa !!!*/
                
                
                if (retcode != OLLP_RETCODE_OK) {
                    retcode = OLLP_RETCODE_INVALID;
                    OLLP_FREE(idx);
                    OLLP_FREE(sortsol);
                    goto exit_cleanup;
                }

                for (int j = inEqMatBeg[i]; j < inEqMatBeg[i + 1]; ++j) {
                    rdx = OLLP_MIN(rdx, (double) floor(resource[j] / (inEqMatElem[j] + 1e-10)));
                }

                rdx = OLLP_MAX(rdx, 0.0);

                if (rdx >= 1) {
                    for (int j = 0; j < m; ++j) {
                        resource[j] -= rdx * inEqMatElem[j];
                    }
                }

                binsol[idx[i]] = rdx;

            }
            
            OLLP_FREE(idx);
            OLLP_FREE(sortsol);
        }
    }
    
exit_cleanup:
    if (retcode != OLLP_RETCODE_OK) {
        OLLP_FREE(p);
        OLLP_FREE(d);
    }

    OLLP_FREE(br);
    
    return retcode;
}

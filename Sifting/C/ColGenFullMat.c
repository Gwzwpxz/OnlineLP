#include "ColGen.h"

int *randperm (int n, int seed)
{
    int *p, k, j, t;
    if (seed == 0) {
        return (NULL);
    }                                   /* return p = NULL (identity) */
    p = malloc(n * sizeof(int));       /* allocate result */
    if (!p) {
        return (NULL);
    }                                   /* out of memory */
    for (k = 0 ; k < n ; k++) {
        p[k] = n - k - 1 ;
    }
    if (seed == -1) {
        return (p);
    }                                   /* return reverse permutation */
    srand (seed);                       /* get new random number seed */
    for (k = 0; k < n; k++) {
        j = k + (rand() % (n - k));     /* j = rand integer in range k to n-1 */
        t = p[j];                       /* swap p[k] and p[j] */
        p[j] = p[k];
        p[k] = t;
    }

    return (p);
}

int ColGenFullMat
(
    int     *matBeg,        // Column pointer in CSC format of constraint matrix
    int     *matIdx,        // Row indices in CSC format of constraint matrix
    double  *matElem,       // Data in CSC format of constraint matrix
    int     nRow,           // Number of rows of constraint matrix
    int     *isRowEq,       // An array that determines constraint type
    double  *rhs,           // Right-hand-side vector of inequality constraints
    int     nCol,           // Dimension of optimization variable
    double  *colObj,        // Linear coefficient of objective function
    int     boostingParam,  // Boosting parameter of online LP algorithm
    double  *sol,           // Approximate solution for evaluation
    double  *y,             // Dual vector
    int     sense           // Objective Sense

/*

Note: we follow the same naming style as in Matlab code

i.e., we have the optimization problem given by
        maximize     c' * x
        subject to
                      A * x <= b
                    Aeq * x  = beq
                     0 <= x <= 1

 and we should ensure that b > 0 (or < 0) for some elements.

*/
) {
    /* Initialization */
    int m = nRow;
    int n = nCol;

    /* Count number of updates and track the last time the update was applied */
    int nUpdate = 0;
    int *lastUpdate = NULL;
    lastUpdate = (int *) malloc(m * sizeof(int));
    for (int i = 0; i < m; ++i) {
        lastUpdate[i] = 0;
    }
    // Define and initialize vector d: d = b / n
    double *d = NULL;
    d = (double*) malloc(m * sizeof(double));

    for (int i = 0; i < m; ++i) {
        d[i] = rhs[i] / n;
    }

    // Initialize stepsize
    double step = 1 / sqrt(boostingParam * n);

    // Initialize solution and price vector
    memset(sol, 0, sizeof(double) * n);
    
    // Initialize overall dual vector y
    if (sense == COLGEN_POSITIVE) {
        for (int i = 0; i < m; ++i) {
            y[i] = 0.0;
        }
    } else {
        for (int i = 0; i < m; ++i) {
            y[i] = 1.0;
        }
    }

    /* Start computation routine */
    /* Start the outer loop */

    // Declare random permutation
    int *p = NULL;

    // Initialize auxiliary variables
    double xk = 0.0;

    for (int k = 0; k < boostingParam; ++k) {

        // Generate random permutation
        srand((unsigned int) time(NULL));
        int seed = rand();
        p = randperm(n, seed + 1);
        srand((unsigned int) seed);
        
        /* Start inner loop */
        for (int i = 0; i < n; ++i) {

            // j here stands for ii
            int j = p[i];

            /* Update any y[s] that is required for computing xk */
            xk = 0.0;
            if (matElem) {
                for (int q = matBeg[j]; q < matBeg[j + 1]; ++q) {
                    int s = matIdx[q];
                    double a = matElem[q];
                    
                    /* Update y[s] till this point */
                    y[s] = y[s] - step * d[s] * (nUpdate - lastUpdate[s]);
                    
                    if (!isRowEq[s]) {
                        y[s] = MAX(y[s], 0.0);
                    }
                    
                    lastUpdate[s] = nUpdate;
                
                    /* Compute xk */
                    xk = xk + a * y[s];
                }
            }

            // Estimate xk
            if (colObj[j] - xk >= 0) {
                xk = 1.0;
            } else {
                xk = 0.0;
            }

            /* Update y[s] by one step with non-zero xk */
            nUpdate++;
            if (matElem && (xk != 0.0)) {
                for (int q = matBeg[j]; q < matBeg[j + 1]; ++q) {
                    int s = matIdx[q];
                    double a = matElem[q];
                    y[s] = y[s] + step * (a - d[s]);
                    
                    if (!isRowEq[s]) {
                        y[s] = MAX(y[s], 0.0);
                    }
                    
                    lastUpdate[s] = nUpdate;
                }
            }
            /* Update sol with xk */
            sol[j] = sol[j] + xk;
        }

        if (p) {
            FREE(p);
        }
    }

    // Output online algorithm estimation
    for (int i = 0; i < n; ++i) {
        sol[i] = sol[i] / boostingParam;
    }

    if (d) {
        FREE(d);
    }

    if (p) {
        FREE(p);
    }

    if (lastUpdate) {
        FREE(lastUpdate);
    }

    return 0;
}

#ifndef ColGen_h
#define ColGen_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int *randperm (int n, int seed);

int ColGen
(
    int     *inEqMatBeg,    // Column pointer in CSC format of inequality constraint matrix
    int     *inEqMatIdx,    // Row indices in CSC format of inequality constraint matrix
    double  *inEqMatElem,   // Data in CSC format of inequality constraint matrix
    int     *eqMatBeg,      // Column pointer in CSC format of equality constraint matrix
    int     *eqMatIdx,      // Row indices in CSC format of equality constraint matrix
    double  *eqMatElem,     // Data in CSC format of equality constraint matrix
    int     nInEqRow,       // Number of rows of inequality matrix
    int     nEqRow,         // Number of rows of equality matrix
    double  *inEqMatRHS,    // Right-hand-side vector of inequality constraints
    double  *eqMatRHS,      // Right-hand-side vector of equality constraints
    int     nCol,           // Dimension of optimization variable
    double  *Colobj,        // Linear coefficient of objective function
    int     boostingParam,  // Boosting parameter of online LP algorithm
    double  *sol,           // Approximate solution for evaluation
    double  *y,             // Price vector
    int     sense           // Objective Sense
);


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define FREE(var) do {free((var)); (var) = NULL;} while (0)

#define COLGEN_POSITIVE 1
#define COLGEN_NEGATIVE 0

#endif



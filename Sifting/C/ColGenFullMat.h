#ifndef ColGen_h
#define ColGen_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int *randperm (int n, int seed);

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
);


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define FREE(var) do {free((var)); (var) = NULL;} while (0)

#define COLGEN_POSITIVE 1
#define COLGEN_NEGATIVE 0

#endif



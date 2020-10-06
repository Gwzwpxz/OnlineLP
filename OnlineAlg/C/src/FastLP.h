//
//  FastLP.h
//  OLLP
//
//  Created by 高文智 on 2020/7/18.
//  Copyright © 2020 高文智. All rights reserved.
//

/* FastLP_h */
#ifndef FastLP_h
#define FastLP_h

#include "Config.h"

/* Online Linear Programming Routine */

/*
 This is a computational routine for Online Linear Programming of form
 
 Maximize         c'x
 Subject to     A * x = b
               0 <= x <= 1
 
 The routine accepts inputs given by
 
 A               : (Sparse) inequality constraint matrix
 b(Non-negative) : Inequality right-hand column vector
 K(Positive)     : Boosting Parameter
 
 The routine relies on sparse matrix computation package SuiteSparse and stores
 data (especially sparse matrix A) in CSC (Column Sparse Compressed) format
 
 In CSC format, the inequality constraint matrix is stored in three arrays
 
 inEqMatBeg      : Column pointer
 inEqMatIdx      : Row indices
 inEqMatElem     : Data stored in matrices
*/

int FastLP
(
    int      *inEqMatBeg,        // Column pointer in CSC format of inequality constraint matrix
    int      *inEqMatIdx,        // Row indices in CSC format of inequality constraint matrix
    double   *inEqMatElem,       // Data in CSC format of inequality constraint matrix
    int      nInEqRow,           // Number of rows of inequality matrix
    double   *inEqMatRHS,        // Right-hand-side vector of inequality constraints
    int      nCol,               // Dimension of optimization variable
    double   *Colobj,            // Linear coefficient of objective function
    int      boostingParam,      // Boosting parameter of online LP algorithm
    double   *sol,               // Solution to OLLP
    double   *y,                 // Price vector
    double   *binsol,            // Solution to Online mixed integer program
    int      CheckInnerFeas,     // Choose whether to check inner feasibility
    double      Xmax
);

OLLP_int FastLP_l
(
    OLLP_int   *inEqMatBeg,      // Column pointer in CSC format of inequality constraint matrix
    OLLP_int   *inEqMatIdx,      // Row indices in CSC format of inequality constraint matrix
    double     *inEqMatElem,     // Data in CSC format of inequality constraint matrix
    OLLP_int   nInEqRow,         // Number of rows of inequality matrix
    double     *inEqMatRHS,      // Right-hand-side vector of inequality constraints
    OLLP_int   nCol,             // Dimension of optimization variable
    double     *Colobj,          // Linear coefficient of objective function
    OLLP_int   boostingParam,    // Boosting parameter of online LP algorithm
    double     *sol,             // Solution to OLLP
    double     *y,               // Price vector
    double     *binsol,          // Solution to Online mixed integer program
    OLLP_int   CheckInnerFeas,   // Choose whether to check inner feasibility
    double   Xmax
);


#endif

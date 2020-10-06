//
//  ollpinterface.h
//  OLLP
//
//  Created by 高文智 on 2020/7/21.
//  Copyright © 2020 高文智. All rights reserved.
//

#ifndef ollpinterface_h
#define ollpinterface_h

#include "Config.h"


#ifndef __OLLP_H__
#define __OLLP_H__

#ifdef _WIN32
#define OLLP_CALL __stdcall
#else
#define OLLP_CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ollp_prob_s ollp_prob;

OLLP_int OLLP_CALL OLLP_CreateProb(ollp_prob **p_prob);
OLLP_int OLLP_CALL OLLP_DeleteProb(ollp_prob **p_prob);

OLLP_int OLLP_LoadLP(ollp_prob *prob,
                     OLLP_int   nCol,
                     OLLP_int   nInEqRow,
                     OLLP_int   boostingParam,
                     OLLP_int   CheckInnerFeas,
                     double     Xmax,
                     double    *colObj,
                     OLLP_int  *inEqMatBeg,
                     OLLP_int  *inEqMatIdx,
                     double    *inEqMatElem,
                     double    *inEqRHS);

int OLLP_CALL OLLP_Solve(ollp_prob *prob);

int OLLP_CALL OLLP_GetSolution(ollp_prob *prob, double *colVal);
int OLLP_CALL OLLP_GetSolutionBin(ollp_prob *prob, double *colVal);
int OLLP_CALL OLLP_GetTime(ollp_prob *prob, double *soltime);
int OLLP_CALL OLLP_GetObjVal(ollp_prob *prob, double *ObjVal);
int OLLP_CALL OLLP_GetObjValbin(ollp_prob *prob, double *ObjVal);

#ifdef __cplusplus
}
#endif

#endif /* __OLLP_H__ */
#endif /* ollpinterface_h */

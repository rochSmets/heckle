
#if defined (FULLP)

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ghosts.h"
#include "hecklebc.h"
#include "bc_constants.h"
#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "smooth.h"
#include "closeModelFullPressureCommon.h"



void subCycledPressure(int it, const STI * const si, const STX * const sx,
    struct st1 *s1, struct st2 *s2, HeckleBC *hbc, Ghosts *ghosts, int ipc) {

    int ijk;
    int n0, n1, n2, n;
    int h, i, j, k, m;
    int idx3D[3];

    double pSub[6];  // P tensor @ [n-1/2] / [n+1/2]
    double cTerm[6];    // explicit part of C (Cyclotron) operator
    double iTerm[6];   // isotropization operator

    double omega;      // electron cycloron frequency

    double vecB[3]; // B field vector
    double vecBnext[3]; // B field vector for next ts
    double unitB[3];   // unit vector along B field
    double modulusB;   // modulus of B
    double *vecBstartAll; // all components of B during subcycling
    double *vecBstepAll; // all components of step for B during subcycling
    double *pSubAll;   // all components of pressure tensor during subcycling
    double *iTermAll;  // all components of isotropization term

    const double subDt = si->ts*si->me;
    const int numOfSubStep = (int)(1.0/si->me);


    /* __ # of g2 grid points __ */
    n0 = sx->n[0] + 2;
    n1 = sx->n[1] + 2;
    n2 = sx->n[2] + 2;
    n = n0*n1*n2;

    vecBstartAll = (double *)malloc(n*3*sizeof(double));
    vecBstepAll = (double *)malloc(n*3*sizeof(double));
    pSubAll = (double *)malloc(n*6*sizeof(double));
    iTermAll = (double *)malloc(n*6*sizeof(double));

    for (i = 1; i < n0 - 1; i++) {
        for (j = 1; j < n1 - 1; j++) {
            for (k = 1; k < n2 - 1; k++) {

                idx3D[0] = i;
                idx3D[1] = j;
                idx3D[2] = k;

                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ set start value of P : [n-1/2] / [n+1/2] __ */
                for (h = 0; h < 6; h++) {
                    switch (ipc) {

                        /* __ predictor __ */
                        case 0: pSub[h] = s2[ijk].pKept[h];
                        break;

                        /* __ corrector __ */
                        case 1: pSub[h] = s2[ijk].ps[0][h];
                        break;

                        /* __ no reason to get there __ */
                        default:
                            IAMDEAD(sx->r);
                    }
                    /* __ keep memory of pressure in a 6 components tensor __ */
                    pSubAll[ijk*6+h] = pSub[h];
                }
                
                
                vectorB(sx, s1, idx3D, BFORE, vecB    );
                vectorB(sx, s1, idx3D, ipc  , vecBnext);
                
                for (h = 0; h < 3; h++) {
                    vecBstartAll[ijk*3+h] = vecB[h];
                    vecBstepAll [ijk*3+h] = (vecBnext[h]-vecB[h])/numOfSubStep;
                }

                /* __ calculation of unitary B field & associated Omega __ */
                miscB(si, sx, vecB, &modulusB, unitB, &omega);
                
                /* __ set isotropization term @ [n-1/2] / [n+1/2] __ */
                setIsotropTerm(pSub, iTerm, unitB, omega);

                /* __ keep memory of isotropization term in a 6 components tensor __ */
                for (h = 0; h < 6; h++) {
                    iTermAll[ijk*6+h] = iTerm[h];
                }

            }
        }
    }

    /* __ loop over the substepping __ */
    for (m = 0; m < numOfSubStep; m++) {
        for (i = 1; i < n0 - 1; i++) {
            for (j = 1; j < n1 - 1; j++) {
                for (k = 1; k < n2 - 1; k++) {

                    idx3D[0] = i;
                    idx3D[1] = j;
                    idx3D[2] = k;

                    /* __ index on g2 __ */
                    ijk = IDX(i, j, k, n0, n1, n2);

                    for (h = 0; h < 6; h++) {
                        pSub[h] = pSubAll[ijk*6+h];
                    }
                    
                    /* __ set current b field __ */
                    for (h = 0; h < 3; h++) {
                        vecB[h] = vecBstartAll[ijk*3+h] + m*vecBstepAll[ijk*3+h];
                    }
                    /* __ calculation of unitary B field & associated Omega __ */
                    miscB(si, sx, vecB, &modulusB, unitB, &omega);
                    
                    /* __ cyclotron term __ */
                    setCyclotronTerm(pSub, cTerm, unitB);

                    /* __ set F := RHS of P(t) eq. from explicit terms __ */
                    for (h = 0; h < 6; h++) {
                        pSubAll[ijk*6+h] += subDt*(s2[ijk].dFull[h]+omega*cTerm[h]+iTerm[h]);
                    }
                }
            }
        }
    }

    /* __ set new pressure values @ [n+1/2] / [n+3/2] __ */
    for (i = 1; i < n0 - 1; i++) {
        for (j = 1; j < n1 - 1; j++) {
            for (k = 1; k < n2 - 1; k++) {
                ijk = IDX(i, j, k, n0, n1, n2);
                for (h = 0; h < 6; h++) {
                    s2[ijk].ps[0][h] = pSubAll[ijk*6+h];
                }
            }
        }
    }

    /* __ free the pointers __ */
    free(vecBstartAll);
    free(vecBstepAll);
    free(pSubAll);
    free(iTermAll);

    /* __ and apply associated boundary conditions __ */
    GhostsSendRecv(ghosts, sx, s2, GHOST_P_FULLP);
    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_P_FULLP);

    /* __ smooth the pressure every PSMOOTH __ */
    if (it%PSMOOTH == 0){
        smoothPressure(sx, s2, ghosts, hbc);
    }

    /* __ keep memory of P @ [n-1/2] __ */
    if (ipc == 0) {
        for (ijk = 0; ijk < n; ijk++) {
            for (h = 0; h < 6; h++) {
                s2[ijk].pKept[h] = s2[ijk].ps[0][h];
            }
        }
    }

    /* __ set driver term for next time step : [n] / [n+1] __ */
    setDriverTerm(si, sx, s2, ipc);

    /* __ then manage boundary conditions for this driver term __ */
    GhostsSendRecv(ghosts, sx, s2, GHOST_DRIVER_FULLP);
    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_DRIVER_FULLP);

    /* __ keep B field @ n-1/2 __ */
    if (ipc == 0) {
        for (i = 0; i < sx->n[0] + 1; i++) {
            for (j = 0; j < sx->n[1] + 1; j++) {
                for (k = 0; k < sx->n[2] + 1; k++) {
                    ijk = IDX(i, j, k, sx->n[0]+1, sx->n[1]+1, sx->n[2]+1);

                    s1[ijk].bKept[0] = s1[ijk].c[0];
                    s1[ijk].bKept[1] = s1[ijk].c[1];
                    s1[ijk].bKept[2] = s1[ijk].c[2];
                }
            }
        }
    }

}

#endif

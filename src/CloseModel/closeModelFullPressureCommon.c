
//#if defined(WINSKE)  && !defined (SUBSTEP_NUMBER)
#if defined(FULLP)
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


// Cyclotron = [ P X b + (P X b)T ]
// Cyclotron(i,j) = P(i,k)*b(n)*E(k,n,j)  +  P(j,k)*b(n)*E(k,n,i)
// E = Levi-Civita symbol
//
//     | P(0,0) = ps[0], P(0,1) = ps[1], P(0,2) = ps[2] |
// P = | P(1,0) = ps[1], P(1,1) = ps[3], P(1,2) = ps[4] |
//     | P(2,0) = ps[2], P(2,1) = ps[4], P(2,2) = ps[5] |
//
// then [*Cyclotron*] term is
//
// | 2(P01b2-P02b1), P02b0-P00b2+P11b2-P12b1, P00b1-P01b0+P21b2-P22b1 |
// | --------------,     2(P12b0-P10b2)     , P10b1-P11b0+P22b0-P20b2 |
// | --------------,------------------------,     2(P20b1-P21b0)      |
//

// P[n+1/2] - (-1)*alpha*dt*Ω*Cyclotron[n+1/2] =
// P[n-1/2] + (-1)*dt*((1-alpha)*Ω*Cyclotron[n-1/2] + Driver[n])
//
// ϰ = alpha*dt*Ω
// P(i,j) - (-1)*ϰ*Cyclotron(i,j) = F(i,j)
//
// in B aligned coordinate system cyclotron term :
//             |   0  ,    P02   ,   -P01   |
// Cyclotron = |  P02 ,   2P12   , -P11+P22 |
//             | -P01 , -P11+P22 ,  -2P21   |
//
// solve the system :
//
// P00 = F00
// P01 + ϰP02 = F01
// P02 - ϰP01 = F02
// P11 + 2ϰP12 = F11
// P12 + ϰ(-P11+P22) = F12
// P22 - 2ϰP21 = F22
//
// so the solution is :
//
// P00 = F00
// P01 = (F01 - ϰF02)/(1+ϰ^2)
// P02 = (F02 + ϰF01)/(1+ϰ^2)
// P11 = ((1+2ϰ^2)*F11 - 2ϰF12 + 2ϰ^2F22)/(1+4ϰ^2)
// P12 = (F12 - ϰ(F11+F22))/(1+4ϰ^2)
// P22 = (2ϰ^2*F11 + 2ϰF12 + (1+2ϰ^2)F22)/(1+4ϰ^2)
//

// __ full electron pressure tensor
// __ P' = - V.∇P - P∇.V - P.∇V - (P.∇V)T - q/M [ P X B + (P X B)T ]
// __ Driver = V.∇P + P∇.V + P.∇V + (P.∇V)T
// __ Cyclotron = [ P X b + (P X b)T ]
// __ P' = (-1)*Driver + (-1)*Ω*Cyclotron
// __ P[n+1/2] - (-1)*alpha*dt*Ω[n+1/2]*Cyclotron[n+1/2] =
// __ P[n-1/2] + (-1)*dt*((1-alpha)*Ω[n-1/2]*Cyclotron[n-1/2] + Driver[n])














/* _____ set driver term for next time step _____ */
void setDriverTerm(const STI * const si, const STX * const sx, struct st2 *s2,
    int ipc){

    int ijk;
    int n0, n1, n2, n;
    int i, j, k, l, m;
    int idx3D[3];
    int h;

    double *pDrive;
    double dTerms[3][3];


    /* __ # of g2 grid points __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    n = n0*n1*n2;

    /* __ driver terms for the full-field tensor  __ */
    pDrive = (double *)malloc(n*6*sizeof(double));

    /* __ electron velocity @ [n+1/2] / [n+3/2] __ */
    for (ijk = 0; ijk < n; ijk++) {
        for (l = 0; l < 3; l++) {
            s2[ijk].vs[0][l] = s2[ijk].vi[l]-s2[ijk].j_smoothed[l]/s2[ijk].ns[0];
        }
    }

    /* __ driver @ [n+1/2] / [n+3/2] __ */
    for (i = 1; i < n0 - 1; i++) {
        for (j = 1; j < n1 - 1; j++) {
            for (k = 1; k < n2 - 1; k++) {

                ijk = IDX(i, j, k, n0, n1, n2);

                idx3D[0] = i;
                idx3D[1] = j;
                idx3D[2] = k;

                /* __ set the 6 driver terms on each grid points __ */
                driverComponents(si, sx, s2, dTerms, idx3D);

                /* __ fill the full-field tensor in a 1d array __ */
                h = 0;
                for (l = 0; l < 3; l++) {
                    for (m = l; m < 3; m++) {
                        pDrive[ijk*6+h] = dTerms[l][m];
                        h++;
                    }
                }
            }
        }
    }

    /* __ extrapolate/interpolate the driver term @ [n+1] __ */
    for (i = 1; i < n0 - 1; i++) {
        for (j = 1; j < n1 - 1; j++) {
            for (k = 1; k < n2 - 1; k++) {

                ijk = IDX(i, j, k, n0, n1, n2);

                for (h = 0; h < 6; h++) {
                    switch (ipc) {
                        case 0: /* __ predictor __ */
                            s2[ijk].dFull[h] = -s2[ijk].dFull[h]+2.0*pDrive[ijk*6+h];

                            /* __ keep memory of driver term @ [n+1/2] __ */
                            s2[ijk].dHalf[h] = pDrive[ijk*6+h];
                        break;

                        case 1: /* __ corrector __ */
                            s2[ijk].dFull[h] = 0.5*(pDrive[ijk*6+h]+s2[ijk].dHalf[h]);
                        break;

                        /* __ no reason to get there __ */
                        default:
                            IAMDEAD(sx->r);
                    }
                }
            }
        }
    }

    free(pDrive);
}


/* _____ calcul unit B vector & cyclotron frequency from local B value _____ */
void vectorB(const STX * const sx, struct st1 *s1, int idx3D[3], int ipc, double componentsB[3]) {


    double sum;
    int g1IDXs[8];
    int n0, n1, n2;
    int h, n;


    /* __ # of g1 grid points __ */
    n0 = sx->n[0]+1;
    n1 = sx->n[1]+1;
    n2 = sx->n[2]+1;

    /* __ index on g1 __ */
    g1IDXs[0] = IDX(idx3D[0]  , idx3D[1]  , idx3D[2]  , n0, n1, n2);
    g1IDXs[1] = IDX(idx3D[0]-1, idx3D[1]  , idx3D[2]  , n0, n1, n2);
    g1IDXs[2] = IDX(idx3D[0]  , idx3D[1]-1, idx3D[2]  , n0, n1, n2);
    g1IDXs[3] = IDX(idx3D[0]-1, idx3D[1]-1, idx3D[2]  , n0, n1, n2);
    g1IDXs[4] = IDX(idx3D[0]  , idx3D[1]  , idx3D[2]-1, n0, n1, n2);
    g1IDXs[5] = IDX(idx3D[0]-1, idx3D[1]  , idx3D[2]-1, n0, n1, n2);
    g1IDXs[6] = IDX(idx3D[0]  , idx3D[1]-1, idx3D[2]-1, n0, n1, n2);
    g1IDXs[7] = IDX(idx3D[0]-1, idx3D[1]-1, idx3D[2]-1, n0, n1, n2);

    for (h = 0; h < 3; h++) {
        sum = 0.0;
        switch (ipc) {
            case 0: /* __ predictor : set b field @ n+1/2 __ */
                for (n = 0; n < 8; n++) {
                    sum += s1[g1IDXs[n]].c[h];
                }
            break;

            case 1: /* __ corrector : set b field @ n+3/2 __ */
                for (n = 0; n < 8; n++) {
                    sum += s1[g1IDXs[n]].b[h];
                }
            break;

            case BFORE: /* __ set b field @ n-1/2 __ */
                for (n = 0; n < 8; n++) {
                    sum += s1[g1IDXs[n]].bKept[h];
                }
            break;

            /* __ no reason to get there __ */
            default:
                IAMDEAD(sx->r);
        }
        componentsB[h] = 0.125*sum;
    }

}


/* _____ calcul unit B vector & cyclotron frequency from local B value _____ */
void miscB(const STI * const si, const STX * const sx, double vectorB[3],
    double *modulusB, double unitB[3], double *omega) {

    if (si->me == 0.0) {
        printf("si.me needs to be nonzero\n");
        printf("file %s func : %s @ line %d on node %d\n", __FILE__, __func__, __LINE__, sx->r);
        exit(-1);
    }

    int h;


    *modulusB = sqrt(vectorB[0]*vectorB[0]
                    +vectorB[1]*vectorB[1]
                    +vectorB[2]*vectorB[2]);

    /* __ if the plasma is unmagnetized __ */
    if (*modulusB < EPS8) {
        for (h = 0; h < 3; h++) {
            unitB[h] = 0.0;
        }

        *omega = 0.0;
    }

    /* __ magnetized case __ */
    else {
        for (h = 0; h < 3; h++) {
            unitB[h] = vectorB[h]/(*modulusB);
        }

        *omega = (*modulusB)/(si->me);
    }
}


/* _____ Driver = P nab.V + V.nab P + P.nab V + (P.nab V)T _______________________*/
void driverComponents(const STI * const si, const STX * const sx, struct st2 *s2,
    double dTerm[3][3], int idx3D[3]) {

    int ijk;
    double pe[3][3];
    double nabV[3][3];
    double nabP[3][3][3];
    int n0, n1, n2;
    int l, m, n;


    /* __ # of g2 grid points __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    ijk = IDX(idx3D[0], idx3D[1], idx3D[2], n0, n1, n2);

    /* __ calculate velocity grad. & pressure grad. __ */
    gradStuff(si, sx, s2, pe, nabV, nabP, idx3D);

    /* __ divergence of velocity __ */
    double divV = nabV[0][0]+nabV[1][1]+nabV[2][2];

    for (l = 0; l < 3; l++) {
        for (m = l; m < 3; m++) {

            /* __ P nabla . V __ */
            dTerm[l][m] = -pe[l][m]*divV;

            for (n = 0; n < 3; n++) {
                /* __ V . nabla P __ */
                dTerm[l][m] -= s2[ijk].vs[0][n]*nabP[n][l][m];
                /* __ P . nabla V __ */
                dTerm[l][m] -= pe[l][n]*nabV[n][m];
                /* __ P . nabla V (transposed) __ */
                dTerm[l][m] -= pe[m][n]*nabV[n][l];
            }
            /* __ & set symmetrical terms __ */
            dTerm[m][l] = dTerm[l][m];
        }
    }
}


/* _____ calculate velocity grad. & pressure grad. _____ */
void gradStuff(const STI * const si, const STX * const sx, struct st2 *s2,
    double pe[3][3], double nabV[3][3], double nabP[3][3][3], int idx3D[3]) {

    int ijk;
    float upwindIDX[3][3];
    float sign;
    int diffIDX[3][2];
    int n0, n1, n2;
    int h, l, m, s;


    /* __ index on g2 __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    ijk = IDX(idx3D[0], idx3D[1], idx3D[2], n0, n1, n2);

    /* __ indices to make difference (on 2 grid mesh) in each dir __ */
    diffIDX[0][0] = IDX(idx3D[0]+1, idx3D[1]  , idx3D[2]  , n0, n1, n2);
    diffIDX[0][1] = IDX(idx3D[0]-1, idx3D[1]  , idx3D[2]  , n0, n1, n2);
    diffIDX[1][0] = IDX(idx3D[0]  , idx3D[1]+1, idx3D[2]  , n0, n1, n2);
    diffIDX[1][1] = IDX(idx3D[0]  , idx3D[1]-1, idx3D[2]  , n0, n1, n2);
    diffIDX[2][0] = IDX(idx3D[0]  , idx3D[1]  , idx3D[2]+1, n0, n1, n2);
    diffIDX[2][1] = IDX(idx3D[0]  , idx3D[1]  , idx3D[2]-1, n0, n1, n2);

    /* __ terms of nabla V, centered in space __ */
    for (l = 0; l < 3; l++) {//index of the spatial component of the derivative
        for (m = 0; m < 3; m++) {//index of the spatial component of the velocity
            nabV[l][m] = 0.5*(s2[diffIDX[l][0]].vs[0][m]
                             -s2[diffIDX[l][1]].vs[0][m])/si->dl[l];
        }
    }

    /* __ set the upwind coefficients, depending on the flow direction __ */
    for (l = 0; l < 3; l++) {//index of the spatial component of the velocity

        /* __ sign of the electron velocity __ */
        sign = copysignf(1, s2[ijk].vs[0][l]);

        upwindIDX[l][0] = +0.5*(1.0-sign)/si->dl[l]; // plus
        upwindIDX[l][1] =           sign /si->dl[l]; //center
        upwindIDX[l][2] = -0.5*(1.0+sign)/si->dl[l]; //minus

    }

    /* __ then use these coefficients for nabla P __ */
    h = 0;
    for (l = 0; l < 3; l++) {//index 1 dimension
        for (m = l; m < 3; m++) {//index 2 dimension

            pe[l][m] = s2[ijk].ps[0][h];
            pe[m][l] = pe[l][m];

            for (s = 0; s < 3; s++) {//index of the spatial component of the derivative
                nabP[s][l][m] = upwindIDX[s][0]*s2[diffIDX[s][0]].ps[0][h]
                               +upwindIDX[s][1]*pe[l][m]
                               +upwindIDX[s][2]*s2[diffIDX[s][1]].ps[0][h];

                /* __ this tensor is symmetrical (2 last indices) __ */
                nabP[s][m][l] = nabP[s][l][m];
            }

            h++;
        }
    }
}

/* _____ calculate the cyclotron term from P tensor & B field _____ */
void setCyclotronTerm(double ps[6], double cyclotron[6], double unitB[3]) {
    int h;
    double bModulus;

    bModulus = sqrt(unitB[0]*unitB[0] + unitB[1]*unitB[1] + unitB[2]*unitB[2]);
    /* __ if the plasma is unmagnetized __ */
    if (bModulus < EPS8) {
        for (h = 0; h < 6; h++) {
            cyclotron[h] = 0.0;
        }
    }
    else {
        cyclotron[0] = -( 2.0*(ps[1]*unitB[2]-ps[2]*unitB[1]) );
        cyclotron[1] = -( ps[2]*unitB[0]-ps[0]*unitB[2]+ps[3]*unitB[2]-ps[4]*unitB[1] );
        cyclotron[2] = -( ps[0]*unitB[1]-ps[1]*unitB[0]+ps[4]*unitB[2]-ps[5]*unitB[1] );
        cyclotron[3] = -( 2.0*(ps[4]*unitB[0]-ps[1]*unitB[2]) );
        cyclotron[4] = -( ps[1]*unitB[1]-ps[3]*unitB[0]+ps[5]*unitB[0]-ps[2]*unitB[2] );
        cyclotron[5] = -( 2.0*(ps[2]*unitB[1]-ps[4]*unitB[0]) );
    }
}


/* _____ calculate isotropization term from P tensor */
void setIsotropTerm(double ps[6], double iTerm[6], double unitB[3],  double omega) {

    double traceP = (ps[0]+ps[3]+ps[5])/3;
    int h;
    double bModulus;


    bModulus = sqrt(unitB[0]*unitB[0] + unitB[1]*unitB[1] + unitB[2]*unitB[2]);

    /* __ if the plasma is unmagnetized __ */
    if (bModulus < EPS8) {
        for (h = 0; h < 6; h++) {
           iTerm[h] = 0.0;
        }
    }
    else {
        double buf = omega/TAU;

        iTerm[0] = -buf*(ps[0] - traceP);
        iTerm[1] = -buf*ps[1];
        iTerm[2] = -buf*ps[2];
        iTerm[3] = -buf*(ps[3] - traceP);
        iTerm[4] = -buf*ps[4];
        iTerm[5] = -buf*(ps[5] - traceP);
    }
}


#endif

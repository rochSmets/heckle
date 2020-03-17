
//#if defined(WINSKE)  && !defined (SUBSTEP_NUMBER)
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
#include "closeModelFullPressureImplicit.h"


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
// P[n-1/2] + (-1)*dt*((1-alpha)*Ω*Cyclotron[n-1/2] + Driver[n])_*/
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















void fromF2P(double[3][3], double[3][3], double);

void transformMatrix(double[3][3], double[3][3], double[3][3], int);


void implicitPressure(int it, const STI * const si, const STX * const sx,
    struct st1 *s1, struct st2 *s2, HeckleBC *hbc, Ghosts *ghosts, int ipc) {

    int ijk;
    int n0, n1, n2, n;
    int h, i, j, k, l, m;
    int idx3D[3];

    double fXYZ[3][3]; // F term : RHS of eq. for P tensor
    double pStart[6];  // P tensor @ [n-1/2] / [n+1/2]
    double cTerm[6];    // explicit part of C (Cyclotron) operator
    double iTerm[6];   // isotropization operator
    double fUVW[3][3]; // F term in B aligned coordinates
    double pUVW[3][3]; // P term in B aligned coordinates
    double pXYZ[3][3]; // P term in lab coordinates

    double omega;      // electron cycloron frequency
    double kappa;      //  := ALPHA . omega . ts

    double vecB[3]; // B field vector
    double unitB[3];   // unit vector along B field
    double modulusB;   // modulus of B
    double uvw[3][3];  // transition matrix from cartesian to B aligned coordinates
    double *pNew;      // pressure tensor after time integration


    /* __ # of g2 grid points __ */
    n0 = sx->n[0] + 2;
    n1 = sx->n[1] + 2;
    n2 = sx->n[2] + 2;
    n = n0*n1*n2;

    pNew = (double *)malloc(n*6*sizeof(double));

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
                        case 0: pStart[h] = s2[ijk].pKept[h];
                        break;

                        /* __ corrector __ */
                        case 1: pStart[h] = s2[ijk].ps[0][h];
                        break;

                        /* __ no reason to get there __ */
                        default:
                            IAMDEAD(sx->r);
                    }
                }

                /* __ calculation of unitary B field & associated Omega __ */
                vectorB(sx, s1, idx3D, BFORE, vecB);
                miscB(si, sx, vecB, &modulusB, unitB, &omega);

                setCyclotronTerm(pStart, cTerm, unitB);

                /* __ set isotropization term @ [n-1/2] / [n+1/2] __ */
                setIsotropTerm(pStart, iTerm, unitB, omega);

                /* __ set F := RHS of P(t) eq. from explicit terms __ */
                h = 0;
                for (l = 0; l < 3; l++) {
                    for (m = l; m < 3; m++) {
                        fXYZ[l][m] = pStart[h]
                            +si->ts*s2[ijk].dFull[h]
                            +si->ts*(1-ALPHA)*omega*cTerm[h]
                            +si->ts*iTerm[h];
                        fXYZ[m][l] = fXYZ[l][m];
                        h++;
                    }
                }

                /* __ then need the implicit part of cyclotron operator __ */
                vectorB(sx, s1, idx3D, ipc, vecB);
                miscB(si, sx, vecB, &modulusB, unitB, &omega);

                /* __ kappa value depends on local value of magnetic field __ */
                kappa = ALPHA*omega*si->ts;

                /* __ get magnetic field aligned basis __ */
                ortho(unitB, uvw);

                /* __ transform to magnetic field aligned coordinates __ */
                transformMatrix(fXYZ, fUVW, uvw, FORWARD);

                /* __ set P terms from F terms __ */
                fromF2P(fUVW, pUVW, kappa);

                /* __ transform back to lab coordinates __ */
                transformMatrix(pUVW, pXYZ, uvw, BACKWARD);

                /* __ convert back 2d to 1d array __ */
                h = 0;
                for (l = 0; l < 3; l++) {
                    for (m = l; m < 3; m++) {
                        pNew[ijk*6+h] = pXYZ[l][m];
                        h++;
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
                    s2[ijk].ps[0][h] = pNew[ijk*6+h];
                }
            }
        }
    }

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

    free(pNew);
}


/* _____ calculate P tensor from F tensor (in B aligned coordinate) _____ */
void fromF2P(double F[3][3], double P[3][3], double k) {

    double k2 = k*k;
    double buff1 = 1.0+k2;
    double buff2 = 1.0+4.0*k2;


    P[0][0] = F[0][0];
    P[0][1] = (F[0][1]-k*F[0][2])/buff1;
    P[0][2] = (F[0][2]+k*F[0][1])/buff1;
    P[1][0] = P[0][1];
    P[1][1] = (F[1][1]*(1.0+2.0*k2)-2.0*k*F[1][2]+2.0*k2*F[2][2])/buff2;
    P[1][2] = (k*F[1][1]+F[1][2]-k*F[2][2])/buff2;
    P[2][0] = P[0][2];
    P[2][1] = P[1][2];
    P[2][2] = (2.0*k2*F[1][1]+2.0*k*F[1][2]+(1.0+2.0*k2)*F[2][2])/buff2;

}


/* _____ transform the tensor "old" 2 "new" using a transition matrix _____ */
void transformMatrix(double old[3][3], double new[3][3], double transit[3][3],
    int way) {

    int k, l, m, n;


    for (k = 0; k < 3; k++) {
        for (l = 0; l < 3; l++) {

            /* __ new tensor after transition __ */
            new[k][l] = 0.0;

            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    /* __ backward transform : from B aligned to XYZ __ */
                    if (way == BACKWARD) {
                        new[k][l] += transit[m][k]*old[m][n]*transit[n][l];

                    /* __ forward transform : from XYZ to B aligned __ */
                    } else if (way == FORWARD) {
                        new[k][l] += transit[k][m]*old[m][n]*transit[l][n];
                    }

                    else {
                        printf("file %s func : %s @ line %d\n", __FILE__, __func__, __LINE__);
                    }
                }
            }
        }
    }
}

#endif

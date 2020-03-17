#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "smooth.h"
//#include "closeModelFullPressureCommon.h"

//#define PSMOOTH 8



/* __ full electron pressure tensor (ag) ______________________________ */
void polytropPressure(int it, const STI * const si, const STX * const sx, struct st1 *s1, struct st2 *s2,
        HeckleBC   *hbc, Ghosts     *ghosts, int ipc)
{

    ST4 *s4a;
    double *pw, *tw4;
    double dlx, dly, dlz;
    double dlx2, dly2, dlz2;
    int nn2, nn4;
    int ijk, ijkw;
    int ijka, ijkb, ijkc, ijkd, ijke, ijkf;
    int n0, n1, n2, o0, o1, o2, p0, p1, p2;
    int i, j, k, l;


    /* __ # of grid points __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;
    o0 = sx->n[0]+3;
    o1 = sx->n[1]+3;
    o2 = sx->n[2]+3;
    p0 = sx->n[0]+4;
    p1 = sx->n[1]+4;
    p2 = sx->n[2]+4;

    /* __ needed for the lax-wendroff scheme __ */
    dlx = 0.5*si->ts/si->dl[0];
    dly = 0.5*si->ts/si->dl[1];
    dlz = 0.5*si->ts/si->dl[2];

    dlx2 = dlx*dlx;
    dly2 = dly*dly;
    dlz2 = dlz*dlz;

    /* __ # of grid points on g2 __ */
    nn2 = n0*n1*n2;

    /* __ # of grid points on g4 __ */
    nn4 = p0*p1*p2;

    s4a = malloc(nn4 * sizeof *s4a);

    tw4 = (double *)malloc(nn4*sizeof(double));
    pw  = (double *)malloc(nn2*sizeof(double));

    for (ijk = 0; ijk < nn2; ijk++) {
        for (l = 0; l < 3; l++) {
            s2[ijk].vs[0][l] = s2[ijk].vi[l]-s2[ijk].j_smoothed[l]/s2[ijk].ns[0];
        }
        pw[ijk] = 0.0;
    }

    for (ijkw = 0; ijkw < nn4; ijkw++) {
        s4a[ijkw].n    = 0.0;
        s4a[ijkw].v[0] = 0.0;
        s4a[ijkw].v[1] = 0.0;
        s4a[ijkw].v[2] = 0.0;

        tw4[ijkw] = 0.0;
    }

    for (i = 0; i < n0; i++) {
        for (j = 0; j < n1; j++) {
            for (k = 0; k < n2; k++) {
                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ index on g4 __ */
                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

                /* __ fill the buffer on g4 __ */
                s4a[ijkw].n    = log(s2[ijk].ps[0][0]*pow(s2[ijk].ms[0], -GAMMA));
                for (l = 0; l < 3; l++) {
                    s4a[ijkw].v[l] = s2[ijk].vs[0][l];
                }
            }
        }
    }


    GhostsSendRecv(ghosts, sx, s4a, GHOST_NV);

    /* __ set tw4 with lax-wendroff : nested loops on the grid points of subdomain __ */
    for (i = 1; i < o0; i++)
    {
        for (j = 1; j < o1; j++)
        {
            for (k = 1; k < o2; k++)
            {
                /* __ index on g2 __ */
                ijkw = IDX(i, j, k, p0, p1, p2);

                /* __ indexes on g2 __ */
                ijka = IDX(i+1, j  , k  , p0, p1, p2);
                ijkb = IDX(i-1, j  , k  , p0, p1, p2);
                ijkc = IDX(i  , j+1, k  , p0, p1, p2);
                ijkd = IDX(i  , j-1, k  , p0, p1, p2);
                ijke = IDX(i  , j  , k+1, p0, p1, p2);
                ijkf = IDX(i  , j  , k-1, p0, p1, p2);

                /* __ new entropy with lw __ */
                tw4[ijkw] = (1.0 - dlx2*s4a[ijkw].v[0]*s4a[ijkw].v[0]
                                 - dly2*s4a[ijkw].v[1]*s4a[ijkw].v[1]
                                 - dlz2*s4a[ijkw].v[2]*s4a[ijkw].v[2])*s4a[ijkw].n
                +0.5*(dlx2 * s4a[ijka].v[0] * s4a[ijka].v[0] - dlx * s4a[ijka].v[0]) * s4a[ijka].n
                +0.5*(dlx2 * s4a[ijkb].v[0] * s4a[ijkb].v[0] + dlx * s4a[ijkb].v[0]) * s4a[ijkb].n
                +0.5*(dly2 * s4a[ijkc].v[1] * s4a[ijkc].v[1] - dly * s4a[ijkc].v[1]) * s4a[ijkc].n
                +0.5*(dly2 * s4a[ijkd].v[1] * s4a[ijkd].v[1] + dly * s4a[ijkd].v[1]) * s4a[ijkd].n
                +0.5*(dlz2 * s4a[ijke].v[2] * s4a[ijke].v[2] - dlz * s4a[ijke].v[2]) * s4a[ijke].n
                +0.5*(dlz2 * s4a[ijkf].v[2] * s4a[ijkf].v[2] + dlz * s4a[ijkf].v[2]) * s4a[ijkf].n;
                /* __ term associated with the electron thermal conductivity __ */
                tw4[ijkw] += 0.0;
            }
        }
    }

    /* __ set pw on g2 fields from g4 fields : nested loops on the grid points of subdomain __ */
    for (i = 0; i < n0; i++)
    {
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n2; k++)
            {
                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ index on g4 __ */
                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

                /* __  __ */
                pw[ijk] = exp(tw4[ijkw]) * pow( s2[ijk].ns[0], GAMMA );
            }
        }
    }

    for (ijk = 0; ijk < nn2; ijk++)
    {
        s2[ijk].ps[0][0] = pw[ijk];
        s2[ijk].ps[0][1] = 0.0;
        s2[ijk].ps[0][2] = 0.0;
        s2[ijk].ps[0][3] = pw[ijk];
        s2[ijk].ps[0][4] = 0.0;
        s2[ijk].ps[0][5] = pw[ijk];
    }

    if (it%PSMOOTH == 0){
        smoothPressure(sx, s2, ghosts, hbc);
    }
    /* __ clean-up the pointers __ */
    free(tw4);
    free(s4a);
    free(pw);
}



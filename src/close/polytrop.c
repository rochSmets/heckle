
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include "structures.h"
#include "defines.h"
#include "fill.h"
#include "smooth.h"


/* __ full electron stress tensor with polytrop hypothesis & therm. cond. ___ */
void stress(struct sti si,
            struct stx sx,
            struct st2 *s2,
            int ipc,
            MPI_Comm com)
{
double *sw, *pw, *uw, *vw, *ww;
double *sw4, *tw4, *uw4, *vw4, *ww4;
double dlx, dly, dlz;
int nn2, nn4;
int ijk, ijkw;
int ijka, ijkb, ijkc, ijkd, ijke, ijkf;
int n0, n1, n2, o0, o1, o2, p0, p1, p2;
int i, j, k;
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ # of grid points __ */
n0 = sx.n[0]+2;
n1 = sx.n[1]+2;
n2 = sx.n[2]+2;
o0 = sx.n[0]+3;
o1 = sx.n[1]+3;
o2 = sx.n[2]+3;
p0 = sx.n[0]+4;
p1 = sx.n[1]+4;
p2 = sx.n[2]+4;

/* __ needed for the lax-wendroff scheme __ */
dlx = 0.5*si.ts/si.dl[0];
dly = 0.5*si.ts/si.dl[1];
dlz = 0.5*si.ts/si.dl[2];

/* __ # of grid points on g2 __ */
nn2 = n0*n1*n2;

/* __ # of grid points on g4 __ */
nn4 = p0*p1*p2;

/* __ memory allocation __ */
sw = (double *)malloc(nn2*sizeof(double));
pw = (double *)malloc(nn2*sizeof(double));
uw = (double *)malloc(nn2*sizeof(double));
vw = (double *)malloc(nn2*sizeof(double));
ww = (double *)malloc(nn2*sizeof(double));

/* __ memory allocation __ */
sw4 = (double *)malloc(nn4*sizeof(double));
tw4 = (double *)malloc(nn4*sizeof(double));
uw4 = (double *)malloc(nn4*sizeof(double));
vw4 = (double *)malloc(nn4*sizeof(double));
ww4 = (double *)malloc(nn4*sizeof(double));

/* __ set electron fluid density @ n+1/2 (ipc=0), n+3/2 (ipc=1)__ */
for (i = 0; i < n0; i++)
    {
    for (j = 0; j < n1; j++)
        {
        for (k = 0; k < n2; k++)
            {
            /* __ index on g2 __ */
            ijk = IDX(i, j, k, n0, n1, n2);

            /* __ set the electron fluid velocity __ */
            uw[ijk] = s2[ijk].vs[0][0];
            vw[ijk] = s2[ijk].vs[0][1];
            ww[ijk] = s2[ijk].vs[0][2];
            }
        }
    }

/* __ smooth electron fluid velocity __ */
//smooth(si, sx, uw, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);
//smooth(si, sx, vw, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);
//smooth(si, sx, ww, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);

/* __ duplicate g2 fields in g4 fields : nested loops on the grid points of subdomain __ */
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

            /* __ set the "entropy" @ n-1/2 (ipc=0), n+1/2 (ipc=1) __ */
            sw[ijk] = log(s2[ijk].ps[0][0]*pow(s2[ijk].ms[0], -GAMMA));

            /* __ fill the buffer on g4 (velocity is @ 1 time step further)__ */
            sw4[ijkw] = sw[ijk];
            uw4[ijkw] = uw[ijk];
            vw4[ijkw] = vw[ijk];
            ww4[ijkw] = ww[ijk];
            }
        }
    }

/* __ fill the buffers on g4 __ */
fillg4(si, sx, sw4, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);
fillg4(si, sx, uw4, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);
fillg4(si, sx, vw4, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);
fillg4(si, sx, ww4, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);

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
            tw4[ijkw] = (1.0-dlx*dlx*uw4[ijkw]*uw4[ijkw]
                            -dly*dly*vw4[ijkw]*vw4[ijkw]
                            -dlz*dlz*ww4[ijkw]*ww4[ijkw])*sw4[ijkw]
                       +0.5*(dlx*dlx*uw4[ijka]*uw4[ijka]
                            -dlx*uw4[ijka])*sw4[ijka]
                       +0.5*(dlx*dlx*uw4[ijkb]*uw4[ijkb]
                            +dlx*uw4[ijkb])*sw4[ijkb]
                       +0.5*(dly*dly*vw4[ijkc]*vw4[ijkc]
                            -dly*vw4[ijkc])*sw4[ijkc]
                       +0.5*(dly*dly*vw4[ijkd]*vw4[ijkd]
                            +dly*vw4[ijkd])*sw4[ijkd]
                       +0.5*(dlz*dlz*ww4[ijke]*ww4[ijke]
                            -dlz*ww4[ijke])*sw4[ijke]
                       +0.5*(dlz*dlz*ww4[ijkf]*ww4[ijkf]
                            +dlz*ww4[ijkf])*sw4[ijkf];
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
            pw[ijk] = exp(tw4[ijkw])*pow(s2[ijk].ns[0],GAMMA);
            }
        }
    }

/* __ smooth stress tensor __ */
//smooth(si, sx, pw, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);

/* __ set ps : nested loops on the grid points of subdomain __ */
for (i = 0; i < n0; i++)
    {
    for (j = 0; j < n1; j++)
        {
        for (k = 0; k < n2; k++)
            {
            /* __ index on g2 __ */
            ijk = IDX(i, j, k, n0, n1, n2);

            /* __ set electron stress on g2 __ */
            s2[ijk].ps[0][0] = pw[ijk];
            s2[ijk].ps[0][1] = 0.0;
            s2[ijk].ps[0][2] = 0.0;
            s2[ijk].ps[0][3] = pw[ijk];
            s2[ijk].ps[0][4] = 0.0;
            s2[ijk].ps[0][5] = pw[ijk];
            }
        }
    }

/* __ clean-up the pointers __ */
free(sw);
free(pw);
free(uw);
free(vw);
free(ww);

/* __ clean-up the pointers __ */
free(sw4);
free(tw4);
free(uw4);
free(vw4);
free(ww4);

/* __ clean up the communicator __ */
MPI_Comm_free(&co);

}

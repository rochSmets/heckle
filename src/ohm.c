
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#include "structures.h"
#include "defines.h"
#include "fill.h"
#include "misc.h"
#include <ghosts.h>
#include <hecklebc.h>
#include <bc_constants.h>


/* _ calculate the e field __________________________________________________ */
void ohm(const STI * const si,
         const STX * const sx,
         ST1 *s1,
         ST2 *s2,
         HeckleBC   *hbc,
         Ghosts     *ghosts,
         int ipc)
{
    double dlx, dly, dlz;
    double dn;
    double dpx, dpy, dpz;
    double ljx, ljy, ljz;
    double bw[3];
    double ilx, ily, ilz;
    int i, j, k, l;
    int ijk;
    int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int ijka, ijkb, ijkc, ijkd, ijke, ijkf, ijkg, ijkh;
    int ijki, ijkj, ijkk, ijkl, ijkm, ijkn, ijko, ijkp;
    int ijkq, ijkr, ijks, ijkt, ijku, ijkv, ijkw, ijkx;
    int ijky, ijkz;
    int m0, m1, m2, n0, n1, n2;
    int nn2;
#if 1

    /* __ # of grid points __ */
    m0 = sx->n[0]+1;
    m1 = sx->n[1]+1;
    m2 = sx->n[2]+1;
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    /* __ # of grid points on g2 __ */
    nn2 = n0*n1*n2;

    /* __ needed for the laplacian __ */
    ilx = 1.0/(si->dl[0]*si->dl[0]);
    ily = 1.0/(si->dl[1]*si->dl[1]);
    ilz = 1.0/(si->dl[2]*si->dl[2]);

    /* __ set e @ n+1/2, n+3/2 : nested loops on subdomain __ */
    for (i = 1; i < m0; i++)
    {
        for (j = 1; j < m1; j++)
        {
            for (k = 1; k < m2; k++)
            {
                /* __ set index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ needed parameters __ */
                dlx = 0.25/(si->dl[0]*s2[ijk].ns[0]);
                dly = 0.25/(si->dl[1]*s2[ijk].ns[0]);
                dlz = 0.25/(si->dl[2]*s2[ijk].ns[0]);

                dn = 1.0/s2[ijk].ns[0];

                /* __ set index on g1 __ */
                ijk1 = IDX(i  , j  , k  , m0, m1, m2);
                ijk2 = IDX(i-1, j  , k  , m0, m1, m2);
                ijk3 = IDX(i  , j-1, k  , m0, m1, m2);
                ijk4 = IDX(i-1, j-1, k  , m0, m1, m2);
                ijk5 = IDX(i  , j  , k-1, m0, m1, m2);
                ijk6 = IDX(i-1, j  , k-1, m0, m1, m2);
                ijk7 = IDX(i  , j-1, k-1, m0, m1, m2);
                ijk8 = IDX(i-1, j-1, k-1, m0, m1, m2);

                switch (ipc)
                {
                case 0 : /* __ loop in each directions __ */
                    for (l = 0; l < 3; l++)
                    {
                        /* __ set the magnetic field value __ */
                        bw[l] = 0.125*(s1[ijk1].c[l]
                                       +s1[ijk2].c[l]
                                       +s1[ijk3].c[l]
                                       +s1[ijk4].c[l]
                                       +s1[ijk5].c[l]
                                       +s1[ijk6].c[l]
                                       +s1[ijk7].c[l]
                                       +s1[ijk8].c[l]);
                    }
                    break;

                case 1 : /* __ loop in each directions __ */
                    for (l = 0; l < 3; l++)
                    {
                        /* __ set the magnetic field value __ */
                        bw[l] = 0.125*(s1[ijk1].b[l]
                                       +s1[ijk2].b[l]
                                       +s1[ijk3].b[l]
                                       +s1[ijk4].b[l]
                                       +s1[ijk5].b[l]
                                       +s1[ijk6].b[l]
                                       +s1[ijk7].b[l]
                                       +s1[ijk8].b[l]);
                    }
                    break;

                    /* __ no reason to get there __ */
                default : IAMDEAD(sx->r);
                }

                /* __ set indexes on g2 __ */
                ijka = IDX(i+1, j+1, k-1, n0, n1, n2);
                ijkb = IDX(i+1, j  , k-1, n0, n1, n2);
                ijkc = IDX(i+1, j-1, k-1, n0, n1, n2);
                ijkd = IDX(i  , j+1, k-1, n0, n1, n2);
                ijke = IDX(i  , j-1, k-1, n0, n1, n2);
                ijkf = IDX(i-1, j+1, k-1, n0, n1, n2);
                ijkg = IDX(i-1, j  , k-1, n0, n1, n2);
                ijkh = IDX(i-1, j-1, k-1, n0, n1, n2);

                ijki = IDX(i+1, j+1, k  , n0, n1, n2);
                ijkj = IDX(i+1, j  , k  , n0, n1, n2);
                ijkk = IDX(i+1, j-1, k  , n0, n1, n2);
                ijkl = IDX(i  , j+1, k  , n0, n1, n2);
                ijkm = IDX(i  , j-1, k  , n0, n1, n2);
                ijkn = IDX(i-1, j+1, k  , n0, n1, n2);
                ijko = IDX(i-1, j  , k  , n0, n1, n2);
                ijkp = IDX(i-1, j-1, k  , n0, n1, n2);

                ijkq = IDX(i+1, j+1, k+1, n0, n1, n2);
                ijkr = IDX(i+1, j  , k+1, n0, n1, n2);
                ijks = IDX(i+1, j-1, k+1, n0, n1, n2);
                ijkt = IDX(i  , j+1, k+1, n0, n1, n2);
                ijku = IDX(i  , j-1, k+1, n0, n1, n2);
                ijkv = IDX(i-1, j+1, k+1, n0, n1, n2);
                ijkw = IDX(i-1, j  , k+1, n0, n1, n2);
                ijkx = IDX(i-1, j-1, k+1, n0, n1, n2);

                ijky = IDX(i  , j  , k+1, n0, n1, n2);
                ijkz = IDX(i  , j  , k-1, n0, n1, n2);

                /* __ grad of electron stress tensor : x comp. of div. ps[0] __ */
                dpx = ((s2[ijka].ps[0][0]-s2[ijkf].ps[0][0]
                        +s2[ijkc].ps[0][0]-s2[ijkh].ps[0][0]
                        +s2[ijkq].ps[0][0]-s2[ijkv].ps[0][0]
                        +s2[ijks].ps[0][0]-s2[ijkx].ps[0][0])*0.125
                        +(s2[ijkb].ps[0][0]-s2[ijkg].ps[0][0]
                        +s2[ijki].ps[0][0]-s2[ijkn].ps[0][0]
                        +s2[ijkk].ps[0][0]-s2[ijkp].ps[0][0]
                        +s2[ijkr].ps[0][0]-s2[ijkw].ps[0][0])*0.250
                        +(s2[ijkj].ps[0][0]-s2[ijko].ps[0][0])*0.500)*dlx

                        +((s2[ijka].ps[0][1]-s2[ijkc].ps[0][1]
                        +s2[ijkf].ps[0][1]-s2[ijkh].ps[0][1]
                        +s2[ijkq].ps[0][1]-s2[ijks].ps[0][1]
                        +s2[ijkv].ps[0][1]-s2[ijkx].ps[0][1])*0.125
                        +(s2[ijkd].ps[0][1]-s2[ijke].ps[0][1]
                        +s2[ijki].ps[0][1]-s2[ijkk].ps[0][1]
                        +s2[ijkn].ps[0][1]-s2[ijkp].ps[0][1]
                        +s2[ijkt].ps[0][1]-s2[ijku].ps[0][1])*0.250
                        +(s2[ijkl].ps[0][1]-s2[ijkm].ps[0][1])*0.500)*dly

                        +((s2[ijkq].ps[0][2]-s2[ijka].ps[0][2]
                        +s2[ijks].ps[0][2]-s2[ijkc].ps[0][2]
                        +s2[ijkv].ps[0][2]-s2[ijkf].ps[0][2]
                        +s2[ijkx].ps[0][2]-s2[ijkh].ps[0][2])*0.125
                        +(s2[ijkr].ps[0][2]-s2[ijkb].ps[0][2]
                        +s2[ijkt].ps[0][2]-s2[ijkd].ps[0][2]
                        +s2[ijku].ps[0][2]-s2[ijke].ps[0][2]
                        +s2[ijkw].ps[0][2]-s2[ijkg].ps[0][2])*0.250
                        +(s2[ijky].ps[0][2]-s2[ijkz].ps[0][2])*0.500)*dlz;

                /* __ y component of div. ps[0] __ */
                dpy = ((s2[ijka].ps[0][1]-s2[ijkf].ps[0][1]
                        +s2[ijkc].ps[0][1]-s2[ijkh].ps[0][1]
                        +s2[ijkq].ps[0][1]-s2[ijkv].ps[0][1]
                        +s2[ijks].ps[0][1]-s2[ijkx].ps[0][1])*0.125
                        +(s2[ijkb].ps[0][1]-s2[ijkg].ps[0][1]
                        +s2[ijki].ps[0][1]-s2[ijkn].ps[0][1]
                        +s2[ijkk].ps[0][1]-s2[ijkp].ps[0][1]
                        +s2[ijkr].ps[0][1]-s2[ijkw].ps[0][1])*0.250
                        +(s2[ijkj].ps[0][1]-s2[ijko].ps[0][1])*0.500)*dlx

                        +((s2[ijka].ps[0][3]-s2[ijkc].ps[0][3]
                        +s2[ijkf].ps[0][3]-s2[ijkh].ps[0][3]
                        +s2[ijkq].ps[0][3]-s2[ijks].ps[0][3]
                        +s2[ijkv].ps[0][3]-s2[ijkx].ps[0][3])*0.125
                        +(s2[ijkd].ps[0][3]-s2[ijke].ps[0][3]
                        +s2[ijki].ps[0][3]-s2[ijkk].ps[0][3]
                        +s2[ijkn].ps[0][3]-s2[ijkp].ps[0][3]
                        +s2[ijkt].ps[0][3]-s2[ijku].ps[0][3])*0.250
                        +(s2[ijkl].ps[0][3]-s2[ijkm].ps[0][3])*0.500)*dly

                        +((s2[ijkq].ps[0][4]-s2[ijka].ps[0][4]
                        +s2[ijks].ps[0][4]-s2[ijkc].ps[0][4]
                        +s2[ijkv].ps[0][4]-s2[ijkf].ps[0][4]
                        +s2[ijkx].ps[0][4]-s2[ijkh].ps[0][4])*0.125
                        +(s2[ijkr].ps[0][4]-s2[ijkb].ps[0][4]
                        +s2[ijkt].ps[0][4]-s2[ijkd].ps[0][4]
                        +s2[ijku].ps[0][4]-s2[ijke].ps[0][4]
                        +s2[ijkw].ps[0][4]-s2[ijkg].ps[0][4])*0.250
                        +(s2[ijky].ps[0][4]-s2[ijkz].ps[0][4])*0.500)*dlz;

                /* __ z component of div. ps[0] __ */
                dpz = ((s2[ijka].ps[0][2]-s2[ijkf].ps[0][2]
                        +s2[ijkc].ps[0][2]-s2[ijkh].ps[0][2]
                        +s2[ijkq].ps[0][2]-s2[ijkv].ps[0][2]
                        +s2[ijks].ps[0][2]-s2[ijkx].ps[0][2])*0.125
                        +(s2[ijkb].ps[0][2]-s2[ijkg].ps[0][2]
                        +s2[ijki].ps[0][2]-s2[ijkn].ps[0][2]
                        +s2[ijkk].ps[0][2]-s2[ijkp].ps[0][2]
                        +s2[ijkr].ps[0][2]-s2[ijkw].ps[0][2])*0.250
                        +(s2[ijkj].ps[0][2]-s2[ijko].ps[0][2])*0.500)*dlx

                        +((s2[ijka].ps[0][4]-s2[ijkc].ps[0][4]
                        +s2[ijkf].ps[0][4]-s2[ijkh].ps[0][4]
                        +s2[ijkq].ps[0][4]-s2[ijks].ps[0][4]
                        +s2[ijkv].ps[0][4]-s2[ijkx].ps[0][4])*0.125
                        +(s2[ijkd].ps[0][4]-s2[ijke].ps[0][4]
                        +s2[ijki].ps[0][4]-s2[ijkk].ps[0][4]
                        +s2[ijkn].ps[0][4]-s2[ijkp].ps[0][4]
                        +s2[ijkt].ps[0][4]-s2[ijku].ps[0][4])*0.250
                        +(s2[ijkl].ps[0][4]-s2[ijkm].ps[0][4])*0.500)*dly

                        +((s2[ijkq].ps[0][5]-s2[ijka].ps[0][5]
                        +s2[ijks].ps[0][5]-s2[ijkc].ps[0][5]
                        +s2[ijkv].ps[0][5]-s2[ijkf].ps[0][5]
                        +s2[ijkx].ps[0][5]-s2[ijkh].ps[0][5])*0.125
                        +(s2[ijkr].ps[0][5]-s2[ijkb].ps[0][5]
                        +s2[ijkt].ps[0][5]-s2[ijkd].ps[0][5]
                        +s2[ijku].ps[0][5]-s2[ijke].ps[0][5]
                        +s2[ijkw].ps[0][5]-s2[ijkg].ps[0][5])*0.250
                        +(s2[ijky].ps[0][5]-s2[ijkz].ps[0][5])*0.500)*dlz;

                /* __ set laplacian of x component of j __ */
                ljx = (s2[ijkj].j[0]-2.0*s2[ijk].j[0]+s2[ijko].j[0])*ilx
                        +(s2[ijkl].j[0]-2.0*s2[ijk].j[0]+s2[ijkm].j[0])*ily
                        +(s2[ijky].j[0]-2.0*s2[ijk].j[0]+s2[ijkz].j[0])*ilz;

                /* __ set laplacian of y component of j __ */
                ljy = (s2[ijkj].j[1]-2.0*s2[ijk].j[1]+s2[ijko].j[1])*ilx
                        +(s2[ijkl].j[1]-2.0*s2[ijk].j[1]+s2[ijkm].j[1])*ily
                        +(s2[ijky].j[1]-2.0*s2[ijk].j[1]+s2[ijkz].j[1])*ilz;

                /* __ set laplacian of z component of j __ */
                ljz = (s2[ijkj].j[2]-2.0*s2[ijk].j[2]+s2[ijko].j[2])*ilx
                        +(s2[ijkl].j[2]-2.0*s2[ijk].j[2]+s2[ijkm].j[2])*ily
                        +(s2[ijky].j[2]-2.0*s2[ijk].j[2]+s2[ijkz].j[2])*ilz;


                /* __ set total electric field __ */
                switch (ipc)
                {
                case 0 : /* __ predictor __ */
                    s2[ijk].f[0] = -(s2[ijk].vi[1]*bw[2]-s2[ijk].vi[2]*bw[1])
                            +(s2[ijk].j[1]*bw[2]-s2[ijk].j[2]*bw[1])*dn
                            -dpx
                            +s2[ijk].r*s2[ijk].j[0]
                            -si->hyvi*ljx;
                    s2[ijk].f[1] = -(s2[ijk].vi[2]*bw[0]-s2[ijk].vi[0]*bw[2])
                            +(s2[ijk].j[2]*bw[0]-s2[ijk].j[0]*bw[2])*dn
                            -dpy
                            +s2[ijk].r*s2[ijk].j[1]
                            -si->hyvi*ljy;
                    s2[ijk].f[2] = -(s2[ijk].vi[0]*bw[1]-s2[ijk].vi[1]*bw[0])
                            +(s2[ijk].j[0]*bw[1]-s2[ijk].j[1]*bw[0])*dn
                            -dpz
                            +s2[ijk].r*s2[ijk].j[2]
                            -si->hyvi*ljz;
                    break;

                case 1 : /* __ corrector __ */
                    s2[ijk].e[0] = -(s2[ijk].vi[1]*bw[2]-s2[ijk].vi[2]*bw[1])
                            +(s2[ijk].j[1]*bw[2]-s2[ijk].j[2]*bw[1])*dn
                            -dpx
                            +s2[ijk].r*s2[ijk].j[0]
                            -si->hyvi*ljx;
                    s2[ijk].e[1] = -(s2[ijk].vi[2]*bw[0]-s2[ijk].vi[0]*bw[2])
                            +(s2[ijk].j[2]*bw[0]-s2[ijk].j[0]*bw[2])*dn
                            -dpy
                            +s2[ijk].r*s2[ijk].j[1]
                            -si->hyvi*ljy;
                    s2[ijk].e[2] = -(s2[ijk].vi[0]*bw[1]-s2[ijk].vi[1]*bw[0])
                            +(s2[ijk].j[0]*bw[1]-s2[ijk].j[1]*bw[0])*dn
                            -dpz
                            +s2[ijk].r*s2[ijk].j[2]
                            -si->hyvi*ljz;
                    break;

                    /* __ no reason to get there __ */
                default : IAMDEAD(sx->r);
                }
            }
        }
    }

#endif

#ifdef OUTDATED
    /* __ memory allocation __ */
    tx = (double *)malloc(nn2*sizeof(double));
    ty = (double *)malloc(nn2*sizeof(double));
    tz = (double *)malloc(nn2*sizeof(double));

    /* __ load the buffers for "fill" : loop on the g2 grid __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* __ set the e field __ */
        switch (ipc)
        {
        case 0 : /* __ predictor __ */
            tx[ijk] = s2[ijk].f[0];
            ty[ijk] = s2[ijk].f[1];
            tz[ijk] = s2[ijk].f[2];
            break;

        case 1 : /* __ corrector __ */
            tx[ijk] = s2[ijk].e[0];
            ty[ijk] = s2[ijk].e[1];
            tz[ijk] = s2[ijk].e[2];
            break;

            /* __ no reason to get there __ */
        default : IAMDEAD(sx.r);
        }
    }

    /* __ fill buffered e field on g2 __ */
    fillg2(si, sx, tx, +1.0, +1.0, -1.0, -1.0, -1.0, -1.0);
    fillg2(si, sx, ty, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0);
    fillg2(si, sx, tz, -1.0, -1.0, -1.0, -1.0, +1.0, +1.0);

    /* __ update the e field : loop on the g2 grid __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* __ set the e field __ */
        switch (ipc)
        {
        case 0 : /* __ predictor __ */
            s2[ijk].f[0] = tx[ijk];
            s2[ijk].f[1] = ty[ijk];
            s2[ijk].f[2] = tz[ijk];
            break;

        case 1 : /* __ corrector __ */
            s2[ijk].e[0] = tx[ijk];
            s2[ijk].e[1] = ty[ijk];
            s2[ijk].e[2] = tz[ijk];
            break;

            /* __ no reason to get there __ */
        default : IAMDEAD(sx.r);
        }
    }
#endif

    /* --------------------------------------------------------------------- */
    /*                     COMMUNICATING GHOST NODES                         */
    /*                      AND BOUNDARY CONDITIONS                          */
    /* --------------------------------------------------------------------- */


    if (ipc == 0)
    {
        GhostsSendRecv(ghosts, sx, s2, GHOST_F);
        HeckleBCFieldApply(hbc,  sx, s2, BC_VAR_F);
    }
    else if (ipc == 1)
    {
        GhostsSendRecv(ghosts, sx, s2, GHOST_E);
        HeckleBCFieldApply(hbc,  sx, s2, BC_VAR_E);
    }





    /* --------------------------------------------------------------------- */
    /*                           EXTRAPOLATION                               */
    /* --------------------------------------------------------------------- */


    /* __ extrapolation of e field : loop on the g2 grid __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* __ loop in each direction __ */
        for (l = 0; l < 3; l++)
        {
            /* __ extrapolate or interpolate the e field __ */
            s2[ijk].e[l] = (ipc == 0) ? -s2[ijk].e[l]+2.0*s2[ijk].f[l]
                                      : 0.5*(s2[ijk].e[l]+s2[ijk].f[l]);
        }
    }

}


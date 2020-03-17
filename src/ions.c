
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"
#include "defines.h"
//#include "plasma.h"
#include "ghosts.h"
#include "close.h"
#include "misc.h"
#include "push.h"
#include "iamdead.h"
#include <ghosts.h>
#include <hecklebc.h>
#include <bc_constants.h>




/* __ set full ion kinetic pressure tensor __________________________________ */
void ions(const STI* const si,
          const STX* const sx,
          ST2 *s2,
          Particle *sp[NS+1],
          Ghosts *ghosts, HeckleBC *hbc,
          int it)
{
    double xw, yw, zw;
    double lx, ly, lz;
    double w1, w2, w3, w4, w5, w6, w7, w8;
    int ijk;
    int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int n0, n1, n2;
    int h, i, j, k, m, s;
    int nn2;
    double odl[3];

    // unused ?
    (void)it;


    /* __ # of grid points __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    /* __ # of grid points on g2 __ */
    nn2 = n0*n1*n2;


    odl[0] = 1./si->dl[0];
    odl[1] = 1./si->dl[1];
    odl[2] = 1./si->dl[2];


    /* __ set specie full pressure : loop on the part __ */
    for (s = 1; s < NS+1; s++)
    {
        /* __ initialization of moments : loop on the g2 indexes __ */
        for (ijk = 0; ijk < nn2; ijk++)
        {
            /* __ loop on the 6 directions __ */
            for (h = 0; h < 6; h++)
            {
                /* __ set specie pressure __ */
                s2[ijk].ps[s][h] = 0.0;
            }
        }
    }


    /* __ set specie charge density & specie velocities : loop on the part __ */
    for (s = 1; s < NS+1; s++)
    {
        //printf("s = %8d  : ws = %8.4f - ms = %8.4f\n", s, si.ws[s], si.ms[s]);
        /* __ loop on the part of specie s __ */
        for (m = 0; m < sx->ns[s]; m++)
        {
            /* __ part "position" in the subdomain @ even # of half ts __ */
            xw = sp[s][m].r[0]*odl[0] - sx->i0[0]+0.5;
            yw = sp[s][m].r[1]*odl[1] - sx->i0[1]+0.5;
            zw = sp[s][m].r[2]*odl[2] - sx->i0[2]+0.5;

            /* __ index for the part. "position" __ */
            i = (int)floor(xw);
            j = (int)floor(yw);
            k = (int)floor(zw);

#ifdef BUG
            if (i < 0 || i >= sx->n[0]+1)
            {
                deadpart(si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
            if (j < 0 || j >= sx->n[1]+1)
            {
                deadpart(si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
            if (k < 0 || k >= sx->n[2]+1)
            {
                deadpart(si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
#endif

            //      #ifdef BUG
            //      if (i < 0 || i >= sx->n[0]+1) shit(sx->r);
            //      if (j < 0 || j >= sx->n[1]+1) shit(sx->r);
            //      if (k < 0 || k >= sx->n[2]+1) shit(sx->r);
            //      #endif

            /* __ part location in the grid __ */
            lx = xw-i;
            ly = yw-j;
            lz = zw-k;

            /* __ weight for each vertices of the rounding grid points __ */
            w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
            w2 = (    lx)*(1.0-ly)*(1.0-lz);
            w3 = (1.0-lx)*(    ly)*(1.0-lz);
            w4 = (    lx)*(    ly)*(1.0-lz);
            w5 = (1.0-lx)*(1.0-ly)*(    lz);
            w6 = (    lx)*(1.0-ly)*(    lz);
            w7 = (1.0-lx)*(    ly)*(    lz);
            w8 = (    lx)*(    ly)*(    lz);
            //if (sx->r == 0) printf("sum w : %8.4f\n", w1+w2+w3+w4+w5+w6+w7+w8);

            /* __ indexes of the rounding grid points on g2 __ */
            ijk1 = IDX(i  , j  , k  , n0, n1, n2);
            ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
            ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
            ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
            ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
            ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
            ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
            ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);



            /* ------------------------------------------------------------- */
            /*                              PXX                              */
            /* ------------------------------------------------------------- */
            s2[ijk1].ps[s][0] += w1*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk1].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk1].vs[s][0]);

            #if 1

            s2[ijk2].ps[s][0] += w2*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk2].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk2].vs[s][0]);

            s2[ijk3].ps[s][0] += w3*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk3].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk3].vs[s][0]);

            s2[ijk4].ps[s][0] += w4*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk4].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk4].vs[s][0]);

            s2[ijk5].ps[s][0] += w5*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk5].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk5].vs[s][0]);

            s2[ijk6].ps[s][0] += w6*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk6].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk6].vs[s][0]);

            s2[ijk7].ps[s][0] += w7*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk7].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk7].vs[s][0]);

            s2[ijk8].ps[s][0] += w8*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk8].vs[s][0])
                    *(sp[s][m].v[0]-s2[ijk8].vs[s][0]);

            /* ------------------------------------------------------------- */
            /*                              PXY                              */
            /* ------------------------------------------------------------- */
            s2[ijk1].ps[s][1] += w1*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk1].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk1].vs[s][1]);

            s2[ijk2].ps[s][1] += w2*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk2].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk2].vs[s][1]);

            s2[ijk3].ps[s][1] += w3*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk3].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk3].vs[s][1]);

            s2[ijk4].ps[s][1] += w4*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk4].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk4].vs[s][1]);

            s2[ijk5].ps[s][1] += w5*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk5].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk5].vs[s][1]);

            s2[ijk6].ps[s][1] += w6*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk6].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk6].vs[s][1]);

            s2[ijk7].ps[s][1] += w7*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk7].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk7].vs[s][1]);

            s2[ijk8].ps[s][1] += w8*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk8].vs[s][0])
                    *(sp[s][m].v[1]-s2[ijk8].vs[s][1]);

            /* ------------------------------------------------------------- */
            /*                              PXZ                              */
            /* ------------------------------------------------------------- */
            s2[ijk1].ps[s][2] += w1*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk1].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk1].vs[s][2]);

            s2[ijk2].ps[s][2] += w2*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk2].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk2].vs[s][2]);

            s2[ijk3].ps[s][2] += w3*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk3].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk3].vs[s][2]);

            s2[ijk4].ps[s][2] += w4*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk4].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk4].vs[s][2]);

            s2[ijk5].ps[s][2] += w5*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk5].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk5].vs[s][2]);

            s2[ijk6].ps[s][2] += w6*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk6].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk6].vs[s][2]);

            s2[ijk7].ps[s][2] += w7*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk7].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk7].vs[s][2]);

            s2[ijk8].ps[s][2] += w8*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[0]-s2[ijk8].vs[s][0])
                    *(sp[s][m].v[2]-s2[ijk8].vs[s][2]);

            /* ------------------------------------------------------------- */
            /*                              PYY                              */
            /* ------------------------------------------------------------- */
            s2[ijk1].ps[s][3] += w1*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk1].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk1].vs[s][1]);

            s2[ijk2].ps[s][3] += w2*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk2].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk2].vs[s][1]);

            s2[ijk3].ps[s][3] += w3*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk3].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk3].vs[s][1]);

            s2[ijk4].ps[s][3] += w4*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk4].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk4].vs[s][1]);

            s2[ijk5].ps[s][3] += w5*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk5].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk5].vs[s][1]);

            s2[ijk6].ps[s][3] += w6*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk6].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk6].vs[s][1]);

            s2[ijk7].ps[s][3] += w7*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk7].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk7].vs[s][1]);

            s2[ijk8].ps[s][3] += w8*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk8].vs[s][1])
                    *(sp[s][m].v[1]-s2[ijk8].vs[s][1]);


            /* ------------------------------------------------------------- */
            /*                              PYZ                              */
            /* ------------------------------------------------------------- */
            s2[ijk1].ps[s][4] += w1*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk1].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk1].vs[s][2]);

            s2[ijk2].ps[s][4] += w2*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk2].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk2].vs[s][2]);

            s2[ijk3].ps[s][4] += w3*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk3].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk3].vs[s][2]);

            s2[ijk4].ps[s][4] += w4*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk4].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk4].vs[s][2]);

            s2[ijk5].ps[s][4] += w5*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk5].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk5].vs[s][2]);

            s2[ijk6].ps[s][4] += w6*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk6].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk6].vs[s][2]);

            s2[ijk7].ps[s][4] += w7*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk7].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk7].vs[s][2]);

            s2[ijk8].ps[s][4] += w8*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[1]-s2[ijk8].vs[s][1])
                    *(sp[s][m].v[2]-s2[ijk8].vs[s][2]);

            /* ------------------------------------------------------------- */
            /*                              PZZ                              */
            /* ------------------------------------------------------------- */
            s2[ijk1].ps[s][5] += w1*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk1].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk1].vs[s][2]);

            s2[ijk2].ps[s][5] += w2*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk2].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk2].vs[s][2]);

            s2[ijk3].ps[s][5] += w3*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk3].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk3].vs[s][2]);

            s2[ijk4].ps[s][5] += w4*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk4].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk4].vs[s][2]);

            s2[ijk5].ps[s][5] += w5*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk5].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk5].vs[s][2]);

            s2[ijk6].ps[s][5] += w6*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk6].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk6].vs[s][2]);

            s2[ijk7].ps[s][5] += w7*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk7].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk7].vs[s][2]);

            s2[ijk8].ps[s][5] += w8*si->ws[s]*si->ms[s]
                    *(sp[s][m].v[2]-s2[ijk8].vs[s][2])
                    *(sp[s][m].v[2]-s2[ijk8].vs[s][2]);

            #endif
        }
    }


    GhostsSendRecv(ghosts, sx, s2, GHOST_P);
    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_P);

}


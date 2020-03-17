
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "structures.h"
#include "iamdead.h"
//#include "plasma.h"
#include "misc.h"
#include "push.h"
#include "wall.h"
#include <particle.h>
#include <particlebc.h>
#include <hecklebc.h>



/* __ push the part. ________________________________________________________ */
void push(STI *si,
          STX *sx,
          Grid0 *s0,
          ST1 *s1,
          ST2 *s2,
          Particle *sp[NS+1],
          PartBC *pbc,
          int ipc,
          int s,
          int it)
{
    double fw, gw;
    double w1, w2, w3, w4, w5, w6, w7, w8;
    double xw, yw, zw;
    double lx, ly, lz;
    double bw[3], ew[3];
    double v[3], u[3];
    double ts;
    int i, j, k, l, m;
    int ijk0, ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int l0, l1, l2, m0, m1, m2, n0, n1, n2;


    /* __ # of grid points __ */
    l0 = sx->n[0];
    l1 = sx->n[1];
    l2 = sx->n[2];
    m0 = sx->n[0]+1;
    m1 = sx->n[1]+1;
    m2 = sx->n[2]+1;
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    /* __ set null time step to initial setting __ */
    ts = (ipc == 0 && it == -1) ? 0.0 : si->ts;


    /* __ needed constant __ */
    fw = 0.5*ts*si->qs[s]/si->ms[s];


    /* __ loop on all the part. __ */
    for (m = 0; m < sx->ns[s]; m++)
    {

        /* __ predictor __ */
        if (ipc == 0)
        {
            /* __ keep memory of older part position __ */
            sp[s][m].s[0] = sp[s][m].r[0];
            sp[s][m].s[1] = sp[s][m].r[1];
            sp[s][m].s[2] = sp[s][m].r[2];
        }

        /* __ interpolate b field __ */
        xw = (sp[s][m].r[0]/si->dl[0]-sx->i0[0] < sx->n[0]) ?
                    sp[s][m].r[0]/si->dl[0]-sx->i0[0] : sx->n[0]-EPS4;

        yw = (sp[s][m].r[1]/si->dl[1]-sx->i0[1] < sx->n[1]) ?
                    sp[s][m].r[1]/si->dl[1]-sx->i0[1] : sx->n[1]-EPS4;

        zw = (sp[s][m].r[2]/si->dl[2]-sx->i0[2] < sx->n[2]) ?
                    sp[s][m].r[2]/si->dl[2]-sx->i0[2] : sx->n[2]-EPS4;

#if BUG
        if (       sp[s][m].r[0] < 0 || sp[s][m].r[0] >= si->l[0]
                   || sp[s][m].r[1] < 0 || sp[s][m].r[1] >= si->l[1]
                   || sp[s][m].r[2] < 0 || sp[s][m].r[2] >= si->l[2])
        {
            printf("MERDE ipc = %d %d %d %f %f %f %f %f %f %f %f %f\n",
                   ipc, m,s,
                   sp[s][m].r[0],sp[s][m].r[1],sp[s][m].r[2],
                    sp[s][m].v[0],sp[s][m].v[1],sp[s][m].v[2],
                    si->l[0], si->l[1], si->l[2]);
        }
#endif
        /* __ index for the part. "position" __ */
        i = (int)xw;
        j = (int)yw;
        k = (int)zw;

#ifdef BUG
        if (i < 0 || i >= sx->n[0])
        {
            deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
        if (j < 0 || j >= sx->n[1])
        {
            deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
        if (k < 0 || k >= sx->n[2])
        {
            deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
#endif

        /* __ part. location in the cell __ */
        lx = xw-i;
        ly = yw-j;
        lz = zw-k;

        /* __ indexes of the rounding grid points on g1 __ */
        ijk0 = IDX(i  , j  , k  , l0, l1, l2);

        /* __ indexes of the rounding grid points on g1 __ */
        ijk1 = IDX(i  , j  , k  , m0, m1, m2);
        ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
        ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
        ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
        ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
        ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
        ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
        ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);

        /* __ weight for each vertices of the rounding grid points __ */
        w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
        w2 = (    lx)*(1.0-ly)*(1.0-lz);
        w3 = (1.0-lx)*(    ly)*(1.0-lz);
        w4 = (    lx)*(    ly)*(1.0-lz);
        w5 = (1.0-lx)*(1.0-ly)*(    lz);
        w6 = (    lx)*(1.0-ly)*(    lz);
        w7 = (1.0-lx)*(    ly)*(    lz);
        w8 = (    lx)*(    ly)*(    lz);

        /* __ loop in the 3 directions __ */
        for (l = 0; l < 3; l++)
        {
            /* __ set b field seen by the part. __ */
            bw[l] = w1*s1[ijk1].b[l]
                    +w2*s1[ijk2].b[l]
                    +w3*s1[ijk3].b[l]
                    +w4*s1[ijk4].b[l]
                    +w5*s1[ijk5].b[l]
                    +w6*s1[ijk6].b[l]
                    +w7*s1[ijk7].b[l]
                    +w8*s1[ijk8].b[l];
        }

        /* __ set # of part per cells __ */
        //s1[ijk1].ppc ++;
        s0[ijk0].ppc ++;

        /* __ interpolate e field __ */
        xw = sp[s][m].r[0]/si->dl[0]-sx->i0[0]+0.5;
        yw = sp[s][m].r[1]/si->dl[1]-sx->i0[1]+0.5;
        zw = sp[s][m].r[2]/si->dl[2]-sx->i0[2]+0.5;

        /* __ index for the part. "position" __ */
        i = (int)xw;
        j = (int)yw;
        k = (int)zw;

#ifdef BUG
        if (i < 0 || i >= sx->n[0]+1)
        {
            deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
        if (j < 0 || j >= sx->n[1]+1)
        {
            deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
        if (k < 0 || k >= sx->n[2]+1)
        {
            deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
#endif

        //  #ifdef BUG
        //  if (i < 0 || i >= sx->n[0]+1) shit(sx->r);
        //  if (j < 0 || j >= sx->n[1]+1) shit(sx->r);
        //  if (k < 0 || k >= sx->n[2]+1) shit(sx->r);
        //  #endif

        /* __ part. location in the cell __ */
        lx = xw-i;
        ly = yw-j;
        lz = zw-k;

        /* __ indexes of the rounding grid points on g2 __ */
        ijk1 = IDX(i  , j  , k  , n0, n1, n2);
        ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
        ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
        ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
        ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
        ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
        ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
        ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

        /* __ weight for each vertices of the rounding grid points __ */
        w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
        w2 = (    lx)*(1.0-ly)*(1.0-lz);
        w3 = (1.0-lx)*(    ly)*(1.0-lz);
        w4 = (    lx)*(    ly)*(1.0-lz);
        w5 = (1.0-lx)*(1.0-ly)*(    lz);
        w6 = (    lx)*(1.0-ly)*(    lz);
        w7 = (1.0-lx)*(    ly)*(    lz);
        w8 = (    lx)*(    ly)*(    lz);

        /* __ loop in the 3 directions __ */
        for (l = 0; l < 3; l++)
        {
            /* __ set e field seen by the part. __ */
            ew[l] = w1*(s2[ijk1].e[l]-s2[ijk1].r*s2[ijk1].j[l])
                    +w2*(s2[ijk2].e[l]-s2[ijk2].r*s2[ijk2].j[l])
                    +w3*(s2[ijk3].e[l]-s2[ijk3].r*s2[ijk3].j[l])
                    +w4*(s2[ijk4].e[l]-s2[ijk4].r*s2[ijk4].j[l])
                    +w5*(s2[ijk5].e[l]-s2[ijk5].r*s2[ijk5].j[l])
                    +w6*(s2[ijk6].e[l]-s2[ijk6].r*s2[ijk6].j[l])
                    +w7*(s2[ijk7].e[l]-s2[ijk7].r*s2[ijk7].j[l])
                    +w8*(s2[ijk8].e[l]-s2[ijk8].r*s2[ijk8].j[l]);
        }

        /* __ needed parameter __ */
        gw = 2.0/(1.0+fw*fw*(bw[0]*bw[0]+bw[1]*bw[1]+bw[2]*bw[2]));

        /* __ new velocities after half acceleration of e field __ */
        v[0] = sp[s][m].v[0]+ew[0]*fw;
        v[1] = sp[s][m].v[1]+ew[1]*fw;
        v[2] = sp[s][m].v[2]+ew[2]*fw;

        /* __ needed for the b field rotation __ */
        u[0] = v[0]+(v[1]*bw[2]-v[2]*bw[1])*fw;
        u[1] = v[1]+(v[2]*bw[0]-v[0]*bw[2])*fw;
        u[2] = v[2]+(v[0]*bw[1]-v[1]*bw[0])*fw;


        /* __ predictor __ */
        if (ipc == 0)
        {
            /* __ set the new velocity __ */
            sp[s][m].v[0] = v[0]+((u[1]*bw[2]-u[2]*bw[1])*gw+ew[0])*fw;
            sp[s][m].v[1] = v[1]+((u[2]*bw[0]-u[0]*bw[2])*gw+ew[1])*fw;
            sp[s][m].v[2] = v[2]+((u[0]*bw[1]-u[1]*bw[0])*gw+ew[2])*fw;

#ifdef BUG
            if (sp[s][m].v[0] > (si->dl[0]/ts))
            {
                deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
            if (sp[s][m].v[1] > (si->dl[1]/ts))
            {
                deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
            if (sp[s][m].v[2] > (si->dl[2]/ts))
            {
                deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
#endif

            //     #ifdef BUG
            //     if (sp[s][m].v[0] > (si->dl[0]/ts)) shit(sx->r);
            //     if (sp[s][m].v[1] > (si->dl[1]/ts)) shit(sx->r);
            //     if (sp[s][m].v[2] > (si->dl[2]/ts)) shit(sx->r);
            //     #endif

            /* __ push the particles __ */
            sp[s][m].r[0] = sp[s][m].r[0]+sp[s][m].v[0]*ts;
            sp[s][m].r[1] = sp[s][m].r[1]+sp[s][m].v[1]*ts;
            sp[s][m].r[2] = sp[s][m].r[2]+sp[s][m].v[2]*ts;


            // now detect wether the particle is out of the domain
            // and if it is store its index and species
            // so that the boundary condition module deals with it later
            if (PartBCisOut(pbc, si, sx, sp[s][m].r) == ISOUT)
            {
                PartBCStore(pbc, m, s);
            }

        }// end predictor


        /* __ corrector __ */
        if (ipc == 1)
        {
            /* __ set the new velocity __ */
            sp[s][m].w[0] = v[0]+((u[1]*bw[2]-u[2]*bw[1])*gw+ew[0])*fw;
            sp[s][m].w[1] = v[1]+((u[2]*bw[0]-u[0]*bw[2])*gw+ew[1])*fw;
            sp[s][m].w[2] = v[2]+((u[0]*bw[1]-u[1]*bw[0])*gw+ew[2])*fw;

#ifdef BUG
            if (sp[s][m].w[0] > (si->dl[0]/ts))
            {
                deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
            if (sp[s][m].w[1] > (si->dl[1]/ts))
            {
                deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
            if (sp[s][m].w[2] > (si->dl[2]/ts))
            {
                deadpart(*si, sx, s0, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
#endif


            //     #ifdef BUG
            //     if (sp[s][m].w[0] > (si->dl[0]/ts)) shit(sx->r);
            //     if (sp[s][m].w[1] > (si->dl[1]/ts)) shit(sx->r);
            //     if (sp[s][m].w[2] > (si->dl[2]/ts)) shit(sx->r);
            //     #endif


            /* __ push the particles __ */
            sp[s][m].s[0] = sp[s][m].r[0]+sp[s][m].w[0]*ts;
            sp[s][m].s[1] = sp[s][m].r[1]+sp[s][m].w[1]*ts;
            sp[s][m].s[2] = sp[s][m].r[2]+sp[s][m].w[2]*ts;



            // now detect wether the particle is out of the domain
            // and if it is store its index and species
            // so that the boundary condition module deals with it later
            if (PartBCisOut(pbc, si, sx, sp[s][m].s) == ISOUT)
            {
                PartBCStore(pbc, m, s);
            }

        } // end corrector

    } // end loop over particles
}


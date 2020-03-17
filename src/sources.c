

// standard C libraries
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// MPI
#include <mpi.h>


// Heckle
#include "structures.h"
#include "defines.h"
#include "close.h"
#include "misc.h"
#include "push.h"
#include "iamdead.h"
#include "smooth.h"
#include <particle.h>           // declaration of the particle type
#include <ghosts.h>             // ghost points, MPI communications
#include <hecklebc.h>           // headers for boundary conditions
#include <bc_constants.h>       // constants for boundary conditions
#include <particlebc.h>         // particle boundary conditions





/*---------------------------------------------------------------------------
  sources()
  ---------------------------------------------------------------------------
  AIM : push particles and accumulate moments on the grid
 ---------------------------------------------------------------------------*/
void sources(STI *si,
             STX *sx,
             Grid0 *s0,
             ST1 *s1,
             ST2 *s2,
             Particle *sp[NS+1],
             HeckleBC *hbc,
             Ghosts *ghosts,
             int ipc,
             int it)
{
    double xw, yw, zw;
    double lx, ly, lz;
    double w1, w2, w3, w4, w5, w6, w7, w8;
    //double *un, *ux, *uy, *uz;
    int ijk;
    int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int n0, n1, n2;
    int i, j, k, l, m, s;
    int nn0, nn2;
    PartBC *pbc;


    // get a particle boundary condition handle
    pbc = HeckleBCGetPartBC(hbc);


    /* __ # of grid pooints __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

    /* __ # of grid points on g0, g1 and g2 __ */
    nn0 = (sx->n[0]  )*(sx->n[1]  )*(sx->n[2]  );
    nn2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);



    /* --------------------------------------------------------------------- */
    /*                   INITIALIZE MOMENT BUFFERS                           */
    /* --------------------------------------------------------------------- */



    /* __ initialization of # of part per cells : loop on the g1 indexes __ */
    for (ijk = 0; ijk < nn0; ijk++)
    {
        /* __ initialization of # of part per cells __ */
        //s1[ijk].ppc = 0;
        s0[ijk].ppc = 0;
    }

    /* __ initialization of moments : loop on the g2 indexes __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* __ loop on species __ */
        for (s = 1; s < NS+1; s++)
        {
            /* __ set specie charge density (even # half ts) __ */
            s2[ijk].os[s] = EPS6;

            /* __ loop on the 3 directions __ */
            for (l = 0; l < 3; l++)
            {
                /* __ set specie velocity __ */
                s2[ijk].vs[s][l] = 0.0;
            }
        }
    }




    /* --------------------------------------------------------------------- */
    /*               PUSH THE PARTICLES + BOUNDARY CONDITION                 */
    /* --------------------------------------------------------------------- */


    /* __ push the part. : loop on the part __ */
    for (s = 1; s < NS+1; s++)
    {
        /* __ push the part of specie s __ */
        push(si, sx, s0, s1, s2, sp, pbc, ipc, s, it);
    }


    // now all particles have been moved we want to apply the boundary condition
    // on those detected as leaving the domain.

    PartBCApply(pbc, si, sx, sp, ipc);
    PartBCReset(pbc);

#if 0
    for (int ispe=1; ispe < NS+1; ispe++)
    {
        //printf("proc %d has %d particles\n", sx->r, sx->ns[ispe]);
        if (sx->ns[ispe] >= si->nm)
        {
            printf("ERROR %d macroparticles on proc %d exceeds max %d\n",
                   sx->ns[ispe], sx->r, si->nm);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
#endif

    /* --------------------------------------------------------------------- */
    /*                      ACCUMULATE MOMENTS ON THE GRID                   */
    /* --------------------------------------------------------------------- */




    /* __ set specie charge density & specie velocities : loop on the part __ */
    for (s = 1; s < NS+1; s++)
    {
        /* __ loop on the part of specie s __ */
        for (m = 0; m < sx->ns[s]; m++)
        {
            /* __ part "position" in the subdomain @ even # of half ts__ */
            xw = 0.5*(sp[s][m].r[0]+sp[s][m].s[0])/si->dl[0]-sx->i0[0]+0.5;
            yw = 0.5*(sp[s][m].r[1]+sp[s][m].s[1])/si->dl[1]-sx->i0[1]+0.5;
            zw = 0.5*(sp[s][m].r[2]+sp[s][m].s[2])/si->dl[2]-sx->i0[2]+0.5;

            /* __ index for the part "position" __ */
            i = (int)xw;
            j = (int)yw;
            k = (int)zw;

            /* __ part location on the grid __ */
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

            switch (ipc)
            {
            case 0 : /* __ predictor : part "position" in the subdomain __ */
                xw = sp[s][m].r[0]/si->dl[0]-sx->i0[0]+0.5;
                yw = sp[s][m].r[1]/si->dl[1]-sx->i0[1]+0.5;
                zw = sp[s][m].r[2]/si->dl[2]-sx->i0[2]+0.5;
                break;

            case 1 : /* __ corrector : part "position" in the subdomain __ */
                xw = sp[s][m].s[0]/si->dl[0]-sx->i0[0]+0.5;
                yw = sp[s][m].s[1]/si->dl[1]-sx->i0[1]+0.5;
                zw = sp[s][m].s[2]/si->dl[2]-sx->i0[2]+0.5;
                break;

                /* __ no reason to get there __ */
            default : IAMDEAD(sx->r);
            }

            /* __ index for the part "position" __ */
            i = (int)xw;
            j = (int)yw;
            k = (int)zw;

#if 1
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

            //      #ifdef BUG
            //      if (i < 0 || i >= sx->n[0]+1) shit(sx->r);
            //      if (j < 0 || j >= sx->n[1]+1) shit(sx->r);
            //      if (k < 0 || k >= sx->n[2]+1) shit(sx->r);
            //      #endif

            /* __ indexes of the rounding grid points on g2 __ */
            ijk1 = IDX(i  , j  , k  , n0, n1, n2);
            ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
            ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
            ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
            ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
            ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
            ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
            ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

            /* __ gather part for species densities __ */
            s2[ijk1].os[s] += w1;
            s2[ijk2].os[s] += w2;
            s2[ijk3].os[s] += w3;
            s2[ijk4].os[s] += w4;
            s2[ijk5].os[s] += w5;
            s2[ijk6].os[s] += w6;
            s2[ijk7].os[s] += w7;
            s2[ijk8].os[s] += w8;

            /* __ loop on the 3 directions __ */
            for (l = 0; l < 3; l++)
            {
                switch (ipc)
                {
                case 0 : /* __ predictor : gather part for species velocities __ */
                    s2[ijk1].vs[s][l] += w1*sp[s][m].v[l];
                    s2[ijk2].vs[s][l] += w2*sp[s][m].v[l];
                    s2[ijk3].vs[s][l] += w3*sp[s][m].v[l];
                    s2[ijk4].vs[s][l] += w4*sp[s][m].v[l];
                    s2[ijk5].vs[s][l] += w5*sp[s][m].v[l];
                    s2[ijk6].vs[s][l] += w6*sp[s][m].v[l];
                    s2[ijk7].vs[s][l] += w7*sp[s][m].v[l];
                    s2[ijk8].vs[s][l] += w8*sp[s][m].v[l];
                    break;

                case 1 : /* __ corrector : gather part for species velocities __ */
                    s2[ijk1].vs[s][l] += w1*sp[s][m].w[l];
                    s2[ijk2].vs[s][l] += w2*sp[s][m].w[l];
                    s2[ijk3].vs[s][l] += w3*sp[s][m].w[l];
                    s2[ijk4].vs[s][l] += w4*sp[s][m].w[l];
                    s2[ijk5].vs[s][l] += w5*sp[s][m].w[l];
                    s2[ijk6].vs[s][l] += w6*sp[s][m].w[l];
                    s2[ijk7].vs[s][l] += w7*sp[s][m].w[l];
                    s2[ijk8].vs[s][l] += w8*sp[s][m].w[l];
                    break;

                    /* __ no reason to get there __ */
                default : IAMDEAD(sx->r);
                }
            }
        }
    }


    /* --------------------------------------------------------------------- */
    /*          COMMUNICATE GHOST NODES AND APPLY BOUNDARY CONDITIONS        */
    /* --------------------------------------------------------------------- */

    GhostsSendRecv(ghosts, sx, s2, GHOST_NS_VS);    // communicate ghost points
    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_OSVS);   // TODO implement this


    /* at this point the density and fluxes are accumulated, ghost nodes have
       been communicated, and boundary condition taken care of for OS and VS */



    /* --------------------------------------------------------------------- */
    /*                    FINISH CALCULATING MOMENTS                         */
    /*                                                                       */
    /*  we still need to :                                                   */
    /*                                                                       */
    /*  - normalization of the density                                       */
    /*  - sum over all ions species                                          */
    /*  - total ion bulk velocity                                            */
    /* --------------------------------------------------------------------- */



    /* __ normalization of specie moments : loop on the part __ */
    for (s = 1; s < NS+1; s++)
    {
        /* __ loop on the g2 indexes __ */
        for (ijk = 0; ijk < nn2; ijk++)
        {
            /* __ loop in the 3 directions __ */
            for (l = 0; l < 3; l++)
            {
                /* __ normalization of specie velocities __ */
                s2[ijk].vs[s][l] = s2[ijk].vs[s][l]/s2[ijk].os[s];
            }

            /* __ normalization of specie density __ */
            s2[ijk].os[s] = s2[ijk].os[s]*si->ws[s];
        }
    }

    /* __ set specie density @ odd # of half si.ts : loop on the g2 indexes __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* __ loop on the part __ */
        for (s = 1; s < NS+1; s++)
        {
            /* __ set specie charge density at t=0 __ */
            if (ipc == 0 && it == -1) s2[ijk].ms[s] = s2[ijk].os[s];

            /* __ set specie charge density (odd # half ts) __ */
            s2[ijk].ns[s] = 0.5*(s2[ijk].ms[s]+s2[ijk].os[s]);

            /* __ keep specie charge density (even # half ts) for next step __ */
            if (ipc == 0) s2[ijk].ms[s] = s2[ijk].os[s];
        }
    }

    /* __ keep memory of s2.ns[0] for the calculation of electron stress __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        s2[ijk].ms[0] = s2[ijk].ns[0];
    }

    /* __ initialization of moments : loop on the g2 indexes __ */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* __ set initial charge density __ */
        s2[ijk].ns[0] = 0.0;

        /* __ loop in each ditection __ */
        for (l = 0; l < 3; l++)
        {
            /* __ set initial fluid velocity __ */
            s2[ijk].vi[l] = 0.0;
        }
    }

    /* __ set electron density : loop on the specie __ */
    for (s = 1; s < NS+1; s++)
    {
        /* __ loop on the g2 indexes __ */
        for (ijk = 0; ijk < nn2; ijk++)
        {
            /* __ integration of charge density __ */
            s2[ijk].ns[0] += s2[ijk].ns[s]*si->qs[s];
        }
    }

    /* __ set total ion fluid velocity : loop on the specie __ */
    for (s = 1; s < NS+1; s++)
    {
        /* __ loop on the g2 indexes __ */
        for (ijk = 0; ijk < nn2; ijk++)
        {
            /* __ loop in each direction __ */
            for (l = 0; l < 3; l++)
            {
                /* __ integration of fluid velocity __ */
                s2[ijk].vi[l] += s2[ijk].ns[s]*s2[ijk].vs[s][l]*si->qs[s]
                        /s2[ijk].ns[0];
            }
        }
    }


    /* --------------------------------------------------------------------- */
    /*                          SMOOTHING MOMENTS                            */
    /* --------------------------------------------------------------------- */

    smooth(sx, s2, ghosts, hbc);

}


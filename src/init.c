
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

#include <structures.h>
#include "defines.h"
#include "maxwell.h"
#include "ohm.h"
#include "misc.h"
#include "closeModel.h"
#include "sources.h"
#include "initModel.h"
#include "particle.h"
#include "ghosts.h"





/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           DEBUG FUNCTIONS                             //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */







/*---------------------------------------------------------------------------
  check_particles()
  ---------------------------------------------------------------------------
  AIM : Check particle arrays. Print a message if not ok.
 ---------------------------------------------------------------------------*/
int check_particles(const STI * const si, const STX * const sx, Particle *sp[NS+1])
{
    int ip;
    int ispe;
    int ok;
    ParticleDBG *dbg;

    dbg = pdbg_init();

    ok = 1;

    for (ispe=1; ispe < NS+1; ispe++)
    {
        for (ip=0; ip < sx->ns[ispe]; ip++)
        {
            if (ParticleDBG_IsInside(si, &sp[ispe][ip]) == 0)
            {
                fprintf(stderr, "ERROR (proc %d)- Particle %d of species %d is out"
                       " : r(%f,%f,%f)   s(%f,%f,%f)\n",
                       sx->r, ip, ispe,
                       sp[ispe][ip].r[0],
                       sp[ispe][ip].r[1],
                       sp[ispe][ip].r[2],
                       sp[ispe][ip].s[0],
                       sp[ispe][ip].s[1],
                       sp[ispe][ip].s[2]);

                ok =0;
            }// end if particle is inside



            if (ParticleDBG_NanInf(&sp[ispe][ip], dbg) == 0)
            {
                fprintf(stderr, "ERROR (proc %d) - Particle %d of species %d has NaNs/Infs"
                       " : r(%f,%f,%f)   s(%f,%f,%f), v(%f,%f,%f),   w(%f,%f,%f)\n",
                       sx->r, ip, ispe,
                       sp[ispe][ip].r[0],
                       sp[ispe][ip].r[1],
                       sp[ispe][ip].r[2],
                       sp[ispe][ip].s[0],
                       sp[ispe][ip].s[1],
                       sp[ispe][ip].s[2],
                       sp[ispe][ip].v[0],
                       sp[ispe][ip].v[1],
                       sp[ispe][ip].v[2],
                       sp[ispe][ip].w[0],
                       sp[ispe][ip].w[1],
                       sp[ispe][ip].w[2]);

                ok =0;
            }
        }
    }

    pdbg_delete(dbg);
    return ok;
}
/*===========================================================================*/












/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */













/*---------------------------------------------------------------------------
  maxwellian();
  ---------------------------------------------------------------------------
  AIM : draw a random velocity from a normal distribution with 0 average and
  standard deviation stddev() in each of the 3 directions.
 ---------------------------------------------------------------------------*/
void maxwellian(double stddev[3], double v[3])
{
    double r1, r2;

    r1 = RNM;
    r2 = RNM;

    r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
    r2   = (fabs(r2 - 1.0) < EPS8) ? r2 - EPS8 : r2;

    r1   = (r1 > EPS8)? r1 : r1 + EPS8;
    r2   = (r2 > EPS8)? r2 : r2 + EPS8;

    v[0] = sqrt(-2*log(r1))*stddev[0] * cos(2*PI*r2);
    v[1] = sqrt(-2*log(r1))*stddev[1] * sin(2*PI*r2);

    r1 = RNM;
    r2 = RNM;

    r1   = (fabs(r1 - 1.0) < EPS8) ? r1 - EPS8 : r1;
    r2   = (fabs(r2 - 1.0) < EPS8) ? r2 - EPS8 : r2;

    r1   = (r1 > EPS8)? r1 : r1 + EPS8;
    r2   = (r2 > EPS8)? r2 : r2 + EPS8;

    v[2] = sqrt(-2*log(r1))*stddev[2] * cos(2*PI*r2);
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  initMagnetic()
  ---------------------------------------------------------------------------
  AIM : initialize the magnetic field profile on g1 according to the
        model the user has chosen.
 ---------------------------------------------------------------------------*/
void initMagnetic(struct sti *si, struct stx *sx, struct st1 *s1)
{
    int npg1[3];    /* number of points on g1 in each direction */
    double pos[3];
    int i,j,k, ijk;
    double B[3];


    npg1[0] = sx->n[0]+1;
    npg1[1] = sx->n[1]+1;
    npg1[2] = sx->n[2]+1;

    for (i = 0; i < npg1[0]; i++)
    {
        for (j = 0; j < npg1[1]; j++)
        {
            for (k = 0; k < npg1[2]; k++)
            {
                /* __ coordinate of grid point __ */
                pos[0] = (i + sx->i0[0]) * (si->dl[0]);
                pos[1] = (j + sx->i0[1]) * (si->dl[1]);
                pos[2] = (k + sx->i0[2]) * (si->dl[2]);


                /* __ index of the grid point on g1 __ */
                ijk = IDX(i, j, k, npg1[0], npg1[1], npg1[2]);

                /* __ set the magnetic field value __ */
                initModelMagnetic(si, sx, pos, B);
                s1[ijk].b[0] = B[0];
                s1[ijk].b[1] = B[1];
                s1[ijk].b[2] = B[2];
                #ifdef FULLP
                s1[ijk].bKept[0] = B[0];
                s1[ijk].bKept[1] = B[1];
                s1[ijk].bKept[2] = B[2];
                #endif
            } // end k loop
        }//end j loop
    } // end i loop
}
/*===========================================================================*/











/*---------------------------------------------------------------------------
  initCurrent()
  ---------------------------------------------------------------------------
  AIM : initialize the current field profile on g2 according to the
        model the user has chosen.
 ---------------------------------------------------------------------------*/
void initCurrent(struct sti *si, struct stx *sx, struct st2 *s2)
{
    int npg2[3];    /* number of points on g2 in each direction */
    double pos[3];
    int i,j,k;
    int ijk;
    double J[3];

    npg2[0] = sx->n[0]+2;
    npg2[1] = sx->n[1]+2;
    npg2[2] = sx->n[2]+2;

    /* __ set e, v, j & r : nested loops on the grid points of subdomain __ */
    for (i = 0; i < npg2[0]; i++)
    {
        for (j = 0; j < npg2[1]; j++)
        {
            for (k = 0; k < npg2[2]; k++)
            {
                /* __ index of the grid point on g2 __ */
                ijk = IDX(i, j, k, npg2[0], npg2[1], npg2[2]);

                /* __ coordinate of grid point __ */
                pos[0] = (i-0.5 + sx->i0[0]) * (si->dl[0]);
                pos[1] = (j-0.5 + sx->i0[1]) * (si->dl[1]);
                pos[2] = (k-0.5 + sx->i0[2]) * (si->dl[2]);

                /* __ loop on the 3 directions __ */
                initModelCurrent(si, sx, pos, J);
                s2[ijk].j[0] = J[0];
                s2[ijk].j[1] = J[1];
                s2[ijk].j[2] = J[2];
            }// end k loop
        }// end j loop
    } // end i loop
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  initElectric()
  ---------------------------------------------------------------------------
  AIM : initialize the electric field profile on g1 according to the
        model the user has chosen.
 ---------------------------------------------------------------------------*/
void initElectric(struct sti *si, struct stx *sx, struct st2 *s2)
{
    int npg2[3];    /* number of points on g2 in each direction */
    double pos[3];
    int i,j,k;
    int ijk;
    double E[3];

    npg2[0] = sx->n[0]+2;
    npg2[1] = sx->n[1]+2;
    npg2[2] = sx->n[2]+2;

    /* __ set e, v, j & r : nested loops on the grid points of subdomain __ */
    for (i = 0; i < npg2[0]; i++)
    {
        for (j = 0; j < npg2[1]; j++)
        {
            for (k = 0; k < npg2[2]; k++)
            {
                /* __ index of the grid point on g2 __ */
                ijk = IDX(i, j, k, npg2[0], npg2[1], npg2[2]);

                /* __ coordinate of grid point __ */
                pos[0] = (i-0.5 + sx->i0[0]) * (si->dl[0]);
                pos[1] = (j-0.5 + sx->i0[1]) * (si->dl[1]);
                pos[2] = (k-0.5 + sx->i0[2]) * (si->dl[2]);

                /* __ loop on the 3 directions __ */
                initModelElectric(si, sx, pos, E);
                s2[ijk].e[0] = E[0];
                s2[ijk].e[1] = E[1];
                s2[ijk].e[2] = E[2];

                //s2[ijk].f[0] = E[0];
                //s2[ijk].f[1] = E[1];
                //s2[ijk].f[2] = E[2];

            }// end k loop
        }// end j loop
    } // end i loop
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  initResistivity()
  ---------------------------------------------------------------------------
  AIM : Initialize the resistivity. Also artificially increases it close
       to non-periodic walls.
 ---------------------------------------------------------------------------*/
void initResistivity(struct sti *si, struct stx *sx, struct st2 *s2)
{
    int npg2[3];    /* number of points on g2 in each direction */
    //double pos[3];
    int i,j,k;
    int ijk;

    npg2[0] = sx->n[0]+2;
    npg2[1] = sx->n[1]+2;
    npg2[2] = sx->n[2]+2;


    /* __ set e, v, j & r : nested loops on the grid points of subdomain __ */
    for (i = 0; i < npg2[0]; i++)
    {
        for (j = 0; j < npg2[1]; j++)
        {
            for (k = 0; k < npg2[2]; k++)
            {
                /* __ index of the grid point on g2 __ */
                ijk = IDX(i, j, k, npg2[0], npg2[1], npg2[2]);

                /* __ coordinate of grid point __ */
                //pos[0] = (i-0.5 + sx->i0[0]) * (si->dl[0]);
                //pos[1] = (j-0.5 + sx->i0[1]) * (si->dl[1]);
                //pos[2] = (k-0.5 + sx->i0[2]) * (si->dl[2]);

                s2[ijk].r = si->rsty;

                /* __ near left wall __ */
                if (sx->nt[ 4] == MPI_PROC_NULL)
                {
                    if (i == 0         ) s2[ijk].r = 125.0*si->rsty;
                    if (i == 1         ) s2[ijk].r =  25.0*si->rsty;
                    if (i == 2         ) s2[ijk].r =   5.0*si->rsty;
                }

                /* __ near right wall __ */
                if (sx->nt[22] == MPI_PROC_NULL)
                {
                    if (i == sx->n[0]+1) s2[ijk].r = 125.0*si->rsty;
                    if (i == sx->n[0]  ) s2[ijk].r =  25.0*si->rsty;
                    if (i == sx->n[0]-1) s2[ijk].r =   5.0*si->rsty;
                }

                /* __ near bottom wall __ */
                if (sx->nt[10] == MPI_PROC_NULL)
                {
                    if (j == 0         ) s2[ijk].r = 125.0*si->rsty;
                    if (j == 1         ) s2[ijk].r =  25.0*si->rsty;
                    if (j == 2         ) s2[ijk].r =   5.0*si->rsty;
                }

                /* __ near top wall __ */
                if (sx->nt[16] == MPI_PROC_NULL)
                {
                    if (j == sx->n[1]+1) s2[ijk].r = 125.0*si->rsty;
                    if (j == sx->n[1]  ) s2[ijk].r =  25.0*si->rsty;
                    if (j == sx->n[1]-1) s2[ijk].r =   5.0*si->rsty;
                }

               /* __ near back wall __ */
                if (sx->nt[12] == MPI_PROC_NULL)
                {
                    if (k == 0         ) s2[ijk].r = 125.0*si->rsty;
                    if (k == 1         ) s2[ijk].r =  25.0*si->rsty;
                    if (k == 2         ) s2[ijk].r =   5.0*si->rsty;
                }

                /* __ near front wall __ */
                if (sx->nt[14] == MPI_PROC_NULL)
                {
                    if (k == sx->n[2]+1) s2[ijk].r = 125.0*si->rsty;
                    if (k == sx->n[2]  ) s2[ijk].r =  25.0*si->rsty;
                    if (k == sx->n[2]-1) s2[ijk].r =   5.0*si->rsty;
                }

            }// end k loop
        }// end j loop
    } // end i loop
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  initElectronTemperature()
  ---------------------------------------------------------------------------
  AIM : Initialize the electron temperature
        This routine asks the model the parallel and perpendicular temperature
        of the electrons, assuming their distribution is gyrotropic.
        Then it initializes the components of the electron pressure tensor
        in the coordinate system used in the simulation (x,y,z).
        To save memory, only 6 out of the 9 components of Pe are kept.
 ---------------------------------------------------------------------------*/
void initElectronTemperature(struct sti *si, struct stx *sx, struct st2 *s2)
{

    int npg2[3];    /* number of points on g2 in each direction */
    double pos[3];
    int i, j, k, ijk;
    double vw[3][3];
    double wte[3][3];
    double ww;
    int e,f,g,h;
    double B[3], T[3];

    npg2[0] = sx->n[0]+2;
    npg2[1] = sx->n[1]+2;
    npg2[2] = sx->n[2]+2;

    for (i = 0; i < npg2[0]; i++)
    {
        for (j = 0; j < npg2[1]; j++)
        {
            for (k = 0; k < npg2[2]; k++)
            {
                /* __ coordinates of grid point __ */
                pos[0] = (i - 0.5 + sx->i0[0]) * (si->dl[0]);
                pos[1] = (j - 0.5 + sx->i0[1]) * (si->dl[1]);
                pos[2] = (k - 0.5 + sx->i0[2]) * (si->dl[2]);

                /* we need the magnetic field at this location
                   to change the pressure tensor into the XYZ basis
                   from the para/perp one. */
                initModelMagnetic(si, sx, pos, B);

                /* parallel & perp. temperature are given
                 * by the model */
                initModelTemperature(si, sx, pos, 0, T);

                /* __ index of the grid point on g2 __ */
                ijk = IDX(i, j, k, npg2[0], npg2[1], npg2[2]);

                /* this returns a change of basis matrix 'vw' */
                ortho(B, vw);

                /* __ nested loops in the 3 directions __ */
                for (e = 0; e < 3; e++)
                {
                    for (f = 0; f < 3; f++)
                    {
                        /* initialize to 0 before incrementation */
                        wte[e][f] = 0.0;

                        /* __ nested loops in the 3 directions __ */
                        for (g = 0; g < 3; g++)
                        {
                            for (h = 0; h < 3; h++)
                            {
                                /* 'ww' is the component Pij in the
                                   field aligned basis. Since we assume
                                   gyrotropy, it is 0 everywhere but on
                                   the diagonal */
                                ww = 0.0;

                                /* take diagonal terms */
                                if (g == 0 && h == 0) ww = T[0];
                                if (g == 1 && h == 1) ww = T[1];
                                if (g == 2 && h == 2) ww = T[1];


                                /* __ cartesian projection of temperature __ */
                                wte[e][f] += vw[g][e]*ww*vw[h][f];
                            }
                        }
                    }
                }


                /* now copy the tensor but keep only 6 components
                   using the symmetry. */
                s2[ijk].te[0] = wte[0][0];
                s2[ijk].te[1] = wte[0][1];
                s2[ijk].te[2] = wte[0][2];
                s2[ijk].te[3] = wte[1][1];
                s2[ijk].te[4] = wte[1][2];
                s2[ijk].te[5] = wte[2][2];

                }
            }
        }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  initPlasma()
  ---------------------------------------------------------------------------
  AIM : This function fills the particle arrays. It generates locally
  Maxwellian distribution functions based on modeled plasma moments.
  It puts particles randomely positionned within cells proportionally to
  the density into that cell.
 ---------------------------------------------------------------------------*/
void initPlasma(struct sti *si, struct stx *sx, struct stp *sp[NS+1])
{
    int ispe;
    int i,j,k;
    int ijk;
    int ncells[3];
    double pos[3];
    int npc;
    int comp;
    double *n, ntot, ntotmpi;
    int ip;
    double vw[3][3];
    double B[3],  curdrift[3], vdrift[3], vavg[3], T[2], vth[3];
    double vpb[3], vps[3];
    struct stp *pcur;
    int64_t *np[NS+1];
    int ts[NS+1];
    int sw[NS+1];
    int rank;

    ncells[0] = sx->n[0];
    ncells[1] = sx->n[1];
    ncells[2] = sx->n[2];

    /* will store the density within each cell */
    n = malloc(ncells[0] * ncells[1] * ncells[2] * sizeof *n);

    sx->ns[0] = 0;
    /* loop over the ION species (start at 1) */
    for (ispe=1; ispe < NS+1; ispe++)
    {
        sw[ispe] = 0;
        sx->ns[ispe] = 0;
        sx->nm[ispe] = 0.;
        /* if there are no particles of this species
           we don't do anything */
        if (si->ns[ispe] == 0)
        {
            continue;
        }

        /* first we need to know the total density in the box */
        ntot = 0;
        for (i=0; i < ncells[0]; i ++)
        {
            for (j=0; j < ncells[1]; j++)
            {
                for (k=0; k < ncells[2];  k++)
                {
                    /* __ coordinates of grid point __ */
                    pos[0] = (i + 0.5 + sx->i0[0]) * (si->dl[0]);
                    pos[1] = (j + 0.5 + sx->i0[1]) * (si->dl[1]);
                    pos[2] = (k + 0.5 + sx->i0[2]) * (si->dl[2]);

                    /* __ index of the grid point on g0 __ */
                    ijk = IDX(i, j, k, ncells[0], ncells[1], ncells[2]);


                    /* get the density value from the model */
                    n[ijk] = initModelDensity(si, sx, pos, ispe);

                    /* __ integrated density value on the subdomain __ */
                    ntot += n[ijk];
                }
            }
        } /* END loop on cells */

        /* __ reduce total density in whole box __ */
        MPI_Allreduce(&ntot, &ntotmpi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        /* now loop over the cells of the MPI domain to place particles */
        for (i = 0; i < ncells[0]; i++)
        {
            for (j=0; j < ncells[1]; j++)
            {
                for (k = 0; k < ncells[2]; k++)
                {
                    /* __ index of the grid point on g0 __ */
                    ijk = IDX(i, j, k, ncells[0], ncells[1], ncells[2]);

                    /* number of particles in this cell */
                    npc = (int) (si->ns[ispe]/ntotmpi * n[ijk]);

                    /* now put npc particles into that cell */
                    for (ip=0; ip < npc; ip++)
                    {
                        /* __ coordinates of grid point __ */
                        pos[0] = (i + RNM) * (si->dl[0]);
                        pos[1] = (j + RNM) * (si->dl[1]);
                        pos[2] = (k + RNM) * (si->dl[2]);

                        /* shortcut for the current particle */
                        pcur = &(sp[ispe][sw[ispe]]);

                        /* set the particle position
                          we dont want the particle to be numerically
                          exactly on the right border of the domain

                          to segmentation faults.
                          so when/if that happens du to bad luck in
                          random numbers, we substract EPS6 */
                       pcur->r[0] = (pos[0] == sx->n[0]*si->dl[0])?
                                     sx->i1[0] * si->dl[0] - EPS6 :
                                     sx->i0[0] * si->dl[0] + pos[0];

                        pcur->r[1] = (pos[1] == sx->n[1]*si->dl[1])?
                                      sx->i1[1] * si->dl[1] - EPS6  :
                                      sx->i0[1] * si->dl[1] + pos[1];

                        pcur->r[2] = (pos[2] == sx->n[2] * si->dl[2]) ?
                                      sx->i1[2] * si->dl[2] - EPS6 :
                                      sx->i0[2] * si->dl[2] + pos[2];

                        pcur->s[0] = pcur->r[0];
                        pcur->s[1] = pcur->r[1];
                        pcur->s[2] = pcur->r[2];


                        /* put 'pos' in the global domain coord. system */
                        pos[0] += sx->i0[0] * si->dl[0];
                        pos[1] += sx->i0[1] * si->dl[1];
                        pos[2] += sx->i0[2] * si->dl[2];

                        /* __ set # of crossing of the boundary domain __ */
                        pcur->b[0] = 0;
                        pcur->b[1] = 0;
                        pcur->b[2] = 0;

                        /* now the particle is positionned we need to set its
                        velocity. For this we need to know at its position :
                            - the ExB drift velocity
                            - the drift associated with the current
                            - the local temperature
                            - the magnetic field (for para/perp defined data)
                        */
                        initModelMagnetic(si, sx, pos, B);
                        initModelTemperature(si, sx, pos, ispe, T);
                        initModelDrift(si, sx, pos, ispe, vdrift);
                        initModelCurDrift(si, sx, pos, ispe, curdrift);

                        vavg[0] = vdrift[0] + curdrift[0];
                        vavg[1] = vdrift[1] + curdrift[1];
                        vavg[2] = vdrift[2] + curdrift[2];

                        vth[0]  = sqrt(T[0]/si->ms[ispe]);
                        vth[1]  = sqrt(T[1]/si->ms[ispe]);
                        vth[2]  = sqrt(T[1]/si->ms[ispe]);

                        /* get a velocity from a maxwellian with
                           the local thermal velocity of this species */
                        maxwellian(vth, vpb);

                        /* 'vp' is in field align basis, we need to put it
                            in the simulation basis. First get a change of
                            basis matrix then find the velocity in the
                            simulation frame */

                        ortho(B, vw);
                        for (comp=0; comp < 3; comp++)
                        {
                            vps[comp] = vw[0][comp] * vpb[0]
                                      + vw[1][comp] * vpb[1]
                                      + vw[2][comp] * vpb[2];
                        }

                        /* assign the velocity to the particle */
                        for (comp = 0; comp < 3; comp++)
                        {
                            /* __ set velocity __ */
                            pcur->v[comp] = vps[comp]+ vavg[comp];

                            /* __ set proj. velocity __ */
                            pcur->w[comp] = pcur->v[comp];
                        }


                        /* __ counter for the part. index __ */
                        sw[ispe]++;
                        if (sw[ispe] >= si->nm) {
                           fprintf(stderr, "you are loading to many particles ! try to increase the si.nm value\n");
                           exit(-1);
                        }

                    } /* End loop particles */
                } /* End loop k */
            } /* End loop j*/
        } /* End loop i */


        /* store the actual number of macroparticle loaded for that species */
        sx->ns[ispe] = sw[ispe];

        /* store the total number of those accross the whole box */
        MPI_Allreduce( &sx->ns[ispe],
                       &si->ns[ispe],
                       1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        /* and also get the max */
        MPI_Allreduce(&sx->ns[ispe],
                      &sx->nm[ispe],
                      1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);


        /* now we give an unique index to all particles. The ID starts
           at 0 on processor 0. For processor r we need to know how many
           particles */

        np[ispe] = malloc(sx->s * sizeof *np[ispe]);

        /* __ set the total # of part. in nodes < r __ */
        MPI_Allgather(&sx->ns[ispe], 1, MPI_LONG_LONG_INT, np[ispe],
                       1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

        /* __ set counter for the total # of part. in nodes < r __ */
        ts[ispe] = 0;

        /* total # of part. in nodes < r */
        for (rank = 0; rank < sx->r; rank++)
        {
            ts[ispe] += np[ispe][rank];
        }

        /* clean-up the pointers */
        free(np[ispe]);

        /* now initialize the index */
        for (ip = 0; ip < sx->ns[ispe]; ip++)
        {
            sp[ispe][ip].i = ts[ispe] + ip;
        }
    }// end species loop


    //sp[1][10].r[0] = log(-10);

    // check that all particles ar OK

    if (! check_particles(si, sx, sp))
    {
        exit(-1);
    }


    free(n);
}
/*===========================================================================*/





/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PUBLIC  FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */



/* __ initialize the part. and the fields ___________________________________ */
void init(STI *si,
          STX *sx,
          Grid0 *s0,
          ST1 *s1,
          ST2 *s2,
          Particle *sp[NS+1],
          HeckleBC *hbc,
          Ghosts *ghosts,
          int *it)
{
    int ijk;
    int npg2[3];    /* number of points on g2 in each direction */
    int nn2;
    int ipc;
    int comp;


    initMagnetic(si, sx, s1);

    initCurrent(si, sx, s2);

    initElectric(si, sx, s2);

    initResistivity(si, sx, s2);

    initElectronTemperature(si, sx, s2);

    initPlasma(si, sx, sp);


    /* set the weight factor */
    normal(si, sx);


    npg2[0] = sx->n[0]+2;
    npg2[1] = sx->n[1]+2;
    npg2[2] = sx->n[2]+2;
    nn2     = npg2[0] * npg2[1] * npg2[2];

    /* --------------------------------------------------------------------- */
    /*                              PREDICTOR                                */
    /* --------------------------------------------------------------------- */
    ipc = 0;
    /* set the iterator to -1 as a flag to NOT push the particles and just
       accumulate moments*/
    *it = -1;

    sources(si, sx, s0, s1, s2, sp, hbc, ghosts, ipc, *it);

    MaxwellFaraday(si, sx, s1, s2, ipc);
    MaxwellAmpere( si, sx, s1, s2, ghosts, hbc, ipc);

    /* initialize the electron bulk velocity and pressure tensor to their
       initial values. Some closure may need it. */
    for (ijk = 0; ijk < nn2; ijk++)
    {
        s2[ijk].ms[0] = s2[ijk].ns[0];
    }

    for (ijk = 0; ijk < nn2; ijk++)
    {
        /* electron pressure tensor */
        for (comp = 0; comp < 6; comp++)
        {
            s2[ijk].ps[0][comp] = s2[ijk].ns[0] * s2[ijk].te[comp];
            #ifdef FULLP
            s2[ijk].pKept[comp] = s2[ijk].ns[0] * s2[ijk].te[comp];
            s2[ijk].dFull[comp] = s2[ijk].ns[0] * s2[ijk].te[comp];
            #endif
        }

        /* electron bulk velocity */
        for (comp = 0; comp < 3; comp++)
        {
            s2[ijk].vs[0][comp] = s2[ijk].vi[comp]
                                - s2[ijk].j[comp] / s2[ijk].ns[0];
        }
    }

    /* __ set electron pressure __ */
    closeModelPressure(*it, si, sx, s1, s2, hbc, ghosts, ipc);

    /* __ set e field with ohm's law __ */
    ohm(si, sx, s1, s2, hbc, ghosts, ipc);

    /* __ tricky : needed to distinguish the first faraday half rotation __ */
    ipc = 2;

    /* __ half rotation with faraday __ */
    MaxwellFaraday(si, sx, s1, s2, ipc);
    MaxwellAmpere(si, sx, s1, s2, ghosts, hbc, ipc);

    /* --------------------------------------------------------------------- */
    /*                              CORRECTOR                                */
    /* --------------------------------------------------------------------- */

    /* __ set corrector __ */
    ipc = 1;
    *it = 0;

    /* now move the plasma, update the magnetic field (faraday)
       the electron closure (stress) and the electric field (ohm) */
    sources(si, sx, s0, s1, s2, sp, hbc, ghosts, ipc, *it);
    MaxwellFaraday(si, sx, s1, s2, ipc);
    MaxwellAmpere(si, sx, s1, s2, ghosts, hbc, ipc);
    closeModelPressure(*it, si, sx, s1, s2, hbc, ghosts, ipc);
    ohm(si, sx, s1, s2, hbc, ghosts, ipc);

    ipc = 2;

    /* __ half rotation with faraday __ */
    MaxwellFaraday(si, sx, s1, s2, ipc);
    MaxwellAmpere( si, sx, s1, s2, ghosts, hbc, ipc);


}






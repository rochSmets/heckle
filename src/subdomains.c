
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "subdomains.h"
#include <mpi.h>
#include <bc_constants.h>


/* __ domain decomposition & initialization of sx tructure ________________ */
void subdomains(struct sti si, struct stx *sx)
{
    double ws[3];
    double wl;
    int wd[3];
    int wp[3];
    int cs[3];
    int ss[3];
    int wn, wr;
    int a, b, c, i, j, k, l;
    MPI_Comm com;
    int nprocs;
    int autodecomp=0;

    /* __ # of subdomains in each direction : arbitrary large double value __ */
    wl = 1.0/EPS6;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* we take the user decomposition unless the number of cores specified
   is not equal to the number of MPI domains required, in which case
   we perform an automatic decomposition */

    if (si.mpidom[0]*si.mpidom[1]*si.mpidom[2] == nprocs)
    {
        sx->d[0] = si.mpidom[0];
        sx->d[1] = si.mpidom[1];
        sx->d[2] = si.mpidom[2];
    }

    else
    {
        autodecomp = 1;
        /* __ nested loops on the total # of nodes __ */
        for (i = 1; i <= sx->s; i++)
        {
            for (j = 1; j <= sx->s; j++)
            {
                for (k = 1; k <= sx->s; k++)
                {
                    /* __ scale ratio for i,j,k subdomains in X,Y,Z directions __ */
                    ws[0] = fabs(si.n[0]/(double)i-si.n[1]/(double)j);
                    ws[1] = fabs(si.n[1]/(double)j-si.n[2]/(double)k);
                    ws[2] = fabs(si.n[2]/(double)k-si.n[0]/(double)i);

                    /* __ test the scale ratio if # of nodes is correct __ */
                    if (i*j*k == sx->s && maxdbl(ws, 3) < wl)
                    {
                        /* __ new scale ratio limit __ */
                        wl = maxdbl(ws, 3);

                        /* __ # of subdomains in each directions __ */
                        sx->d[0] = i;
                        sx->d[1] = j;
                        sx->d[2] = k;
                    }
                }
            }
        }
    } // end of automatic decomposition




#ifdef BUG
    if (sx->d[0] > si.n[0]) IAMDEAD(sx->r);
    if (sx->d[1] > si.n[1]) IAMDEAD(sx->r);
    if (sx->d[2] > si.n[2]) IAMDEAD(sx->r);
#endif

    /* __ # of subdomains for the topology __ */
    for (l = 0; l < 3; l++)
    {
        wd[l] = sx->d[l];
    }

    /* __ bc for the topology __ */
    for (l = 0; l < 3; l++)
    {
        wp[l] = (si.bc[l] == BC_TYPE_PERIODIC) ? 1 : 0;
    }

    /* __ create a cartesian topology __ */
    MPI_Cart_create(MPI_COMM_WORLD, 3, wd, wp, 1, &com);

    /* __ get the coordinates of the subdomain __ */
    MPI_Cart_coords(com, sx->r, 3, cs);

    /* __ min & max indexes for the subdomains in each directions : loop in the 3 directions __ */
    for (l = 0; l < 3; l++)
    {
        /* __ approximate # of cells & remainder __ */
        wn = (int)(si.n[l]/((double)sx->d[l]));
        wr = si.n[l]-wn*sx->d[l];

        /* __ set indexes if remainder is null __ */
        if (wr == 0)
        {
            /* __ set min & max indexes (g1) in each directions __ */
            sx->i0[l] = cs[l]*wn;
            sx->i1[l] = cs[l]*wn+wn;
        }

        /* __ if remainder is nonzero __ */
        if (wr != 0)
        {
            /* __ larger # of cells for small indexes of subdomains ( < wr) __ */
            if (cs[l] < wr)
            {
                /* __ set min & max indexes (g1) in each directions __ */
                sx->i0[l] = cs[l]*(wn+1);
                sx->i1[l] = cs[l]*(wn+1)+wn+1;
            }
            /* __ smaller # of cells for large indexes of subdomains ( >= wr) __ */
            if (cs[l] >= wr)
            {
                /* __ set min & max indexes (g1) in each directions __ */
                sx->i0[l] = wr*(wn+1)+(cs[l]-wr)*wn;
                sx->i1[l] = wr*(wn+1)+(cs[l]-wr)*wn+wn;
            }
        }
    }

    /* __ # of cells & size of subdomain __ */
    for (l = 0; l < 3; l++)
    {
        sx->n[l] = sx->i1[l]-sx->i0[l];
        sx->l[l] = sx->n[l]*si.dl[l];
    }

    /* __ find the neighbors : nested loops on the neighbors __ */
    for (a = -1; a <= +1; a++)
    {
        for (b = -1; b <= +1; b++)
        {
            for (c = -1; c <= +1; c++)
            {
                /* __ guess the coordinate of the neighbor __ */
                ss[0] = cs[0]+a;
                ss[1] = cs[1]+b;
                ss[2] = cs[2]+c;

                /* __ if nonperiodic & coordinate out of range : no neighbor __ */
                if ((si.bc[0] !=  BC_TYPE_PERIODIC && (ss[0] < 0 || ss[0] >= wd[0])) ||
                        (si.bc[1] !=  BC_TYPE_PERIODIC && (ss[1] < 0 || ss[1] >= wd[1])) ||
                        (si.bc[2] !=  BC_TYPE_PERIODIC && (ss[2] < 0 || ss[2] >= wd[2])))
                {
                    wr = MPI_PROC_NULL;
                }
                /* __ get the rank of the neighbor @ the specified coordinate __ */
                else
                {
                    MPI_Cart_rank(com, ss, &wr);
                }

                /* __ rank to receive (from) & to send (to) __ */
                sx->nf[(1-c)+3*((1-b)+3*(1-a))] = wr;
                sx->nt[(1+c)+3*((1+b)+3*(1+a))] = wr;
            }
        }
    }

    /* __ set irun index __ */
    sx->irun = 0;

    /* __ set initial seed for the random numbers __ */
    srand(1+sx->r);


    if (sx->r == 0)
    {
        printf("________________ domain decomposition ________________\n");
        printf("decomposition    : %s\n", (autodecomp)? "automatic":"user");
        printf("\n");
        printf("domains (x)      : %d\n", sx->d[0]);
        printf("domains (y)      : %d\n", sx->d[1]);
        printf("domains (z)      : %d\n", sx->d[2]);
    }


}


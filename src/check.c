
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>

#include "structures.h"
#include "check.h"
#include "defines.h"
#include "misc.h"




/* __ check time evolution of div b, total energy & electric fluctuations ___ */
/*---------------------------------------------------------------------------
  check()
  ---------------------------------------------------------------------------
  AIM : print some information about the simulation on screen
 ---------------------------------------------------------------------------*/
void check(struct sti si,
           struct stx sx,
           struct st0 *s0,
           struct st1 *s1,
           struct st2 *s2,
           struct stp *sp[NS+1],
           struct std *sd,
           int it,
           clock_t ticks[2],
           time_t ttime[2])
{
    const int nenergy = 1+3*(NS+1);
    double wdb[2], *dbw[2];
    double wet[1+3*(NS+1)], etw[1+3*(NS+1)];
    double wle[3], lew[3];
    double wfa[3], *faw[3];
    double wfi[3], *fiw[3];
    double dv;
    double bte;
    int wpc, *pcw;
    int nsplus1;
    int h, l, s;


    for (s = 0; s < NS+1; s++)
    {
        /* __ total # of part. __ */
        MPI_Allreduce(&(sx.ns[s]), &(si.ns[s]), 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

        /* __ min # of part. in domain __ */
        MPI_Allreduce(&(sx.ns[s]), &(sx.no[s]), 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);

        /* __ max # of part. in domain __ */
        MPI_Allreduce(&(sx.ns[s]), &(sx.nm[s]), 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
//      printf("s = %1d   r = %1d ... sx.ns[%1d] = %zd   sx.nm[%1d] = %zd\n", s, sx.r, s, sx.ns[s], s, sx.nm[s]);
    }


    /* __ volume element __ */
    dv = si.l[0]*si.l[1]*si.l[2];

    /* __ get local maximum of div b __ */
    db(si, sx, s0, s1, wdb, &wpc);

    /* __ div. b & pseudo-div. b __ */
    for (h = 0; h < 2; h++)
    {
        /* __ memory allocation __ */
        dbw[h] = (double *)malloc(sx.s*sizeof(double));


    }

    pcw = (int *)malloc(sx.s*sizeof(int));


    /* __ gather all local maximum of div b __ */
    for (h = 0; h < 2; h++)
    {
        MPI_Allgather(&wdb[h], 1, MPI_DOUBLE, dbw[h],
                      1, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    /* __ gather all local minimum values of parts per cells __ */
    MPI_Allgather(&wpc, 1, MPI_INT, pcw, 1, MPI_INT, MPI_COMM_WORLD);

    /* __ get local values of energies __ */
    et(si, sx, s1, s2, sp, wet, wle);

    /* __ reduce total energy __ */
    MPI_Allreduce(wet, etw, nenergy, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    /* __ reduce CFL condition __ */
    MPI_Allreduce(wle, lew, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    /* __ electric fluctuations __ */
    ef(sx, s2, wfa, wfi);

    /* __ loop in each direction : buffer for min & max electric field __ */
    for (l = 0; l < 3; l++)
    {
        /* __ memory allocation __ */
        faw[l] = (double *)malloc(sx.s*sizeof(double));
        fiw[l] = (double *)malloc(sx.s*sizeof(double));
    }

    /* __ loop in each direction __ */
    for (l = 0; l < 3; l++)
    {
        /* __ gather all local max & min of electric field components __ */
        MPI_Allgather(&wfa[l], 1, MPI_DOUBLE, faw[l], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(&wfi[l], 1, MPI_DOUBLE, fiw[l], 1, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    /* __ maximum of div b __ */
    sd->db = maxdbl(dbw[0], sx.s);
    sd->dw = maxdbl(dbw[1], sx.s);

    /* __ minimum # of parts per cells __ */
    sd->pc = minint(pcw, sx.s);

    /* __ min & maxi # of parts in subdomain __ */
    nsplus1 = NS+1;
    sd->no = minint64(sx.nm, nsplus1);
    sd->nm = maxint64(sx.nm, nsplus1);

    /* __ total energy __ */
    sd->ma = etw[0]/dv;
    bte = sd->ma;
    for (s = 0; s < NS+1; s++) {
        sd->pb[s] = etw[1+3*s]/dv;
        sd->ta[s] = etw[2+3*s]/dv;
        sd->te[s] = etw[3+3*s]/dv;
        bte += sd->pb[s]+sd->ta[s]+sd->te[s];
    }

    /* __ electric fluctuations in each directions __ */
    sd->fx = maxdbl(faw[0], sx.s)-mindbl(fiw[0], sx.s);
    sd->fy = maxdbl(faw[1], sx.s)-mindbl(fiw[1], sx.s);
    sd->fz = maxdbl(faw[2], sx.s)-mindbl(fiw[2], sx.s);

    /* __ set initial value of total energy __ */
    if (it == 0) {
        sd->e0 = sd->ma;
        for (s = 0; s < NS+1; s++) {
            sd->e0 += sd->pb[s]+sd->ta[s]+sd->te[s];
        }
    }

    /* __ only on node 0 __ */
    if (sx.r == 0)
    {
        /* __ print the checked parameters __ */
        printf("\n");
        printf("divergence b     :%12.6lf\n", sd->db);
        printf("pseudo div b     :%12.6lf\n", sd->dw);
        printf("cfl (x)          :%12.6lf\n", lew[0]*si.ts/si.dl[0]);
        printf("cfl (y)          :%12.6lf\n", lew[1]*si.ts/si.dl[2]);
        printf("cfl (z)          :%12.6lf\n", lew[2]*si.ts/si.dl[2]);
        printf("min # part/cell  :%12d\n", sd->pc);
        printf("min # part/sub   :%12zu\n", sd->no);
        printf("max # part/sub   :%12zu\n", sd->nm);
        printf("magnetic energy  :%12.6lf\n", sd->ma);

        for (s = 0; s < NS+1; s++) {
            printf("bulk energy [%1d]  :%12.6lf\n", s, sd->pb[s]);
            printf("para energy [%1d]  :%12.6lf\n", s, sd->ta[s]);
            printf("perp energy [%1d]  :%12.6lf\n", s, sd->te[s]);
        }
        printf("total energy     :%12.6lf\n", bte);
        printf("delta e / e0 (%%) :%12.6lf\n", 100.0*fabs(bte-sd->e0)/sd->e0);
        printf("e. (x) fluct.    :%12.6lf\n", sd->fx);
        printf("e. (y) fluct.    :%12.6lf\n", sd->fy);
        printf("e. (z) fluct.    :%12.6lf\n", sd->fz);
        for (s = 1; s < NS+1; s++) {
//          printf("# of part [%1d]    :%12" PRId64 "\n", s, si.ns[s]);
            printf("# of part [%1d]    :%12zu\n", s, si.ns[s]);
        }

        if (si.coll == 1 || si.coll == 3 || si.coll ==  5 || si.coll == 7) {
            printf("coulomb log i-i  :%12.6lf\n", sd->lii);
            printf("min nu i-i       :%12.6lf\n", sd->niimin);
        }

        printf("\n");

        /* __  __ */
        ticks[1] = clock();
        ttime[1] = time(NULL);
        printf("proc. time (s)   :%12ld\n", (ticks[1]-ticks[0])/CLOCKS_PER_SEC);
        printf("cpu time (s)     :%12d\n", (int)difftime(ttime[1], ttime[0]));
        printf("\n");
    }

    /* __ div. b & pseudo div. b __ */
    for (h = 0; h < 2; h++)
    {
        /* __ clean-up the pointers __ */
        free(dbw[h]);
    }

    /* __ clean-up the pointers __ */
    free(pcw);

    /* __ electric fluctuations __ */
    for (l = 0; l < 3; l++)
    {
        /* __ clean-up the pointers __ */
        free(faw[l]);
        free(fiw[l]);
    }
}
/*===========================================================================*/






/* __ calculation of divergence b ___________________________________________ */
void db(struct sti si, struct stx sx, struct st0 *s0, struct st1 *s1, double wdb[2], int *ppc)
{
double *bw, *cw;
double w1, w2, w3;
int *pw;
int i, j, k;
int ijk;
int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
int l0, l1, l2, m0, m1, m2;


/* __ # of grid points __ */
l0 = sx.n[0];
l1 = sx.n[1];
l2 = sx.n[2];
m0 = sx.n[0]+1;
m1 = sx.n[1]+1;
m2 = sx.n[2]+1;

/* __ memory allocation __ */
bw = (double *)calloc(sx.n[0]*sx.n[1]*sx.n[2], sizeof(double));
cw = (double *)calloc(sx.n[0]*sx.n[1]*sx.n[2], sizeof(double));

pw = (int *)calloc(sx.n[0]*sx.n[1]*sx.n[2], sizeof(int));

/* __ calculation of div b : nested loops on the subdomain __ */
for (i = 0; i < sx.n[0]; i++)
    {
    for (j = 0; j < sx.n[1]; j++)
        {
        for (k = 0; k < sx.n[2]; k++)
            {
            /* __ set index on g0 __ */
            ijk = IDX(i, j, k, l0, l1, l2);

            /* __ set indexes on g1 __ */
            ijk1 = IDX(i  , j  , k  , m0, m1, m2);
            ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
            ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
            ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
            ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
            ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
            ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
            ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);

            /* __ # of parts per cells __ */
            //pw[ijk] = s1[ijk1].ppc;
            pw[ijk] = s0[ijk].ppc;

            /* __ x component of div. b __ */
            w1 = 0.25*(s1[ijk2].b[0]-s1[ijk1].b[0]
                      +s1[ijk4].b[0]-s1[ijk3].b[0]
                      +s1[ijk6].b[0]-s1[ijk5].b[0]
                      +s1[ijk8].b[0]-s1[ijk7].b[0])/si.dl[0];


            /* __ y component of div. b __ */
            w2 = 0.25*(s1[ijk3].b[1]-s1[ijk1].b[1]
                      +s1[ijk4].b[1]-s1[ijk2].b[1]
                      +s1[ijk7].b[1]-s1[ijk5].b[1]
                      +s1[ijk8].b[1]-s1[ijk6].b[1])/si.dl[1];


            /* __ z component of div. b __ */
            w3 = 0.25*(s1[ijk5].b[2]-s1[ijk1].b[2]
                      +s1[ijk6].b[2]-s1[ijk2].b[2]
                      +s1[ijk7].b[2]-s1[ijk3].b[2]
                      +s1[ijk8].b[2]-s1[ijk4].b[2])/si.dl[2];

            /* __ div. b __ */
            bw[ijk] = fabs(w1+w2+w3);

            /* __ pseudo-div. b __ */
            cw[ijk] = ((fabs(w1)+fabs(w2)+fabs(w3)) == 0.0) ? 0.0 : bw[ijk]/(fabs(w1)+fabs(w2)+fabs(w3));
            }
        }
    }

/* __ calculation of the maximum value of div. b & pseudo-div. b __ */
wdb[0] = maxdbl(bw, sx.n[0]*sx.n[1]*sx.n[2]);
wdb[1] = maxdbl(cw, sx.n[0]*sx.n[1]*sx.n[2]);

/* __ min # of parts per cells __ */
*ppc = minint(pw, sx.n[0]*sx.n[1]*sx.n[2]);

/* __ clean-up the pointers __ */
free(bw);
free(cw);
free(pw);

}









/* __ calculation of total energy ___________________________________________ */
void et(struct sti si, struct stx sx, struct st1 *s1, struct st2 *s2, struct stp *sp[1+3*(NS+1)], double wet[9], double wle[3])
{
double bw[3];
double xw, yw, zw;
double lx, ly, lz;
double w1, w2, w3, w4, w5, w6, w7, w8;
double wsv[3], wsb, wst[3], wsa[3];
double aw[3][3];
double wep[3][3], wpe[3][3];
int d, e, f, g, h, i, j, k, l, m, s;
int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
int ijk;
int m0, m1, m2, n0, n1, n2;


/* __ # of grid points __ */
m0 = sx.n[0]+1;
m1 = sx.n[1]+1;
m2 = sx.n[2]+1;
n0 = sx.n[0]+2;
n1 = sx.n[1]+2;
n2 = sx.n[2]+2;

/* __ loop on all components of energy __ */
for (d = 0; d < 1+3*(NS+1); d++)
    {
    /* __ reset energy buffers __ */
    wet[d] = 0.0;
    }

/* __ magnetic energy : nested loops on subdomain indexes __ */
for (i = 0; i < sx.n[0]; i++)
    {
    for (j = 0; j < sx.n[1]; j++)
        {
        for (k = 0; k < sx.n[2]; k++)
            {
            /* __ set indexes on g1 __ */
            ijk1 = IDX(i  , j  , k  , m0, m1, m2);
            ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
            ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
            ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
            ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
            ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
            ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
            ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);

            /* __ loop in the 3 directions __ */
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

                /* __ magnetic energy value __ */
                wet[0] += bw[l]*bw[l];
                }
            }
        }
    }

/* __ initi cfl condition __ */
for (l = 0; l < 3; l++)
    {
    wle[l] = 0.0;
    }

/* __ protons & alphas energies : loop on the species __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ loop on the part. of specie s __ */
    for (m = 0; m < sx.ns[s]; m++)
        {
        /* __ cfl condition __ */
        if (fabs(sp[s][m].v[0]) > wle[0]) wle[0] = fabs(sp[s][m].v[0]);
        if (fabs(sp[s][m].v[1]) > wle[1]) wle[1] = fabs(sp[s][m].v[1]);
        if (fabs(sp[s][m].v[2]) > wle[2]) wle[2] = fabs(sp[s][m].v[2]);

        /* __ part. "position" in the subdomain __ */
        xw = sp[s][m].r[0]/si.dl[0]-sx.i0[0]+0.5;
        yw = sp[s][m].r[1]/si.dl[1]-sx.i0[1]+0.5;
        zw = sp[s][m].r[2]/si.dl[2]-sx.i0[2]+0.5;

        /* __ index for the part. "position" __ */
        i = (int)floor(xw);
        j = (int)floor(yw);
        k = (int)floor(zw);

        /* __ part. location in the grid __ */
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

        /* __ indexes of the rounding grid points on g2 __ */
        ijk1 = IDX(i  , j  , k  , n0, n1, n2);
        ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
        ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
        ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
        ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
        ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
        ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
        ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

        /* __ loop in each directions __ */
        for (l = 0; l < 3; l++)
            {
            /* __ fluid velocity of specie s __ */
            wsv[l] = w1*s2[ijk1].vs[s][l]
                    +w2*s2[ijk2].vs[s][l]
                    +w3*s2[ijk3].vs[s][l]
                    +w4*s2[ijk4].vs[s][l]
                    +w5*s2[ijk5].vs[s][l]
                    +w6*s2[ijk6].vs[s][l]
                    +w7*s2[ijk7].vs[s][l]
                    +w8*s2[ijk8].vs[s][l];
            }

        /* __ kinetic energy of specie s __ */
        wsb = wsv[0]*wsv[0]+wsv[1]*wsv[1]+wsv[2]*wsv[2];

        /* __ loop in each direction __ */
        for (l = 0; l < 3; l++)
            {
            /* __ set thermal velocity __ */
            wst[l] = sp[s][m].v[l]-wsv[l];
            }

        /* __ part. "position" in the subdomain __ */
        xw = (sp[s][m].r[0]/si.dl[0]-sx.i0[0] < sx.n[0]) ?
              sp[s][m].r[0]/si.dl[0]-sx.i0[0] : sx.n[0]-EPS4;
        yw = (sp[s][m].r[1]/si.dl[1]-sx.i0[1] < sx.n[1]) ?
              sp[s][m].r[1]/si.dl[1]-sx.i0[1] : sx.n[1]-EPS4;
        zw = (sp[s][m].r[2]/si.dl[2]-sx.i0[2] < sx.n[2]) ?
              sp[s][m].r[2]/si.dl[2]-sx.i0[2] : sx.n[2]-EPS4;

        /* __ index for the part. "position" __ */
        i = (int)floor(xw);
        j = (int)floor(yw);
        k = (int)floor(zw);

        /* __ part. location in the grid __ */
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

        /* __ indexes of the rounding grid points on g2 __ */
        ijk1 = IDX(i  , j  , k  , m0, m1, m2);
        ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
        ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
        ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
        ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
        ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
        ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
        ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);

        /* __ loop in each direction __ */
        for (l = 0; l < 3; l++)
            {
            /* __ set the magnetic field value __ */
            bw[l] = w1*s1[ijk1].b[l]
                   +w2*s1[ijk2].b[l]
                   +w3*s1[ijk3].b[l]
                   +w4*s1[ijk4].b[l]
                   +w5*s1[ijk5].b[l]
                   +w6*s1[ijk6].b[l]
                   +w7*s1[ijk7].b[l]
                   +w8*s1[ijk8].b[l];
            }

        /* __ base of orthonormal vectors __ */
        ortho(bw, aw);

        /* __ loop in each directions __ */
        for (l = 0; l < 3; l++)
            {
            /* __ set value of field-aligned pressure __ */
            wsa[l] = 0.0;

            /* __ nested loop in each direction __ */
            for (e = 0; e < 3; e++)
                {
                for (f = 0; f < 3; f++)
                    {
                    /* __ thermal pressure in field-aligned coordinates __ */
                    wsa[l] += aw[l][e]*wst[e]*wst[f]*aw[l][f];
                    }
                }
            }

        // then cumulative sum for "s specie" energies
        wet[1+3*s] += si.ws[s]*si.ms[s]*wsb;             // bulk energy
        wet[2+3*s] += si.ws[s]*si.ms[s]*wsa[0];          // parallel thermal energy
        wet[3+3*s] += si.ws[s]*si.ms[s]*(wsa[1]+wsa[2]); // perp thermal energy

        }
    }

/* __ electron energies : nested loops in each direction __ */
for (i = 1; i < sx.n[0]+1; i++)
    {
    for (j = 1; j < sx.n[1]+1; j++)
        {
        for (k = 1; k < sx.n[2]+1; k++)
            {
            /* __ set index on g2 __ */
            ijk = IDX(i  , j  , k  , n0, n1, n2);

            /* __ set indexes on g1 __ */
            ijk1 = IDX(i  , j  , k  , m0, m1, m2);
            ijk2 = IDX(i-1, j  , k  , m0, m1, m2);
            ijk3 = IDX(i  , j-1, k  , m0, m1, m2);
            ijk4 = IDX(i-1, j-1, k  , m0, m1, m2);
            ijk5 = IDX(i  , j  , k-1, m0, m1, m2);
            ijk6 = IDX(i-1, j  , k-1, m0, m1, m2);
            ijk7 = IDX(i  , j-1, k-1, m0, m1, m2);
            ijk8 = IDX(i-1, j-1, k-1, m0, m1, m2);

            /* __ loop in each direction __ */
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

            /* __ base of orthonormal vectors __ */
            ortho(bw, aw);

            /* __ set electron pressure tensor __ */
            wep[0][0] = s2[ijk].ps[0][0];
            wep[0][1] = s2[ijk].ps[0][1];
            wep[0][2] = s2[ijk].ps[0][2];
            wep[1][0] = s2[ijk].ps[0][1];
            wep[1][1] = s2[ijk].ps[0][3];
            wep[1][2] = s2[ijk].ps[0][4];
            wep[2][0] = s2[ijk].ps[0][2];
            wep[2][1] = s2[ijk].ps[0][4];
            wep[2][2] = s2[ijk].ps[0][5];

            /* __ nested loops in each direction __ */
            for (e = 0; e < 3; e++)
                {
                for (f = 0; f < 3; f++)
                    {
                    /* __ set value of field-aligned pressure __ */
                    wpe[e][f] = 0.0;

                    /* __ nested loops in each direction __ */
                    for (g = 0; g < 3; g++)
                        {
                        for (h = 0; h < 3; h++)
                            {
                            /* __ set value of field-aligned pressure __ */
                            wpe[e][f] += aw[e][g] * wep[h][g] * aw[f][h];
                            }
                        }
                    }
                }

            /* __ set parallel & perpendicular thermal energies __ */
            wet[2] +=           wpe[0][0];
            wet[3] += wpe[1][1]+wpe[2][2];
            }
        }
    }

/* __ loop on all components of energy __ */
for (d = 0; d < 1+3*(NS+1); d++)
    {
    /* _____ normalization _____ */
    wet[d] *= 0.5*si.dl[0]*si.dl[1]*si.dl[2];
    }

}


/* __ calculation of the level of electric fluctuations _____________________ */
void ef(struct stx sx, struct st2 *s2, double wa[3], double wi[3])
{
double *ew[3];
int l;
int ijk;
int n2;


/* __ # of grid points on g2 __ */
n2 = (sx.n[0]+2)*(sx.n[1]+2)*(sx.n[2]+2);

/* __ loop in each direction __ */
for (l = 0; l < 3; l++)
    {
    /* __ memory allocation __ */
    ew[l] = (double *)malloc(n2*sizeof(double));
    }

/* __ loop on g2 grid __ */
for (ijk = 0; ijk < n2; ijk++)
    {
    /* __ loop in each direction __ */
    for (l = 0; l < 3; l++)
        {
        /* __ fill buffer with electric field components __ */
        ew[l][ijk] = s2[ijk].e[l];
        }
    }

/* __ loop in each direction __ */
for (l = 0; l < 3; l++)
    {
    /* __ find max and min value of electric field __ */
    wa[l] = maxdbl(ew[l], n2);
    wi[l] = mindbl(ew[l], n2);
    }

/*___ loop in each direction __ */
for (l = 0; l < 3; l++)
    {
    /*___ clean-up the pointers __ */
    free(ew[l]);
    }

}


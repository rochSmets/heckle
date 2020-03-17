

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defines.h"
#include "structures.h"
#include "misc.h"
#include "plasma.h"


/* __ read "munster.txt" to set stt structure (magnetic fluctuations) _____ */
void topo(struct sti si, struct stx sx, struct stt *st, MPI_Comm com)
{
double w1b, w1e;
double wb, wbx, wby;
double we, wex, wey, wez;
double kz;
double kzb, kze;
int nm;
int f, z;
MPI_Comm co;
char junk[80];
FILE *fp;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ open the file __ */
fp = fopen("munster.txt", "r");
if (fp == NULL) printf("problem in opening file munster.txt\n");

/* __ read the topology parameters __ */
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", st->b, st->b+1, st->b+2); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", st->n, st->n+1, st->n+2); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->k)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->e)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->p)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", st->z, st->z+1); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->gb)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", st->t, st->t+1, st->t+2); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->ge)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%d %d", &(st->m[0]), &(st->m[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%d", &(st->td)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%d", &(st->sd)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);

/* __ close the file __ */
fclose(fp);

/* __ # of modes __ */
nm = st->m[1]-st->m[0]+1;

/* __ memory allocation for magnetic & electric energy density & phase __ */
st->bx = (double *)malloc(nm*sizeof(double));
st->by = (double *)malloc(nm*sizeof(double));
st->ex = (double *)malloc(nm*sizeof(double));
st->ey = (double *)malloc(nm*sizeof(double));
st->ez = (double *)malloc(nm*sizeof(double));
st->px = (double *)malloc(nm*sizeof(double));
st->py = (double *)malloc(nm*sizeof(double));
st->fx = (double *)malloc(nm*sizeof(double));
st->fy = (double *)malloc(nm*sizeof(double));
st->fz = (double *)malloc(nm*sizeof(double));

/* __ set initial sum __ */
w1b = 0.0;
w1e = 0.0;

/* __ summation over the wave numbers : loop on the wave numbers __ */
for (f = st->m[0]; f <= st->m[1]; f++)
    {
    /* __ wave number in z dir. __ */
    kz = 2.0*f*PI/(st->sd*si.dl[2]);

    /* __ wave number at -gamma power __ */
    kzb = pow(kz,-st->gb);
    kze = pow(kz,-st->ge);

    /* __ sum for magnetic field __ */
    w1b += kzb;

    /* __ sum for electric field __ */
    w1e += kze;
    }

/* __ magnetic energy __ */
wb = 0.5;

/* __ electric energy __ */
we = 0.5;

/* __ magnetic energy density in x & y dir. __ */
wbx = sqrt(st->z[0]*wb/w1b)*(si.n[2]/st->sd);
wby = sqrt(st->z[1]*wb/w1b)*(si.n[2]/st->sd);

/* __ electric energy density in x, y or z dir. __ */
wex = sqrt(st->t[0]*we/w1e)*(si.n[2]/st->sd);
wey = sqrt(st->t[1]*we/w1e)*(si.n[2]/st->sd);
wez = sqrt(st->t[2]*we/w1e)*(si.n[2]/st->sd);

/* __ weighted energy for each modes : loop on the wave numbers __ */
for (f = st->m[0]; f <= st->m[1]; f++)
    {
    /* __ wave number in z dir. __ */
    kz = 2.0*f*PI/(st->sd*si.dl[2]);

    /* __ wave number at -gamma/2 power __ */
    kzb = pow(kz,-0.5*st->gb);
    kze = pow(kz,-0.5*st->ge);

    /* __ weighted energy for each modes __ */
    st->bx[f] = wbx*kzb;
    st->by[f] = wby*kzb;
    st->ex[f] = wex*kze;
    st->ey[f] = wey*kze;
    st->ez[f] = wez*kze;
    }

/* __ set the modes' phases : only on node 0 __ */
if (sx.r == 0)
   {
   /* __ nested loops on the wave numbers __ */
   for (f = st->m[0]; f <= st->m[1]; f++)
       {
       /* __ set randomly the wave phase __ */
       st->px[f] = RNM*2.0*PI;
       st->py[f] = RNM*2.0*PI;
       st->fx[f] = RNM*2.0*PI;
       st->fy[f] = RNM*2.0*PI;
       st->fz[f] = RNM*2.0*PI;
       }
   }

/* __ broadcast the phases of each modes __ */
MPI_Bcast(st->px, nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->py, nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->fx, nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->fy, nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->fz, nm, MPI_DOUBLE, 0, co);

/* __ display the topology parameters __ */
if (sx.r == 0)
   {
   printf("________________ magnetic topology : b & e fluct. ____\n");
   printf("b              : %8.4f    %8.4f    %8.4f\n", st->b[0], st->b[1], st->b[2]);
   printf("proton density : %8.4f\n", st->n[1]);
   printf("alpha density  : %8.4f\n", st->n[2]);
   printf("beta protons   : %8.4f\n", st->p*(st->k/wb-1.0));
   printf("beta electrons : %8.4f\n", (1.0-st->p)*(st->k/wb-1.0));
   printf("beta total     : %8.4f\n", st->k/wb-1.0);
   printf("ener. b (x)    : %8.4f\n", st->z[0]);
   printf("ener. b (y)    : %8.4f\n", st->z[1]);
   printf("gamma (b)      : %8.4f\n", st->gb);
   printf("ener. e (x)    : %8.4f\n", st->t[0]);
   printf("ener. e (y)    : %8.4f\n", st->t[1]);
   printf("ener. e (z)    : %8.4f\n", st->t[2]);
   printf("gamma (v)      : %8.4f\n", st->ge);
   printf("min & max modes: %8d    %8d\n", st->m[0], st->m[1]);
   printf("# ts for drive : %8d\n", st->td);
   printf("\n");
   printf("______________________________________________________\n");
   printf("\n");
   }

/* __ clean up the communicator __ */
MPI_Comm_free(&co);
}


/* __ density _______________________________________________________________ */
void density(struct sti si, struct stx sx, struct stt *st, double q[3], double n[NS+1], MPI_Comm com)
{


/* __ electron density __ */
n[0] = st->n[0];

/* __ proton density __ */
n[1] = st->n[1];

/* __ alpha density __ */
n[2] = st->n[2];

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{


/* __ dc magnetic field components __ */
b[0] = st->b[0];
b[1] = st->b[1];
b[2] = st->b[2];

}


/* __ drive the magnetic & electric field ___________________________________ */
void drive(struct sti si, struct stx sx, struct stt *st, struct st1 *s1, struct st2 *s2, MPI_Comm com)
{
double q[3];
double kz;
double kq;
double wt;
int ijk;
int f, i, j, k;
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ factor associated to time forcing __ */
wt = 1.0/(st->td+1);

/* __ superposition of each modes for mag. field : nested loops on subdomains __ */
for (i = 0; i < sx.n[0]+1; i++)
    {
    for (j = 0; j < sx.n[1]+1; j++)
        {
        for (k = 0; k < sx.n[2]+1; k++)
            {
            /* __ set index on g2 __ */
            ijk = idx(i, j, k, sx.n[0]+1, sx.n[1]+1, sx.n[2]+1);

            /* __ driving only if in the appropriate space location __ */
            if (sx.i0[2]+k <= st->sd)
               {
               /* __ coordinates of grid point __ */
               q[0] = (sx.i0[0]+i)*si.dl[0];
               q[1] = (sx.i0[1]+j)*si.dl[1];
               q[2] = (sx.i0[2]+k)*si.dl[2];

               /* __ loops on the modes __ */
               for (f = st->m[0]; f <= st->m[1]; f++)
                   {
                   /* __ wave number in z dir. __ */
                   kz = 2.0*f*PI/(st->sd*si.dl[2]);

                   /* __ set k dot q __ */
                   kq = kz*q[2];

                   /* __ set the mag. field components __ */
                   s1[ijk].b[0] += wt*st->bx[f]*cos(kq+st->px[f]);
                   s1[ijk].b[1] += wt*st->by[f]*cos(kq+st->py[f]);
                   }
               }
            }
        }
    }

/* __ superposition of each modes for elec field : nested loops on subdomains __ */
for (i = 0; i < sx.n[0]+2; i++)
    {
    for (j = 0; j < sx.n[1]+2; j++)
        {
        for (k = 0; k < sx.n[2]+2; k++)
            {
            /* __ set index on g2 __ */
            ijk = idx(i, j, k, sx.n[0]+2, sx.n[1]+2, sx.n[2]+2);

            /* __ driving only if in the appropriate space location __ */
            if (sx.i0[2]+k <= st->sd)
               {
               /* __ coordinates of grid point __ */
               q[0] = (sx.i0[0]+i-0.5)*si.dl[0];
               q[1] = (sx.i0[1]+j-0.5)*si.dl[1];
               q[2] = (sx.i0[2]+k-0.5)*si.dl[2];

               /* __ loops on the modes __ */
               for (f = st->m[0]; f <= st->m[1]; f++)
                   {
                   /* __ wave number in z dir. __ */
                   kz = 2.0*f*PI/(st->sd*si.dl[2]);

                   /* __ set k dot q __ */
                   kq = kz*q[2];

                   /* __ set the elec. field components __ */
                   s2[ijk].e[0] += wt*st->ex[f]*cos(kq+st->fx[f]);
                   s2[ijk].e[1] += wt*st->ey[f]*cos(kq+st->fy[f]);
                   s2[ijk].e[2] += wt*st->ez[f]*cos(kq+st->fz[f]);
                   }
               }
            }
        }
    }

/* __ clean up the communicator __ */
MPI_Comm_free(&co);

}


/* __ electric field ________________________________________________________ */
void electric(struct sti si, struct stx sx, struct stt *st, double q[3], double e[3], MPI_Comm com)
{


/* __ electric field components __ */
e[0] = 0.0;
e[1] = 0.0;
e[2] = 0.0;

}


/* __ rotational of magnetic field __________________________________________ */
void current(struct sti si, struct stx sx, struct stt *st, double q[3], double j[3], MPI_Comm com)
{


/* __ current components __ */
j[0] = 0.0;
j[1] = 0.0;
j[2] = 0.0;

}


/* __ drift velocity ________________________________________________________ */
void drift(struct sti si, struct stx sx, struct stt *st, double q[3], double v[NS+1][3], MPI_Comm com)
{


/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = 0.0;

/* __ proton velocity __ */
v[1][0] = 0.0;
v[1][1] = 0.0;
v[1][2] = 0.0;

/* __ alpha velocity __ */
v[2][0] = 0.0;
v[2][1] = 0.0;
v[2][2] = 0.0;

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{
double wt;
double we, wp, wa;
double wi0, wi1;


/* __ total kinetic pressure __ */
wt = st->k-0.5*(st->b[0]*st->b[0]+st->b[1]*st->b[1]+st->b[2]*st->b[2]);

/* __ kinetic pressure of electrons, protons & alphas __ */
we = st->e*wt;
wp = st->p*wt;
wa = (1.0-st->e-st->p)*wt;

/* __ electron temperatures __ */
ta[0] = we/st->n[0];
te[0] = we/st->n[0];

/* __  proton temperatures __ */
ta[1] = wp/st->n[1];
te[1] = wp/st->n[1];

/* __  alpha temperatures __ */
ta[2] = (st->n[2] != 0.0) ? wa/st->n[2] : 0.0;
te[2] = (st->n[2] != 0.0) ? wa/st->n[2] : 0.0;

}


/* __ calculation of the normalization coefficients _________________________ */
void normal(struct sti *si, struct stx sx, struct stt *st, MPI_Comm com)
{
double wq[3], wn[3];
double wd, wt;
int i, j, k, s;
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ volume of a cell __ */
wd = si->dl[0]*si->dl[1]*si->dl[2];

/* __ loop on the part. __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ set initial value of counter __ */
    wt = 0.0;

    /* __ nested loops on the cells in the total domain __ */
    for (i = 0; i < si->n[0]; i++)
        {
        for (j = 0; j < si->n[1]; j++)
            {
            for (k = 0; k < si->n[2]; k++)
                {
                /* __ crux of the cell __ */
                wq[0] = ((double)i+0.5)*si->dl[0];
                wq[1] = ((double)j+0.5)*si->dl[1];
                wq[2] = ((double)k+0.5)*si->dl[2];

                /* __ local density for part. s __ */
                density(*si, sx, st, wq, wn, co);

                /* __ integrated local "mass" of part. s __ */
                wt += wn[s]*wd;
                }
            }
        }

    /* __ weight of part. "s" __ */
    si->ws[s] = (si->ns[s] == 0) ? EPS6 : wt/(si->ns[s]*wd);
    }

/* __ clean up the communicator __ */
MPI_Comm_free(&co);

}


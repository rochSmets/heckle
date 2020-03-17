

#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defines.h"
#include "structures.h"
#include "misc.h"
#include "plasma.h"


/* __ read "munster.txt" to set stt structure (magnetic fluctuations) _____ */
void topo(struct sti si, struct stx *sx, struct stt *st, MPI_Comm com)
{
double w1b;
double wb, wbx, wby, wbz;
double kx, ky;
double kxb, kyb;
int nm;
int f, h, fh;
MPI_Comm co;
char junk[80];
FILE *fp;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ open the file __ */
fp = fopen("electrosmog.txt", "r");
if (fp == NULL) printf("problem in opening file munster.txt\n");

/* __ read the topology parameters __ */
fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->m[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->m[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->l[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->l[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->te), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->th1), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->d[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->d[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->a[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->a[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->eb[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->eb[1]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->eb[2]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->slb[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->slb[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->ix[0]), &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->ix[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->iy[0]), &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->iy[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->td), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->sd[0]), &(sx->irun));
fscanint(__FILE__, __LINE__, fp, sx->r, &(st->sd[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->c[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->c[1]), &(sx->irun));


/* __ close the file __ */
fclose(fp);

/* __ # of modes __ */
nm = (st->ix[1]-st->ix[0]+1)*(st->iy[1]-st->iy[0]+1);

/* __ memory allocation for magnetic energy density & phase __ */
st->bx = (double *)malloc(nm*sizeof(double));
st->by = (double *)malloc(nm*sizeof(double));
st->bz = (double *)malloc(nm*sizeof(double));
st->px = (double *)malloc(nm*sizeof(double));
st->py = (double *)malloc(nm*sizeof(double));
st->pz = (double *)malloc(nm*sizeof(double));

/* __ set initial sum __ */
w1b = 0.0;

/* __ summation over the wave numbers x direction: loop on the wave numbers __ */
for (f = st->ix[0]; f <= st->ix[1]; f++)
    {
    /* __ summation over the wave numbers y direction : loop on the wave numbers __ */
    for (h = st->iy[0]; h <= st->iy[1]; h++)
        {
        /* __ wave number in x dir. __ */
        kx = 2.0*f*PI/(st->sd[0]*si.dl[0]);

        /* __ wave number at -gamma power __ */
        kxb = (kx != 0.0) ? pow(kx, -st->slb[0]) : 0.0;

        /* __ wave number in z dir. __ */
        ky = 2.0*h*PI/(st->sd[1]*si.dl[1]);

        /* __ wave number at -gamma power __ */
        kyb = (ky != 0.0) ? pow(ky, -st->slb[1]) : 0.0;

        /* __ sum over the energies __*/
        w1b += (ky != 0.0) ? kxb*kyb*(1+(pow(kx, 2)/pow(ky, 2))) : 0.0;
        }
    }

/* __ magnetic energy __ */
wb = 0.5;

/* __ magnetic energy density in x & y dir(C in  the latex document). __ */
wbx = (st->eb[0]*wb/w1b)*((si.n[0]*si.n[1])/(st->sd[0]*st->sd[1]) >= 0.0) ? sqrt(st->eb[0]*wb/w1b)*((si.n[0]*si.n[1])/(st->sd[0]*st->sd[1])) : 0.0;
wby = (st->eb[1]*wb/w1b)*((si.n[0]*si.n[1])/(st->sd[0]*st->sd[1]) >= 0.0) ? sqrt(st->eb[1]*wb/w1b)*((si.n[0]*si.n[1])/(st->sd[0]*st->sd[1])) : 0.0;
wbz = 0.0;

/* __ weighted energy for each modes : loop on the wave numbers __ */
for (f = st->ix[0]; f <= st->ix[1]; f++)
    {
    for (h = st->iy[0]; h <= st->iy[1]; h++)
        {
//	fh = ID(f, h, st->ix[1]-st->ix[0]+1, st->iy[1]-st->iy[0]+1);
        fh = IDX(f, h, 0, st->ix[1]-st->ix[0]+1, st->iy[1]-st->iy[0]+1, 1);

        /* __ wave number in z dir. __ */
        kx = 2.0*f*PI/(st->sd[0]*si.dl[0]);
        ky = 2.0*h*PI/(st->sd[1]*si.dl[1]);

        /* __ wave number at -gamma/2 power __ */
        kxb = (kx != 0.0) ? pow(kx, -0.5*st->slb[0]) : 0.0;
        kyb = (ky != 0.0) ? pow(ky, -0.5*st->slb[1]) : 0.0;

        /* __ weighted energy for each modes __ */
        st->bx[fh] = wbx*kxb*kyb;
        st->by[fh] = (ky !=0.0) ? wby*kxb*kyb*(-kx/ky) : 0.0;
        st->bz[fh] = 0.0;
        }
     }

/* __ set the modes' phases : only on node 0 __ */
if (sx->r == 0)
   {
   /* __ nested loops on the wave numbers __ */
   for (f = st->ix[0]; f <= st->ix[1]; f++)
       {
       for (h = st->iy[0]; h <= st->iy[1]; h++)
           {
//         fh = ID(f, h, st->ix[1]-st->ix[0]+1, st->iy[1]-st->iy[0]+1);
           fh = IDX(f, h, 0, st->ix[1]-st->ix[0]+1, st->iy[1]-st->iy[0]+1, 1);

           /* __ set randomly the wave phase __ */
           st->px[fh] = RNM*2.0*PI;
           st->py[fh] = RNM*2.0*PI;
           st->pz[fh] = RNM*2.0*PI;
           }
       }
   }

/* __ broadcast the phases of each modes __ */
MPI_Bcast(st->px, nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->py, nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->pz, nm, MPI_DOUBLE, 0, co);

/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : b & e fluct. ____\n");
   printf("Magnetic field                         : %8.4f    %8.4f\n", st->b[0], st->b[1]);
   printf("Proton and alpha density               : %8.4f    %8.4f\n", st->m[0], st->m[1]);
   printf("sheet thickness 1 and 2                : %8.4f    %8.4f\n", st->l[0], st->l[1]);
   printf("temp(e)                                : %8.4f\n", st->te);
   printf("Temperature in the first zone (protons): %8.4f\n", st->th1); 
   printf("1 & 2 discontinuity positions          : %8.4f    %8.4f\n", st->d[0], st->d[1]);
   printf("Angle of magnetic field 1 & 2          : %8.4f    %8.4f\n", st->a[0], st->a[1]);  
   printf("spec. mag. en. : %8.4f    %8.4f\n", st->eb[0], st->eb[1]);
   printf("slope b spec.  : %8.4f    %8.4f\n", st->slb[0], st->slb[1]);
   printf("min & max modes x durection: %8d    %8d\n", st->ix[0], st->ix[1]);
   printf("min & max modes y direction: %8d    %8d\n", st->iy[0], st->iy[1]);
   printf("time steps 4 dr: %8d\n", st->td);
   printf("grid cells 4 dr: %8d    %8d\n", st->sd[0], st->sd[1]);
   printf("antenne center position x & y  : %8f    %8f\n", st->c[0], st->c[1]);
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

double y1, y2;


/* __ location of the 2 current sheets __ */
y1 = (q[1]-st->d[0]*si.l[1])/st->l[0];
y2 = (q[1]-st->d[1]*si.l[1])/st->l[1];

/* __ electron density __ */
n[0] = 0.5*(st->m[0])*(tanh(y2)-tanh(y1))
      +st->m[0]
      +0.5*(st->m[1])*(tanh(y1)-tanh(y2));

/* __ proton density __ */
n[1] = 0.5*(st->m[0])*(tanh(y2)-tanh(y1))
      +st->m[0];

/* __ alpha density __ */
n[2] = 0.5*(st->m[1])*(tanh(y1)-tanh(y2));

}

/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{

double y1, y2;
double btot, alpha;
double w1;
double w2;

w1 = 0.5*(st->b[0]*st->b[0])+(st->th1+st->te)*st->m[0];
w2 = (q[1] < st->d[0]*si.l[1] && q[1] > st->d[1]*si.l[1]) ? st->th1 : ((w1-0.5*(st->b[1]*st->b[1]))/st->m[1])-st->te;

/* __ location of the 2 current sheets __ */
y1 = (q[1]-st->d[0]*si.l[1])/st->l[0];
y2 = (q[1]-st->d[1]*si.l[1])/st->l[1];

/* __ Magnetic field module __ */
btot = (2*w1-2*(w2+st->te)*((0.5*(st->m[1]-st->m[0]))*(tanh(y1)-tanh(y2))+st->m[0]) >= 0.0) ? pow(2*w1-2*(w2+st->te)*((0.5*(st->m[1]-st->m[0]))*(tanh(y1)-tanh(y2))+st->m[0]) , 0.5) : 0.0;

/* __ Angle __ */
alpha = 0.5*(st->a[1]-st->a[0])*(tanh(y1)-tanh(y2))+st->a[0];

printf("alpha= %8.4f\n", alpha);

/* __ magnetic field components __ */
b[0] = btot*cos(alpha);
b[1] = 0.0;
b[2] = btot*sin(alpha);
}


/* __ drive the magnetic & electric field ___________________________________ */
void drive(struct sti si, struct stx sx, struct stt *st, struct st1 *s1, struct st2 *s2, MPI_Comm com)
{
double q[3], x0, y0;
double kx, ky;
double kqx, kqy, kq;
double wt;
int ijk;
int f, i, j, k, h, fh;
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ factor associated to time forcing __ */
wt = 1.0/(st->td+1);

/* __ Position of the center of the waves source __ */
x0 = st->c[0]*si.n[0];  /* There is a # of grid points */
y0 = st->c[1]*si.n[1];  /* There is a # of grid points */

/* __ superposition of each modes for mag. field : nested loops on subdomains __ */
for (i = 0; i < sx.n[0]+1; i++)
    {
    for (j = 0; j < sx.n[1]+1; j++)
        {
        for (k = 0; k < sx.n[2]+1; k++)
            {
            /* __ set index on g2 __ */
            ijk = IDX(i, j, k, sx.n[0]+1, sx.n[1]+1, sx.n[2]+1);

            /* __ driving only if in the appropriate space location __ */
            if ((x0 - 0.5*st->sd[0]) <= (i+sx.i0[0])
             && (i+sx.i0[0]) <= (x0 + 0.5*st->sd[0])
             && (y0 - 0.5*st->sd[1]) <= (j+sx.i0[1])
             && (j+sx.i0[1]) <= (y0 + 0.5*st->sd[1]))
               {
               /* __ coordinates of grid point __ */
               q[0] = (sx.i0[0]+i)*si.dl[0];
               q[1] = (sx.i0[1]+j)*si.dl[1];
               q[2] = (sx.i0[2]+k)*si.dl[2];


               /* __ loops on the modes x __ */
               for (f = st->ix[0]; f <= st->ix[1]; f++)
                   {
                   /* __ loops on the modes y __ */
                   for (h = st->iy[0]; h <= st->iy[1]; h++)
                       {
                       fh = IDX(f, h, 0, st->ix[1]-st->ix[0]+1, st->iy[1]-st->iy[0]+1, 1);
//                     fh = ID(f, h, st->ix[1]-st->ix[0]+1, st->iy[1]-st->iy[0]+1);
                       /* __ wave number in x dir. __ */
                       kx = 2.0*f*PI/(st->sd[0]*si.dl[0]);
                       /* __ set k dot q __ */
                       kqx = kx*q[0];
                       /* __ wave number in y dir. __ */
                       ky = 2.0*h*PI/(st->sd[1]*si.dl[1]);
                       /* __ set k dot q __ */
                       kqy = ky*q[1];
                       kq = kqy + kqx;
                       /* __ set the mag. field components __ */
                       s1[ijk].b[0] += wt*st->bx[fh]*cos(kq+st->px[fh]);
                       s1[ijk].b[1] += wt*st->by[fh]*cos(kq+st->py[fh]);
                       s1[ijk].b[2] = 0.0; //+= wt*st->bz[fh]*cos(kq+st->pz[fh]);
                       }
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

double w1 = 0.5*(st->b[0]*st->b[0])+(st->th1+st->te)*st->m[0]; /* w1 is the total pressure*/

/* __ electron temperatures __ */
ta[0] = st->te;
te[0] = st->te;

/* __  proton temperatures __ */
ta[1] = st->th1;
te[1] = st->th1;

/* __  alpha temperatures __ */
ta[2] = ((w1-0.5*(st->b[1]*st->b[1]))/st->m[1])-st->te;
te[2] = ((w1-0.5*(st->b[1]*st->b[1]))/st->m[1])-st->te;
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


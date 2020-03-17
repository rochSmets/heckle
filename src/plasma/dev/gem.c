
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "plasma.h"


/* __ read the initial parameters to set stt structure (uniform plasma) _____ */
void topo(struct sti si, struct stx sx, struct stt *st, MPI_Comm com)
{
int z;
char junk[80];
FILE *fp;


/* __ read the topology parameters __ */
fp = fopen("gem.txt", "r");
if (fp == NULL) printf("problem in opening file uniform.txt\n");

z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->b)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->n)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->m)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->k)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->l)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->a)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->p)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);

fclose(fp);

/* __ display the topology parameters __ */
if (sx.r == 0)
   {
   printf("________________ magnetic topology : uniform _________\n");
   printf("b              : %8.4f\n", st->b);
   printf("n protons      : %8.4f\n", st->n);
   printf("n alphas       : %8.4f\n", st->m);
   printf("beta total     : %8.4f\n", 2.0*st->k/(st->b*st->b)-1.0);
   printf("sheet thickness: %8.4f\n", st->l);
   printf("beta(e)/beta(p): %8.4f\n", st->a);
   printf("mag. perturb.  : %8.4f\n", st->p);
   printf("\n");
   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ density _______________________________________________________________ */
void density(struct sti si, struct stx sx, struct stt *st, double q[3], double n[NS+1], MPI_Comm com)
{


/* __ shift y component __ */
q[1] -= 0.5*si.l[1];

/* __ electron density __ */
n[0] = st->n/pow(cosh(q[1]/st->l),2)+st->m;

/* __ proton density __ */
n[1] = st->n/pow(cosh(q[1]/st->l),2);

/* __ alpha density __ */
n[2] = st->m;

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{
double wx, wy;


/* __ needed constants __ */
wx = 2.0*PI/si.l[0];
wy = PI/si.l[1];

/* __ shift x & y component __ */
q[0] -= 0.5*si.l[0];
q[1] -= 0.5*si.l[1];

/* __ magnetic field components __ */
b[0] = -wy*st->p*cos(wx*q[0])*sin(wy*q[1])+st->b*tanh(q[1]/st->l);
b[1] = +wx*st->p*sin(wx*q[0])*cos(wy*q[1]);
b[2] = 0.0;

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
double beta;


/* __ "beta" value __ */
beta = 2.0*st->k/(st->b*st->b)-1.0;

/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = 0.0;

/* __ proton velocity __ */
v[1][0] = 0.0;
v[1][1] = 0.0;
v[1][2] = sqrt(2.0*beta)/st->l;

/* __ alpha velocity __ */
v[2][0] = 0.0;
v[2][1] = 0.0;
v[2][2] = 0.0;

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{
double wt;


/* __ total kinetic pressure __ */
wt = st->k-0.5*(st->b*st->b);

/* __ electron temperatures __ */
ta[0] = (1.0-st->a)*wt;
te[0] = (1.0-st->a)*wt;

/* __  proton temperatures __ */
ta[1] = st->a*wt;
te[1] = st->a*wt;

/* __  alpha temperatures __ */
ta[2] = 0.0;
te[2] = 0.0;

}


/* __ calculation of the normalization coefficients _________________________ */
void normal(struct sti *si, struct stx sx, struct stt *st, MPI_Comm com)
{
double wq[3], wn[NS+1];
double wd, wt;
int i, j, k, s;


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
                density(*si, sx, st, wq, wn, com);

                /* __ integrated local "mass" of part. s __ */
                wt += wn[s]*wd;
                }
            }
        }

    /* __ weight of part. "s" __ */
    si->ws[s] = (si->ns[s] == 0) ? EPS6 : wt/(si->ns[s]*wd);
    }

}


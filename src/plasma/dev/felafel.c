
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "structures.h"
#include "plasma.h"


/* __ read the initial parameters to set stt structure _____ */
void topo(struct sti si, struct stx sx, struct stt *st, MPI_Comm com)
{
int z;
char junk[80];
FILE *fp;


/* __ read the topology parameters __ */
fp = fopen("felafel.txt", "r");
if (fp == NULL) printf("problem in opening file felafel.txt\n");

z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", &(st->b[0]), &(st->b[1]), &(st->b[2])); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", &(st->n[0]), &(st->n[1]), &(st->n[2])); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", &(st->v[0]), &(st->v[1]), &(st->v[2])); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->l)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->t)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->k)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->e)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->p)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->i)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->j)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);

fclose(fp);

/* __ display the topology parameters __ */
if (sx.r == 0)
   {
   printf("________________ magnetic topology : felafel _________\n");
   printf("b              : %8.4f    %8.4f    %8.4f\n", st->b[0], st->b[1], st->b[2]);
   printf("n protons      : %8.4f\n", st->n[1]);
   printf("n alphas       : %8.4f\n", st->n[2]);
   printf("v species      : %8.4f    %8.4f    %8.4f\n", st->v[0], st->v[1], st->v[2]);
   printf("i (|/#) protons: %8.4f\n", st->i);
   printf("i (|/#) alphas : %8.4f\n", st->j);
   printf("beta protons   : %8.4f\n", st->p);
   printf("beta alphas    : %8.4f\n", 1.0-st->p-st->e);
   printf("beta electrons : %8.4f\n", st->e);
   printf("beta total     : %8.4f\n", 2.0*st->k/(st->b[0]*st->b[0]+st->b[1]*st->b[1]+st->b[2]*st->b[2])-1.0);
   printf("step position  : %8.4f\n", st->l);
   printf("step thickness : %8.4f\n", st->t);
   printf("\n");
   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ density _______________________________________________________________ */
void density(struct sti si, struct stx sx, struct stt *st, double q[3], double n[NS+1], MPI_Comm com)
{


/* __ specie density __ */
n[0] = st->n[0];
n[1] = st->n[1];
n[2] = st->n[2];

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{


/* __ magnetic field components __ */
b[0] = st->b[0];
b[1] = st->b[1];
b[2] = st->b[2];

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
double ww;


/* __ step profile __ */
ww = -0.5*(tanh((q[0]-st->l)/st->t)-1.0);

/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = 0.0;

/* __ proton velocity __ */
v[1][0] = st->v[1]*ww;
v[1][1] = 0.0;
v[1][2] = 0.0;

/* __ alpha velocity __ */
v[2][0] = st->v[2]*ww;
v[2][1] = 0.0;
v[2][2] = 0.0;

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{
double wt;
double wp, wq, we;
double wi0, wi1, wj0, wj1;


/* __ total kinetic pressure __ */
wt = st->k-0.5*(st->b[0]*st->b[0]+st->b[1]*st->b[1]+st->b[2]*st->b[2]);

/* __ kinetic pressure of electrons, protons & alphas __ */
we = st->e*wt;
wp = st->p*wt;
wq = (1.0-st->p-st->e)*wt;

/* __ anisotropy factors __ */
wi0 = 3.0*st->i/(st->i+2.0);
wi1 = 3.0      /(st->i+2.0);
wj0 = 3.0*st->j/(st->j+2.0);
wj1 = 3.0      /(st->j+2.0);

/* __ electron temperatures __ */
ta[0] = we/(st->n[0]);
te[0] = we/(st->n[0]);

/* __  proton temperatures __ */
ta[1] = (st->n[1] == 0.0) ? 0.0 : wp*wi0/st->n[1];
te[1] = (st->n[1] == 0.0) ? 0.0 : wp*wi1/st->n[1];

/* __  alpha temperatures __ */
ta[2] = (st->n[2] == 0.0) ? 0.0 : wq*wj0/st->n[2];
te[2] = (st->n[2] == 0.0) ? 0.0 : wq*wj1/st->n[2];

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


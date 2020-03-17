
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "plasma.h"


/* __ read the initial parameters to set stt structure (uniform plasma) _____ */
void topo(struct sti si, struct stx *sx, struct stt *st, MPI_Comm com)
{
char junk[80];
FILE *fp;


/* __ read the topology parameters __ */
fp = fopen("2gem.txt", "r");
if (fp == NULL) printf("problem in opening file uniform.txt\n");

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->n), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->m), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->l), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->a), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->p), &(sx->irun));

fclose(fp);

/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : 2 gem conf ______\n");
   printf("b              : %8.4f\n", st->b);
   printf("n protons      : %8.4f\n", st->n);
   printf("n alphas       : %8.4f\n", st->m);
   printf("sheet thickness: %8.4f\n", st->l);
   printf("temp(e)/temp(p): %8.4f\n", st->a);
   printf("mag. perturb.  : %8.4f\n", st->p);
   printf("\n");
   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ density _______________________________________________________________ */
void density(struct sti si, struct stx sx, struct stt *st, double q[3], double n[NS+1], MPI_Comm com)
{
double y1, y2;


/* __ location of the 2 current sheets __ */
y1 = (q[1]-0.25*si.l[1])/st->l;
y2 = (q[1]-0.75*si.l[1])/st->l;

/* __ electron density __ */
n[0] = st->n/pow(cosh(y1), 2)+st->n/pow(cosh(y2), 2)+st->m;

/* __ proton density __ */
n[1] = st->n/pow(cosh(y1), 2)+st->n/pow(cosh(y2), 2);

/* __ alpha density __ */
n[2] = st->m;

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{
double x0, y1, y2;
const double w1 = 0.1, w2 = 2.0;
double w3, w4, w5;


/* __ location of the 2 current sheets __ */
y1 = (q[1]-0.25*si.l[1])/st->l;
y2 = (q[1]-0.75*si.l[1])/st->l;

/* __ middle of the box in x direction __ */
x0 = (q[0]-0.5*si.l[0])/st->l;

/* __ 2 needed constants... from seiji __ */
w3 = exp(-(x0*x0+y1*y1)/(w2*w2));
w4 = exp(-(x0*x0+y2*y2)/(w2*w2));
w5 = 2.0*w1/w2;

/* __ magnetic field components __ */
b[0] = st->b*(tanh(y1)-tanh(y2)-1.0+(-w5*y1*w3)+(+w5*y2*w4));
b[1] = st->b*((+w5*x0*w3)+(-w5*x0*w4));
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
double w1, w2;


/* __ needed constants __ */
w1 = st->b/(si.qs[1]*st->n*st->l);
w2 = 1.0/(1.0+st->a);

/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = q[1] < 0.5*si.l[1] ? -w1*w2*st->a : +w1*w2*st->a;

/* __ proton velocity __ */
v[1][0] = 0.0;
v[1][1] = 0.0;
v[1][2] = q[1] < 0.5*si.l[1] ? -w1*w2 : +w1*w2;

/* __ alpha velocity __ */
v[2][0] = 0.0;
v[2][1] = 0.0;
v[2][2] = 0.0;

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{
double wt;


/* __ assymptotic pressure diveded by density __ */
wt = 0.5*(st->b*st->b);

/* __ electron temperatures __ */
ta[0] = wt*st->a/(1+st->a);
te[0] = wt*st->a/(1+st->a);

/* __  proton temperatures __ */
ta[1] = wt/(1+st->a);
te[1] = wt/(1+st->a);

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


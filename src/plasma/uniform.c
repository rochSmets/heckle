
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "defines.h"
#include "structures.h"
#include "misc.h"
#include "plasma.h"


/* __ read the initial parameters to set stt structure (uniform plasma) _____ */
void topo(struct sti si, struct stx *sx, struct stt *st, MPI_Comm com)
{
//int z;
char junk[80];
FILE *fp;


/* __ read the topology parameters __ */
fp = fopen("uniform.txt", "r");
if (fp == NULL) printf("problem in opening file uniform.txt\n");

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b[1]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b[2]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->n), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->m), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->v[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->v[1]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->v[2]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->w[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->w[1]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->w[2]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->te), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->tp), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->ta), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));

fclose(fp);

/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : uniform _________\n");
   printf("b[0]             :%12.6lf\n", st->b[0]);
   printf("  1              :%12.6lf\n", st->b[1]);
   printf("  2              :%12.6lf\n", st->b[2]);
   printf("electron temp.   :%12.6lf\n", st->te);
   printf("n protons        :%12.6lf\n", st->n);
   printf("v[0] protons     :%12.6lf\n", st->v[0]);
   printf("  1              :%12.6lf\n", st->v[1]);
   printf("  2              :%12.6lf\n", st->v[2]);
   printf("proton temp.     :%12.6lf\n", st->tp);
   printf("n alphas         :%12.6lf\n", st->m);
   printf("v[0] alphas      :%12.6lf\n", st->w[0]);
   printf("  1              :%12.6lf\n", st->w[1]);
   printf("  2              :%12.6lf\n", st->w[2]);
   printf("alpha temp.      :%12.6lf\n", st->ta);
   printf("\n");
   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ density _______________________________________________________________ */
void density(struct sti si, struct stx sx, struct stt *st, double q[3], double n[NS+1], MPI_Comm com)
{


/* __ electron density __ */
n[0] = st->n+st->m;

/* __ proton density __ */
n[1] = st->n;

/* __ alpha density __ */
n[2] = st->m;

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


/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = 0.0;

/* __ proton velocity __ */
v[1][0] = st->v[0];
v[1][1] = st->v[1];
v[1][2] = st->v[2];

/* __ alpha velocity __ */
v[2][0] = st->w[0];
v[2][1] = st->w[1];
v[2][2] = st->w[2];

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{


/* __ electron temperatures __ */
ta[0] = st->te;
te[0] = st->te;

/* __  proton temperatures __ */
ta[1] = st->tp;
te[1] = st->tp;

/* __  alpha temperatures __ */
ta[2] = st->ta;
te[2] = st->ta;

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


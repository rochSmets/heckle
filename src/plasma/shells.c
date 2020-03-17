
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "structures.h"
#include "misc.h"
#include "plasma.h"


/* __ read the initial parameters to set stt structure (omega plasma) _____ */
void topo(struct sti si, struct stx *sx, struct stt *st, MPI_Comm com)
{
char junk[80];
FILE *fp;


/* __ read the topology parameters __ */
fp = fopen("shells.txt", "r");
if (fp == NULL) printf("problem in opening file shells.txt\n");

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->fi[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->fi[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->psi[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->psi[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b0[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->b0[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->nb[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->nb[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->n0), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->v0[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->v0[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->te), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->tp[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->tp[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->ta), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->ls[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->ls[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->lw[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->lw[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->lz[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->lz[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->xs[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->xs[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->ys[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->ys[1]), &(sx->irun));

fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->zs[0]), &(sx->irun));
fscandbl(__FILE__, __LINE__, fp, sx->r, &(st->zs[1]), &(sx->irun));

fclose(fp);

/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : 2 shells ________\n");
   printf("fi tilts         : %10.6f    %10.6lf\n", st->fi[0], st->fi[1]);
   printf("psi tilts        : %10.6f    %10.6lf\n", st->psi[0], st->psi[1]);
   printf("b                : %10.6f    %10.6f\n", st->b0[0], st->b0[1]);
   printf("n protons        : %10.6f    %10.6lf\n", st->nb[0], st->nb[1]);
   printf("n alphas         : %10.6f\n", st->n0);
   printf("v protons        : %10.6f    %10.6lf\n", st->v0[0], st->v0[1]);
   printf("electrons temp   : %10.6f\n", st->te);
   printf("protons temp     : %10.6f    %10.6lf\n", st->tp[0], st->tp[1]);
   printf("alpha temp       : %10.6f\n", st->ta);
   printf("\"beta\"           : %10.6f    %10.6lf\n", 2.0*(st->nb[0]*st->te
                                                          +st->n0*st->te
                                                          +st->nb[0]*st->tp[0]
                                                          +st->n0*st->ta)/pow(st->b0[0], 2),
                                                      2.0*(st->nb[1]*st->te
                                                          +st->n0*st->te
                                                          +st->nb[1]*st->tp[1]
                                                          +st->n0*st->ta)/pow(st->b0[1], 2));
   printf("shells radius    : %10.6f    %10.6lf\n", st->ls[0], st->ls[1]);
   printf("shells 1/2 width : %10.6f    %10.6lf\n", st->lw[0], st->lw[1]);
   printf("shells z extens. : %10.6f    %10.6lf\n", st->lz[0], st->lz[1]);
   printf("shells x loc.    : %10.6f    %10.6lf\n", st->xs[0], st->xs[1]);
   printf("shells y loc.    : %10.6f    %10.6lf\n", st->ys[0], st->ys[1]);
   printf("shells z loc.    : %10.6f    %10.6lf\n", st->zs[0], st->zs[1]);
   printf("\n");
   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ polynom used for the interpolation __ */
double polynom(double x)
{
double w, y;


w = -6.0*fabs(x)*fabs(x)*fabs(x)*fabs(x)*fabs(x)
    +15.0*x*x*x*x
    -10.0*fabs(x)*fabs(x)*fabs(x)
    +1;

y = (fabs(x) <= 1.0) ? w : 0.0;

return y;

}


/* __ density _______________________________________________________________ */
void density(struct sti si, struct stx sx, struct stt *st, double q[3], double n[NS+1], MPI_Comm com)
{
double wf, wp;
double wx, wy;
double tx, ty;
double wr;
double wa[2];
int h;


/* __ loop on the 2 shells __ */
for (h = 0; h < 2; h++)
    {
    /* __ set fi & psi angles __ */
    wf = -st->fi[h]*PI/180.0;
    wp = +st->psi[h]*PI/180.0;

    /* __ backward rotation with fi angle __ */
    wx = (q[0]-st->xs[h])*cos(wf)-(q[1]-st->ys[h])*sin(wf);
    wy = (q[0]-st->xs[h])*sin(wf)+(q[1]-st->ys[h])*cos(wf);

    /* __ rotation with psi angle __ */
    tx = wx;
    ty = wy/cos(wp);

    /* __ modulus of arg __ */
    wr = sqrt(pow(tx, 2)+pow(ty, 2));

    /* __ set arg of the polynom __ */
    wa[h] = wr/st->ls[h];
    }

/* __ electron density __ */
n[0] = st->nb[0]*polynom(wa[0])
      +st->nb[1]*polynom(wa[1])
      +st->n0;

/* __ proton density __ */
n[1] = st->nb[0]*polynom(wa[0])
      +st->nb[1]*polynom(wa[1]);

/* __ alpha density __ */
n[2] = st->n0;

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{
double wf, wp;
double wx, wy;
double tx, ty;
double wr;
double wa[2];
double ux[2], uy[2], uz[2];
int h;


/* __ loop on the 2 shells __ */
for (h = 0; h < 2; h++)
    {
    /* __ set fi & psi angles __ */
    wf = -st->fi[h]*PI/180;
    wp = +st->psi[h]*PI/180;

    /* __ backward rotation with fi angle __ */
    wx = (q[0]-st->xs[h])*cos(wf)-(q[1]-st->ys[h])*sin(wf);
    wy = (q[0]-st->xs[h])*sin(wf)+(q[1]-st->ys[h])*cos(wf);

    /* __ rotation with psi angle __ */
    tx = wx;
    ty = wy/cos(wp);

    /* __ modulus of arg __ */
    wr = sqrt(pow(tx, 2)+pow(ty, 2));

    /* __ set arg of the polynom __ */
    wa[h] = (wr-st->ls[h])/st->lw[h];

    /* __ unit vector along magnetic field lines __ */
    ux[h] = (wr == 0) ? 0.0 : +(+ty/wr)*cos(wf)+(-tx/wr)*sin(wf)*cos(wp);
    uy[h] = (wr == 0) ? 0.0 : -(+ty/wr)*sin(wf)+(-tx/wr)*cos(wf)*cos(wp);
    uz[h] = (wr == 0) ? 0.0 : +(-tx/wr)*sin(wp);
    }

/* __ magnetic field __ */
b[0] = st->b0[0]*polynom(wa[0])*ux[0]
      +st->b0[1]*polynom(wa[1])*ux[1];
b[1] = st->b0[0]*polynom(wa[0])*uy[0]
      +st->b0[1]*polynom(wa[1])*uy[1];
b[2] = st->b0[0]*polynom(wa[0])*uz[0]
      +st->b0[1]*polynom(wa[1])*uz[1];
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
double wf, wp;
double wx, wy;
double tx, ty;
double wr;
double wa[2];
double ux[2], uy[2], uz[2];
int h;


/* __ loop on the 2 shells __ */
for (h = 0; h < 2; h++)
    {
    /* __ set fi & psi angles __ */
    wf = -st->fi[h]*PI/180.0;
    wp = +st->psi[h]*PI/180.0;

    /* __ backward rotation with fi angle __ */
    wx = +(q[0]-st->xs[h])*cos(wf)-(q[1]-st->ys[h])*sin(wf);
    wy = +(q[0]-st->xs[h])*sin(wf)+(q[1]-st->ys[h])*cos(wf);

    /* __ rotation with psi angle __ */
    tx = wx;
    ty = wy/cos(wp);

    /* __ modulus of arg __ */
    wr = sqrt(pow(tx, 2)+pow(ty, 2));

    /* __ set arg of the polynom __ */
    wa[h] = (2.0*wr-st->ls[h])/st->ls[h];

    /* __ unit vector in radial direction __ */
    ux[h] = (wr == 0) ? 0.0 : +(+tx/wr)*cos(wf)+(+ty/wr)*sin(wf)*cos(wp);
    uy[h] = (wr == 0) ? 0.0 : -(+tx/wr)*sin(wf)+(+ty/wr)*cos(wf)*cos(wp);
    uz[h] = (wr == 0) ? 0.0 : +(+ty/wr)*sin(wp);
    }

/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = 0.0;

/* __ proton velocity __ */
v[1][0] = st->v0[0]*polynom(wa[0])*ux[0]
         +st->v0[1]*polynom(wa[1])*ux[1];
v[1][1] = st->v0[0]*polynom(wa[0])*uy[0]
         +st->v0[1]*polynom(wa[1])*uy[1];
v[1][2] = st->v0[0]*polynom(wa[0])*uz[0]
         +st->v0[1]*polynom(wa[1])*uz[1];

/* __ alpha velocity __ */
v[2][0] = 0.0;
v[2][1] = 0.0;
v[2][2] = 0.0;

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{
double wf, wp;
double wx, wy, wa[2];
double n0, n1;
int h;


/* __ loop on the 2 shells __ */
for (h = 0; h < 2; h++)
    {
    /* __ set fi & psi angles __ */
    wf = -st->fi[h]*PI/180.0;
    wp = +st->psi[h]*PI/180.0;

    /* __ set rotation with fi angle __ */
    wx = (q[0]-st->xs[h])*cos(wf)+(q[1]-st->ys[h])*sin(wf);
    wy = (q[0]-st->xs[h])*sin(wf)-(q[1]-st->ys[h])*cos(wf);

    /* __ set arg of polynom & rotation with psi angle __ */
    wa[h] = sqrt(pow(wx, 2)+pow(wy/cos(wp), 2))/st->ls[h];
    }

/* __ proton density __ */
n0 = st->nb[0]*polynom(wa[0]);
n1 = st->nb[1]*polynom(wa[1]);

/* __ electron temperatures __ */
ta[0] = st->te;
te[0] = st->te;

/* __  proton temperatures __ */
ta[1] = (n0+n1 != 0.0) ? (n0*st->tp[0]+n1*st->tp[1])/(n0+n1) : 0.0;
te[1] = (n0+n1 != 0.0) ? (n0*st->tp[0]+n1*st->tp[1])/(n0+n1) : 0.0;

/* __  alpha temperatures __ */
ta[2] = st->ta;
te[2] = st->ta;

}


/* __ calculation of the normalization coefficients _________________________ */
void normal(struct sti *si, struct stx sx, struct stt *st, MPI_Comm com)
{
double wq[3], wn[3];
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

                /* __ local density for macro-particles s __ */
                density(*si, sx, st, wq, wn, com);

                /* __ integrated local "mass" of macro-particles s __ */
                wt += wn[s]*wd;
                }
            }
        }

    /* __ weight of part. "s" __ */
    si->ws[s] = (si->ns[s] == 0) ? EPS6 : wt/(si->ns[s]*wd);
    }

}


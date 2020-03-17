
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "structures.h"
#include "plasma.h"


/* __ read the initial parameters to set stt structure (omega plasma) _____ */
void topo(struct sti si, struct stx sx, struct stt *st, MPI_Comm com)
{
int z;
char junk[80];
FILE *fp;


/* __ read the topology parameters __ */
fp = fopen("yaourt.txt", "r");
if (fp == NULL) printf("problem in opening file yaourt.txt\n");

z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->b0[0]),  &(st->b0[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->fi), &(st->psi)); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->nb[0]), &(st->nb[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->n0)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->v0[0]), &(st->v0[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->k[0]), &(st->k[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->p[0]), &(st->p[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->e[0]), &(st->e[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->lb[0]), &(st->lb[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", &(st->lr[0]), &(st->lr[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);

fclose(fp);

/* __ display the topology parameters __ */
if (sx.r == 0)
   {
   printf("________________ magnetic topology : yaourt __________\n");
   printf("b                : %10.6f    %10.6f\n", st->b0[0], st->b0[1]);
   printf("targets tilts    : %10.6f    %10.6lf\n", st->fi, st->psi);
   printf("n protons        : %10.6f    %10.6lf\n", st->nb[0], st->nb[1]);
   printf("n alphas         : %10.6f\n", st->n0);
   printf("v protons        : %10.6f    %10.6lf\n", st->v0[0], st->v0[1]);
   printf("l bubble         : %10.6f    %10.6lf\n", st->lb[0], st->lb[1]);
   printf("l ribbon         : %10.6f    %10.6lf\n", st->lr[0], st->lr[1]);
   printf("beta protons     : %10.6f    %10.6lf\n",
          st->p[0]*(2.0*st->k[0]/(st->b0[0]*st->b0[0])-1.0),
          st->p[1]*(2.0*st->k[1]/(st->b0[1]*st->b0[1])-1.0));
   printf("beta alphas      : %10.6f    %10.6lf\n",
          (1.0-st->p[0]-st->e[0])*(2.0*st->k[0]/(st->b0[0]*st->b0[0])-1.0),
          (1.0-st->p[1]-st->e[1])*(2.0*st->k[1]/(st->b0[1]*st->b0[1])-1.0));
   printf("beta electrons   : %10.6f    %10.6lf\n",
          st->e[0]*(2.0*st->k[0]/(st->b0[0]*st->b0[0])-1.0),
          st->e[1]*(2.0*st->k[1]/(st->b0[1]*st->b0[1])-1.0));
   printf("beta total       : %10.6f    %10.6lf\n",
          2.0*st->k[0]/(st->b0[0]*st->b0[0])-1.0,
          2.0*st->k[1]/(st->b0[1]*st->b0[1])-1.0);
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
double wx, wy, wf, wp;
double xf, yf;
double xp, yp;
double r0, r1, r2;
double a0, a1, a2;


/* __ fi angle __ */
if (st->fi != 0.0 && st->psi == 0.0)
   {
   /* __ bubble 0 __ */
   if (q[1] < 0.25*si.l[1])
      {
      wy = 0.0*si.l[1];
      wf = +st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r0 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-0.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a0 = r0/st->lb[0];

      /* __ electron density __ */
      n[0] = st->n0+st->nb[0]*polynom(a0);

      /* __ proton density __ */
      n[1] = st->nb[0]*polynom(a0);

      /* __ alpha density __ */
      n[2] = st->n0;
      }

   /* __ bubble 1 __ */
   if (q[1] >= 0.25*si.l[1] && q[1] <= 0.75*si.l[1])
      {
      wy = 0.5*si.l[1];
      wf = -st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r1 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-0.5*si.l[1], 2));

      /* __ arg of the polynom __ */
      a1 = r1/st->lb[1];

      /* __ electron density __ */
      n[0] = st->n0+st->nb[1]*polynom(a1);

      /* __ proton density __ */
      n[1] = st->nb[1]*polynom(a1);

      /* __ alpha density __ */
      n[2] = st->n0;
      }

   /* __ bubble 2 __ */
   if (q[1] > 0.75*si.l[1])
      {
      wy = 1.0*si.l[1];
      wf = +st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r2 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-1.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a2 = r2/st->lb[0];

      /* __ electron density __ */
      n[0] = st->n0+st->nb[0]*polynom(a2);

      /* __ proton density __ */
      n[1] = st->nb[0]*polynom(a2);

      /* __ alpha density __ */
      n[2] = st->n0;
      }
   }

/* __ psi angle __ */
if (st->psi != 0.0 && st->fi == 0.0)
   {
   /* __ bubble 0 __ */
   if (q[1] < 0.25*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = +st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r0 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-0.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a0 = r0/st->lb[0];

      /* __ electron density __ */
      n[0] = st->n0+st->nb[0]*polynom(a0);

      /* __ proton density __ */
      n[1] = st->nb[0]*polynom(a0);

      /* __ alpha density __ */
      n[2] = st->n0;
      }

   /* __ bubble 1 __ */
   if (q[1] >= 0.25*si.l[1] && q[1] <= 0.75*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = -st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r1 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-0.5*si.l[1], 2));

      /* __ arg of the polynom __ */
      a1 = r1/st->lb[1];

      /* __ electron density __ */
      n[0] = st->n0+st->nb[1]*polynom(a1);

      /* __ proton density __ */
      n[1] = st->nb[1]*polynom(a1);

      /* __ alpha density __ */
      n[2] = st->n0;
      }

   /* __ bubble 2 __ */
   if (q[1] > 0.75*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = +st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r2 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-1.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a2 = r2/st->lb[0];

      /* __ electron density __ */
      n[0] = st->n0+st->nb[0]*polynom(a2);

      /* __ proton density __ */
      n[1] = st->nb[0]*polynom(a2);

      /* __ alpha density __ */
      n[2] = st->n0;
      }
   }

/* __ coplanar targets __ */
if (st->fi == 0.0 && st->psi == 0.0)
   {
   /* __ radius in the xy plan __ */
   r0 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-0.0*si.l[1], 2));
   r1 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-0.5*si.l[1], 2));
   r2 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-1.0*si.l[1], 2));

   /* __ arg of the polynom __ */
   a0 = r0/st->lb[0];
   a1 = r1/st->lb[1];
   a2 = r2/st->lb[0];

   /* __ electron density __ */
   n[0] = st->n0+st->nb[0]*polynom(a0)
                +st->nb[1]*polynom(a1)
                +st->nb[0]*polynom(a2);

   /* __ proton density __ */
   n[1] = st->nb[0]*polynom(a0)
         +st->nb[1]*polynom(a1)
         +st->nb[0]*polynom(a2);

   /* __ alpha density __ */
   n[2] = st->n0;
   }

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{
double wx, wy, wf, wp;
double xf, yf;
double xp, yp;
double wbx, wby;
double r0, r1, r2;
double a0, a1, a2;
double b0, b1, b2;
double xx, y0, y1, y2;
double u0[2], u1[2], u2[2];


/* __ fi angle __ */
if (st->fi != 0.0 && st->psi == 0.0)
   {
   /* __ bubble 0 __ */
   if (q[1] < 0.25*si.l[1])
      {
      wy = 0.0*si.l[1];
      wf = +st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r0 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-0.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a0 = (st->lb[0]-r0)/st->lr[0];

      /* __ modulus of the magnetic field __ */
      b0 = st->b0[0]*polynom(a0);

      /* __ particle coordinate __ */
      xx = xf-0.5*si.l[0];
      y0 = yf-0.0*si.l[1];

      /* __ vector tangential to the ribbon __ */
      u0[0] = (r0 == 0.0) ? 0.0 : +y0/r0;
      u0[1] = (r0 == 0.0) ? 0.0 : -xx/r0;

      /* __ magnetic field components (in the backward plane) __ */
      wbx = b0*u0[0];
      wby = b0*u0[1];

      /* __ set magnetic field components __ */
      b[0] = wbx;
      b[1] = wby*cos(wf);
      b[2] = wby*sin(wf);
      }

   /* __ bubble 1 __ */
   if (q[1] >= 0.25*si.l[1] && q[1] <= 0.75*si.l[1])
      {
      wy = 0.5*si.l[1];
      wf = -st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r1 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-0.5*si.l[1], 2));

      /* __ arg of the polynom __ */
      a1 = (st->lb[1]-r1)/st->lr[1];

      /* __ modulus of the magnetic field __ */
      b1 = st->b0[1]*polynom(a1);

      /* __ particle coordinate __ */
      xx = xf-0.5*si.l[0];
      y1 = yf-0.5*si.l[1];

      /* __ vector tangential to the ribbon __ */
      u1[0] = (r1 == 0.0) ? 0.0 : +y1/r1;
      u1[1] = (r1 == 0.0) ? 0.0 : -xx/r1;

      /* __ magnetic field components (in the backward plane) __ */
      wbx = b1*u1[0];
      wby = b1*u1[1];

      /* __ set magnetic field components __ */
      b[0] = wbx;
      b[1] = wby*cos(wf);
      b[2] = wby*sin(wf);
      }

   /* __ bubble 2 __ */
   if (q[1] > 0.75*si.l[1])
      {
      wy = 1.0*si.l[1];
      wf = +st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r2 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-1.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a2 = (st->lb[0]-r2)/st->lr[0];

      /* __ modulus of the magnetic field __ */
      b2 = st->b0[0]*polynom(a2);

      /* __ particle coordinate __ */
      xx = xf-0.5*si.l[0];
      y2 = yf-1.0*si.l[1];

      /* __ vector tangential to the ribbon __ */
      u2[0] = (r2 == 0.0) ? 0.0 : +y2/r2;
      u2[1] = (r2 == 0.0) ? 0.0 : -xx/r2;

      /* __ magnetic field components (in the backward plane) __ */
      wbx = b2*u2[0];
      wby = b2*u2[1];

      /* __ set magnetic field components __ */
      b[0] = wbx;
      b[1] = wby*cos(wf);
      b[2] = wby*sin(wf);
      }
   }

/* __ psi angle __ */
if (st->psi != 0.0 && st->fi == 0.0)
   {
   /* __ bubble 0 __ */
   if (q[1] < 0.25*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = +st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r0 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-0.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a0 = (st->lb[0]-r0)/st->lr[0];

      /* __ modulus of the magnetic field __ */
      b0 = st->b0[0]*polynom(a0);

      /* __ particle coordinate __ */
      xx = xp-0.5*si.l[0];
      y0 = yp-0.0*si.l[1];

      /* __ vector tangential to the ribbon __ */
      u0[0] = (r0 == 0.0) ? 0.0 : +y0/r0;
      u0[1] = (r0 == 0.0) ? 0.0 : -xx/r0;

      /* __ magnetic field components (in the backward plane) __ */
      wbx = b0*u0[0];
      wby = b0*u0[1];

      /* __ set magnetic field components __ */
      b[0] = wbx*cos(wp);
      b[1] = wby;
      b[2] = wbx*sin(wp);
      }

   /* __ bubble 1 __ */
   if (q[1] >= 0.25*si.l[1] && q[1] <= 0.75*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = -st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r1 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-0.5*si.l[1], 2));

      /* __ arg of the polynom __ */
      a1 = (st->lb[1]-r1)/st->lr[1];

      /* __ modulus of the magnetic field __ */
      b1 = st->b0[1]*polynom(a1);

      /* __ particle coordinate __ */
      xx = xp-0.5*si.l[0];
      y1 = yp-0.5*si.l[1];

      /* __ vector tangential to the ribbon __ */
      u1[0] = (r1 == 0.0) ? 0.0 : +y1/r1;
      u1[1] = (r1 == 0.0) ? 0.0 : -xx/r1;

      /* __ magnetic field components (in the backward plane) __ */
      wbx = b1*u1[0];
      wby = b1*u1[1];

      /* __ set magnetic field components __ */
      b[0] = wbx*cos(wp);
      b[1] = wby;
      b[2] = wbx*sin(wp);
      }

   /* __ bubble 2 __ */
   if (q[1] > 0.75*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = +st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r2 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-1.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a2 = (st->lb[0]-r2)/st->lr[0];

      /* __ modulus of the magnetic field __ */
      b2 = st->b0[0]*polynom(a2);

      /* __ particle coordinate __ */
      xx = xp-0.5*si.l[0];
      y2 = yp-1.0*si.l[1];

      /* __ vector tangential to the ribbon __ */
      u2[0] = (r2 == 0.0) ? 0.0 : +y2/r2;
      u2[1] = (r2 == 0.0) ? 0.0 : -xx/r2;

      /* __ magnetic field components (in the backward plane) __ */
      wbx = b2*u2[0];
      wby = b2*u2[1];

      /* __ set magnetic field components __ */
      b[0] = wbx*cos(wp);
      b[1] = wby;
      b[2] = wbx*sin(wp);
      }
   }

/* __ coplanar targets __ */
if (st->fi == 0.0 && st->psi == 0.0)
   {
   /* __ radius in the xy plan __ */
   r0 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-0.0*si.l[1], 2));
   r1 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-0.5*si.l[1], 2));
   r2 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-1.0*si.l[1], 2));

   /* __ arg of the polynom __ */
   a0 = (st->lb[0]-r0)/st->lr[0];
   a1 = (st->lb[1]-r1)/st->lr[1];
   a2 = (st->lb[0]-r2)/st->lr[0];

   /* __ modulus of the magnetic field __ */
   b0 = st->b0[0]*polynom(a0);
   b1 = st->b0[1]*polynom(a1);
   b2 = st->b0[0]*polynom(a2);

   /* __ particle coordinate __ */
   xx = q[0]-0.5*si.l[0];
   y0 = q[1]-0.0*si.l[1];
   y1 = q[1]-0.5*si.l[1];
   y2 = q[1]-1.0*si.l[1];

   /* __ vector tangential to the ribbon __ */
   u0[0] = (r0 == 0.0) ? 0.0 : +y0/r0;
   u0[1] = (r0 == 0.0) ? 0.0 : -xx/r0;
   u1[0] = (r1 == 0.0) ? 0.0 : +y1/r1;
   u1[1] = (r1 == 0.0) ? 0.0 : -xx/r1;
   u2[0] = (r2 == 0.0) ? 0.0 : +y2/r2;
   u2[1] = (r2 == 0.0) ? 0.0 : -xx/r2;

   /* __ magnetic field components (in the backward plane) __ */
   b[0] = b0*u0[0]+b1*u1[0]+b2*u2[0];
   b[1] = b0*u0[1]+b1*u1[1]+b2*u2[1];
   b[2] = 0.0;
   }

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
double wx, wy, wf, wp;
double xf, yf;
double xp, yp;
double wvx, wvy;
double r0, r1, r2;
double a0, a1, a2;
double v0, v1, v2;
double xx, y0, y1, y2;
double u0[2], u1[2], u2[2];


/* __ fi angle __ */
if (st->fi != 0.0 && st->psi == 0.0)
   {
   /* __ bubble 0 __ */
   if (q[1] < 0.25*si.l[1])
      {
      wy = 0.0*si.l[1];
      wf = -st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r0 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-0.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a0 = (2.0*r0-st->lb[0])/st->lb[0];

      /* __ modulus of the proton velocity __ */
      v0 = st->v0[0]*polynom(a0);

      /* __ particle coordinate __ */
      xx = xf-0.5*si.l[0];
      y0 = yf-0.0*si.l[1];

      /* __ vector radial within the ribbon __ */
      u0[0] = (r0 == 0.0) ? 0.0 : +xx/r0;
      u0[1] = (r0 == 0.0) ? 0.0 : +y0/r0;

      /* __ magnetic field components (in the backward plane) __ */
      wvx = v0*u0[0];
      wvy = v0*u0[1];

      /* __ proton velocity __ */
      v[1][0] = wvx;
      v[1][1] = wvy*cos(wf);
      v[1][2] = wvy*sin(wf);
      }

   /* __ bubble 1 __ */
   if (q[1] >= 0.25*si.l[1] && q[1] <= 0.75*si.l[1])
      {
      wy = 0.5*si.l[1];
      wf = -st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r1 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-0.5*si.l[1], 2));

      /* __ arg of the polynom __ */
      a1 = (2.0*r1-st->lb[1])/st->lb[1];

      /* __ modulus of the proton velocity __ */
      v1 = st->v0[1]*polynom(a1);

      /* __ particle coordinate __ */
      xx = xf-0.5*si.l[0];
      y1 = yf-0.5*si.l[1];

      /* __ vector radial within the ribbon __ */
      u1[0] = (r1 == 0.0) ? 0.0 : +xx/r1;
      u1[1] = (r1 == 0.0) ? 0.0 : +y1/r1;

      /* __ magnetic field components (in the backward plane) __ */
      wvx = v1*u1[0];
      wvy = v1*u1[1];

      /* __ proton velocity __ */
      v[1][0] = wvx;
      v[1][1] = wvy*cos(wf);
      v[1][2] = wvy*sin(wf);
      }

   /* __ bubble 2 __ */
   if (q[1] > 0.75*si.l[1])
      {
      wy = 1.0*si.l[1];
      wf = -st->fi*PI/180.0;

      /* __ set position with backward rotation __ */
      xf = q[0];
      yf = (q[1]-wy)/cos(wf)+wy;

      /* __ radius in the xy plan __ */
      r2 = sqrt(pow(xf-0.5*si.l[0], 2)+pow(yf-1.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a2 = (2.0*r2-st->lb[0])/st->lb[0];

      /* __ modulus of the proton velocity __ */
      v2 = st->v0[0]*polynom(a2);

      /* __ particle coordinate __ */
      xx = xf-0.5*si.l[0];
      y2 = yf-1.0*si.l[1];

      /* __ vector radial within the ribbon __ */
      u2[0] = (r2 == 0.0) ? 0.0 : +xx/r2;
      u2[1] = (r2 == 0.0) ? 0.0 : +y2/r2;

      /* __ magnetic field components (in the backward plane) __ */
      wvx = v2*u2[0];
      wvy = v2*u2[1];

      /* __ proton velocity __ */
      v[1][0] = wvx;
      v[1][1] = wvy*cos(wf);
      v[1][2] = wvy*sin(wf);
      }
   }

/* __ psi angle __ */
if (st->psi != 0.0 && st->fi == 0.0)
   {
   /* __ bubble 0 __ */
   if (q[1] < 0.25*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = -st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r0 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-0.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a0 = (2.0*r0-st->lb[0])/st->lb[0];

      /* __ modulus of the proton velocity __ */
      v0 = st->v0[0]*polynom(a0);

      /* __ particle coordinate __ */
      xx = xp-0.5*si.l[0];
      y0 = yp-0.0*si.l[1];

      /* __ vector radial within the ribbon __ */
      u0[0] = (r0 == 0.0) ? 0.0 : +xx/r0;
      u0[1] = (r0 == 0.0) ? 0.0 : +y0/r0;

      /* __ magnetic field components (in the backward plane) __ */
      wvx = v0*u0[0];
      wvy = v0*u0[1];

      /* __ proton velocity __ */
      v[1][0] = wvx*cos(wp);
      v[1][1] = wvy;
      v[1][2] = wvx*sin(wp);
      }

   /* __ bubble 1 __ */
   if (q[1] >= 0.25*si.l[1] && q[1] <= 0.75*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = -st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r1 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-0.5*si.l[1], 2));

      /* __ arg of the polynom __ */
      a1 = (2.0*r1-st->lb[1])/st->lb[1];

      /* __ modulus of the proton velocity __ */
      v1 = st->v0[1]*polynom(a1);

      /* __ particle coordinate __ */
      xx = xp-0.5*si.l[0];
      y1 = yp-0.5*si.l[1];

      /* __ vector radial within the ribbon __ */
      u1[0] = (r1 == 0.0) ? 0.0 : +xx/r1;
      u1[1] = (r1 == 0.0) ? 0.0 : +y1/r1;

      /* __ magnetic field components (in the backward plane) __ */
      wvx = v1*u1[0];
      wvy = v1*u1[1];

      /* __ proton velocity __ */
      v[1][0] = wvx*cos(wp);
      v[1][1] = wvy;
      v[1][2] = wvx*sin(wp);
      }

   /* __ bubble 2 __ */
   if (q[1] > 0.75*si.l[1])
      {
      wx = 0.5*si.l[0];
      wp = -st->psi*PI/180.0;

      /* __ set position with backward rotation __ */
      xp = (q[0]-wx)/cos(wp)+wx;
      yp = q[1];

      /* __ radius in the xy plan __ */
      r2 = sqrt(pow(xp-0.5*si.l[0], 2)+pow(yp-1.0*si.l[1], 2));

      /* __ arg of the polynom __ */
      a2 = (2.0*r2-st->lb[0])/st->lb[0];

      /* __ modulus of the proton velocity __ */
      v2 = st->v0[0]*polynom(a2);

      /* __ particle coordinate __ */
      xx = xp-0.5*si.l[0];
      y2 = yp-1.0*si.l[1];

      /* __ vector radial within the ribbon __ */
      u2[0] = (r2 == 0.0) ? 0.0 : +xx/r2;
      u2[1] = (r2 == 0.0) ? 0.0 : +y2/r2;

      /* __ magnetic field components (in the backward plane) __ */
      wvx = v2*u2[0];
      wvy = v2*u2[1];

      /* __ proton velocity __ */
      v[1][0] = wvx*cos(wp);
      v[1][1] = wvy;
      v[1][2] = wvx*sin(wp);
      }
   }

/* __ coplanar targets __ */
if (st->fi == 0.0 && st->psi == 0.0)
   {
   /* __ radius in the xy plan __ */
   r0 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-0.0*si.l[1], 2));
   r1 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-0.5*si.l[1], 2));
   r2 = sqrt(pow(q[0]-0.5*si.l[0], 2)+pow(q[1]-1.0*si.l[1], 2));

   /* __ arg of the polynom __ */
   a0 = (2.0*r0-st->lb[0])/st->lb[0];
   a1 = (2.0*r1-st->lb[1])/st->lb[1];
   a2 = (2.0*r2-st->lb[0])/st->lb[0];

   /* __ modulus of the proton velocity __ */
   v0 = st->v0[0]*polynom(a0);
   v1 = st->v0[1]*polynom(a1);
   v2 = st->v0[0]*polynom(a2);

   /* __ particle coordinate __ */
   xx = q[0]-0.5*si.l[0];
   y0 = q[1]-0.0*si.l[1];
   y1 = q[1]-0.5*si.l[1];
   y2 = q[1]-1.0*si.l[1];

   /* __ vector radial within the ribbon __ */
   u0[0] = (r0 == 0.0) ? 0.0 : +xx/r0;
   u0[1] = (r0 == 0.0) ? 0.0 : +y0/r0;
   u1[0] = (r1 == 0.0) ? 0.0 : +xx/r1;
   u1[1] = (r1 == 0.0) ? 0.0 : +y1/r1;
   u2[0] = (r2 == 0.0) ? 0.0 : +xx/r2;
   u2[1] = (r2 == 0.0) ? 0.0 : +y2/r2;
// u0[0] = (r0 == 0.0) ? 0.0 : +y0/r0;
// u0[1] = (r0 == 0.0) ? 0.0 : +xx/r0;
// u1[0] = (r1 == 0.0) ? 0.0 : +y1/r1;
// u1[1] = (r1 == 0.0) ? 0.0 : +xx/r1;
// u2[0] = (r2 == 0.0) ? 0.0 : +y2/r2;
// u2[1] = (r2 == 0.0) ? 0.0 : +xx/r2;

   /* __ proton velocity __ */
   v[1][0] = v0*u0[0]+v1*u1[0]+v2*u2[0];
   v[1][1] = v0*u0[1]+v1*u1[1]+v2*u2[1];
   v[1][2] = 0.0;
   }

/* __ electron velocity __ */
v[0][0] = 0.0;
v[0][1] = 0.0;
v[0][2] = 0.0;

/* __ alpha velocity __ */
v[2][0] = 0.0;
v[2][1] = 0.0;
v[2][2] = 0.0;

}


/* __ temperature & thermal velocity ________________________________________ */
void kinetic(struct sti si, struct stx sx, struct stt *st, double q[3], double ta[NS+1], double te[NS+1], MPI_Comm com)
{
double wt[2];
double wp[2], wq[2], we[2];


/* __ total kinetic pressure __ */
wt[0] = st->k[0]-0.5*(st->b0[0]*st->b0[0]);
wt[1] = st->k[1]-0.5*(st->b0[1]*st->b0[1]);

/* __ electron pressure __ */
we[0] = st->e[0]*wt[0];
we[1] = st->e[1]*wt[1];

/* __ proton pressure __ */
wp[0] = st->p[0]*wt[0];
wp[1] = st->p[1]*wt[1];

/* __ alpha pressure __ */
wq[0] = (1.0-st->p[0]-st->e[0])*wt[0];
wq[1] = (1.0-st->p[1]-st->e[1])*wt[1];

/* __ in the first bubble __ */
if (q[1] <= 0.25*si.l[1] || q[1] >= 0.75*si.l[1])
   {
   /* __ electron temperatures __ */
   ta[0] = we[0]/(st->nb[0]+st->n0);
   te[0] = we[0]/(st->nb[0]+st->n0);

   /* __  proton temperatures __ */
   ta[1] = wp[0]/(st->nb[0]);
   te[1] = wp[0]/(st->nb[0]);

   /* __  alpha temperatures __ */
   ta[2] = wq[0]/(st->n0);
   te[2] = wq[0]/(st->n0);
   }

/* __ in the second bubble __ */
else
   {
   /* __ electron temperatures __ */
   ta[0] = we[1]/(st->nb[1]+st->n0);
   te[0] = we[1]/(st->nb[1]+st->n0);

   /* __  proton temperatures __ */
   ta[1] = wp[1]/(st->nb[1]);
   te[1] = wp[1]/(st->nb[1]);

   /* __  alpha temperatures __ */
   ta[2] = wq[1]/(st->n0);
   te[2] = wq[1]/(st->n0);
   }

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


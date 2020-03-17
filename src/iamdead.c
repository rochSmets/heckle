
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>

#include "structures.h"
#include <particle.h>
#include "defines.h"


/* __ report problem if part is not in the subdomain ________________________ */
void deadpart(struct sti si,
                 struct stx *sx,
                 struct st0 *s0,
                 struct st1 *s1,
                 struct st2 *s2,
                 struct stp *sp[NS+1],
                 int s,
                 int m,
                 int it,
                 char *file,
                 int line)
{
double xw, yw, zw;
int i,j, k;
int ijk0, ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
int l0, l1, l2, m0, m1, m2, n0, n1, n2;


/* __ # of grid points __ */
l0 = sx->n[0];
l1 = sx->n[1];
l2 = sx->n[2];
m0 = sx->n[0]+1;
m1 = sx->n[1]+1;
m2 = sx->n[2]+1;
n0 = sx->n[0]+2;
n1 = sx->n[1]+2;
n2 = sx->n[2]+2;

printf("________________ iamdead : a part seems fucked-up ! ___\n");
printf("\n");
printf("________________ time step : %10d ______________\n", it);
printf("\n");
printf("________________ part on node %6d _________________\n", sx->r);
printf("\n");
printf("________________ subdomain size ______________________\n");
printf("in x direction   :%12.6lf  %12.6lf\n", sx->i0[0]*si.dl[0],
                                               sx->i1[0]*si.dl[0]);
printf("in y direction   :%12.6lf  %12.6lf\n", sx->i0[1]*si.dl[1],
                                               sx->i1[1]*si.dl[1]);
printf("in z direction   :%12.6lf  %12.6lf\n", sx->i0[2]*si.dl[2],
                                               sx->i1[2]*si.dl[2]);
printf("\n");
printf("________________ part identity _______________________\n");
printf("part specie      :%12d\n", s);
printf("part index       :%12d\n", m);
printf("part id #        :%12d\n", sp[s][m].i);
printf("\n");
printf("________________ max # of part in subdomain __________\n");
//printf("nm               :%12" PRId64 "\n", sx->nm[s]);
printf("nm               :%12zu\n", sx->nm[s]);
printf("\n");
printf("________________ part position _______________________\n");
printf("r[0]             :%12.6lf\n", sp[s][m].r[0]);
printf("r[1]             :%12.6lf\n", sp[s][m].r[1]);
printf("r[2]             :%12.6lf\n", sp[s][m].r[2]);
printf("s[0]             :%12.6lf\n", sp[s][m].s[0]);
printf("s[1]             :%12.6lf\n", sp[s][m].s[1]);
printf("s[2]             :%12.6lf\n", sp[s][m].s[2]);
printf("\n");
printf("________________ part velocity _______________________\n");
printf("v[0]             :%12.6lf\n", sp[s][m].v[0]);
printf("v[1]             :%12.6lf\n", sp[s][m].v[1]);
printf("v[2]             :%12.6lf\n", sp[s][m].v[2]);
printf("w[0]             :%12.6lf\n", sp[s][m].w[0]);
printf("w[1]             :%12.6lf\n", sp[s][m].w[1]);
printf("w[2]             :%12.6lf\n", sp[s][m].w[2]);
printf("\n");
printf("________________ cfl for part ________________________\n");
printf("cfl (ipc=0) x    :%12.6lf\n", fabs(sp[s][m].v[0]*si.ts/si.dl[0]));
printf("cfl (ipc=0) y    :%12.6lf\n", fabs(sp[s][m].v[1]*si.ts/si.dl[1]));
printf("cfl (ipc=0) z    :%12.6lf\n", fabs(sp[s][m].v[2]*si.ts/si.dl[2]));
printf("cfl (ipc=1) x    :%12.6lf\n", fabs(sp[s][m].w[0]*si.ts/si.dl[0]));
printf("cfl (ipc=1) y    :%12.6lf\n", fabs(sp[s][m].w[1]*si.ts/si.dl[1]));
printf("cfl (ipc=1) z    :%12.6lf\n", fabs(sp[s][m].w[2]*si.ts/si.dl[2]));
printf("\n");

/* __ part position on g1 __ */
xw = (sp[s][m].r[0]/si.dl[0]-sx->i0[0] < sx->n[0]) ?
      sp[s][m].r[0]/si.dl[0]-sx->i0[0] : sx->n[0]-EPS4;
yw = (sp[s][m].r[1]/si.dl[1]-sx->i0[1] < sx->n[1]) ?
      sp[s][m].r[1]/si.dl[1]-sx->i0[1] : sx->n[1]-EPS4;
zw = (sp[s][m].r[2]/si.dl[2]-sx->i0[2] < sx->n[2]) ?
      sp[s][m].r[2]/si.dl[2]-sx->i0[2] : sx->n[2]-EPS4;

/* __ index for the part. "position" __ */
i = (int)floor(xw);
j = (int)floor(yw);
k = (int)floor(zw);

/* __ indexes of the rounding grid points on g1 __ */
ijk0 = IDX(i  , j  , k  , l0, l1, l2);
ijk1 = IDX(i  , j  , k  , m0, m1, m2);
ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);
printf("________________ coordinates of the cell in g1 _______\n");
printf("i                :%12d\n", i);
printf("j                :%12d\n", j);
printf("k                :%12d\n", k);
printf("\n");
//printf("________________ # of part per cell : %6d _________\n", s1[ijk1].ppc);
printf("________________ # of part per cell : %6d _________\n", s0[ijk0].ppc);
printf("\n");

/* __ magnetic field components on the rounding grid points __ */
printf("________________ magnetic field components ___________\n");
printf("s1[ijk1].b[0]    :%12.6lf\n", s1[ijk1].b[0]);
printf("           1     :%12.6lf\n", s1[ijk1].b[1]);
printf("           2     :%12.6lf\n", s1[ijk1].b[2]);
printf("s1[ijk2].b[0]    :%12.6lf\n", s1[ijk2].b[0]);
printf("           1     :%12.6lf\n", s1[ijk2].b[1]);
printf("           2     :%12.6lf\n", s1[ijk2].b[2]);
printf("s1[ijk3].b[0]    :%12.6lf\n", s1[ijk3].b[0]);
printf("           1     :%12.6lf\n", s1[ijk3].b[1]);
printf("           2     :%12.6lf\n", s1[ijk3].b[2]);
printf("s1[ijk4].b[0]    :%12.6lf\n", s1[ijk4].b[0]);
printf("           1     :%12.6lf\n", s1[ijk4].b[1]);
printf("           2     :%12.6lf\n", s1[ijk4].b[2]);
printf("s1[ijk5].b[0]    :%12.6lf\n", s1[ijk5].b[0]);
printf("           1     :%12.6lf\n", s1[ijk5].b[1]);
printf("           2     :%12.6lf\n", s1[ijk5].b[2]);
printf("s1[ijk6].b[0]    :%12.6lf\n", s1[ijk6].b[0]);
printf("           1     :%12.6lf\n", s1[ijk6].b[1]);
printf("           2     :%12.6lf\n", s1[ijk6].b[2]);
printf("s1[ijk7].b[0]    :%12.6lf\n", s1[ijk7].b[0]);
printf("           1     :%12.6lf\n", s1[ijk7].b[1]);
printf("           2     :%12.6lf\n", s1[ijk7].b[2]);
printf("s1[ijk8].b[0]    :%12.6lf\n", s1[ijk8].b[0]);
printf("           1     :%12.6lf\n", s1[ijk8].b[1]);
printf("           2     :%12.6lf\n", s1[ijk8].b[2]);
printf("\n");

/* __ part position on g2 __ */
xw = sp[s][m].r[0]/si.dl[0]-sx->i0[0]+0.5;
yw = sp[s][m].r[1]/si.dl[1]-sx->i0[1]+0.5;
zw = sp[s][m].r[2]/si.dl[2]-sx->i0[2]+0.5;

/* __ index for the part. "position" __ */
i = (int)floor(xw);
j = (int)floor(yw);
k = (int)floor(zw);

/* __ indexes of the rounding grid points on g2 __ */
ijk1 = IDX(i  , j  , k  , n0, n1, n2);
ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

printf("________________ coordinates of the cell in g2 _______\n");
printf("i                :%12d\n", i);
printf("j                :%12d\n", j);
printf("k                :%12d\n", k);
printf("\n");

/* __ electric field components on the rounding grid points __ */
printf("________________ electric field components ___________\n");
printf("s2[ijk1].e[0]    :%12.6lf\n", s2[ijk1].e[0]);
printf("           1     :%12.6lf\n", s2[ijk1].e[1]);
printf("           2     :%12.6lf\n", s2[ijk1].e[2]);
printf("s2[ijk2].e[0]    :%12.6lf\n", s2[ijk2].e[0]);
printf("           1     :%12.6lf\n", s2[ijk2].e[1]);
printf("           2     :%12.6lf\n", s2[ijk2].e[2]);
printf("s2[ijk3].e[0]    :%12.6lf\n", s2[ijk3].e[0]);
printf("           1     :%12.6lf\n", s2[ijk3].e[1]);
printf("           2     :%12.6lf\n", s2[ijk3].e[2]);
printf("s2[ijk4].e[0]    :%12.6lf\n", s2[ijk4].e[0]);
printf("           1     :%12.6lf\n", s2[ijk4].e[1]);
printf("           2     :%12.6lf\n", s2[ijk4].e[2]);
printf("s2[ijk5].e[0]    :%12.6lf\n", s2[ijk5].e[0]);
printf("           1     :%12.6lf\n", s2[ijk5].e[1]);
printf("           2     :%12.6lf\n", s2[ijk5].e[2]);
printf("s2[ijk6].e[0]    :%12.6lf\n", s2[ijk6].e[0]);
printf("           1     :%12.6lf\n", s2[ijk6].e[1]);
printf("           2     :%12.6lf\n", s2[ijk6].e[2]);
printf("s2[ijk7].e[0]    :%12.6lf\n", s2[ijk7].e[0]);
printf("           1     :%12.6lf\n", s2[ijk7].e[1]);
printf("           2     :%12.6lf\n", s2[ijk7].e[2]);
printf("s2[ijk8].e[0]    :%12.6lf\n", s2[ijk8].e[0]);
printf("           1     :%12.6lf\n", s2[ijk8].e[1]);
printf("           2     :%12.6lf\n", s2[ijk8].e[2]);
printf("\n");

/* __ electric field components on the rounding grid points __ */
printf("________________ ion velocity components _____________\n");
printf("s2[ijk1].vi[0]   :%12.6lf\n", s2[ijk1].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk1].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk1].vi[2]);
printf("s2[ijk2].vi[0]   :%12.6lf\n", s2[ijk2].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk2].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk2].vi[2]);
printf("s2[ijk3].vi[0]   :%12.6lf\n", s2[ijk3].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk3].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk3].vi[2]);
printf("s2[ijk4].vi[0]   :%12.6lf\n", s2[ijk4].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk4].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk4].vi[2]);
printf("s2[ijk5].vi[0]   :%12.6lf\n", s2[ijk5].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk5].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk5].vi[2]);
printf("s2[ijk6].vi[0]   :%12.6lf\n", s2[ijk6].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk6].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk6].vi[2]);
printf("s2[ijk7].vi[0]   :%12.6lf\n", s2[ijk7].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk7].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk7].vi[2]);
printf("s2[ijk8].vi[0]   :%12.6lf\n", s2[ijk8].vi[0]);
printf("            1    :%12.6lf\n", s2[ijk8].vi[1]);
printf("            2    :%12.6lf\n", s2[ijk8].vi[2]);
printf("\n");

/* __ electron density on the rounding grid points ___________ */
printf("________________ electron density ____________________\n");
printf("s2[ijk1].ns[0]   :%12.6lf\n", s2[ijk1].ns[0]);
printf("s2[ijk2].ns[0]   :%12.6lf\n", s2[ijk2].ns[0]);
printf("s2[ijk3].ns[0]   :%12.6lf\n", s2[ijk3].ns[0]);
printf("s2[ijk4].ns[0]   :%12.6lf\n", s2[ijk4].ns[0]);
printf("s2[ijk5].ns[0]   :%12.6lf\n", s2[ijk5].ns[0]);
printf("s2[ijk6].ns[0]   :%12.6lf\n", s2[ijk6].ns[0]);
printf("s2[ijk7].ns[0]   :%12.6lf\n", s2[ijk7].ns[0]);
printf("s2[ijk8].ns[0]   :%12.6lf\n", s2[ijk8].ns[0]);
printf("\n");

/* __ current density components on the rounding grid points _ */
printf("________________ current density components __________\n");
printf("s2[ijk1].j[0]    :%12.6lf\n", s2[ijk1].j[0]);
printf("           1     :%12.6lf\n", s2[ijk1].j[1]);
printf("           2     :%12.6lf\n", s2[ijk1].j[2]);
printf("s2[ijk2].j[0]    :%12.6lf\n", s2[ijk2].j[0]);
printf("           1     :%12.6lf\n", s2[ijk2].j[1]);
printf("           2     :%12.6lf\n", s2[ijk2].j[2]);
printf("s2[ijk3].j[0]    :%12.6lf\n", s2[ijk3].j[0]);
printf("           1     :%12.6lf\n", s2[ijk3].j[1]);
printf("           2     :%12.6lf\n", s2[ijk3].j[2]);
printf("s2[ijk4].j[0]    :%12.6lf\n", s2[ijk4].j[0]);
printf("           1     :%12.6lf\n", s2[ijk4].j[1]);
printf("           2     :%12.6lf\n", s2[ijk4].j[2]);
printf("s2[ijk5].j[0]    :%12.6lf\n", s2[ijk5].j[0]);
printf("           1     :%12.6lf\n", s2[ijk5].j[1]);
printf("           2     :%12.6lf\n", s2[ijk5].j[2]);
printf("s2[ijk6].j[0]    :%12.6lf\n", s2[ijk6].j[0]);
printf("           1     :%12.6lf\n", s2[ijk6].j[1]);
printf("           2     :%12.6lf\n", s2[ijk6].j[2]);
printf("s2[ijk7].j[0]    :%12.6lf\n", s2[ijk7].j[0]);
printf("           1     :%12.6lf\n", s2[ijk7].j[1]);
printf("           2     :%12.6lf\n", s2[ijk7].j[2]);
printf("s2[ijk8].j[0]    :%12.6lf\n", s2[ijk8].j[0]);
printf("           1     :%12.6lf\n", s2[ijk8].j[1]);
printf("           2     :%12.6lf\n", s2[ijk8].j[2]);
printf("\n");

/* __ current density components on the rounding grid points _ */
printf("________________ current density components __________\n");
printf("s2[ijk1].ps[0][0]:%12.6lf\n", s2[ijk1].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk1].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk1].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk1].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk1].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk1].ps[0][5]);
printf("s2[ijk2].ps[0][0]:%12.6lf\n", s2[ijk2].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk2].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk2].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk2].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk2].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk2].ps[0][5]);
printf("s2[ijk3].ps[0][0]:%12.6lf\n", s2[ijk3].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk3].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk3].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk3].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk3].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk3].ps[0][5]);
printf("s2[ijk4].ps[0][0]:%12.6lf\n", s2[ijk4].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk4].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk4].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk4].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk4].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk4].ps[0][5]);
printf("s2[ijk5].ps[0][0]:%12.6lf\n", s2[ijk5].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk5].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk5].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk5].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk5].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk5].ps[0][5]);
printf("s2[ijk6].ps[0][0]:%12.6lf\n", s2[ijk6].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk6].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk6].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk6].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk6].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk6].ps[0][5]);
printf("s2[ijk7].ps[0][0]:%12.6lf\n", s2[ijk7].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk7].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk7].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk7].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk7].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk7].ps[0][5]);
printf("s2[ijk8].ps[0][0]:%12.6lf\n", s2[ijk8].ps[0][0]);
printf("               1 :%12.6lf\n", s2[ijk8].ps[0][1]);
printf("               2 :%12.6lf\n", s2[ijk8].ps[0][2]);
printf("               3 :%12.6lf\n", s2[ijk8].ps[0][3]);
printf("               4 :%12.6lf\n", s2[ijk8].ps[0][4]);
printf("               5 :%12.6lf\n", s2[ijk8].ps[0][5]);
printf("\n");

/* __ misc informations ______________________________________ */
printf("________________ iamdead : file ______________________\n");
printf("file             : %s\n", file);
printf("line             :%12d\n", line);
printf("date             :%12s\n", __DATE__);
printf("time             :%12s\n", __TIME__);
printf("\n");

/* __ if too many part. in the domain __ */
if (m >= sx->nm[s])
   {
   printf("________________ too many part. in the subdomain _____\n\n");
   sx->irun = 1;
   }

/* __ if the part is out of domain __ */
if (sp[s][m].r[0] < sx->i0[0]*si.dl[0])
   {
   printf("________________ part out of domain : r[0] < i0*dl ___\n\n");
   sx->irun = 2;
sp[s][m].r[0] = 0.0;
   }
if (sp[s][m].r[0] > sx->i1[0]*si.dl[0])
   {
   printf("________________ part out of domain : r[0] > i1*dl ___\n\n");
   sx->irun = 3;
sp[s][m].r[0] = 0.0;
   }
if (sp[s][m].r[1] < sx->i0[1]*si.dl[1])
   {
   printf("________________ part out of domain : r[1] < i0*dl ___\n\n");
   sx->irun = 4;
sp[s][m].r[1] = 0.0;
   }
if (sp[s][m].r[1] > sx->i1[1]*si.dl[1])
   {
   printf("________________ part out of domain : r[1] > i1*dl ___\n\n");
   sx->irun = 5;
sp[s][m].r[1] = 0.0;
   }
if (sp[s][m].r[2] < sx->i0[2]*si.dl[2])
   {
   printf("________________ part out of domain : r[2] < i0*dl ___\n\n");
   sx->irun = 6;
sp[s][m].r[2] = 0.0;
   }
if (sp[s][m].r[2] > sx->i1[2]*si.dl[2])
   {
   printf("________________ part out of domain : r[2] > i1*dl ___\n\n");
   sx->irun = 7;
sp[s][m].r[2] = 0.0;
   }
if (sp[s][m].s[0] < sx->i0[0]*si.dl[0])
   {
   printf("________________ part out of domain : s[0] < i0*dl ___\n\n");
   sx->irun = 8;
sp[s][m].s[0] = 0.0;
   }
if (sp[s][m].s[0] > sx->i1[0]*si.dl[0])
   {
   printf("________________ part out of domain : s[0] > i1*dl ___\n\n");
   sx->irun = 9;
sp[s][m].s[0] = 0.0;
   }
if (sp[s][m].s[1] < sx->i0[1]*si.dl[1])
   {
   printf("________________ part out of domain : s[1] < i0*dl ___\n\n");
   sx->irun = 10;
sp[s][m].s[1] = 0.0;
   }
if (sp[s][m].s[1] > sx->i1[1]*si.dl[1])
   {
   printf("________________ part out of domain : s[1] > i1*dl ___\n\n");
   sx->irun = 11;
sp[s][m].s[1] = 0.0;
   }
if (sp[s][m].s[2] < sx->i0[2]*si.dl[2])
   {
   printf("________________ part out of domain : s[2] < i0*dl ___\n\n");
   sx->irun = 12;
sp[s][m].s[2] = 0.0;
   }
if (sp[s][m].s[2] > sx->i1[2]*si.dl[2])
   {
   printf("________________ part out of domain : s[2] > i1*dl ___\n\n");
   sx->irun = 13;
sp[s][m].s[2] = 0.0;
   }

/* __ if cfl condition is not satisfied __ */
if (fabs(sp[s][m].v[0]*si.ts) > si.dl[0])
   {
   printf("________________ cfl not satisfied : v[0] > dl/ts ____\n\n");
   sx->irun = 14;
sp[s][m].v[0] = 0.0;
   }
if (fabs(sp[s][m].v[1]*si.ts) > si.dl[1])
   {
   printf("________________ cfl not satisfied : v[1] > dl/ts ____\n\n");
   sx->irun = 15;
sp[s][m].v[1] = 0.0;
   }
if (fabs(sp[s][m].v[2]*si.ts) > si.dl[2])
   {
   printf("________________ cfl not satisfied : v[2] > dl/ts ____\n\n");
   sx->irun = 16;
sp[s][m].v[2] = 0.0;
   }
if (fabs(sp[s][m].w[0]*si.ts) > si.dl[0])
   {
   printf("________________ cfl not satisfied : w[0] > dl/ts ____\n\n");
   sx->irun = 17;
sp[s][m].w[0] = 0.0;
   }
if (fabs(sp[s][m].w[1]*si.ts) > si.dl[1])
   {
   printf("________________ cfl not satisfied : w[1] > dl/ts ____\n\n");
   sx->irun = 18;
sp[s][m].w[1] = 0.0;
   }
if (fabs(sp[s][m].w[2]*si.ts) > si.dl[2])
   {
   printf("________________ cfl not satisfied : w[2] > dl/ts ____\n\n");
   sx->irun = 19;
sp[s][m].w[2] = 0.0;
   }

printf("\n");

}


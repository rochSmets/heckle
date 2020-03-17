
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//     attention : le forcage est sur toute la boite...                       //
//                 adapter si on ne le veux que sur un bout de la boite       //
//                 voir munster...                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "defines.h"
#include "structures.h"
#include "plasma.h"


/* __ read "munster.txt" to set stt structure (magnetic fluctuations) _____ */
void topo(struct sti si, struct stx sx, struct stt *st, MPI_Comm com)
{
double w1b, w2b, w1v;
double wb, wb1, wb2;
double we, wvx, wvy, wvz;
double kx, ky;
double kxb, kyb, kxv, kyv;
int fg;
int nm;
int f, g, z;
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
z = fscanf(fp, "%lf", &(st->n)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->k)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->p)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->a)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf", st->z, st->z+1); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->gb)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf %lf %lf", st->e, st->e+1, st->e+2); if (z != 3) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%lf", &(st->gv)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%d %d", &(st->m[0]), &(st->m[1])); if (z != 2) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);
z = fscanf(fp, "%d", &(st->id)); if (z != 1) shit(sx.r);
z = fscanf(fp, "%s", junk); if (z != 1) shit(sx.r);

/* __ close the file __ */
fclose(fp);

/* __ # of modes __ */
nm = st->m[1]-st->m[0]+1;

/* __ memory allocation for magnetic energy density & phase __ */
st->bxy = (double *)malloc(nm*nm*sizeof(double));
st->bzz = (double *)malloc(nm*nm*sizeof(double));
st->vxx = (double *)malloc(nm*nm*sizeof(double));
st->vyy = (double *)malloc(nm*nm*sizeof(double));
st->vzz = (double *)malloc(nm*nm*sizeof(double));
st->pxy = (double *)malloc(nm*nm*sizeof(double));
st->pzz = (double *)malloc(nm*nm*sizeof(double));
st->fxx = (double *)malloc(nm*nm*sizeof(double));
st->fyy = (double *)malloc(nm*nm*sizeof(double));
st->fzz = (double *)malloc(nm*nm*sizeof(double));

/* __ set initial sum for 1d & 2d __ */
w1b = 0.0;
w2b = 0.0;
w1v = 0.0;

/* __ summation over the wave numbers : nested loops on the wave numbers __ */
for (f = st->m[0]; f <= st->m[1]; f++)
    {
    for (g = st->m[0]; g <= st->m[1]; g++)
        {
        /* __ wave number in x & y dir. __ */
        kx = 2.0*f*PI/si.l[0];
        ky = 2.0*g*PI/si.l[1];

        /* __ wave number at -gamma power __ */
        kxb = pow(kx,-st->gb);
        kyb = pow(ky,-st->gb);
        kxv = pow(kx,-st->gv);
        kyv = pow(ky,-st->gv);

        /* __ sum for 1d __ */
        w1b += kxb*kyb;

        /* __ sum for 2d __ */
        w2b += kxb*kyb*(1.0+(kxb*kxb)/(kyb*kyb));

        /* __ sum for 1d __ */
        w1v += kxv*kyv;
        }
    }

/* __ magnetic energy __ */
wb = 0.5*(st->b[0]*st->b[0]+st->b[1]*st->b[1]+st->b[2]*st->b[2]);

/* __ kinetic thermal energy (1d) __ */
we = 0.5*(st->k-wb);

/* __ magnetic energy density in z dir. __ */
wb1 = sqrt(st->z[1]*wb/w1b);

/* __ magnetic energy density in xy plan __ */
wb2 = sqrt(st->z[0]*wb/w2b);

/* __ kinetic energy density in x, y or z dir. __ */
wvx = sqrt(st->e[0]*we/w1v);
wvy = sqrt(st->e[1]*we/w1v);
wvz = sqrt(st->e[2]*we/w1v);

/* __ weighted energy for each modes : nested loops on the wave numbers __ */
for (f = st->m[0]; f <= st->m[1]; f++)
    {
    for (g = st->m[0]; g <= st->m[1]; g++)
        {
        /* __ set index __ */
        fg = idx(f-st->m[0], g-st->m[0], 0, nm, nm, 1);

        /* __ wave number in x & y dir. __ */
        kx = 2.0*f*PI/si.l[0];
        ky = 2.0*g*PI/si.l[1];

        /* __ wave number at -gamma/2 power __ */
        kxb = pow(kx,-0.5*st->gb);
        kyb = pow(ky,-0.5*st->gb);
        kxv = pow(kx,-0.5*st->gv);
        kyv = pow(ky,-0.5*st->gv);

        /* __ weighted energy for each modes __ */
        st->bxy[fg] = wb2*kxb*kyb;
        st->bzz[fg] = wb1*kxb*kyb;
        st->vxx[fg] = wvx*kxv*kyv;
        st->vyy[fg] = wvy*kxv*kyv;
        st->vzz[fg] = wvz*kxv*kyv;
        }
    }

/* __ set the modes' phases : only on node 0 __ */
if (sx.r == 0)
   {
   /* __ nested loops on the wave numbers __ */
   for (f = st->m[0]; f <= st->m[1]; f++)
       {
       for (g = st->m[0]; g <= st->m[1]; g++)
           {
           /* __ set index __ */
           fg = idx(f-st->m[0], g-st->m[0], 0, nm, nm, 1);

           /* __ set randomly the wave phase __ */
           st->pxy[fg] = RNM*2.0*PI;
           st->pzz[fg] = RNM*2.0*PI;
           st->fxx[fg] = RNM*2.0*PI;
           st->fyy[fg] = RNM*2.0*PI;
           st->fzz[fg] = RNM*2.0*PI;
           }
       }
   }

/* __ broadcast the phases of each modes __ */
MPI_Bcast(st->pxy, nm*nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->pzz, nm*nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->fxx, nm*nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->fyy, nm*nm, MPI_DOUBLE, 0, co);
MPI_Bcast(st->fzz, nm*nm, MPI_DOUBLE, 0, co);

/* __ display the topology parameters __ */
if (sx.r == 0)
   {
   printf("________________ magnetic topology : b & v fluct. ____\n");
   printf("b              : %8.4f    %8.4f    %8.4f\n", st->b[0], st->b[1], st->b[2]);
   printf("n protons      : %8.4f\n", st->n);
   printf("anisotropy     : %8.4f\n", st->a);
   printf("beta protons   : %8.4f\n", st->p*(st->k/wb-1.0));
   printf("beta electrons : %8.4f\n", (1.0-st->p)*(st->k/wb-1.0));
   printf("beta total     : %8.4f\n", st->k/wb-1.0);
   printf("ener. b (xz)   : %8.4f\n", st->z[0]);
   printf("ener. b (z)    : %8.4f\n", st->z[1]);
   printf("gamma (b)      : %8.4f\n", st->gb);
   printf("ener. v (x)    : %8.4f\n", st->e[0]);
   printf("ener. v (y)    : %8.4f\n", st->e[1]);
   printf("ener. v (z)    : %8.4f\n", st->e[2]);
   printf("gamma (v)      : %8.4f\n", st->gv);
   printf("min & max modes: %8d    %8d\n", st->m[0], st->m[1]);
   printf("# ts for drive : %8d\n", st->id);
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
n[0] = st->n;

/* __ proton density __ */
n[1] = st->n;

/* __ alpha density __ */
n[2] = 0.0;

}


/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti si, struct stx sx, struct stt *st, double q[3], double b[3], MPI_Comm com)
{


/* __ dc magnetic field components __ */
b[0] = st->b[0];
b[1] = st->b[1];
b[2] = st->b[2];

}


/* __ drive the magnetic field ______________________________________________ */
void drive(struct sti si, struct stx sx, struct stt *st, struct st1 *s1, MPI_Comm com)
{
double q[3];
double kx, ky;
double kq;
double wt;
int fg;
int nm;
int ijk;
int f, g, i, j, k;
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ # of modes __ */
nm = st->m[1]-st->m[0]+1;

/* __ factor associated to time forcing __ */
wt = 1.0/(st->id+1);

/* __ superposition of each modes : nested loops on subdomains __ */
for (i = 0; i < sx.n[0]+1; i++)
    {
    for (j = 0; j < sx.n[1]+1; j++)
        {
        for (k = 0; k < sx.n[2]+1; k++)
            {
            /* __ set index on g2 __ */
            ijk = idx(i, j, k, sx.n[0]+1, sx.n[1]+1, sx.n[2]+1);

            /* __ coordinates of grid point __ */
            q[0] = (sx.i0[0]+i)*si.dl[0];
            q[1] = (sx.i0[1]+j)*si.dl[1];
            q[2] = (sx.i0[2]+k)*si.dl[2];

            /* __ set the vel. field components __ */
            s1[ijk].v[0] = 0.0;
            s1[ijk].v[1] = 0.0;
            s1[ijk].v[2] = 0.0;

            /* __ nested loops on the modes __ */
            for (f = st->m[0]; f <= st->m[1]; f++)
                {
                for (g = st->m[0]; g <= st->m[1]; g++)
                    {
                    /* __ set index __ */
                    fg = idx(f-st->m[0], g-st->m[0], 0, nm, nm, 1);

                    /* __ wave number in x & y dir. __ */
                    kx = 2.0*f*PI/si.l[0];
                    ky = 2.0*g*PI/si.l[1];

                    /* __ set k dot q __ */
                    kq = kx*q[0]+ky*q[1];

                    /* __ set the mag. field components __ */
                    s1[ijk].b[0] += wt*st->bxy[fg]*cos(kq+st->pxy[fg]);
                    s1[ijk].b[1] += wt*st->bxy[fg]*cos(kq+st->pxy[fg])*(-kx/ky);
                    s1[ijk].b[2] += wt*st->bzz[fg]*cos(kq+st->pzz[fg]);

                    /* __ set the vel. field components __ */
                    s1[ijk].v[0] += wt*st->vxx[fg]*cos(kq+st->fxx[fg]);
                    s1[ijk].v[1] += wt*st->vyy[fg]*cos(kq+st->fyy[fg]);
                    s1[ijk].v[2] += wt*st->vzz[fg]*cos(kq+st->fzz[fg]);
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
double wp, we;
double wi0, wi1;
double ww;


/* __ total beta __ */
wt = 2.0*st->k/(st->b[0]*st->b[0]+st->b[1]*st->b[1]+st->b[2]*st->b[2])-1.0;

/* __ beta of protons, alphas & electrons __ */
wp = st->p*wt;
we = (1.0-st->p)*wt;

/* __ anisotropy factors __ */
wi0 = 3.0*st->a/(st->a+2.0);
wi1 = 3.0      /(st->a+2.0);

/* __ magnetic pressure __ */
ww = 0.5*(st->b[0]*st->b[0]+st->b[1]*st->b[1]+st->b[2]*st->b[2]);

/* __ electron temperatures __ */
ta[0] = we*ww/st->n;
te[0] = we*ww/st->n;

/* __  proton temperatures __ */
ta[1] = wp*wi0*ww/st->n;
te[1] = wp*wi1*ww/st->n;

/* __  alpha temperatures __ */
ta[2] = 0.0;
te[2] = 0.0;

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


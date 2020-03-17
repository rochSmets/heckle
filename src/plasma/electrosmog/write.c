
#ifndef PYP
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structures.h"
#include "plasma.h"
#include "defines.h"
#include "misc.h"
#include "iamdead.h"


/* __ write the xplosed fields ______________________________________________ */
void writefield(struct sti si, struct stx sx, struct st1 *s1, struct st2 *s2, int it)
{
float *bw[3], *ew[3], *jw[3], *iw[3];
float *nw[NS+1], *vw[NS+1][3], *pw[NS+1][6];
int i, j, k, l, s;
int ijk;
int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
int m0, m1, m2, n0, n1, n2;
int wi;
int nn1;
char nfile[16];
FILE *fp;


/* __ #of grid points __ */
m0 = sx.n[0]+1;
m1 = sx.n[1]+1;
m2 = sx.n[2]+1;
n0 = sx.n[0]+2;
n1 = sx.n[1]+2;
n2 = sx.n[2]+2;

/* __ set current ts __ */
wi = it/si.tf;

/* __ # of grid points on g1 __ */
nn1 = m0*m1*m2;

/* __ print informations __ */
if (sx.r == 0) {
   printf("________________ write f-dump # %4i for node %3i ____\n", wi, sx.r);
   printf("\n");
}

/* __ make the file name __ */
sprintf(nfile, "xf%03i-%04i.dat", sx.r, wi);

/* __ open the file __ */
fp = fopen(nfile, "wb");
if (fp == NULL) printf("problem in opening file %s\n", nfile);

/* __ memory allocation __ */
for (l = 0; l < 3; l++) {
   bw[l] = (float *) malloc(nn1*sizeof(float));
   ew[l] = (float *) malloc(nn1*sizeof(float));
   jw[l] = (float *) malloc(nn1*sizeof(float));
   iw[l] = (float *) malloc(nn1*sizeof(float));
}

for (s = 0; s < NS+1; s++) {
   nw[s] = (float *) malloc(nn1*sizeof(float));

   for (l = 0; l < 3; l++) {
      vw[s][l] = (float *) malloc(nn1*sizeof(float));
   }

   for (l = 0; l < 6; l++) {
      pw[s][l] = (float *) malloc(nn1*sizeof(float));
   }
}

/* __ fill the arrays : nested loops on the subdomain __ */
for (i = 0; i < m0; i++) {
   for (j = 0; j < m1; j++) {
      for (k = 0; k < m2; k++) {
         /* __ set index on g1 __ */
         ijk = IDX(i, j, k, m0, m1, m2);

         /* __ set indexes on g2 __ */
         ijk1 = IDX(i  , j  , k  , n0, n1, n2);
         ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
         ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
         ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
         ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
         ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
         ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
         ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

         /* __ for g1 grid __ */
         for (l = 0; l < 3; l++) {
            bw[l][ijk] = (float)s1[ijk].b[l];
         }

         /* __ for g2 grid __ */
         for (l = 0; l < 3; l++) {
            ew[l][ijk] = (float)0.125*(s2[ijk1].e[l]
                                      +s2[ijk2].e[l]
                                      +s2[ijk3].e[l]
                                      +s2[ijk4].e[l]
                                      +s2[ijk5].e[l]
                                      +s2[ijk6].e[l]
                                      +s2[ijk7].e[l]
                                      +s2[ijk8].e[l]);

            jw[l][ijk] = (float)0.125*(s2[ijk1].j[l]
                                      +s2[ijk2].j[l]
                                      +s2[ijk3].j[l]
                                      +s2[ijk4].j[l]
                                      +s2[ijk5].j[l]
                                      +s2[ijk6].j[l]
                                      +s2[ijk7].j[l]
                                      +s2[ijk8].j[l]);

            iw[l][ijk] = (float)0.125*(s2[ijk1].vi[l]
                                      +s2[ijk2].vi[l]
                                      +s2[ijk3].vi[l]
                                      +s2[ijk4].vi[l]
                                      +s2[ijk5].vi[l]
                                      +s2[ijk6].vi[l]
                                      +s2[ijk7].vi[l]
                                      +s2[ijk8].vi[l]);
         }

         for (s = 0; s < NS+1; s++) {
            nw[s][ijk] = (float)0.125*(s2[ijk1].ns[s]
                                      +s2[ijk2].ns[s]
                                      +s2[ijk3].ns[s]
                                      +s2[ijk4].ns[s]
                                      +s2[ijk5].ns[s]
                                      +s2[ijk6].ns[s]
                                      +s2[ijk7].ns[s]
                                      +s2[ijk8].ns[s]);

            for (l = 0; l < 3; l++) {
               vw[s][l][ijk] = (float)0.125*(s2[ijk1].vs[s][l]
                                            +s2[ijk2].vs[s][l]
                                            +s2[ijk3].vs[s][l]
                                            +s2[ijk4].vs[s][l]
                                            +s2[ijk5].vs[s][l]
                                            +s2[ijk6].vs[s][l]
                                            +s2[ijk7].vs[s][l]
                                            +s2[ijk8].vs[s][l]);
            }

            for (l = 0; l < 6; l++) {
               pw[s][l][ijk] = (float)0.125*(s2[ijk1].ps[s][l]
                                            +s2[ijk2].ps[s][l]
                                            +s2[ijk3].ps[s][l]
                                            +s2[ijk4].ps[s][l]
                                            +s2[ijk5].ps[s][l]
                                            +s2[ijk6].ps[s][l]
                                            +s2[ijk7].ps[s][l]
                                            +s2[ijk8].ps[s][l]);
            }
         }
      }
   }
}

/* __ write the file __ */
fwrite(sx.i0, sizeof(int), 3, fp);
fwrite(sx.n, sizeof(int), 3, fp);

for (l = 0; l < 3; l++) {
   fwrite(bw[l], sizeof(float), nn1, fp);
}

for (l = 0; l < 3; l++) {
   fwrite(ew[l], sizeof(float), nn1, fp);
}

for (l = 0; l < 3; l++) {
   fwrite(jw[l], sizeof(float), nn1, fp);
}

for (l = 0; l < 3; l++) {
   fwrite(iw[l], sizeof(float), nn1, fp);
}

for (s = 0; s < NS+1; s++) {
   fwrite(nw[s], sizeof(float), nn1, fp);
}

for (s = 0; s < NS+1; s++) {
   for (l = 0; l < 3; l++) {
      fwrite(vw[s][l], sizeof(float), nn1, fp);
   }
}

for (s = 0; s < NS+1; s++) {
   for (l = 0; l < 6; l++) {
      fwrite(pw[s][l], sizeof(float), nn1, fp);
   }
}

/* __ clean-up the pointers __ */
for (l = 0; l < 3; l++) {
   free(bw[l]);
   free(ew[l]);
   free(jw[l]);
   free(iw[l]);
}

for (s = 0; s < NS+1; s++) {
   free(nw[s]);

   for (l = 0; l < 3; l++) {
      free(vw[s][l]);
   }

   for (l = 0; l < 6; l++) {
      free(pw[s][l]);
   }
}

/* __ close the file __ */
fclose(fp);

}


/* __ write the xplosed part. _______________________________________________ */
void writepart(struct sti si, struct stx sx, struct stp *sp[NS+1], int it)
{
float *xw[NS+1], *yw[NS+1], *zw[NS+1];
float *uw[NS+1], *vw[NS+1], *ww[NS+1];
int *iw[NS+1];
int m, s;
int wi;
char nfile[16];
FILE *fp;


/* __ set current ts __ */
wi = it/si.tp;

/* __ print informations __ */
if (sx.r == 0)
   {
   printf("________________ write p-dump # %4i for node %3i ____\n", wi, sx.r);
   printf("\n");
   }

/* __ make the file name __ */
sprintf(nfile, "xp%03i-%04i.dat", sx.r, wi);

/* __ open the file __ */
fp = fopen(nfile, "wb");
if (fp == NULL) printf("problem in opening file %s\n", nfile);

/* __ loop on the part. __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ memory allocation __ */
    if (sx.ns[s] != 0)
       {
       iw[s] = (int *)malloc(sx.ns[s]*sizeof(int));

       xw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
       yw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
       zw[s] = (float *)malloc(sx.ns[s]*sizeof(float));

       uw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
       vw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
       ww[s] = (float *)malloc(sx.ns[s]*sizeof(float));
       }
    }

/* __ fill the arrays : loop on the species __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ loop on the part. of specie "s" __ */
    for (m = 0; m < sx.ns[s]; m++)
        {
        /* __ fill the buffers __ */
        iw[s][m] = sp[s][m].i;

        xw[s][m] = (float)sp[s][m].r[0];
        yw[s][m] = (float)sp[s][m].r[1];
        zw[s][m] = (float)sp[s][m].r[2];

        uw[s][m] = (float)sp[s][m].v[0];
        vw[s][m] = (float)sp[s][m].v[1];
        ww[s][m] = (float)sp[s][m].v[2];
        }
    }

/* __ write the file __ */
fwrite(sx.ns, sizeof(int), NS+1, fp);
fwrite(si.ns, sizeof(int), NS+1, fp);

/* __ loop on the part. __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ write the buffers in the file __ */
    if (sx.ns[s] != 0)
       {
       fwrite(iw[s], sizeof(int), sx.ns[s], fp);

       fwrite(xw[s], sizeof(float), sx.ns[s], fp);
       fwrite(yw[s], sizeof(float), sx.ns[s], fp);
       fwrite(zw[s], sizeof(float), sx.ns[s], fp);

       fwrite(uw[s], sizeof(float), sx.ns[s], fp);
       fwrite(vw[s], sizeof(float), sx.ns[s], fp);
       fwrite(ww[s], sizeof(float), sx.ns[s], fp);
       }
    }

/* __ loop on the species __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ clean-up the pointers __ */
    if (sx.ns[s] != 0)
       {
       free(iw[s]);

       free(xw[s]);
       free(yw[s]);
       free(zw[s]);

       free(uw[s]);
       free(vw[s]);
       free(ww[s]);
       }
    }

/* __ close the file __ */
fclose(fp);

}


/* __ write the time dump ___________________________________________________ */
void writedump(struct sti si, struct stx sx, struct std sd, int it)
{
float tt;
float db, dw, ma, pb,pa, pe, ab, aa, ae, ea, ee, fx, fy, fz;
FILE *fp;


/* __ only on node 0 __ */
if (sx.r == 0)
   {
   tt = (float)it*si.ts;
   db = (float)sd.db;
   dw = (float)sd.dw;
   ma = (float)sd.ma;
   pb = (float)sd.pb;
   pa = (float)sd.pa;
   pe = (float)sd.pe;
   ab = (float)sd.ab;
   aa = (float)sd.aa;
   ae = (float)sd.ae;
   ea = (float)sd.ea;
   ee = (float)sd.ee;
   fx = (float)sd.fx;
   fy = (float)sd.fy;
   fz = (float)sd.fz;

   /* __ print informations __ */
   printf("________________ time step : %8d ________________\n", it);
   printf("\n");

   /* __ remove old dump __ */
   if (it == 0) remove("hdump.dat");

   /* __ open the file __ */
   fp = fopen("hdump.dat", "a+b");
   if (fp == NULL) printf("problem in opening file hdump.dat\n");

   /* __ write the file __ */
   fwrite(&tt, sizeof(float), 1, fp);
   fwrite(&db, sizeof(float), 1, fp);
   fwrite(&dw, sizeof(float), 1, fp);
   fwrite(&ma, sizeof(float), 1, fp);
   fwrite(&pb, sizeof(float), 1, fp);
   fwrite(&pa, sizeof(float), 1, fp);
   fwrite(&pe, sizeof(float), 1, fp);
   fwrite(&ab, sizeof(float), 1, fp);
   fwrite(&aa, sizeof(float), 1, fp);
   fwrite(&ae, sizeof(float), 1, fp);
   fwrite(&ea, sizeof(float), 1, fp);
   fwrite(&ee, sizeof(float), 1, fp);
   fwrite(&fx, sizeof(float), 1, fp);
   fwrite(&fy, sizeof(float), 1, fp);
   fwrite(&fz, sizeof(float), 1, fp);

   /* __ close the file __ */
   fclose(fp);
   }

}


/* __ write the xplosed orbits ______________________________________________ */
void writeorbit(struct sti si, struct stx *sx, struct st1 *s1, struct st2 *s2, struct stp *sp[NS+1], struct sto *so, int it)
{
float xw, yw, zw;
float ro[3], wo[3], bo[3], eo[3], jo[3], io[3];
float no[NS+1], vo[NS+1][3], po[NS+1][6];
float lx, ly, lz;
float w1, w2, w3, w4, w5, w6, w7, w8;
int sw, iw;
int wi;
int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
int m0, m1, m2, n0, n1, n2;
int i, j, k, l, m, n, s;
char nfile[16];
FILE *fp;


/* __ # of grid points __ */
m0 = sx->n[0]+1;
m1 = sx->n[1]+1;
m2 = sx->n[2]+1;
n0 = sx->n[0]+2;
n1 = sx->n[1]+2;
n2 = sx->n[2]+2;

/* __ set the so.no value __ */
so->no = 0;

/* __ identify the orbits to follow (with s & m) : loop on the species __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ loop on all the part. __ */
    for (m = 0; m < sx->ns[s]; m++)
        {
        /* __ loop on all the orbits to follow __ */
        for (n = 0; n < so->wo; n++)
            {
            /* __ if the part. is an orbit __ */
            if (s == so->so[n] && sp[s][m].i == so->io[n])
               {
               /* __ keep memory of s & m __ */
               so->s[so->no] = s;
               so->m[so->no] = m;

               /* __ increase the # of orbits on the node __ */
               so->no++;
               }
            }
        }
    }

/* __ set current ts __ */
wi = it/si.tt;

/* __ print informations __ */
if (sx->r == 0)
   {
   printf("________________ write o-dump # %4i for node %3i ____\n", wi, sx->r);
   printf("\n");
   }

/* __ make the file name __ */
sprintf(nfile, "xo%03i-%04i.dat", sx->r, wi);

/* __ open the file __ */
fp = fopen(nfile, "a+b");
if (fp == NULL) printf("problem in opening file %s\n", nfile);

/* __ write the # of orbits & the current ts __ */
fwrite(&so->no, sizeof(int), 1, fp);
fwrite(&wi, sizeof(int), 1, fp);

/* __ loop on the particles __ */
for (n = 0; n < so->no; n++)
    {
    /* __ part orbit __ */
    sw = so->s[n];
    iw = sp[so->s[n]][so->m[n]].i;

    /* __ interpolate b field __ */
    xw = (sp[so->s[n]][so->m[n]].r[0]/si.dl[0]-sx->i0[0] < sx->n[0]) ?
          sp[so->s[n]][so->m[n]].r[0]/si.dl[0]-sx->i0[0] : sx->n[0]-EPS4;
    yw = (sp[so->s[n]][so->m[n]].r[1]/si.dl[1]-sx->i0[1] < sx->n[1]) ?
          sp[so->s[n]][so->m[n]].r[1]/si.dl[1]-sx->i0[1] : sx->n[1]-EPS4;
    zw = (sp[so->s[n]][so->m[n]].r[2]/si.dl[2]-sx->i0[2] < sx->n[2]) ?
          sp[so->s[n]][so->m[n]].r[2]/si.dl[2]-sx->i0[2] : sx->n[2]-EPS4;

    /* __ index for the part. "position" __ */
    i = (int)floor(xw);
    j = (int)floor(yw);
    k = (int)floor(zw);

    #ifdef BUG
    if (i < 0 || i >= sx->n[0])
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (j < 0 || j >= sx->n[1])
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (k < 0 || k >= sx->n[2])
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    #endif

//  #ifdef BUG
//  if (i < 0 || i >= sx->n[0]) shit(sx->r);
//  if (j < 0 || j >= sx->n[1]) shit(sx->r);
//  if (k < 0 || k >= sx->n[2]) shit(sx->r);
//  #endif

    /* __ part. location in the cell __ */
    lx = xw-i;
    ly = yw-j;
    lz = zw-k;

    /* __ indexes of the rounding grid points on g1 __ */
    ijk1 = IDX(i  , j  , k  , m0, m1, m2);
    ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
    ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
    ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
    ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
    ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
    ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
    ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);

    /* __ weight for each vertices of the rounding grid points __ */
    w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
    w2 = (    lx)*(1.0-ly)*(1.0-lz);
    w3 = (1.0-lx)*(    ly)*(1.0-lz);
    w4 = (    lx)*(    ly)*(1.0-lz);
    w5 = (1.0-lx)*(1.0-ly)*(    lz);
    w6 = (    lx)*(1.0-ly)*(    lz);
    w7 = (1.0-lx)*(    ly)*(    lz);
    w8 = (    lx)*(    ly)*(    lz);

    /* __ loop in the 3 directions __ */
    for (l = 0; l < 3; l++)
        {
        /* __ set b field seen by the part. __ */
        bo[l] = w1*s1[ijk1].b[l]
               +w2*s1[ijk2].b[l]
               +w3*s1[ijk3].b[l]
               +w4*s1[ijk4].b[l]
               +w5*s1[ijk5].b[l]
               +w6*s1[ijk6].b[l]
               +w7*s1[ijk7].b[l]
               +w8*s1[ijk8].b[l];
        }

    /* __ interpolate e, n, j, v & pe fields __ */
    xw = sp[so->s[n]][so->m[n]].r[0]/si.dl[0]-sx->i0[0]+0.5;
    yw = sp[so->s[n]][so->m[n]].r[1]/si.dl[1]-sx->i0[1]+0.5;
    zw = sp[so->s[n]][so->m[n]].r[2]/si.dl[2]-sx->i0[2]+0.5;

    /* __ index for the part. "position" __ */
    i = (int)floor(xw);
    j = (int)floor(yw);
    k = (int)floor(zw);

    #ifdef BUG
    if (i < 0 || i >= sx->n[0]+1)
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (j < 0 || j >= sx->n[1]+1)
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (k < 0 || k >= sx->n[2]+1)
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    #endif

//  #ifdef BUG
//  if (i < 0 || i >= sx->n[0]+1) shit(sx->r);
//  if (j < 0 || j >= sx->n[1]+1) shit(sx->r);
//  if (k < 0 || k >= sx->n[2]+1) shit(sx->r);
//  #endif

    /* __ part. location in the cell __ */
    lx = xw-i;
    ly = yw-j;
    lz = zw-k;

    /* __ indexes of the rounding grid points on g2 __ */
    ijk1 = IDX(i  , j  , k  , n0, n1, n2);
    ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
    ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
    ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
    ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
    ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
    ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
    ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

    /* __ weight for each vertices of the rounding grid points __ */
    w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
    w2 = (    lx)*(1.0-ly)*(1.0-lz);
    w3 = (1.0-lx)*(    ly)*(1.0-lz);
    w4 = (    lx)*(    ly)*(1.0-lz);
    w5 = (1.0-lx)*(1.0-ly)*(    lz);
    w6 = (    lx)*(1.0-ly)*(    lz);
    w7 = (1.0-lx)*(    ly)*(    lz);
    w8 = (    lx)*(    ly)*(    lz);

    /* __ loop in the 3 directions __ */
    for (l = 0; l < 3; l++)
        {
        /* __ set e field seen by the part. __ */
        eo[l] = w1*s2[ijk1].e[l]
               +w2*s2[ijk2].e[l]
               +w3*s2[ijk3].e[l]
               +w4*s2[ijk4].e[l]
               +w5*s2[ijk5].e[l]
               +w6*s2[ijk6].e[l]
               +w7*s2[ijk7].e[l]
               +w8*s2[ijk8].e[l];

        /* __ set e field seen by the part. __ */
        jo[l] = w1*s2[ijk1].j[l]
               +w2*s2[ijk2].j[l]
               +w3*s2[ijk3].j[l]
               +w4*s2[ijk4].j[l]
               +w5*s2[ijk5].j[l]
               +w6*s2[ijk6].j[l]
               +w7*s2[ijk7].j[l]
               +w8*s2[ijk8].j[l];

        /* __ set e field seen by the part. __ */
        io[l] = w1*s2[ijk1].vi[l]
               +w2*s2[ijk2].vi[l]
               +w3*s2[ijk3].vi[l]
               +w4*s2[ijk4].vi[l]
               +w5*s2[ijk5].vi[l]
               +w6*s2[ijk6].vi[l]
               +w7*s2[ijk7].vi[l]
               +w8*s2[ijk8].vi[l];

        /* __ set position of the part. __ */
        ro[l] = (float) sp[so->s[n]][so->m[n]].r[l]+sp[so->s[n]][so->m[n]].b[l]*si.l[l];

        /* __ set velocity of the part. __ */
        wo[l] = (float) sp[so->s[n]][so->m[n]].v[l];
        }

    /* __ loop on the species __ */
    for (s = 0; s < NS+1; s++)
        {
        no[s] = w1*s2[ijk1].ns[s]
               +w2*s2[ijk2].ns[s]
               +w3*s2[ijk3].ns[s]
               +w4*s2[ijk4].ns[s]
               +w5*s2[ijk5].ns[s]
               +w6*s2[ijk6].ns[s]
               +w7*s2[ijk7].ns[s]
               +w8*s2[ijk8].ns[s];

        for (l = 0; l < 3; l++)
            {
            vo[s][l] = w1*s2[ijk1].vs[s][l]
                      +w2*s2[ijk2].vs[s][l]
                      +w3*s2[ijk3].vs[s][l]
                      +w4*s2[ijk4].vs[s][l]
                      +w5*s2[ijk5].vs[s][l]
                      +w6*s2[ijk6].vs[s][l]
                      +w7*s2[ijk7].vs[s][l]
                      +w8*s2[ijk8].vs[s][l];
            }

        for (l = 0; l < 6; l++)
            {
            po[s][l] = w1*s2[ijk1].ps[s][l]
                      +w2*s2[ijk2].ps[s][l]
                      +w3*s2[ijk3].ps[s][l]
                      +w4*s2[ijk4].ps[s][l]
                      +w5*s2[ijk5].ps[s][l]
                      +w6*s2[ijk6].ps[s][l]
                      +w7*s2[ijk7].ps[s][l]
                      +w8*s2[ijk8].ps[s][l];
            }
        }


    /* __ write the file __ */
    fwrite(&sw, sizeof(int), 1, fp);
    fwrite(&iw, sizeof(int), 1, fp);
    fwrite(ro, sizeof(float), 3, fp);
    fwrite(wo, sizeof(float), 3, fp);
    fwrite(bo, sizeof(float), 3, fp);
    fwrite(eo, sizeof(float), 3, fp);
    fwrite(jo, sizeof(float), 3, fp);
    fwrite(io, sizeof(float), 3, fp);
    for (s = 0; s < NS+1; s++) fwrite(&no[s], sizeof(float), 1, fp);
    for (s = 0; s < NS+1; s++) fwrite(vo[s], sizeof(float), 3, fp);
    for (s = 0; s < NS+1; s++) fwrite(po[s], sizeof(float), 6, fp);
    }

/* __ close the file __ */
fclose(fp);

/* __ last time step __ */
if (it == si.nt)
   {
   /* __ clean-up the pointers __ */
   free(so->so);
   free(so->io);
   free(so->s);
   free(so->m);
   }

}


/* _____ write the xplosed restart __________________________________________ */
void writerestart(struct sti si, struct stx sx, struct st1 *s1, struct st2 *s2, struct stp *sp[NS+1], struct std sd, struct stt st, int it)
{
int nn1, nn2;
int nm;
int s;
char nfile[10];
FILE *fp;


/* __ # of grid points on g1 & g2 __ */
nn1 = (sx.n[0]+1)*(sx.n[1]+1)*(sx.n[2]+1);
nn2 = (sx.n[0]+2)*(sx.n[1]+2)*(sx.n[2]+2);

/* __ # of modes __ */
nm = st.m[1]-st.m[0]+1;

/* __ print informations __ */
if (sx.r == 0)
   {
   printf("________________ write restart file for node  %3i ____\n", sx.r);
   printf("\n");
   }

/* __ make the file name __ */
if ((it/si.tr) % 2 == 0) sprintf(nfile, "hr%03i.dat", sx.r);
if ((it/si.tr) % 2 == 1) sprintf(nfile, "hR%03i.dat", sx.r);

/* __ open the file __ */
fp = fopen(nfile, "wb");
if (fp == NULL) printf("problem in opening file %s\n", nfile);

/* __ write the time step __ */
fwrite(&it, sizeof(int), 1, fp);

///* __ write the sti structure __ */
//fwrite(&si, sizeof(struct sti), 1, fp);

/* __ write ns & ws __ */
for (s = 1; s < NS+1; s++) {
   fwrite(&(si.ns[s]), sizeof(int), 1, fp);
   fwrite(&(si.ws[s]), sizeof(double), 1, fp);
}

/* __ write the stx structure __ */
fwrite(&sx, sizeof(struct stx), 1, fp);

/* __ write the st1 structure __ */
fwrite(s1, sizeof(struct st1), nn1, fp);

/* __ write the st2 structure __ */
fwrite(s2, sizeof(struct st2), nn2, fp);

/* __ write the stp structure __ */
for (s = 1; s < NS+1; s++) {
   fwrite(sp[s], sizeof(struct stp), sx.ns[s], fp);
}

/* __ write the initial total energy __ */
fwrite(&sd.e0, sizeof(double), 1, fp);

/* __ write the stt structure __ */
fwrite(st.bx, sizeof(double), nm, fp);
fwrite(st.by, sizeof(double), nm, fp);
//fwrite(st.ex, sizeof(double), nm, fp);
//fwrite(st.ey, sizeof(double), nm, fp);
//fwrite(st.ez, sizeof(double), nm, fp);
fwrite(st.px, sizeof(double), nm, fp);
fwrite(st.py, sizeof(double), nm, fp);
//fwrite(st.fx, sizeof(double), nm, fp);
//fwrite(st.fy, sizeof(double), nm, fp);
//fwrite(st.fz, sizeof(double), nm, fp);

/* __ close the file __ */
fclose(fp);

}



#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "structures.h"
#include "misc.h"
#include "memory.h"
#include "initModel.h"
#include "closeModel.h"




/*---------------------------------------------------------------------------
    readheckle()
  ---------------------------------------------------------------------------
    AIM : read the heckle.txt file and fills the parameter structure
 ---------------------------------------------------------------------------*/
void readheckle(struct sti *si, char *dir, int icpu)
{
    int ir;
    int g, h, l;
    char sfh[80];
    char junk[80];
    char modelname[80];
    char closename[80];
    FILE *fp;
    char *harrisfp    = "harrisfp";
    char *harris      = "harris";
    char *uniform     = "uniform";
    char *asymangle   = "recoasymangle";
    char *nbeams      = "Nbeams";
    char *shearb      = "shearB";
    char *nlasers     = "Nlasers";
    char *isotherm    = "isotherm";
    char *polytrop    = "polytrop";
    char *fullimpl    = "fullimpl";
    char *fullsub     = "fullsub";


    // unsused
    (void)icpu;

    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "heckle.txt");

    /* __ open the heckle.txt file __ */
    fp = fopen(sfh, "r+");
    if (fp == NULL) printf("\n\nproblem in opening file %s\n", sfh);

    /* reading the number of cells */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    for (h = 0; h < 3; h++) {
        fscanint(__FILE__, __LINE__, fp, 0, &(si->n[h]), &ir);
    }
    /* domain size in physical units */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    for (h = 0; h < 3; h++) {
        fscandbl(__FILE__, __LINE__, fp, 0, &(si->l[h]), &ir);
    }

    /* boundary conditions */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    for (h = 0; h < 3; h++) {
        fscanint(__FILE__, __LINE__, fp, 0, &(si->bc[h]), &ir);
    }

    /* time step and number of iterations */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(si->ts), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->nt), &ir);

    /* diagnostic time steps (fields, particles, info) */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->tf), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->tp), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->tt), &ir);

    /* max number of macroparticles per subdomain */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanint64(__FILE__, __LINE__, fp, 0, &(si->nm), &ir);

    /* open particle boundary and collisions */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->feed), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->coll), &ir);

    /* electron mass */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(si->me), &ir);

    /* electron thermal conduction */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(si->kapp), &ir);

    /* resistivity and hyper-resistivity */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(si->rsty), &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(si->hyvi), &ir);

    /* restart flag and time */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->rst), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->tr), &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(si->time4rst), &ir);

    /* number of particle per species (protons, alphas) */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    for (g = 1; g < NS+1; g++) {
        fscanint64(__FILE__, __LINE__, fp, 0, &(si->ns[g]), &ir);
    }

    /* mass of the ion species */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    for (g = 1; g < NS+1; g++) {
        fscandbl(__FILE__, __LINE__, fp, 0, &(si->ms[g]), &ir);
    }

    /* charge of ion species */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    for (g = 1; g < NS+1; g++) {
        fscandbl(__FILE__, __LINE__, fp, 0, &(si->qs[g]), &ir);
    }


    /* initial model */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanstr(__FILE__, __LINE__, fp, 0, modelname, &ir);

    if (strcmp(harris,modelname) == 0)
    {
        si->InitModelID = INITMODEL_HARRIS;
        si->drive = 0;
    }
    else if (strcmp(harrisfp,modelname) == 0)
    {
        si->InitModelID = INITMODEL_DOUBLEHARRIS;
        si->drive = 0;
    }
    else if (strcmp(uniform,modelname) == 0)
    {
        si->InitModelID = INITMODEL_UNIFORM;
        si->drive = 0;
    }
    else if (strcmp(asymangle,modelname) == 0)
    {
        si->InitModelID = INITMODEL_ASYMANGLE;
        si->drive = 0;
    }
    else if (strcmp(nbeams,modelname) == 0)
    {
        si->InitModelID = INITMODEL_NBEAMS;
        si->drive = 0;
    }
    else if (strcmp(shearb,modelname) == 0)
    {
        si->InitModelID = INITMODEL_SHEARB;
        si->drive = 0;
    }
    else if (strcmp(nlasers,modelname) == 0)
    {
        si->InitModelID = INITMODEL_NLASERS;
        si->drive = 1;
    }
    else
    {
        printf("Error - Init Model not implemented (%s)\n", modelname);
    }

    /* closure model */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanstr(__FILE__, __LINE__, fp, 0, closename, &ir);

    if (strcmp(isotherm, closename) == 0)
    {
        si->CloseModelID = CLOSE_ISOTHERM;
    }
    else if (strcmp(polytrop, closename) == 0)
    {
        si->CloseModelID = CLOSE_POLYTROP;
    }
    else if (strcmp(fullimpl, closename) == 0)
    {
        si->CloseModelID = CLOSE_FULLIMPL;
    }
    else if (strcmp(fullsub, closename) == 0)
    {
        si->CloseModelID = CLOSE_FULLSUB;
    }
    else
    {
        printf("Error - Closure Model not implemented (%s)\n", closename);
    }


    /* number of MPI domains in each direction */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->mpidom[0]), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->mpidom[1]), &ir);
    fscanint(__FILE__, __LINE__, fp, 0, &(si->mpidom[2]), &ir);

    /* __ close the file __ */
    fclose(fp);

    /* __ set grid steps __ */
    for (l = 0; l < 3; l++)
    {
       si->dl[l] = si->l[l]/(double)si->n[l];
    }

}
/*===========================================================================*/









/* _____ read the horbi.txt file ____________________________________________ */
void readhorbi(struct sto *so)
{
int ir;
int m;
char sfo[80];
char junk[80];
FILE *fp;


/* __ build the file nemae __ */
strcpy(sfo, "horbi.txt");

/* __ open the horbi.txt file __ */
fp = fopen(sfo, "rb");
if (fp == NULL) printf("\n\nproblem in opening file %s\n", sfo);

//if ((z = fscanf(fp, "%s", junk)) != 1) CRASH(0);
//if ((z = fscanf(fp, "%d", &(so->wo))) != 1) CRASH(0);
//if ((z = fscanf(fp, "%s", junk)) != 1) CRASH(0);
fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
fscanint(__FILE__, __LINE__, fp, 0, &(so->wo), &ir);
fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);

/* __ memory allocation __ */
so->so = (int *)malloc(so->wo*sizeof(int));
so->io = (int *)malloc(so->wo*sizeof(int));
so->s = (int *)malloc(so->wo*sizeof(int));
so->m = (int *)malloc(so->wo*sizeof(int));

/* __ loop on the part. to track __ */
for (m = 0; m < so->wo; m++) {
   /* __ read the part. i index & s specie to track __ */
//  if ((z = fscanf(fp, "%d", (so->so)+m)) != 1) CRASH(0);
//  if ((z = fscanf(fp, "%d", (so->io)+m)) != 1) CRASH(0);
   fscanint(__FILE__, __LINE__, fp, 0, (so->so)+m, &ir);
   fscanint(__FILE__, __LINE__, fp, 0, (so->io)+m, &ir);
   }

/* __ close the file __ */
fclose(fp);

}


/* __ read the xplosed field files __________________________________________ */
void readfield(struct sti si, struct stf *sf, int r, int f)
{
float *bw[3], *ew[3], *jw[3], *iw[3];
float *nw[NS+1], *vw[NS+1][3], *pw[NS+1][6];
int ijk, ijkw;
int nn2;
int i, j, k, l, s;
int n0, n1, n2;
int ir;
struct stx sx;
char sff[80];
char nfile[16];
FILE *fp;


/* __ build the name of the field file __ */
sprintf(nfile, "xf%03i-%04i.dat", r, f);

/* __ build the file nemae __ */
strcpy(sff,nfile);

/* __ open the file __ */
fp = fopen(sff, "rb");
if (fp == NULL) printf("\n\nproblem in opening file %s\n", sff);

/* __ read the data __ */
//z = fread(sx.i0, sizeof(int), 3, fp); if (z != 3) shit(0);
//z = fread(sx.n, sizeof(int), 3, fp); if (z != 3) shit(0);
freadint(__FILE__, __LINE__, fp, r, 3, sx.i0, &ir);
freadint(__FILE__, __LINE__, fp, r, 3, sx.n, &ir);

/* __ #of grid points __ */
n0 = sx.n[0]+1;
n1 = sx.n[1]+1;
n2 = sx.n[2]+1;

/* __ # of points for arrays __ */
nn2 = (sx.n[0]+1)*(sx.n[1]+1)*(sx.n[2]+1);

/* __ memory allocation for electromagnetic & fluid fields __ */
for (l = 0; l < 3; l++) {
   bw[l] = (float *) malloc(nn2*sizeof(float));
   ew[l] = (float *) malloc(nn2*sizeof(float));
   jw[l] = (float *) malloc(nn2*sizeof(float));
   iw[l] = (float *) malloc(nn2*sizeof(float));
   }

/* __ memory allocation for specie moments __ */
for (s = 0; s < NS+1; s++) {
   /* __ specie density __ */
   nw[s] = (float *) malloc(nn2*sizeof(float));

   /* __ specie velocity __ */
   for (l = 0; l < 3; l++) {
      vw[s][l] = (float *) malloc(nn2*sizeof(float));
      }

   /* __ specie pressure __ */
   for (l = 0; l < 6; l++) {
      pw[s][l] = (float *) malloc(nn2*sizeof(float));
      }
   }

/* __ read the data : e-m & fluid fields __ */
for (l = 0; l < 3; l++) {
//  z = fread(bw[l], sizeof(float), nn2, fp); if (z != nn2) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn2, bw[l], &ir);
   }

/* __ read the data : e-m & fluid fields __ */
for (l = 0; l < 3; l++) {
//  z = fread(ew[l], sizeof(float), nn2, fp); if (z != nn2) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn2, ew[l], &ir);
   }

/* __ read the data : e-m & fluid fields __ */
for (l = 0; l < 3; l++) {
//  z = fread(jw[l], sizeof(float), nn2, fp); if (z != nn2) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn2, jw[l], &ir);
   }

/* __ read the data : e-m & fluid fields __ */
for (l = 0; l < 3; l++) {
//  z = fread(iw[l], sizeof(float), nn2, fp); if (z != nn2) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn2, iw[l], &ir);
   }

/* __ read the data : specie density __ */
for (s = 0; s < NS+1; s++) {
//  z = fread(nw[s], sizeof(float), nn2, fp); if (z != nn2) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn2, nw[s], &ir);
   }

/* __ read the data : specie velocity __ */
for (s = 0; s < NS+1; s++) {
    for (l = 0; l < 3; l++) {
//      z = fread(vw[s][l], sizeof(float), nn2, fp); if (z != nn2) shit(0);
      freadflt(__FILE__, __LINE__, fp, 0, nn2, vw[s][l], &ir);
      }
   }

/* __ read the data : specie pressure __ */
for (s = 0; s < NS+1; s++) {
    for (l = 0; l < 6; l++) {
//      z = fread(pw[s][l], sizeof(float), nn2, fp); if (z != nn2) shit(0);
      freadflt(__FILE__, __LINE__, fp, 0, nn2, pw[s][l], &ir);
      }
   }

/* __ close the file __ */
fclose(fp);

/* __ nested loops on the subdomain __ */
for (i = 0; i < sx.n[0]+1; i++) {
   for (j = 0; j < sx.n[1]+1; j++) {
      for (k = 0; k < sx.n[2]+1; k++) {
         /* __ set the index on large domain __ */
         ijk = IDX(i+sx.i0[0], j+sx.i0[1], k+sx.i0[2], si.n[0]+1, si.n[1]+1, si.n[2]+1);

         /* __ set the index on the subdomain __ */
         ijkw = IDX(i, j, k, n0, n1, n2);

         /* __ fill the main arrays __ */
         for (l = 0; l < 3; l++) {
            sf->b[l][ijk] = bw[l][ijkw];
            sf->e[l][ijk] = ew[l][ijkw];
            sf->j[l][ijk] = jw[l][ijkw];
            sf->vi[l][ijk] = iw[l][ijkw];
            }

         for (s = 0; s < NS+1; s++) {
            sf->n[s][ijk] = nw[s][ijkw];

            for (l = 0; l < 3; l++) {
               sf->v[s][l][ijk] = vw[s][l][ijkw];
               }

             for (l = 0; l < 6; l++) {
                sf->p[s][l][ijk] = pw[s][l][ijkw];
                }
            }
         }
      }
   }

/* __ clean-up the pointers : e-m & fluid fields __ */
for (l = 0; l < 3; l++) {
   free(bw[l]);
   free(ew[l]);
   free(jw[l]);
   free(iw[l]);
   }

/* __ clean-up the pointers : specie moments __ */
for (s = 0; s < NS+1; s++) {
   free(nw[s]);

   for (l = 0; l < 3; l++) {
      free(vw[s][l]);
      }

   for (l = 0; l < 6; l++) {
       free(pw[s][l]);
       }
   }

}


///* __ read the xplosed particles files ______________________________________ */
//void readpart(struct sti *si, struct stq *sq, int r, int p)
//{
//float *xw[NS+1], *yw[NS+1], *zw[NS+1];
//float *uw[NS+1], *vw[NS+1], *ww[NS+1];
//int *iw[NS+1];
//int ir;
//int m, s;
//struct stx sx;
//char sfp[80];
//char nfile[16];
//FILE *fp;
//
//
///* __ build the name of the particles file __ */
//sprintf(nfile, "xp%03i-%04i.dat", r, p);
//
///* __ build the file nemae __ */
//strcpy(sfp,nfile);
//
///* __ open the file __ */
//fp = fopen(sfp, "rb");
//if (fp == NULL) printf("\n\nproblem in opening file %s\n", sfp);
//
///* __ read the data __ */
////z = fread(sx.ns, sizeof(int), NS+1, fp); if (z != NS+1) shit(0);
////z = fread(&(si->ns), sizeof(int), NS+1, fp); if (z != NS+1) shit(0);
//freadint(__FILE__, __LINE__, fp, 0, NS+1, sx.ns, &ir);
//freadint(__FILE__, __LINE__, fp, 0, NS+1, si->ns, &ir);
//
///* __ loop on the part. __ */
//for (s = 1; s < NS+1; s++) {
//   /* __ memory allocation __ */
//   iw[s] = (int *) malloc(sx.ns[s]*sizeof(int));
//
//   xw[s] = (float *) malloc(sx.ns[s]*sizeof(float));
//   yw[s] = (float *) malloc(sx.ns[s]*sizeof(float));
//   zw[s] = (float *) malloc(sx.ns[s]*sizeof(float));
//
//   uw[s] = (float *) malloc(sx.ns[s]*sizeof(float));
//   vw[s] = (float *) malloc(sx.ns[s]*sizeof(float));
//   ww[s] = (float *) malloc(sx.ns[s]*sizeof(float));
//
//   /* __ read the data __ */
////  z = fread(iw[s], sizeof(int), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
//   freadint(__FILE__, __LINE__, fp, 0, sx.ns[s], iw[s], &ir);
//
////  z = fread(xw[s], sizeof(float), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
////  z = fread(yw[s], sizeof(float), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
////  z = fread(zw[s], sizeof(float), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
//   freadflt(__FILE__, __LINE__, fp, 0, sx.ns[s], xw[s], &ir);
//   freadflt(__FILE__, __LINE__, fp, 0, sx.ns[s], yw[s], &ir);
//   freadflt(__FILE__, __LINE__, fp, 0, sx.ns[s], zw[s], &ir);
//
////  z = fread(uw[s], sizeof(float), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
////  z = fread(vw[s], sizeof(float), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
////  z = fread(ww[s], sizeof(float), sx.ns[s], fp); if (z != sx.ns[s]) shit(0);
//   freadflt(__FILE__, __LINE__, fp, 0, sx.ns[s], uw[s], &ir);
//   freadflt(__FILE__, __LINE__, fp, 0, sx.ns[s], vw[s], &ir);
//   freadflt(__FILE__, __LINE__, fp, 0, sx.ns[s], ww[s], &ir);
//
//   /* __ loop on the part : fill the stq struct __ */
//   for (m = 0; m < sx.ns[s]; m++) {
//      /* __ fill the main arrays __ */
//      sq->r[s][0][iw[s][m]] = xw[s][m];
//      sq->r[s][1][iw[s][m]] = yw[s][m];
//      sq->r[s][2][iw[s][m]] = zw[s][m];
//
//      sq->v[s][0][iw[s][m]] = uw[s][m];
//      sq->v[s][1][iw[s][m]] = vw[s][m];
//      sq->v[s][2][iw[s][m]] = ww[s][m];
//      }
//
//   /* __ clean-up the pointers __ */
//   free(iw[s]);
//
//   free(xw[s]);
//   free(yw[s]);
//   free(zw[s]);
//
//   free(uw[s]);
//   free(vw[s]);
//   free(ww[s]);
//   }
//
///* __ close the file __ */
//fclose(fp);
//
//}


/* __ read the xplosed orbits files _________________________________________ */
void readorbit(struct sto so, struct stm *sm, int r, int o)
{
float ro[3], wo[3], bo[3], eo[3], jo[3], io[3];
float no[NS+1], vo[NS+1][3], po[NS+1][6];
int ir;
int sw, iw;
int nmap, wi;
int l, m, n, s;
char sfm[80];
char nfile[16];
FILE *fp;


/* __ build the name of the particles file __ */
sprintf(nfile, "xo%03i-%04i.dat", r, o);

/* __ build the file nemae __ */
strcpy(sfm, nfile);

/* __ open the file __ */
fp = fopen(sfm, "rb");
if (fp == NULL) printf("\n\nproblem in opening file %s\n", sfm);

/* __ read the data __ */
//z = fread(&nmap, sizeof(int), 1, fp); if (z != 1) shit(0);
//z = fread(&wi, sizeof(int), 1, fp); if (z != 1) shit(0);
freadint(__FILE__, __LINE__, fp, 0, 1, &nmap, &ir);
freadint(__FILE__, __LINE__, fp, 0, 1, &wi, &ir);

/* __ loop on the # of orbits written in the file __ */
for (m = 0; m < nmap; m++) {
   /* __ read the data __ */
//  z = fread(&sw, sizeof(int), 1, fp); if (z != 1) shit(0);
//  z = fread(&iw, sizeof(int), 1, fp); if (z != 1) shit(0);
   freadint(__FILE__, __LINE__, fp, 0, 1, &sw, &ir);
   freadint(__FILE__, __LINE__, fp, 0, 1, &iw, &ir);

    /* __ read the data __ */
//  z = fread(ro, sizeof(float), 3, fp); if (z != 3) shit(0);
//  z = fread(wo, sizeof(float), 3, fp); if (z != 3) shit(0);
//  z = fread(bo, sizeof(float), 3, fp); if (z != 3) shit(0);
//  z = fread(eo, sizeof(float), 3, fp); if (z != 3) shit(0);
//  z = fread(jo, sizeof(float), 3, fp); if (z != 3) shit(0);
//  z = fread(io, sizeof(float), 3, fp); if (z != 3) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, 3, ro, &ir);
   freadflt(__FILE__, __LINE__, fp, 0, 3, wo, &ir);
   freadflt(__FILE__, __LINE__, fp, 0, 3, bo, &ir);
   freadflt(__FILE__, __LINE__, fp, 0, 3, eo, &ir);
   freadflt(__FILE__, __LINE__, fp, 0, 3, jo, &ir);
   freadflt(__FILE__, __LINE__, fp, 0, 3, io, &ir);

   for (s = 0; s < NS+1; s++) {
//      z = fread(&no[s], sizeof(float), 1, fp); if (z != 1) shit(0);
      freadflt(__FILE__, __LINE__, fp, 0, 1, &no[s], &ir);
   }

   for (s = 0; s < NS+1; s++) {
//      z = fread(vo[s], sizeof(float), 3, fp); if (z != 3) shit(0);
      freadflt(__FILE__, __LINE__, fp, 0, 3, vo[s], &ir);
   }

   for (s = 0; s < NS+1; s++) {
//      z = fread(po[s], sizeof(float), 6, fp); if (z != 6) shit(0);
      freadflt(__FILE__, __LINE__, fp, 0, 6, po[s], &ir);
   }

   /* __ find the # of the map : loop on the total # of orbits __ */
   for (n = 0; n < so.wo; n++) {
      /* __ select the right part. __ */
      if (sw == so.so[n] && iw == so.io[n]) break;
   }

    /* __ fill the sm structure __ */
   for (l = 0; l < 3; l++) {
      sm[n].r[l][wi] = ro[l];
      sm[n].w[l][wi] = wo[l];
      sm[n].b[l][wi] = bo[l];
      sm[n].e[l][wi] = eo[l];
      sm[n].j[l][wi] = jo[l];
      sm[n].vi[l][wi] = io[l];
   }

   /* __ fill the sm structure __ */
   for (s = 0; s < NS+1; s++) {
      sm[n].n[s][wi] = no[s];

      for (l = 0; l < 3; l++) {
         sm[n].v[s][l][wi] = vo[s][l];
         }

      for (l = 0; l < 6; l++) {
         sm[n].p[s][l][wi] = po[s][l];
      }
   }
}

/* __ close the file __ */
fclose(fp);

}


/* __ read the full field file ______________________________________________ */
void readfullfield(struct sti si, struct stf *sf, char *sff)
{
int nn1;
int ir;
int l, s;
FILE *fp;


/* __ # of grid points __ */
nn1 = (si.n[0]+1)*(si.n[1]+1)*(si.n[2]+1);

/* __ open the file __ */
fp = fopen(sff, "rb");
if (fp == NULL) printf("\n\nproblem in opening file %s\n", sff);

/* __ read the hf___.dat file __ */
for (l = 0; l < 3; l++) {
//  z = fread(sf->b[l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->b[l], &ir);
}

for (l = 0; l < 3; l++) {
//  z = fread(sf->a[l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->a[l], &ir);
}

for (l = 0; l < 3; l++) {
//  z = fread(sf->e[l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->e[l], &ir);
}

for (l = 0; l < 3; l++) {
//  z = fread(sf->j[l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->j[l], &ir);
}

for (l = 0; l < 3; l++) {
//  z = fread(sf->vi[l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->vi[l], &ir);
}

for (s = 0; s < NS+1; s++) {
//  z = fread(sf->n[s], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->n[s], &ir);
}

for (s = 0; s < NS+1; s++) {
   for (l = 0; l < 3; l++) {
//    z = fread(sf->v[s][l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->v[s][l], &ir);
   }
}

for (s = 0; s < NS+1; s++) {
    for (l = 0; l < 6; l++) {
//      z = fread(sf->p[s][l], sizeof(float), nn1, fp); if (z != nn1) shit(0);
   freadflt(__FILE__, __LINE__, fp, 0, nn1, sf->p[s][l], &ir);
    }
}

/* __ close the file __ */
fclose(fp);

}


///* __ read the full part. file ______________________________________________ */
//void readfullpart(struct sti *si, struct stq *sq, char *sfp)
//{
//int ir;
//int s;
//FILE *fp;
//
//
///* __ open the file __ */
//fp = fopen(sfp, "rb");
//if (fp == NULL) printf("\n\nproblem in opening file %s\n", sfp);
//
///* __ read the part. in the hp file : loop on the specie __ */
//for (s = 1; s < NS+1; s++) {
//   /* __ read the part. in the hp files __ */
////  z = fread(&(si->ns[s]), sizeof(int), 1, fp); if (z != 1) shit(0);
//   freadint(__FILE__, __LINE__, fp, 0, 1, &(si->ns[s]), &ir);
//
//   /* __ read if there are part. __ */
//   if (si->ns[s] != 0) {
////     z = fread(sq->r[s][0], sizeof(float), si->ns[s], fp); if (z != si->ns[s]) shit(0);
////     z = fread(sq->r[s][1], sizeof(float), si->ns[s], fp); if (z != si->ns[s]) shit(0);
////     z = fread(sq->r[s][2], sizeof(float), si->ns[s], fp); if (z != si->ns[s]) shit(0);
//      freadflt(__FILE__, __LINE__, fp, 0, si->ns[s], sq->r[s][0], &ir);
//      freadflt(__FILE__, __LINE__, fp, 0, si->ns[s], sq->r[s][1], &ir);
//      freadflt(__FILE__, __LINE__, fp, 0, si->ns[s], sq->r[s][2], &ir);
//
////     z = fread(sq->v[s][0], sizeof(float), si->ns[s], fp); if (z != si->ns[s]) shit(0);
////     z = fread(sq->v[s][1], sizeof(float), si->ns[s], fp); if (z != si->ns[s]) shit(0);
////     z = fread(sq->v[s][2], sizeof(float), si->ns[s], fp); if (z != si->ns[s]) shit(0);
//      freadflt(__FILE__, __LINE__, fp, 0, si->ns[s], sq->v[s][0], &ir);
//      freadflt(__FILE__, __LINE__, fp, 0, si->ns[s], sq->v[s][1], &ir);
//      freadflt(__FILE__, __LINE__, fp, 0, si->ns[s], sq->v[s][2], &ir);
//   }
//}
//
///* __ close the file __ */
//fclose(fp);
//
//}


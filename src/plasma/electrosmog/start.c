
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "plasma.h"
#include "check.h"
#include "read.h"
#include "init.h"
#include "ions.h"
#include "subdomains.h"
#include "write.h"


/* __ use a restart file as initial condition _______________________________ */
void start(struct sti *si, struct stt *st, struct stx *sx, struct st1 **s1, struct st2 **s2, struct stp *(*sp)[NS+1], struct std *sd, struct sto *so, int *it, MPI_Comm *com, clock_t ticks[2], time_t ttime[2])
{
int nn1, nn2;
int nm;
int s;
char nfile[10];
FILE *fp;


/* __  __ */
ticks[0] = clock();
ttime[0] = time(NULL);

/* __ if no need of a restart __ */
if (si->rst == 0) {
   /* __ define the subdomains __ */
   subdomains(*si, sx, com);

   /* __ # of grid points __ */
   nn1 = (sx->n[0]+1)*(sx->n[1]+1)*(sx->n[2]+1);
   nn2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);

   /* __ memory allocation __ */
   *s1 = (struct st1 *)malloc(nn1*sizeof(struct st1));
   *s2 = (struct st2 *)malloc(nn2*sizeof(struct st2));

   /* __ memory allocation : loop on the species __ */
   for (s = 1; s < NS+1; s++) {
      (*sp)[s] = (struct stp *)malloc(si->nm*sizeof(struct stp));
   }

   /* __ set the initial drive __ */
   drive(*si, *sx, st, *s1, *s2, *com);

   /* __ initialization __ */
   init(si, sx, st, *s1, *s2, *sp, it, *com);
}

/* __ if start from a restart __ */
if (si->rst != 0) {
   /* __ print informations __ */
   printf("________________ read restart file for node   %3i ____\n", sx->r);
   printf("\n");

   /* __ make the file name __ */
   if (si->rst == 1) sprintf(nfile, "hr%03i.dat", sx->r);
   if (si->rst == 2) sprintf(nfile, "hR%03i.dat", sx->r);

   /* __ open the file __ */
   fp = fopen(nfile, "rb");
   if (fp == NULL) printf("problem in opening file %s\n", nfile);

   /* __ set it value __ */
   freadint(__FILE__, __LINE__, fp, sx->r, 1, it, &(sx->irun));

   /* __ read ns & ws __ */
   for (s = 1; s < NS+1; s++) {
      freadint(__FILE__, __LINE__, fp, sx->r, 1, &(si->ns[s]), &(sx->irun));
      freaddbl(__FILE__, __LINE__, fp, sx->r, 1, &(si->ws[s]), &(sx->irun));
   }

   /* __ read the stx structure __ */
   freadstx(__FILE__, __LINE__, fp, sx->r, 1, sx, &(sx->irun));

   /* __ # of points for s1 & s2 __ */
   nn1 = (sx->n[0]+1)*(sx->n[1]+1)*(sx->n[2]+1);
   nn2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);

   /* __ # of modes __ */
   nm = st->m[1]-st->m[0]+1;

   /* __ memory allocation __ */
   *s1 = (struct st1 *)malloc(nn1*sizeof(struct st1));
   *s2 = (struct st2 *)malloc(nn2*sizeof(struct st2));

   /* __ memory allocation : loop on the species __ */
   for (s = 1; s < NS+1; s++) {
      (*sp)[s] = (struct stp *)malloc(si->nm*sizeof(struct stp));
   }

   /* __ read the st1 structure __ */
   freadst1(__FILE__, __LINE__, fp, sx->r, nn1, *s1, &(sx->irun));

   /* __ read the st2 structure __ */
   freadst2(__FILE__, __LINE__, fp, sx->r, nn2, *s2, &(sx->irun));

   /* __ read the stp structure __ */
   for (s = 1; s < NS+1; s++) {
      freadstp(__FILE__, __LINE__, fp, sx->r, sx->ns[s], (*sp)[s], &(sx->irun));
   }

   /* __ set initial total energy __ */
   freaddbl(__FILE__, __LINE__, fp, sx->r, 1, &((*sd).e0), &(sx->irun));

   /* __ read the stt structure : amplitude of b & e __ */
   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->bx, &(sx->irun));
   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->by, &(sx->irun));
//   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->ex, &(sx->irun));
//   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->ey, &(sx->irun));
//   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->ez, &(sx->irun));
   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->px, &(sx->irun));
   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->py, &(sx->irun));
//   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->fx, &(sx->irun));
//   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->fy, &(sx->irun));
//   freaddbl(__FILE__, __LINE__, fp, sx->r, nm, st->fz, &(sx->irun));

   /* __ close the file __ */
   fclose(fp);
}

/* __ set the full ion pressure tensor __ */
ions(*si, sx, *s1, *s2, *sp, *it, *com);

/* __ diagnostic __ */
check(*si, *sx, *s1, *s2, *sp, sd, *it, *com, ticks, ttime);
writedump(*si, *sx, *sd, *it);
writefield(*si, *sx, *s1, *s2, *it);
writepart(*si, *sx, *sp, *it);
writeorbit(*si, sx, *s1, *s2, *sp, so, *it);

}


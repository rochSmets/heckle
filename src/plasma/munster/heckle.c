
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "structures.h"
#include "plasma.h"
#include "add.h"
#include "check.h"
#include "close.h"
#include "defines.h"
#include "maxwell.h"
#include "init.h"
#include "ions.h"
#include "ohm.h"
#include "read.h"
#include "setup.h"
#include "sort.h"
#include "start.h"
#include "subdomains.h"
#include "sources.h"
#include "HeckleIO.h"


/* __ main __________________________________________________________________ */
int main(int ac, char **av)
{
struct sti si;
struct stt st;
struct st1 *s1;
struct st2 *s2;
struct stp *sp[NS+1];
struct sto so;
struct std sd;
struct stx sx;
MPI_Comm com;
int it;
int ipc;
int s;
int irun;
clock_t ticks[2];
time_t ttime[2], timer;
char fulltime[24];


/* __ mpi set-up __ */
MPI_Init(&ac, &av);
MPI_Comm_size(MPI_COMM_WORLD, &(sx.s));
MPI_Comm_rank(MPI_COMM_WORLD, &(sx.r));
MPI_Comm_dup(MPI_COMM_WORLD, &com);

/* __ set the run parameters __ */
si = setup(sx, av[1]);

/* __ time & date of the begining of the run __ */
if (sx.r == 0)
   {
   timer = time(NULL);
   strncpy(fulltime, asctime(localtime(&timer)), 24);
   printf("________________ %s ____________\n",  fulltime);
   printf("\n");
   }

/* __ set the initial topology __ */
topo(si, &sx, &st, com);

/* __ setup the so struct. __ */
readhorbi(&so);

/* __ heckle running from a restart ? __ */
start(&si, &st, &sx, &s1, &s2, &sp, &sd, &so, &it, &com, ticks, ttime);

/* __ loop on the total # of time step __ */
while (it < si.nt)
      {
      /* __ predictor __ */
      ipc = 0;

      /* __ drive the magnetic fluctuations __ */
      if (it < st.td) drive(si, sx, &st, s1, s2, com);

      /* __ integration of part. moments __ */
      sources(&si, &sx, &st, s1, s2, sp, ipc, it, com);

      /* __ build the full ion kinetic pressure tensor __ */
      if (it%si.tf == 0) ions(si, &sx, s1, s2, sp, it, com);

      /* __ half rotation with faraday __ */
      faraday(si, sx, s1, s2, ipc, it, com);

      /* __ set electron stress tensor __ */
      stress(si, sx, s2, ipc, com);

      /* __ set e field with ohm's law __ */
      ohm(si, sx, s1, s2, ipc, it, com);

      /* __ tricky : needed to distinguish the first faraday half rotation __ */
      ipc = 2;

      /* __ half rotation with faraday __ */
      faraday(si, sx, s1, s2, ipc, it, com);

      /* __ corrector __ */
      ipc = 1;

      /* __ integration of part. moments __ */
      sources(&si, &sx, &st, s1, s2, sp, ipc, it, com);

      /* __ half rotation with faraday __ */
      faraday(si, sx, s1, s2, ipc, it, com);

      /* __ set electron stress tensor __ */
      stress(si, sx, s2, ipc, com);

      /* __ set e field with ohm's law __ */
      ohm(si, sx, s1, s2, ipc, it, com);

      /* __ tricky : needed to distinguish the first faraday half rotation __ */
      ipc = 2;

      /* __ half rotation with faraday __ */
      faraday(si, sx, s1, s2, ipc, it, com);

      /* __ add part. if the density is too low __ */
      if (si.feed == 1) add(&si, &sx, &st, s2, sp, com);

      /* __ sort the part __ */
      if (it%NSORT == 0) sort(&si, &sx, sp);

      /* __ increase time step __ */
      it ++;

      /* __ diagnostic __ */
      if (it%si.tt == 0) check(si, sx, s1, s2, sp, &sd, it, com, ticks, ttime);
      if (it%si.tt == 0) writedump(si, sx, sd, it);
      if (it%si.tf == 0) writefield(si, sx, s1, s2, it);
      if (it%si.tp == 0) writepart(si, sx, sp, it);
      if (it%si.tt == 0) writeorbit(si, &sx, s1, s2, sp, &so, it);
      if (it%si.tr == 0) writerestart(si, sx, s1, s2, sp, sd, st, it);

      /* __ reduce the irun index __ */
      MPI_Allreduce(&(sx.irun), &irun, 1, MPI_INT, MPI_SUM, com);

      /* __ if the run has to be stopped __ */
      if (irun != 0)
         {
         /* __ print fancy stuff __ */
         if (sx.r == 0)
            {
            timer = time(NULL);
            strncpy(fulltime, asctime(localtime(&timer)), 24);
            printf("________________ %s ____________\n",  fulltime);
            printf("\n");
            printf("________________ this is the end... __________________\n");
            printf("\n");
            printf("________________ run is dead : err. # %4d ___________\n", irun);
            printf("\n");
            }

         /* __ get out of the while __ */
         break;

         /* __ mpi termination __ */
         MPI_Comm_free(&com);
         MPI_Finalize();

         return 0;
         }
      }

/* __ time & date of the end of the run __ */
if (sx.r == 0)
   {
   timer = time(NULL);
   strncpy(fulltime, asctime(localtime(&timer)), 24);
   printf("________________ %s ____________\n",  fulltime);
   printf("\n");
   printf("________________ this is the end... __________________\n");
   printf("\n");
   }

/* __ clean-up the pointers __ */
free(s1);
free(s2);

for (s = 1; s < NS+1; s++)
    {
    free(sp[s]);
    }

/* __ mpi termination __ */
MPI_Comm_free(&com);
MPI_Finalize();

return 0;

}


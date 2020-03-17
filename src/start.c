
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
//#include "plasma.h"
#include "check.h"
#include "read.h"
#include "init.h"
#include "ions.h"
#include "subdomains.h"

#include <HeckleIOFields.h>
#include <HeckleIOSpecies.h>
#include <HeckleIORestarts.h>
#include <HeckleIOTime.h>
#include <particle.h>
#include <hecklebc.h>
#include <ghosts.h>
#include <HeckleMem.h>
#include <string.h>




/*---------------------------------------------------------------------------
  start()
  ---------------------------------------------------------------------------
  AIM : start a simulation. Given the input parameters, the routine will

          - do the domain decomposition
          - allocate the memory for the fields and particles
          - initialize electromagnetic fields and arrays

        the routine can start a simulation from scratch or restart a previous
        one.
 ---------------------------------------------------------------------------*/

void start(HeckleIOFields **hiof,
           HeckleIOSpecies **hios,
           HeckleIORestart **hior,
           HeckleIOTime **hiot,
           STI *si,
           STX *sx,
           Collision *sc,
           Grid0 **s0,
           struct st1 **s1,
           ST2 **s2,
           Particle *(*sp)[NS+1],
           HeckleBC **hbc,
           Ghosts **ghost,
           struct std *sd,
           struct sto *so,
           int *it,
           clock_t ticks[2],
           time_t ttime[2],
           int argc, char **argv)
{
    int nn0, nn1, nn2;
    int s;
    int ijk0;
//  char nfile[10];
//  FILE *fp;
    int iarg;


    ticks[0] = clock();
    ttime[0] = time(NULL);

    memory_start();

    /* --------------------------------------------------------------- */
    /*                  START RUNNING FROM SCRATCH                     */
    /* --------------------------------------------------------------- */


    if (si->rst == 0)
    {

        /* domain decomposition */
        subdomains(*si, sx);



        /* _____ if starting just for memory check  ______ */

        for (iarg =0; iarg < argc; iarg++)
        {
            if (strcmp(argv[iarg], "memcheck") == 0)
            {
                memory_estimate(si, sx);
                exit(-1);
                break;
            }
        }
        /* ________________________________________________ */




        /*# of points and memory allocation for the grids g1 and g2 */
        nn0 = (sx->n[0]  )*(sx->n[1]  )*(sx->n[2]  );
        nn1 = (sx->n[0]+1)*(sx->n[1]+1)*(sx->n[2]+1);
        nn2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);

        *s0 = memory_allocate(nn0 * sizeof **s0);
        *s1 = memory_allocate(nn1 * sizeof **s1);
        *s2 = memory_allocate(nn2 * sizeof **s2);

        /* memory allocation : loop on the species */
        for (s = 1; s < NS+1; s++) {
            //(*sp)[s] = malloc(si->nm * sizeof (Particle));
            (*sp)[s] = memory_allocate(sizeof(Particle) * si->nm);

            if ((*sp)[s] == NULL) {
                printf("Error memory allocation on proc %d at %s %d\n", sx->r, __FILE__, __LINE__);
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }

        }

        /* set to zero the # of part in each cell */
        if (si->coll != 0) {
            for (ijk0 = 0; ijk0 < nn0; ijk0++) {
                for (s = 1; s < NS+1; s++) {
                    (*s0)[ijk0].npart[s] = 0;
                }
            }
        }

        /* initialize the boundary condition module */
        *hbc = HeckleBCInit(si);

        // initialize the ghost point communication module
        *ghost  = GhostsInit(si, sx);

        memory_status();

        /* now we set initial values for fields and particles */
        init(si, sx, *s0, *s1, *s2, *sp, *hbc, *ghost, it);

        if(sx->r == 0) {
            printf("\n________________ particles weights ___________________\n");
            printf("\n");
            for (s = 1; s < NS+1; s++) {
                printf("... weight [%1d]   :%12.6lf\n", s, si->ws[s]);
            }
        printf("\n______________________________________________________\n");
        printf("\n");
        }

        /* start the IO module */
        *hiof = HeckleIOInitFields(si, sx);
        *hios = HeckleIOInitSpecies(si, sx);
        *hior = HeckleIOInitRestarts(si, sx, argv[1]);
        *hiot = HeckleIOInitTime(si, sx);
    }






    /* --------------------------------------------------------------- */
    /*              START RUNNING FROM A PREVIOUS RUN                  */
    /* --------------------------------------------------------------- */

    // TODO re-write the restart module

    if (si->rst != 0)
    {
//      printf("________________ read restart file for node   %3i ____\n\n", sx->r);

//      /* the filename depends on the type of restart file */
//      if (si->rst == 1) sprintf(nfile, "hr%03i.dat", sx->r);
//      if (si->rst == 2) sprintf(nfile, "hR%03i.dat", sx->r);
//
//      fp = fopen(nfile, "rb");
//      if (fp == NULL) printf("problem in opening file %s\n", nfile);
//
//      /* __ get the iteration 'it' of the restart  __ */
//      freadint(__FILE__, __LINE__, fp, 0, 1, it, &(sx->irun));
//
//        /* now for each specieitreadsxe # of particles and their weight */
//        for (s = 1; s < NS+1; s++)
//        {
//            freadint(__FILE__, __LINE__, fp, 0, 1, &(si->ns[s]), &(sx->irun));
//            freaddbl(__FILE__, __LINE__, fp, 0, 1, &(si->ws[s]), &(sx->irun));
//        }
//
//        /* get Parallel MPIsxnformation */
//        freadstx(__FILE__, __LINE__, fp, 0, 1, sx, &(sx->irun));
//
//
//        /* # of grid points ansxmemosx allocation for the grid */
//        nn1 = (sx->n[0]+1)*(sx->n[1]+1)*(sx->n[2]+1);
//        nn2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);
//
//        *s1 = (ST1 *)malloc(nn1*sizeof(ST1));
//        *s2 = (ST2 *)malloc(nn2*sizeof(ST2));
//
//        /* allocate memory for each species' particle array */
//        for (s = 1; s < NS+1; s++)
//        {
//            (*sp)[s] = malloc(si->nm*sizeof (Particle));
//        }
//
//        /* read the st1 and st2 structures */
//        freadst1(__FILE__, __LINE__, fp, 0, nn1, *s1, &(sx->irun));
//        freadst2(__FILE__, __LINE__, fp, 0, nn2, *s2, &(sx->irun));
//
//        /* read the stp structus2  */
//        for (s = 1; s < NS+1; s++)
//        {
//            freadstp(__FILE__, __LINE__, fp, 0, sx->ns[s], (*sp)[s], &(sx->irun));
//        }
//
//        /* __ set initial tsxal energy __ */
//        freaddbl(__FILE__, __LINE__, fp, 0, 1, &((*sd).e0), &(sx->irun));
//
//        /* __ close the fisd __ */
//        fclose(fp);

        // start the IO module just for restart
        *hior = HeckleIOInitRestarts(si, sx, argv[1]);

        // read the appropriate group in the restart file
        readRestarts(*hior, si, sx, s0, s1, s2, sp, sd, it, argv[1]);

        // start the other IO modules (now all struct are correctly filled, so do sx)
        *hiof = HeckleIOInitFields(si, sx);
        *hios = HeckleIOInitSpecies(si, sx);
        *hiot = HeckleIOInitTime(si, sx);

        /* initialize the boundary condition module */
        *hbc = HeckleBCInit(si);

        // initialize the ghost point communication module
        *ghost  = GhostsInit(si, sx);

        /* set to zero the # of part in each cell */
        if (si->coll != 0) {
            /*# of points and memory allocation for the grids g1 and g2 */
            nn0 = (sx->n[0])*(sx->n[1])*(sx->n[2]);

            for (ijk0 = 0; ijk0 < nn0; ijk0++) {
                for (s = 1; s < NS+1; s++) {
                    (*s0)[ijk0].npart[s] = 0;
                }
            }
        }


    } /* end if RESTART */



    /* calculate the full ion pressure tensor */
    ions(si, sx, *s2, *sp, *ghost, *hbc, *it);


    /* __ diagnostic __ */
    check(*si, *sx, *s0, *s1, *s2, *sp, sd, *it, ticks, ttime);
    writeTime(*hiot, si, sx, sd, *it);
    writeFields(*hiof, sx, *s1, *s2, *it*si->ts);
    writeSpecies(*hios, si, sx, *sp, *it*si->ts);
    //writeRestarts(*hior, si, sx, *s0, *s1, *s2, *sp, *sd, *it, argv[1]);
  //writeorbit(*si, sx, *s1, *s2, *sp, so, *it);

}





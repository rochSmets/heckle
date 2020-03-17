
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "structures.h"
#include "check.h"
#include "initModel.h"
#include "closeModel.h"
#include "defines.h"
#include "maxwell.h"
#include "init.h"
#include "ions.h"
#include "ohm.h"
#include "read.h"
#include "setup.h"
#include "start.h"
#include "subdomains.h"
#include "sources.h"
#include "collision.h"
#include <HeckleIOFields.h>
#include <HeckleIOSpecies.h>
#include <HeckleIORestarts.h>
#include <HeckleIOTime.h>
#include <ghosts.h>
#include <hecklebc.h>
#include "Particle/particle.h"



/* __ main __________________________________________________________________ */
int main(int ac, char **av)
{
    STI si;
    STX sx;
    Collision sc;
    Grid0 *s0;
    ST1 *s1;
    ST2 *s2;
    Particle *sp[NS+1];
    struct sto so;
    struct std sd;
    HeckleBC *hbc;
    Ghosts   *ghosts;
    HeckleIOFields *hiof = NULL;
    HeckleIOSpecies *hios = NULL;
    HeckleIORestart *hior = NULL;
    HeckleIOTime *hiot = NULL;
    int it;
    int ipc;
    int s;
    int irun;
    clock_t ticks[2];
    time_t ttime[2];//, timer;
//  char fulltime[24];


    /* __ mpi set-up __ */
    MPI_Init(&ac, &av);
    MPI_Comm_size(MPI_COMM_WORLD, &(sx.s));
    MPI_Comm_rank(MPI_COMM_WORLD, &(sx.r));

    /* __ set the run parameters __ */
    si = setup(sx, av[1]);
    sc = setCollision(si, sx, &sd, av[1]);

    /* __ time & date of the begining of the run __ */
//  if (sx.r == 0) {
//      timer = time(NULL);
//      strncpy(fulltime, asctime(localtime(&timer)), 24);
//      printf("________________ %s ____________\n",  fulltime);
//      printf("\n");
//  }

    /* __ set the initial topology & closure equation __ */
    initModelStart(&si, &sx, av[1]);
    closeModelStart(&si, &sx, av[1]);

    /* __ setup the so struct. __ */
//  readhorbi(&so);
//  no need to read horbi.txt... while writeorbits is not ok for now,

    /* __ heckle running from a restart ? __ */
    start(&hiof, &hios, &hior, &hiot, &si, &sx, &sc, &s0, &s1, &s2, &sp, &hbc, &ghosts, &sd, &so, &it, ticks, ttime, ac, av);
    /* __ loop on the total # of time step __ */
    while (it < si.nt) {
        // .....................................................................
        //                              PREDICTOR
        // .....................................................................

        ipc = 0;

        /* __ call the collision operator __ */
        if (si.coll != 0) {
            collision(sc, si, sx, &sd, s0, s1, s2, sp, ghosts, hbc, ipc);
        }

        /* __ move particles and accumulate density and bulk vel. __ */
        sources(&si, &sx, s0, s1, s2, sp, hbc, ghosts, ipc, it);

        /* __ build the full ion kinetic pressure tensor __ */
        if (it%si.tf == 0) ions(&si, &sx, s2, sp, ghosts, hbc, it);

        /* __ Get B^{n+1/2} from B^n __ */
        /* __ and J^{n+1/2} from B^{n+1/2} __ */
        MaxwellFaraday(&si, &sx, s1, s2, ipc);
        MaxwellAmpere(&si, &sx, s1, s2, ghosts, hbc, ipc);


        /* __ Calculate electron pressure tensor from the electron closure __ */
        closeModelPressure(it, &si, &sx, s1, s2, hbc, ghosts, ipc);

        /* __ drive the system if needed __ */
        if (si.drive != 0) {
            initModelDrive(si, sx, s0, s1, s2, sp, ipc, it);
        }


        /* __ Get the electric field E^{n+1/2} from Ohm's law __ */
        /* __ then predict E^{n+1} using E^{n+1/2} and E^n __ */
        ohm(&si, &sx, s1, s2, hbc, ghosts, ipc);


        /* __ now get B^{n+1}  from B^{n+1/2} and the predicted E^{n+1} __ */
        /* __ then get J^{n+1} from B^{n+1} __ */
        ipc = 2;
        MaxwellFaraday(&si, &sx, s1, s2, ipc);
        MaxwellAmpere(&si, &sx, s1, s2, ghosts, hbc, ipc);


        // .....................................................................
        //                              CORRECTOR
        // .....................................................................

        ipc = 1;

        if (si.coll != 0) {
            collision(sc, si, sx, &sd, s0, s1, s2, sp, ghosts, hbc, ipc);
        }

        /* __ move particles and accumulate density and bulk vel. __ */
        sources(&si, &sx, s0, s1, s2, sp, hbc, ghosts, ipc, it);

        /* __ Get B^{n+3/2} from B^{n+1} __ */
        /* __ and J^{n+3/2} from B^{n+3/2} __ */
        MaxwellFaraday(&si, &sx, s1, s2, ipc);
        MaxwellAmpere(&si, &sx, s1, s2, ghosts, hbc, ipc);

        /* __ Calculate electron pressure tensor from the electron closure __ */
        closeModelPressure(it, &si, &sx, s1, s2, hbc, ghosts, ipc);

        /* __ drive the system if needed __ */
        if (si.drive != 0) {
              initModelDrive(si, sx, s0, s1, s2, sp, ipc, it);
        }


        /* __ Get the electric field E^{n+3/2} from Ohm's law __ */
        /* __ then correct E^{n+1} using E^{n+1/2} and E^{n+3/2} __ */
        ohm(&si, &sx, s1, s2, hbc, ghosts, ipc);


        /* __ now get B^{n+1}  from B^{n+1/2} and the corrected E^{n+1} __ */
        /* __ then get J^{n+1} from B^{n+1} __ */
        ipc = 2;
        MaxwellFaraday(&si, &sx, s1, s2, ipc);
        MaxwellAmpere(&si, &sx, s1, s2, ghosts, hbc, ipc);



        // .....................................................................
        //                              DIAGNOSTICS
        // .....................................................................

        /* __ add part. if the density is too low __ */
        //if (si.feed == 1) add(&si, &sx, &st, s2, sp, com);

        /* __ increase time step __ */
        it ++;

        /* __ diagnostic __ */
       if (it%si.tt == 0) check(si, sx, s0, s1, s2, sp, &sd, it, ticks, ttime);
       if (it%si.tt == 0) writeTime(hiot, &si, &sx, &sd, it);
       if (it%si.tf == 0) writeFields(hiof, &sx, s1, s2, it*si.ts);
       if (it%si.tp == 0) writeSpecies(hios, &si, &sx, sp, it*si.ts);
       if (it%si.tr == 0) writeRestarts(hior, &si, &sx, s0, s1, s2, sp, sd, it, av[1]);

        /* __ reduce the irun index __ */
        MPI_Allreduce(&(sx.irun), &irun, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        /* __ if the run has to be stopped __ */
        if (irun != 0) {
            /* __ print fancy stuff __ */
            if (sx.r == 0) {
//              timer = time(NULL);
//              strncpy(fulltime, asctime(localtime(&timer)), 24);
//              printf("________________ %s ____________\n",  fulltime);
//              printf("\n");
                printf("________________ this is the end... __________________\n");
                printf("\n");
                printf("________________ run is dead : err. # %4d ___________\n", irun);
                printf("\n");
            }

            /* __ get out of the while __ */
            break;

            MPI_Finalize();

            return 0;
        }
    }

    /* __ time & date of the end of the run __ */
    if (sx.r == 0) {
//      timer = time(NULL);
//      strncpy(fulltime, asctime(localtime(&timer)), 24);
//      printf("________________ %s ____________\n",  fulltime);
//      printf("\n");
        printf("________________ this is the end... __________________\n");
        printf("\n");
    }

    /* __ clean-up the pointers __ */
    free(s1);
    free(s2);

    for (s = 1; s < NS+1; s++) {
        free(sp[s]);
    }


    /* __ close the BC module __ */
    HeckleBCDelete(hbc);


    /* __ close the IO module __ */
    HeckleIODeleteFields(hiof);
    HeckleIODeleteSpecies(hios);
    HeckleIODeleteRestarts(hior);
    HeckleIODeleteTime(hiot);

    MPI_Finalize();

    return 0;

}


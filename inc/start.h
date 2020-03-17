
#ifndef START
#define START

#include <HeckleIOFields.h>
#include <HeckleIOSpecies.h>
#include <HeckleIORestarts.h>
#include <HeckleIOTime.h>
#include <structures.h>
#include <mpi.h>
#include <particle.h>
#include <particlebc.h>
#include <hecklebc.h>
#include <ghosts.h>


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
           Particle *(*)[NS+1],
           HeckleBC **hbc,
           Ghosts **gc,
           struct std *sd,
           struct sto *so,
           int *it,
           clock_t ticks[2],
           time_t ttime[2],
           int argc,
           char **argv);

#endif


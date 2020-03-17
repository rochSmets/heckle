
#ifndef INIT
#define INIT

#include<structures.h>
#include <particle.h>
#include <particlebc.h>
#include <hecklebc.h>
#include <ghosts.h>



/* _____ initialize the particles and the field _____________________________ */
void init(STI *si,
          STX *sx,
          Grid0 *s0,
          ST1 *s1,
          ST2 * s2,
          Particle *sp[NS+1],
          HeckleBC *hbc,
          Ghosts *ghosts,
          int *it);

#endif


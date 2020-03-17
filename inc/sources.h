
#ifndef SOURCES
#define SOURCES


#include<structures.h>
#include<particle.h>
#include<particlebc.h>
#include <ghosts.h>
#include <hecklebc.h>



/*---------------------------------------------------------------------------
  sources()
  ---------------------------------------------------------------------------
  AIM : push particles and accumulate moments on the grid
 ---------------------------------------------------------------------------*/
void sources(STI *si,
             STX *sx,
             Grid0 *s0,
             ST1 *s1,
             ST2 *s2,
             Particle *[NS+1],
             HeckleBC *hbc,
             Ghosts *ghosts,
             int ipc,
             int it);
#endif

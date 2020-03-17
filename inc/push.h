
#ifndef PUSH
#define PUSH


#include <structures.h>
#include<particle.h>
#include<particlebc.h>


/* ____ push the particles __________________________________________________ */
void push(STI *si,
          STX *sx,
          Grid0 *s0,
          ST1 *s1,
          ST2 *s2,
          Particle *sp[NS+1],
          PartBC *pbc,
          int ipc,
          int s,
          int it);

#endif


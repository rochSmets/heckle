
#ifndef IONS
#define IONS

#include <structures.h>
#include <particle.h>
#include <hecklebc.h>
#include <ghosts.h>

/* ____ update moments from particles dynamics ______________________________ */
void ions(const STI* const si,
          const STX* const sx,
          ST2 *,
          Particle *stp[NS+1],
          Ghosts *ghosts, HeckleBC *hbc,
          int it);

#endif

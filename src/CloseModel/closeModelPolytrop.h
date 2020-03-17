
#ifndef CLOSE_POLYTROP_H
#define CLOSE_POLYTROP_H

#include "structures.h"



void polytropPressure(int it,
                      const STI * const si,
                      const STX * const sx,
                      struct st1 *s1,
                      struct st2 *s2,
                      HeckleBC *hbc,
                      Ghosts *ghosts,
                      int ipc);

#endif // endif CLOSE_POLYTROP_H


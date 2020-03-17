
#ifndef CLOSE_FULLIMPL_H
#define CLOSE_FULLIMPL_H

#include "structures.h"

#define ALPHA 0.5 // coeff for the fraction of implicit C operator


void implicitPressure(int it,
                    const STI * const si,
                    const STX * const sx,
                    struct st1 *s1,
                    struct st2 *s2,
                    HeckleBC *hbc,
                    Ghosts *ghosts,
                    int ipc);

#endif // endif CLOSE_FULLIMPL_H

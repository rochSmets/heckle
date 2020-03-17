
#ifndef SMOOTH
#define SMOOTH

#include <ghosts.h>
#include <structures.h>
#include <hecklebc.h>


void smooth(const STX * const sx,
            ST2 *s2, Ghosts *ghosts,
            HeckleBC *hbc);

void smoothPressure(const STX * const sx,
            ST2 *s2, Ghosts *ghosts,
            HeckleBC *hbc);

void smoothCurrent(const STX * const sx,
            ST2 *s2, Ghosts *ghosts,
            HeckleBC *hbc);
#endif


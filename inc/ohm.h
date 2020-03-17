
#ifndef OHM
#define OHM


#include <structures.h>
#include <ghosts.h>
#include <hecklebc.h>



/* _____ calculate the e fields _____________________________________________ */
void ohm(const STI * const si,
         const STX * const sx,
         ST1 *s1, ST2 *s2,
         HeckleBC *hbc,
         Ghosts *ghosts,
         int ipc);

#endif


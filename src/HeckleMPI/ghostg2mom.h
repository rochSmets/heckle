#ifndef GHOSTMOMG2_H
#define GHOSTMOMG2_H


#include <misc.h>



typedef struct ghostg2mom_s GhostG2Mom;



GhostG2Mom* GhostG2MomInit(const STX * const sx, const int n[3]);



void GhostG2MomDelete(GhostG2Mom *self);



void GhostG2MomSendRecv(GhostG2Mom *self, const STX* const sx, ST2 * s2);



// this function resets the internal buffers
// should be private and called once data has been received
// and unpacked in g2.
// the buffer we reset depend on ghostID (mrecv/send or vrecv/send)
void reset(GhostG2Mom *self);






#endif // GHOSTMOMG2_H


#ifndef GhostG4Winske_H
#define GhostG4Winske_H



#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>


typedef struct ghostg4winske_s GhostG4Winske;






/*---------------------------------------------------------------------------
  GhostG4WinskeInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
GhostG4Winske* GhostG4WinskeInit(const STX * const sx);









/*---------------------------------------------------------------------------
  GhostG4WinskeDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG4WinskeDelete(GhostG4Winske* self);









/*---------------------------------------------------------------------------
  GhostG4WinskeSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG4WinskeSendRecv(GhostG4Winske* self, const STX * const sx, ST4 * s4);













#endif // GhostG4Winske_H


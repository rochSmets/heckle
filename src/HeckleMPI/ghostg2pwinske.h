#ifndef GhostG2Winske_H
#define GhostG2Winske_H



#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>

#include <mpi.h>



#define GHOSTFIELD_P 0
#define GHOSTFIELD_DRIVER 1

typedef struct GhostG2winske_s GhostG2Winske;








/*---------------------------------------------------------------------------
  GhostG2WinskeInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (P)
 ---------------------------------------------------------------------------*/
GhostG2Winske* GhostG2WinskeInit(const STX * const sx, int fieldID);






/*---------------------------------------------------------------------------
  GhostG2WinskeDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG2WinskeDelete(GhostG2Winske* self);






/*---------------------------------------------------------------------------
  GhostG2WinskeSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG2WinskeSendRecv(GhostG2Winske* self, const STX * const sx, ST2 * s2);












#endif // GhostG2Winske_H


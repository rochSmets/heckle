#ifndef GHOSTG2P_H
#define GHOSTG2P_H



#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>

#include <mpi.h>





typedef struct ghostg2mom_s GhostG2P;








/*---------------------------------------------------------------------------
  GhostG2PInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (P)
 ---------------------------------------------------------------------------*/
GhostG2P* GhostG2PInit(const STX * const sx);






/*---------------------------------------------------------------------------
  GhostG2PDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG2PDelete(GhostG2P* self);






/*---------------------------------------------------------------------------
  GhostG2PSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG2PSendRecv(GhostG2P* self, const STX * const sx, ST2 * s2);












#endif // GHOSTG2P_H


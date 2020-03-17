#ifndef GhostG4Current_H
#define GhostG4Current_H



#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>


typedef struct GhostG4Current_s GhostG4Current;






/*---------------------------------------------------------------------------
  GhostG4CurrentInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
GhostG4Current* GhostG4CurrentInit(const STX * const sx);









/*---------------------------------------------------------------------------
  GhostG4CurrentDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG4CurrentDelete(GhostG4Current* self);









/*---------------------------------------------------------------------------
  GhostG4CurrentSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG4CurrentSendRecv(GhostG4Current* self, const STX * const sx, ST4 * s4);













#endif // GhostG4Current_H


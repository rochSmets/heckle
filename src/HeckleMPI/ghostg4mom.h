#ifndef GHOSTG4MOM_H
#define GHOSTG4MOM_H



#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>


typedef struct ghostg4mom_s GhostG4Mom;






/*---------------------------------------------------------------------------
  GhostG4MomInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
GhostG4Mom* GhostG4MomInit(const STX * const sx);









/*---------------------------------------------------------------------------
  GhostG4MomDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG4MomDelete(GhostG4Mom* self);









/*---------------------------------------------------------------------------
  GhostG4MomSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG4MomSendRecv(GhostG4Mom* self, const STX * const sx, ST4 * s4);













#endif // GHOSTG4MOM_H


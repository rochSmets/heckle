#ifndef GHOSTS_H
#define GHOSTS_H


#include <misc.h>
#include <ghostg2mom.h>
#include <ghostg2fields.h>




#define GHOST_NS_VS            0
#define GHOST_NV               1
#define GHOST_P                2
#define GHOST_E                3
#define GHOST_F                4
#define GHOST_J                5
#define GHOST_P_FULLP          6
#define GHOST_DRIVER_FULLP     7
#define GHOST_P_FULLP_G4       8
#define GHOST_J_G4             9



// public object
typedef struct ghosts_s Ghosts;




/*---------------------------------------------------------------------------
  GhostsInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
Ghosts* GhostsInit(const STI*   const si, const STX * const sx);





/*---------------------------------------------------------------------------
  GhostsDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostsDelete(Ghosts* self);



/*---------------------------------------------------------------------------
  GhostsSendRecv()
  ---------------------------------------------------------------------------
  AIM : Communicate ghost points for the specified quantity
 ---------------------------------------------------------------------------*/
void GhostsSendRecv(Ghosts *self, const STX* const sx, void *data, int qtyID);







#endif // GHOSTS_H







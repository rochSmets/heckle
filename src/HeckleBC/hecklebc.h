#ifndef HECKLEBC_H
#define HECKLEBC_H


#include <structures.h>
#include <particlebc.h>





// public object for heckle boundary condition
typedef struct s_hecklebc HeckleBC;





/*---------------------------------------------------------------------------
  HeckleBCInit()
  ---------------------------------------------------------------------------
  AIM : Initialize the Heckle Boundary Conditions module
 ---------------------------------------------------------------------------*/
HeckleBC* HeckleBCInit(const STI* const si);





/*---------------------------------------------------------------------------
  HeckleBCDelete()
  ---------------------------------------------------------------------------
  AIM : delete a heckleBC object
 ---------------------------------------------------------------------------*/
void HeckleBCDelete(HeckleBC *self);






/*---------------------------------------------------------------------------
  HeckleBCFieldApply()
  ---------------------------------------------------------------------------
  AIM : apply boundary condition for the selected quantity
 ---------------------------------------------------------------------------*/
void HeckleBCFieldApply(HeckleBC *self, const STX* const sx, ST2 *s2, int qtyID);




/*---------------------------------------------------------------------------
  HeckleBCGetPartBC()
  ---------------------------------------------------------------------------
  AIM : returns a handle on the particle boundary condition
 ---------------------------------------------------------------------------*/
PartBC* HeckleBCGetPartBC(HeckleBC *self);





#endif // HECKLEBC_H






#ifndef WRITEFIELDS_H
#define WRITEFIELDS_H


#include "structures.h"


typedef struct s_HeckleIOFields HeckleIOFields;




/*---------------------------------------------------------------------------
  HeckleIOInitFields()
  ---------------------------------------------------------------------------
  AIM : Initialize the heckleIOFields module
 ---------------------------------------------------------------------------*/
HeckleIOFields *HeckleIOInitFields(const STI* const si, const STX* const sx);



/*---------------------------------------------------------------------------
  HeckleIODeleteFields()
  ---------------------------------------------------------------------------
  AIM : Delete the heckleIOFields module handle
 ---------------------------------------------------------------------------*/
void HeckleIODeleteFields(HeckleIOFields *iof);



/*---------------------------------------------------------------------------
  writeFields()
  ---------------------------------------------------------------------------
  AIM : this routine writes the fields (electromagnetic and fluid moments)
  in a HDF5 file.
 ---------------------------------------------------------------------------*/

void writeFields(HeckleIOFields *hiof,  /* IO module handle  */
                 const STX* const sx,   /* MPI parameters    */
                 const ST1* const s1,   /* g1 grid           */
                 const ST2* const s2,   /* g2 grid           */
                 double time);          /* current time      */


#endif /* WRITEFIELDS_H */


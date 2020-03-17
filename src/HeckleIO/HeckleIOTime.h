
#ifndef WRITETIME_H
#define WRITETIME_H



typedef struct s_HeckleIOTime HeckleIOTime;




/*---------------------------------------------------------------------------
  HeckleIOInitTime()
  ---------------------------------------------------------------------------
  AIM : Initialize the heckleIOTime module
 ---------------------------------------------------------------------------*/
HeckleIOTime *HeckleIOInitTime(const STI* const si,
                               const STX* const sx);



/*---------------------------------------------------------------------------
  HeckleIODeleteTime()
  ---------------------------------------------------------------------------
  AIM : Delete the heckleIOTime module handle
 ---------------------------------------------------------------------------*/
void HeckleIODeleteTime(HeckleIOTime *hiot);



/*---------------------------------------------------------------------------
  writeTime()
  ---------------------------------------------------------------------------
  AIM : this routine writes the time dumps (div.b, energy components...
  in a HDF5 file.
 ---------------------------------------------------------------------------*/

void writeTime(HeckleIOTime *hiot,           /* IO module handle  */
               const STI* const si,          /* sti struct        */
               const STX* const sx,          /* MPI parameters    */
               struct std *sd,               /* std struct        */
               int it);                      /* current time      */


#endif /* WRITETIME_H */


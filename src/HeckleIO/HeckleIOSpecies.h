
#ifndef WRITESPECIES
#define WRITESPECIES

#include "structures.h"

typedef struct s_HeckleIOSpecies HeckleIOSpecies;




/*---------------------------------------------------------------------------
  HeckleIOInitSpecies()
  ---------------------------------------------------------------------------
  AIM : Initialize the heckleIOSpecies module
 ---------------------------------------------------------------------------*/
HeckleIOSpecies *HeckleIOInitSpecies(const STI* const si,
                                     const STX* const sx);



/*---------------------------------------------------------------------------
  HeckleIODeleteSpecies()
  ---------------------------------------------------------------------------
  AIM : Delete the heckleIOSpecies module handle
 ---------------------------------------------------------------------------*/
void HeckleIODeleteSpecies(HeckleIOSpecies *ios);



/*---------------------------------------------------------------------------
  writeSpecies()
  ---------------------------------------------------------------------------
  AIM : this routine writes the species (position, velocity & id)
  in a HDF5 file.
 ---------------------------------------------------------------------------*/
void writeSpecies(HeckleIOSpecies *hios,    /* IO module handle  */
                  const STI* const si,      /* run parameters    */
                  const STX* const sx,      /* MPI parameters    */
                  struct stp *sp[NS+1],     /* part positions, velocities & id */
                  double time);             /* current time      */

#endif


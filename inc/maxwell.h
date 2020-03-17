
#ifndef MAXWELL_H
#define MAXWELL_H


#include <structures.h>
#include <ghosts.h>
#include <hecklebc.h>



/*---------------------------------------------------------------------------
  MaxwellFaraday()
  ---------------------------------------------------------------------------
  AIM : update the magnetic field from the curl of the electric field

  PREDICTOR STEP :

  ipc = 0 :  B^{n+1/2} = B^n - dt/2 * curl(E^n)
  (ohm gets E^{n+1/2} and extrapolate to E^{n+1}
  ipc = 2 : B^{n+1} = B^{n+1/2} - dt/2 * curl(E^n+1)


  CORRECTOR STEP :

  ipc = 0 : B^{n+3/2} = B^{n+1} - curl(E^{n+1})*dt/2
  (ohm gets E^{n+3/2} and interpolate E^{n+1})
  ipc = 2 : B^{n+1} = B^{n+1/2} - dt/2 curl(E^{n+1})


  ipc=0 and ipc=2 lead to a centered scheme :
  (B^n+1 - B^n)/dt = -0.5 curl(E^n+1 + E^n)

 ---------------------------------------------------------------------------*/
void MaxwellFaraday(const STI * const si,
                    const STX * const sx,
                    ST1 *s1, ST2 *s2, int ipc);




/*---------------------------------------------------------------------------
  MaxwellAmpere()
  ---------------------------------------------------------------------------
  AIM : Calculate the electric current density J given the magnetic field B
        from Maxwell Ampere's law, assuming no displacement current
 ---------------------------------------------------------------------------*/
void MaxwellAmpere(const STI * const si,
                   const STX * const sx,
                   const ST1 * const s1, ST2 *s2,
                   Ghosts *ghosts, HeckleBC *hbc, int ipc);




#endif // MAXWELL_H


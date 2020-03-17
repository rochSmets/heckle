
#ifndef INIT_RECOASYM_H
#define INIT_RECOASYM_H

#include "structures.h"


/*---------------------------------------------------------------------------
    AsymAngleStart()
  ---------------------------------------------------------------------------
    AIM : read the parameters for the AsymAngle model
 ---------------------------------------------------------------------------*/
void asymangle_start(struct sti *si, struct stx *sx, char *dir);




/*---------------------------------------------------------------------------
    asymangleDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for all species at the given position (pos)
 ---------------------------------------------------------------------------*/
double asymangleDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe);







/*---------------------------------------------------------------------------
    asymangleMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleMagnetic(struct sti *si, struct stx *sx,
                    double pos[3], double B[3]);





/*---------------------------------------------------------------------------
    asymangleElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleElectric(struct sti *si, struct stx *sx,
                       double pos[3], double E[3]);






/*---------------------------------------------------------------------------
    asymangleCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleCurrent(struct sti *si, struct stx *sx,
                      double pos[3], double J[3]);






/*---------------------------------------------------------------------------
    asymangleTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void asymangleTemperature(struct sti *si, struct stx *sx,
                          double pos[3], int ispe, double T[2]);





/*---------------------------------------------------------------------------
    asymangleCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void asymangleCurDrift(struct sti *si, struct stx *sx, double pos[3],
                       int ispe, double curdrift[3]);








/*---------------------------------------------------------------------------
    asymangleDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleDrift(struct sti *si, struct stx *sx,
                    double pos[3], int ispe, double vdrift[3]);


void asymangleDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);





#endif // INIT_RECOASYM_H

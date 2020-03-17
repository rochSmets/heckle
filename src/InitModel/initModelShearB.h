
#ifndef INIT_SHEARB_H
#define INIT_SHEARB_H

#include "structures.h"

void shearB_start(struct sti *si, struct stx *sx, char *dir);

double shearBDensity(struct sti *si, struct stx *sx, double pos[3], int ispe);

void shearBMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3]);

void shearBElectric(struct sti *si, struct stx *sx, double pos[3], double E[3]);

void shearBCurrent(struct sti *si, struct stx *sx, double pos[3], double J[3]);

void shearBCurDrift(struct sti *si, struct stx *sx, double pos[3], int ispe,
                    double curdrift[3]);

void shearBDrift(struct sti *si, struct stx *sx, double pos[3], int ispe,
                 double vdrift[3]);

void shearBTemperature(struct sti *si, struct stx *sx, double pos[3], int ispe,
                       double T[2]);

void shearBDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);

#endif // endif INIT_SHEARB_H

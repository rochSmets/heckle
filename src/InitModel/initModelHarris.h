
#ifndef INIT_HARRIS_H
#define INIT_HARRIS_H

#include "structures.h"


void harris_start(struct sti *si, struct stx *sx, char *dir);

double harrisDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe);


void harrisMagnetic(struct sti *si, struct stx *sx,
                      double pos[3], double B[3]);


void harrisElectric(struct sti *si, struct stx *sx,
                      double pos[3], double E[3]);


void harrisCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3]);


void harrisCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3]);


void harrisDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3]);


void harrisTemperature(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double T[2]);

void harrisDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);

#endif // endif INIT_HARRIS_H

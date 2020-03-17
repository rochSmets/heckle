
#ifndef INIT_UNIFORM_H
#define INIT_UNIFORM_H

#include "structures.h"


void uniform_start(struct sti *si, struct stx *sx, char *dir);

double uniformDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe);


void uniformMagnetic(struct sti *si, struct stx *sx,
                      double pos[3], double B[3]);


void uniformElectric(struct sti *si, struct stx *sx,
                      double pos[3], double E[3]);


void uniformCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3]);


void uniformCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3]);


void uniformDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3]);


void uniformTemperature(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double T[2]);




void uniformDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);



#endif // INIT_UNIFORM_H

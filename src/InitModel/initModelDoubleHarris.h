
#ifndef INIT_DOUBLEHARRIS_H
#define INIT_DOUBLEHARRIS_H

#include "structures.h"


void doubleharris_start(struct sti *si, struct stx *sx, char *dir);

double doubleharrisDensity(struct sti *si, struct stx *sx,
                                 double pos[3], int ispe);


void doubleharrisMagnetic(struct sti *si, struct stx *sx,
                                double pos[3], double B[3]);


void doubleharrisElectric(struct sti *si, struct stx *sx,
                                double pos[3], double E[3]);


void doubleharrisCurrent(struct sti *si, struct stx *sx,
                               double pos[3], double J[3]);


void doubleharrisCurDrift(struct sti *si, struct stx *sx, double pos[3],
                                int ispe, double curdrift[3]);


void doubleharrisDrift(struct sti *si, struct stx *sx,
                             double pos[3], int ispe, double vdrift[3]);


void doubleharrisTemperature(struct sti *si, struct stx *sx,
                                   double pos[3], int ispe, double T[2]);

void doubleharrisDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);



#endif // INIT_DOUBLEHARRIS_H

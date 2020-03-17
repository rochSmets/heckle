#ifndef MODEL_LASERS_H
#define MODEL_LASERS_H



#include <structures.h>



void nlasers_start(struct sti *si, struct stx *sx, char *dir);

double nlasersDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe);


void nlasersMagnetic(struct sti *si, struct stx *sx,
                      double pos[3], double B[3]);


void nlasersElectric(struct sti *si, struct stx *sx,
                      double pos[3], double E[3]);


void nlasersCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3]);


void nlasersCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3]);


void nlasersDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3]);


void nlasersTemperature(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double T[2]);

void nlasersDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);







#endif // MODEL_LASERS_H


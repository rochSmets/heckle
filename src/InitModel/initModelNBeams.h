#ifndef MODEL_NBEAMS_H
#define MODEL_NBEAMS_H



#include <structures.h>



void nbeams_start(struct sti *si, struct stx *sx, char *dir);

double nbeamsDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe);


void nbeamsMagnetic(struct sti *si, struct stx *sx,
                      double pos[3], double B[3]);


void nbeamsElectric(struct sti *si, struct stx *sx,
                      double pos[3], double E[3]);


void nbeamsCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3]);


void nbeamsCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3]);


void nbeamsDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3]);


void nbeamsTemperature(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double T[2]);

void nbeamsDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);







#endif // MODEL_NBEAMS_H


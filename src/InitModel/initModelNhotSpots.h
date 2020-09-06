#ifndef MODEL_NHOTSPOTS_H
#define MODEL_NHOTSPOTS_H



#include <structures.h>



void nhotspots_start(struct sti *si, struct stx *sx, char *dir);

double nhotspotsDensity(struct sti *si, struct stx *sx, double pos[3], int ispe);

void nhotspotsMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3]);

void nhotspotsElectric(struct sti *si, struct stx *sx, double pos[3], double E[3]);

void nhotspotsCurrent(struct sti *si, struct stx *sx, double pos[3], double J[3]);

void nhotspotsCurDrift(struct sti *si, struct stx *sx, double pos[3], int ispe, double curdrift[3]);

void nhotspotsDrift(struct sti *si, struct stx *sx, double pos[3], int ispe, double vdrift[3]);

void nhotspotsTemperature(struct sti *si, struct stx *sx, double pos[3], int ispe, double T[2]);

void nhotspotsDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);



#endif // MODEL_NHOTSPOTS_H


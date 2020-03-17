

#ifndef INIT_MODEL_H
#define INIT_MODEL_H

#include <defines.h>
#include <structures.h>
#include <Particle/particle.h>

typedef enum {
  INITMODEL_DOUBLEHARRIS,
  INITMODEL_HARRIS,
  INITMODEL_UNIFORM,
  INITMODEL_ASYMANGLE,
  INITMODEL_NBEAMS,
  INITMODEL_SHEARB,
  INITMODEL_NLASERS
} kind_initModel;

void initModelStart(struct sti *si, struct stx *sx, char *dir);

void initModelMagnetic(struct sti *si, struct stx *sx, double pos[3],
                       double B[3]);

void initModelElectric(struct sti *si, struct stx *sx, double pos[3],
                       double B[3]);

void initModelCurrent(struct sti *si, struct stx *sx, double pos[3],
                      double J[3]);

double initModelDensity(struct sti *si, struct stx *sx, double pos[3],
                        int speciesID);

void initModelTemperature(struct sti *si, struct stx *sx, double pos[3],
                          int ispe, double T[2]);

void initModelCurDrift(struct sti *si, struct stx *sx, double pos[3], int ispe,
                       double curdrift[3]);

void initModelDrift(struct sti *si, struct stx *sx, double pos[3], int ispe,
                    double vdrift[3]);

void initModelDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);

void normal(struct sti *si, struct stx *sx);

#endif // endif INIT_MODEL_H

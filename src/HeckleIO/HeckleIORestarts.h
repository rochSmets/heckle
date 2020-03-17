
#ifndef WRITERESTARTS_H
#define WRITERESTARTS_H

#include "structures.h"

#define RUN_FROM_SCRATCH 0
#define RUN_RESTARTING 1

typedef struct s_HeckleIORestart HeckleIORestart;

/*---------------------------------------------------------------------------
  HeckleIOInitRestarts()
  ---------------------------------------------------------------------------
  AIM : Initialize the HeckleIORestarts module
 ---------------------------------------------------------------------------*/
HeckleIORestart *HeckleIOInitRestarts(const STI *const si, const STX *const sx,
                                      char *dir);

/*---------------------------------------------------------------------------
  HeckleIODeleteRestarts()
  ---------------------------------------------------------------------------
  AIM :
 ---------------------------------------------------------------------------*/
void HeckleIODeleteRestarts(HeckleIORestart *hior);

/*---------------------------------------------------------------------------
  writeRestarts()
  ---------------------------------------------------------------------------
  AIM : this routine writes the restarts in a HDF5 file.
 ---------------------------------------------------------------------------*/
void writeRestarts(HeckleIORestart *hior, STI *si, STX *sx, Grid0 *s0,
                   const ST1 *const s1, const ST2 *const s2,
                   struct stp *sp[NS + 1], struct std sd, int it, char *dir);

void readRestarts(HeckleIORestart *hior, struct sti *si, struct stx *sx,
                  Grid0 **s0, struct st1 **s1, struct st2 **s2,
                  struct stp *(*sp)[NS + 1], struct std *sd, int *it,
                  char *dir);

#endif /* WRITERESTARTS_H */

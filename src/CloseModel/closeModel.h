

#ifndef CLOSEMODEL_H
#define CLOSEMODEL_H

#include "structures.h"
#include "defines.h"
#include "ghosts.h"
#include "hecklebc.h"


typedef enum{CLOSE_ISOTHERM,
             CLOSE_POLYTROP,
             CLOSE_FULLIMPL,
             CLOSE_FULLSUB} kind_closeModel;



void closeModelStart(struct sti *si,
                     struct stx *sxi,
                     char *dir);


void closeModelPressure(int it,
                        const STI * const si,
                        const STX * const sx,
                        struct st1 *s1,
                        struct st2 *s2,
                        HeckleBC *hbc,
                        Ghosts *ghosts,
                        int ipc);

#endif // CLOSEMODEL_H

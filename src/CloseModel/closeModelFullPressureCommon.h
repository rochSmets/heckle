
//#if defined(WINSKE)  && !defined (SUBSTEP_NUMBER)
#if defined(FULLP)
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ghosts.h"
#include "hecklebc.h"
#include "bc_constants.h"
#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "smooth.h"

#define FORWARD 0
#define BACKWARD 1

#define TAU 100.0 // isotropisation time for I operator
//#define PSMOOTH 8 // stride for smoothing the pressure
#define BFORE 2 // get memory of magnetic field at [n-1/2]



void setDriverTerm(const STI * const si, const STX * const sx, struct st2 *s2,
    int ipc);

void driverComponents(const STI * const si, const STX * const sx, struct st2 *,
    double[3][3], int[3]);

void setCyclotronTerm(double[6], double[6], double[3]);

void setIsotropTerm(double[6], double[6], double[3], double);

void vectorB(const STX * const sx, struct st1 *, int[3], int, double[3]);

void miscB(const STI *, const STX * const sx, double [3], double *, double [3], double *);

void gradStuff(const STI * const si, const STX * const sx, struct st2 *,
    double[3][3], double[3][3], double[3][3][3], int[3]);

#endif

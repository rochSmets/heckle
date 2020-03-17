
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "smooth.h"




/* __ full electron pressure tensor (isotherm) ______________________________ */
void isothermPressure(int it, const STI * const si, const STX * const sx, struct st1 *s1, struct st2 *s2,
        HeckleBC   *hbc,
        Ghosts     *ghosts, int ipc)
{
double *pw;
int ijk;
int n2;


// unused variables but may be used in other 'stress' routines
(void)si;
(void)ipc;

/* __ # of grid points on g2 __ */
n2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);

/* __ memory allocation __ */
pw = (double *)malloc(n2*sizeof(double));

/* __ set the scalar stress with isotherm hypothesis : loop on the subdomain __ */
for (ijk = 0; ijk < n2; ijk++)
    {
    pw[ijk] = s2[ijk].ns[0]*s2[ijk].te[0];
    }

/* __ smooth stress tensor __ */
//smooth(si, sx, pw, +1.0, +1.0, +1.0, +1.0, +1.0, +1.0, co);

/* __ loop on the subdomain __ */
for (ijk = 0; ijk < n2; ijk++)
    {
    s2[ijk].ps[0][0] = pw[ijk];
    s2[ijk].ps[0][1] = 0.0;
    s2[ijk].ps[0][2] = 0.0;
    s2[ijk].ps[0][3] = pw[ijk];
    s2[ijk].ps[0][4] = 0.0;
    s2[ijk].ps[0][5] = pw[ijk];
    }

/* __ clean-up the pointers __ */
free(pw);

}


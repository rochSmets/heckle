
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"
#include "defines.h"
#include "smooth.h"


/* __ full electron pressure tensor (isotherm) ______________________________ */
void stress(struct sti si, struct stx sx, struct st2 *s2, int ipc, MPI_Comm com)
{
double *pw;
int ijk;
int n2;
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ # of grid points on g2 __ */
n2 = (sx.n[0]+2)*(sx.n[1]+2)*(sx.n[2]+2);

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

/* __ clean up the communicator __ */
MPI_Comm_free(&co);

}

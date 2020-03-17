
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#include "structures.h"
#include "defines.h"
#include "fill.h"
#include "misc.h"
#include "smooth.h"
#include <ghosts.h>
#include <hecklebc.h>
#include <bc_constants.h>


#define bx(i,j,k) s1[k + (j)*nzg1 + (i)*(nzg1*nyg1)].b[0]
#define by(i,j,k) s1[k + (j)*nzg1 + (i)*(nzg1*nyg1)].b[1]
#define bz(i,j,k) s1[k + (j)*nzg1 + (i)*(nzg1*nyg1)].b[2]

#define cx(i,j,k) s1[k + (j)*nzg1 + (i)*(nzg1*nyg1)].c[0]
#define cy(i,j,k) s1[k + (j)*nzg1 + (i)*(nzg1*nyg1)].c[1]
#define cz(i,j,k) s1[k + (j)*nzg1 + (i)*(nzg1*nyg1)].c[2]

#define jx(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].j[0]
#define jy(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].j[1]
#define jz(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].j[2]


#define ex(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[0]
#define ey(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[1]
#define ez(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[2]





/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */






/*---------------------------------------------------------------------------
  MaxwellFaraday()
  ---------------------------------------------------------------------------
  AIM : update the magnetic field from the curl of the electric field

  PREDICTOR STEP :

  ipc = 0 :  B^{n+1/2} = B^n - dt/2 * curl(E^n)
  (ohm gets E^{n+1/2} and extrapolate to E^{n+1}
  ipc = 2 : B^{n+1} = B^{n+1/2} - dt/2 * curl(E^n+1)


  CORRECTOR STEP :

  ipc = 0 : B^{n+3/2} = B^{n+1} - curl(E^{n+1})*dt/2
  (ohm gets E^{n+3/2} and interpolate E^{n+1})
  ipc = 2 : B^{n+1} = B^{n+1/2} - dt/2 curl(E^{n+1})


  ipc=0 and ipc=2 lead to a centered scheme :
  (B^n+1 - B^n)/dt = -0.5 curl(E^n+1 + E^n)

 ---------------------------------------------------------------------------*/
void MaxwellFaraday(const STI * const si,
                    const STX * const sx,
                    ST1 *s1, ST2 *s2, int ipc)
{
    double dtx, dty, dtz;
    double dexdy, dexdz, deydx, deydz, dezdx, dezdy;
    int i, j, k;
    int nxg1, nyg1, nzg1;
    int nyg2, nzg2;


    /* __ # of grid points __ */

    nxg1 = sx->n[0]+1;
    nyg1 = sx->n[1]+1;
    nzg1 = sx->n[2]+1;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;


    /* __ set some needed parameters __ */
    dtx = 0.125*si->ts/si->dl[0];
    dty = 0.125*si->ts/si->dl[1];
    dtz = 0.125*si->ts/si->dl[2];


    /* __ set b @ n+1/2 (ipc=0), n+3/2(ipc=1), n+1(ipc=2) : nested loops on all indexes __ */
    for (i = 0; i < nxg1; i++)
    {
        for (j = 0; j < nyg1; j++)
        {
            for (k = 0; k < nzg1; k++)
            {
                 /* __ curl(E) __ */

                dexdy =   ex(i,   j+1, k  ) - ex(i  , j  , k  )
                        + ex(i+1, j+1, k  ) - ex(i+1, j  , k  )
                        + ex(i  , j+1, k+1) - ex(i  , j  , k+1)
                        + ex(i+1, j+1, k+1) - ex(i+1, j  , k+1);

                dexdz =   ex(i  , j  , k+1) - ex(i  , j  , k  )
                        + ex(i+1, j  , k+1) - ex(i+1, j  , k  )
                        + ex(i  , j+1, k+1) - ex(i  , j+1, k  )
                        + ex(i+1, j+1, k+1) - ex(i+1, j+1, k  );

                deydx =   ey(i+1, j  , k  ) - ey(i  , j  , k  )
                        + ey(i+1, j+1, k  ) - ey(i  , j+1, k  )
                        + ey(i+1, j  , k+1) - ey(i  , j  , k+1)
                        + ey(i+1, j+1, k+1) - ey(i  , j+1, k+1);

                deydz =   ey(i  , j  , k+1) - ey(i  , j  , k  )
                        + ey(i+1, j  , k+1) - ey(i+1, j  , k  )
                        + ey(i  , j+1, k+1) - ey(i  , j+1, k  )
                        + ey(i+1, j+1, k+1) - ey(i+1, j+1, k  );

                dezdx =   ez(i+1, j  , k  ) - ez(i  , j  , k  )
                        + ez(i+1, j+1, k  ) - ez(i  , j+1, k  )
                        + ez(i+1, j  , k+1) - ez(i  , j  , k+1)
                        + ez(i+1, j+1, k+1) - ez(i  , j+1, k+1);

                dezdy =   ez(i  , j+1, k  ) - ez(i  , j  , k  )
                        + ez(i+1, j+1, k  ) - ez(i+1, j  , k  )
                        + ez(i  , j+1, k+1) - ez(i  , j  , k+1)
                        + ez(i+1, j+1, k+1) - ez(i+1, j  , k+1);


                //dexdy = s2[ijk3].e[0]-s2[ijk1].e[0]+s2[ijk4].e[0]-s2[ijk2].e[0]
                //        +s2[ijk7].e[0]-s2[ijk5].e[0]+s2[ijk8].e[0]-s2[ijk6].e[0];

                //dexdz = s2[ijk5].e[0]-s2[ijk1].e[0]+s2[ijk6].e[0]-s2[ijk2].e[0]
                //        +s2[ijk7].e[0]-s2[ijk3].e[0]+s2[ijk8].e[0]-s2[ijk4].e[0];

                /* __ set components of y divergence __ */
                //deydx = s2[ijk2].e[1]-s2[ijk1].e[1]+s2[ijk4].e[1]-s2[ijk3].e[1]
                //        +s2[ijk6].e[1]-s2[ijk5].e[1]+s2[ijk8].e[1]-s2[ijk7].e[1];


                //deydz = s2[ijk5].e[1]-s2[ijk1].e[1]+s2[ijk6].e[1]-s2[ijk2].e[1]
//                        +s2[ijk7].e[1]-s2[ijk3].e[1]+s2[ijk8].e[1]-s2[ijk4].e[1];

                /* __ set components of z divergence __ */
                //dezdx = s2[ijk2].e[2]-s2[ijk1].e[2]+s2[ijk4].e[2]-s2[ijk3].e[2]
                //        +s2[ijk6].e[2]-s2[ijk5].e[2]+s2[ijk8].e[2]-s2[ijk7].e[2];
               // dezdy = s2[ijk3].e[2]-s2[ijk1].e[2]+s2[ijk4].e[2]-s2[ijk2].e[2]
               //         +s2[ijk7].e[2]-s2[ijk5].e[2]+s2[ijk8].e[2]-s2[ijk6].e[2];

                switch (ipc)
                {
                case 0 : /* __ predictor : set b field @ n+1/2 __ */
                    cx(i,j,k) = bx(i,j,k) + deydz*dtz - dezdy*dty;
                    cy(i,j,k) = by(i,j,k) + dezdx*dtx - dexdz*dtz;
                    cz(i,j,k) = bz(i,j,k) + dexdy*dty - deydx*dtx;
                    break;

                case 1 : /* __ corrector : set b field @ n+3/2 __ */
                    bx(i,j,k) = bx(i,j,k) + deydz*dtz - dezdy*dty;
                    by(i,j,k) = by(i,j,k) + dezdx*dtx - dexdz*dtz;
                    bz(i,j,k) = bz(i,j,k) + dexdy*dty - deydx*dtx;
                    break ;

                case 2 : /* __ set b field @ n+1 __ */
                    bx(i,j,k) = cx(i,j,k) + deydz*dtz - dezdy*dty;
                    by(i,j,k) = cy(i,j,k) + dezdx*dtx - dexdz*dtz;
                    bz(i,j,k) = cz(i,j,k) + dexdy*dty - deydx*dtx;
                    break ;

                    /* __ no reason to get there __ */
                default : IAMDEAD(sx->r);
                } // end switch ipc
            } // end loop on k
        } // end loop on j
    } // end loop on i
}
/*===========================================================================*/

















/*---------------------------------------------------------------------------
  MaxwellAmpere()
  ---------------------------------------------------------------------------
  AIM : Calculate the electric current density J given the magnetic field B
        from Maxwell Ampere's law, assuming no displacement current
 ---------------------------------------------------------------------------*/
void MaxwellAmpere(const STI * const si,
                   const STX * const sx,
                   const ST1 * const s1, ST2 *s2,
                   Ghosts *ghosts, HeckleBC *hbc, int ipc)
{
    int i, j,k;
    double dbxdy, dbxdz,dbydx,dbydz,dbzdx,dbzdy;
    int nxg1, nyg1, nzg1;
    int nyg2, nzg2;
    double odl[3];


    odl[0] = 1./si->dl[0];
    odl[1] = 1./si->dl[1];
    odl[2] = 1./si->dl[2];


    /* __ # of grid points __ */

    nxg1 = sx->n[0]+1;
    nyg1 = sx->n[1]+1;
    nzg1 = sx->n[2]+1;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;

    /* __ set j @ n+1/2 (ipc=0), n+3/2(ipc=1), n+1(ipc=2) : nested loops on all indexes __ */
    for (i = 1; i < nxg1; i++)
    {
        for (j = 1; j < nyg1; j++)
        {
            for (k = 1; k < nzg1; k++)
            {
                switch (ipc)
                {
                    case 0 : /* __ predictor : set rot(b) field @ n+1/2 __ */


                        dbxdy =   cx(i  , j  , k  ) - cx(i  , j-1, k  )
                                + cx(i-1, j  , k  ) - cx(i-1, j-1, k  )
                                + cx(i  , j  , k-1) - cx(i  , j-1, k-1)
                                + cx(i-1, j  , k-1) - cx(i-1, j-1, k-1);


                        dbxdz =   cx(i  , j  , k  ) - cx(i  , j  , k-1)
                                + cx(i-1, j  , k  ) - cx(i-1, j  , k-1)
                                + cx(i  , j-1, k  ) - cx(i  , j-1, k-1)
                                + cx(i-1, j-1, k  ) - cx(i-1, j-1, k-1);


                        dbydx =   cy(i  , j  , k  ) - cy(i-1, j  , k  )
                                + cy(i  , j-1, k  ) - cy(i-1, j-1, k  )
                                + cy(i  , j  , k-1) - cy(i-1, j  , k-1)
                                + cy(i  , j-1, k-1) - cy(i-1, j-1, k-1);


                        dbydz =   cy(i  , j  , k  ) - cy(i  , j  , k-1)
                                + cy(i-1, j  , k  ) - cy(i-1, j  , k-1)
                                + cy(i  , j-1, k  ) - cy(i  , j-1, k-1)
                                + cy(i-1, j-1, k  ) - cy(i-1, j-1, k-1);


                        dbzdx =   cz(i  , j  , k  ) - cz(i-1, j  , k  )
                                + cz(i  , j-1, k  ) - cz(i-1, j-1, k  )
                                + cz(i  , j  , k-1) - cz(i-1, j  , k-1)
                                + cz(i  , j-1, k-1) - cz(i-1, j-1, k-1);


                        dbzdy =   cz(i  , j  , k  ) - cz(i  , j-1, k  )
                                + cz(i-1, j  , k  ) - cz(i-1, j-1, k  )
                                + cz(i  , j  , k-1) - cz(i  , j-1, k-1)
                                + cz(i-1, j  , k-1) - cz(i-1, j-1, k-1);

                    break;



                 case 1 : /* __ predictor : set rot(b) field @ n+3/2 __ */


                    // ipc == 1 and 2 both get J from 'b'

                 case 2 : /* __ set rot(b) field @ n+1 __ */

                        dbxdy =   bx(i  , j  , k  ) - bx(i  , j-1, k  )
                                + bx(i-1, j  , k  ) - bx(i-1, j-1, k  )
                                + bx(i  , j  , k-1) - bx(i  , j-1, k-1)
                                + bx(i-1, j  , k-1) - bx(i-1, j-1, k-1);


                        dbxdz =   bx(i  , j  , k  ) - bx(i  , j  , k-1)
                                + bx(i-1, j  , k  ) - bx(i-1, j  , k-1)
                                + bx(i  , j-1, k  ) - bx(i  , j-1, k-1)
                                + bx(i-1, j-1, k  ) - bx(i-1, j-1, k-1);


                        dbydx =   by(i  , j  , k  ) - by(i-1, j  , k  )
                                + by(i  , j-1, k  ) - by(i-1, j-1, k  )
                                + by(i  , j  , k-1) - by(i-1, j  , k-1)
                                + by(i  , j-1, k-1) - by(i-1, j-1, k-1);


                        dbydz =   by(i  , j  , k  ) - by(i  , j  , k-1)
                                + by(i-1, j  , k  ) - by(i-1, j  , k-1)
                                + by(i  , j-1, k  ) - by(i  , j-1, k-1)
                                + by(i-1, j-1, k  ) - by(i-1, j-1, k-1);


                        dbzdx =   bz(i  , j  , k  ) - bz(i-1, j  , k  )
                                + bz(i  , j-1, k  ) - bz(i-1, j-1, k  )
                                + bz(i  , j  , k-1) - bz(i-1, j  , k-1)
                                + bz(i  , j-1, k-1) - bz(i-1, j-1, k-1);


                        dbzdy =   bz(i  , j  , k  ) - bz(i  , j-1, k  )
                                + bz(i-1, j  , k  ) - bz(i-1, j-1, k  )
                                + bz(i  , j  , k-1) - bz(i  , j-1, k-1)
                                + bz(i-1, j  , k-1) - bz(i-1, j-1, k-1);

                    break;

                    /* __ no reason to get there __ */
                default : IAMDEAD(sx->r);
                }

                /* __ set j __ */
                    jx(i,j,k) =  0.25*(dbzdy * odl[1] - dbydz * odl[2]);
                    jy(i,j,k) =  0.25*(dbxdz * odl[2] - dbzdx * odl[0]);
                    jz(i,j,k) =  0.25*(dbydx * odl[0] - dbxdy * odl[1]);

            }
        }
    }


    /* --------------------------------------------------------------------- */
    /*                 GHOST POINTS AND BOUNDARY CONDITION                   */
    /* --------------------------------------------------------------------- */

   GhostsSendRecv(ghosts,  sx, s2, GHOST_J);

   HeckleBCFieldApply(hbc, sx, s2, BC_VAR_J);

   //set smoothed current
	int l, ijk;
	int n0, n1, n2, n; // # of g2 grid points
	n0 = sx->n[0] + 2;
	n1 = sx->n[1] + 2;
	n2 = sx->n[2] + 2;
	n = n0 * n1 * n2;
    for (ijk = 0; ijk < n; ijk++) {
   		/* __ loop on the direction __ */
		for (l = 0; l < 3; l++) {
			s2[ijk].j_smoothed[l] = s2[ijk].j[l];
   		}
   	}
    smoothCurrent(sx, s2, ghosts, hbc);
}
/*===========================================================================*/
















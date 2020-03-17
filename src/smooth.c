
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"
#include "defines.h"
#include "fill.h"
#include "misc.h"

#include <ghosts.h>
#include <hecklebc.h>
#include <bc_constants.h>


/* __ smoothing a field defined on g2 _______________________________________ */
void smooth(const STX * const sx,
            ST2 *s2, Ghosts *ghosts,
            HeckleBC *hbc)
{
    ST4 *s4a, *s4b;
    int i, j, k;
    int ijk;
    int ijka, ijkb, ijkc, ijkd, ijke, ijkf, ijkg, ijkh;
    int ijki, ijkj, ijkk, ijkl, ijkm, ijkn, ijko, ijkp;
    int ijkq, ijkr, ijks, ijkt, ijku, ijkv, ijkw, ijkx;
    int ijky, ijkz;
    int n0, n1, n2, p0, p1, p2;
    int nn4;
    int c;
    const double k2 = 0.125;
    const double k3 = 0.0625;
    const double k4 = 0.03125;
    const double k5 = 0.015625;


    /* __ # of grid points __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;
    p0 = sx->n[0]+4;
    p1 = sx->n[1]+4;
    p2 = sx->n[2]+4;

    /* __ # of grid points on g4 __ */
    nn4 = p0*p1*p2;

    /* __ memory allocation __ */
    s4a = malloc(nn4 * sizeof *s4a);
    s4b = malloc(nn4 * sizeof *s4b);


    for (ijkw = 0; ijkw < p0*p1*p2; ijkw++)
    {
        s4a[ijkw].n    = 0;
        s4a[ijkw].v[0] = 0;
        s4a[ijkw].v[1] = 0;
        s4a[ijkw].v[2] = 0;
    }

    /* __ duplicate tw in uw : nested loops on the grid points of subdomain __ */
    for (i = 0; i < n0; i++)
    {
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n2; k++)
            {
                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ index on g4 __ */
                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

                /* __ fill the buffer on g4 __ */
                s4a[ijkw].n    = s2[ijk].ns[0];
                s4a[ijkw].v[0] = s2[ijk].vi[0];
                s4a[ijkw].v[1] = s2[ijk].vi[1];
                s4a[ijkw].v[2] = s2[ijk].vi[2];
            }
        }
    }

    /* __ fill the buffer uw on g4 __ */
    GhostsSendRecv(ghosts, sx, s4a, GHOST_NV);


    /* __ calculate the smoothing : nested loops on the grid points of subdomain __ */
    for (i = 1; i < sx->n[0]+3; i++)
    {
        for (j = 1; j < sx->n[1]+3; j++)
        {
            for (k = 1; k < sx->n[2]+3; k++)
            {
                /* __ index on g4 __ */
                ijk = IDX(i, j, k, p0, p1, p2);

                /* __ index on g4 __ */
                ijka = IDX(i-1, j  , k  , p0, p1, p2);
                ijkb = IDX(i+1, j  , k  , p0, p1, p2);
                ijkc = IDX(i  , j-1, k  , p0, p1, p2);
                ijkd = IDX(i  , j+1, k  , p0, p1, p2);
                ijke = IDX(i  , j  , k-1, p0, p1, p2);
                ijkf = IDX(i  , j  , k+1, p0, p1, p2);

                /* __ index on g4 __ */
                ijkg = IDX(i-1, j-1, k  , p0, p1, p2);
                ijkh = IDX(i-1, j+1, k  , p0, p1, p2);
                ijki = IDX(i-1, j  , k-1, p0, p1, p2);
                ijkj = IDX(i-1, j  , k+1, p0, p1, p2);
                ijkk = IDX(i  , j-1, k-1, p0, p1, p2);
                ijkl = IDX(i  , j-1, k+1, p0, p1, p2);
                ijkm = IDX(i  , j+1, k-1, p0, p1, p2);
                ijkn = IDX(i  , j+1, k+1, p0, p1, p2);
                ijko = IDX(i+1, j  , k-1, p0, p1, p2);
                ijkp = IDX(i+1, j  , k+1, p0, p1, p2);
                ijkq = IDX(i+1, j-1, k  , p0, p1, p2);
                ijkr = IDX(i+1, j+1, k  , p0, p1, p2);

                /* __ index on g4 __ */
                ijks = IDX(i-1, j-1, k-1, p0, p1, p2);
                ijkt = IDX(i-1, j-1, k+1, p0, p1, p2);
                ijku = IDX(i-1, j+1, k-1, p0, p1, p2);
                ijkv = IDX(i-1, j+1, k+1, p0, p1, p2);
                ijkw = IDX(i+1, j-1, k-1, p0, p1, p2);
                ijkx = IDX(i+1, j-1, k+1, p0, p1, p2);
                ijky = IDX(i+1, j+1, k-1, p0, p1, p2);
                ijkz = IDX(i+1, j+1, k+1, p0, p1, p2);

                /* ---------------------------------------------------------- */
                /*                      SMOOTHING DENSITY                     */
                /* ---------------------------------------------------------- */
                s4b[ijk].n = k2*(s4a[ijk].n)

                            +k3*(  s4a[ijka].n
                                 + s4a[ijkb].n
                                 + s4a[ijkc].n
                                 + s4a[ijkd].n
                                 + s4a[ijke].n
                                 + s4a[ijkf].n
                                 )

                            +k4*(  s4a[ijkg].n
                                 + s4a[ijkh].n
                                 + s4a[ijki].n
                                 + s4a[ijkj].n
                                 + s4a[ijkk].n
                                 + s4a[ijkl].n
                                 + s4a[ijkm].n
                                 + s4a[ijkn].n
                                 + s4a[ijko].n
                                 + s4a[ijkp].n
                                 + s4a[ijkq].n
                                 + s4a[ijkr].n
                                )

                            +k5*(  s4a[ijks].n
                                 + s4a[ijkt].n
                                 + s4a[ijku].n
                                 + s4a[ijkv].n
                                 + s4a[ijkw].n
                                 + s4a[ijkx].n
                                 + s4a[ijky].n
                                 + s4a[ijkz].n
                                );


                /* ---------------------------------------------------------- */
                /*                      SMOOTHING VELOCITY                    */
                /* ---------------------------------------------------------- */

                for (c=0; c < 3; c++)
                {
                    s4b[ijk].v[c] = k2*(s4a[ijk].v[c])

                                +k3*(  s4a[ijka].v[c]
                                     + s4a[ijkb].v[c]
                                     + s4a[ijkc].v[c]
                                     + s4a[ijkd].v[c]
                                     + s4a[ijke].v[c]
                                     + s4a[ijkf].v[c]
                                     )

                                +k4*(  s4a[ijkg].v[c]
                                     + s4a[ijkh].v[c]
                                     + s4a[ijki].v[c]
                                     + s4a[ijkj].v[c]
                                     + s4a[ijkk].v[c]
                                     + s4a[ijkl].v[c]
                                     + s4a[ijkm].v[c]
                                     + s4a[ijkn].v[c]
                                     + s4a[ijko].v[c]
                                     + s4a[ijkp].v[c]
                                     + s4a[ijkq].v[c]
                                     + s4a[ijkr].v[c]
                                    )

                                +k5*(  s4a[ijks].v[c]
                                     + s4a[ijkt].v[c]
                                     + s4a[ijku].v[c]
                                     + s4a[ijkv].v[c]
                                     + s4a[ijkw].v[c]
                                     + s4a[ijkx].v[c]
                                     + s4a[ijky].v[c]
                                     + s4a[ijkz].v[c]
                                    );
                } // end of velocity smoothing


            }
        }
    }

    /* __ reduce sw in tw : nested loops on the grid points of subdomain __ */
    for (i = 0; i < n0; i++)
    {
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n2; k++)
            {
                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ index on g4 __ */
                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

                // unpack s4 into s2
                s2[ijk].ns[0] = s4b[ijkw].n;
                s2[ijk].vi[0] = s4b[ijkw].v[0];
                s2[ijk].vi[1] = s4b[ijkw].v[1];
                s2[ijk].vi[2] = s4b[ijkw].v[2];
            }
        }
    }


    /* --------------------------------------------------------------------- */
    /*                       Applying boundary conditions                    */
    /* --------------------------------------------------------------------- */

    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_NV);
    free(s4a);
    free(s4b);
}

void smoothPressure(const STX * const sx,
	            ST2 *s2, Ghosts *ghosts,
	            HeckleBC *hbc)
	{
	    ST4 *s4a, *s4b;
	    int i, j, k;
	    int ijk;
	    int ijka, ijkb, ijkc, ijkd, ijke, ijkf, ijkg, ijkh;
	    int ijki, ijkj, ijkk, ijkl, ijkm, ijkn, ijko, ijkp;
	    int ijkq, ijkr, ijks, ijkt, ijku, ijkv, ijkw, ijkx;
	    int ijky, ijkz;
	    int n0, n1, n2, p0, p1, p2;
	    int nn4;
	    int c;
	    const double k2 = 0.125;
	    const double k3 = 0.0625;
	    const double k4 = 0.03125;
	    const double k5 = 0.015625;


	    /* __ # of grid points __ */
	    n0 = sx->n[0]+2;
	    n1 = sx->n[1]+2;
	    n2 = sx->n[2]+2;
	    p0 = sx->n[0]+4;
	    p1 = sx->n[1]+4;
	    p2 = sx->n[2]+4;

	    /* __ # of grid points on g4 __ */
	    nn4 = p0*p1*p2;

	    /* __ memory allocation __ */
	    s4a = malloc(nn4 * sizeof *s4a);
	    s4b = malloc(nn4 * sizeof *s4b);


	    for (ijkw = 0; ijkw < p0*p1*p2; ijkw++)
	    {
	        for (c = 0; c < 6; c++) {
	        	s4a[ijkw].ps[c]    = 0.0;
	        }
	    }

	    /* __ duplicate tw in uw : nested loops on the grid points of subdomain __ */
	    for (i = 0; i < n0; i++)
	    {
	        for (j = 0; j < n1; j++)
	        {
	            for (k = 0; k < n2; k++)
	            {
	                /* __ index on g2 __ */
	                ijk = IDX(i, j, k, n0, n1, n2);

	                /* __ index on g4 __ */
	                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

	                /* __ fill the buffer on g4 __ */
	                for (c = 0; c < 6; c++) {
	                	s4a[ijkw].ps[c] = s2[ijk].ps[0][c];
	                }
	            }
	        }
	    }

	    /* __ fill the buffer uw on g4 __ */
	    GhostsSendRecv(ghosts, sx, s4a, GHOST_P_FULLP_G4);


	    /* __ calculate the smoothing : nested loops on the grid points of subdomain __ */
	    for (i = 1; i < sx->n[0]+3; i++)
	    {
	        for (j = 1; j < sx->n[1]+3; j++)
	        {
	            for (k = 1; k < sx->n[2]+3; k++)
	            {
	                /* __ index on g4 __ */
	                ijk = IDX(i, j, k, p0, p1, p2);

	                /* __ index on g4 __ */
	                ijka = IDX(i-1, j  , k  , p0, p1, p2);
	                ijkb = IDX(i+1, j  , k  , p0, p1, p2);
	                ijkc = IDX(i  , j-1, k  , p0, p1, p2);
	                ijkd = IDX(i  , j+1, k  , p0, p1, p2);
	                ijke = IDX(i  , j  , k-1, p0, p1, p2);
	                ijkf = IDX(i  , j  , k+1, p0, p1, p2);

	                /* __ index on g4 __ */
	                ijkg = IDX(i-1, j-1, k  , p0, p1, p2);
	                ijkh = IDX(i-1, j+1, k  , p0, p1, p2);
	                ijki = IDX(i-1, j  , k-1, p0, p1, p2);
	                ijkj = IDX(i-1, j  , k+1, p0, p1, p2);
	                ijkk = IDX(i  , j-1, k-1, p0, p1, p2);
	                ijkl = IDX(i  , j-1, k+1, p0, p1, p2);
	                ijkm = IDX(i  , j+1, k-1, p0, p1, p2);
	                ijkn = IDX(i  , j+1, k+1, p0, p1, p2);
	                ijko = IDX(i+1, j  , k-1, p0, p1, p2);
	                ijkp = IDX(i+1, j  , k+1, p0, p1, p2);
	                ijkq = IDX(i+1, j-1, k  , p0, p1, p2);
	                ijkr = IDX(i+1, j+1, k  , p0, p1, p2);

	                /* __ index on g4 __ */
	                ijks = IDX(i-1, j-1, k-1, p0, p1, p2);
	                ijkt = IDX(i-1, j-1, k+1, p0, p1, p2);
	                ijku = IDX(i-1, j+1, k-1, p0, p1, p2);
	                ijkv = IDX(i-1, j+1, k+1, p0, p1, p2);
	                ijkw = IDX(i+1, j-1, k-1, p0, p1, p2);
	                ijkx = IDX(i+1, j-1, k+1, p0, p1, p2);
	                ijky = IDX(i+1, j+1, k-1, p0, p1, p2);
	                ijkz = IDX(i+1, j+1, k+1, p0, p1, p2);

					for (c = 0; c < 6; c++) {
						s4b[ijk].ps[c] = k2 * (s4a[ijk].ps[c])
						                                        + k3*(s4a[ijka].ps[c] + s4a[ijkb].ps[c]
						                                            + s4a[ijkc].ps[c] + s4a[ijkd].ps[c]
						                                            + s4a[ijke].ps[c] + s4a[ijkf].ps[c])

						                                        + k4*(s4a[ijkg].ps[c] + s4a[ijkh].ps[c]
						                                            + s4a[ijki].ps[c] + s4a[ijkj].ps[c]
						                                            + s4a[ijkk].ps[c] + s4a[ijkl].ps[c]
						                                            + s4a[ijkm].ps[c] + s4a[ijkn].ps[c]
						                                            + s4a[ijko].ps[c] + s4a[ijkp].ps[c]
						                                            + s4a[ijkq].ps[c] + s4a[ijkr].ps[c])

						                                        + k5*(s4a[ijks].ps[c] + s4a[ijkt].ps[c]
						                                            + s4a[ijku].ps[c] + s4a[ijkv].ps[c]
						                                            + s4a[ijkw].ps[c] + s4a[ijkx].ps[c]
						                                            + s4a[ijky].ps[c] + s4a[ijkz].ps[c]);
					}


	            }
	        }
	    }

	    /* __ reduce sw in tw : nested loops on the grid points of subdomain __ */
	    for (i = 0; i < n0; i++)
	    {
	        for (j = 0; j < n1; j++)
	        {
	            for (k = 0; k < n2; k++)
	            {
	                /* __ index on g2 __ */
	                ijk = IDX(i, j, k, n0, n1, n2);

	                /* __ index on g4 __ */
	                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

	                // unpack s4 into s2
	                for (c = 0; c < 6; c++) {
	                	s2[ijk].ps[0][c] = s4b[ijkw].ps[c];
	                }
	            }
	        }
	    }


	    /* --------------------------------------------------------------------- */
	    /*                       Applying boundary conditions                    */
	    /* --------------------------------------------------------------------- */

	    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_P_FULLP);
	    free(s4a);
	    free(s4b);
	}

void smoothCurrent(const STX * const sx,
            ST2 *s2, Ghosts *ghosts,
            HeckleBC *hbc)
{
    ST4 *s4a, *s4b;
    int i, j, k;
    int ijk;
    int ijka, ijkb, ijkc, ijkd, ijke, ijkf, ijkg, ijkh;
    int ijki, ijkj, ijkk, ijkl, ijkm, ijkn, ijko, ijkp;
    int ijkq, ijkr, ijks, ijkt, ijku, ijkv, ijkw, ijkx;
    int ijky, ijkz;
    int n0, n1, n2, p0, p1, p2;
    int nn4;
    int c;
    const double k2 = 0.125;
    const double k3 = 0.0625;
    const double k4 = 0.03125;
    const double k5 = 0.015625;


    /* __ # of grid points __ */
    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;
    p0 = sx->n[0]+4;
    p1 = sx->n[1]+4;
    p2 = sx->n[2]+4;

    /* __ # of grid points on g4 __ */
    nn4 = p0*p1*p2;

    /* __ memory allocation __ */
    s4a = malloc(nn4 * sizeof *s4a);
    s4b = malloc(nn4 * sizeof *s4b);


    for (ijkw = 0; ijkw < p0*p1*p2; ijkw++)
    {
        s4a[ijkw].j_smoothed[0] = 0;
        s4a[ijkw].j_smoothed[1] = 0;
        s4a[ijkw].j_smoothed[2] = 0;
    }

    /* __ duplicate tw in uw : nested loops on the grid points of subdomain __ */
    for (i = 0; i < n0; i++)
    {
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n2; k++)
            {
                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ index on g4 __ */
                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

                /* __ fill the buffer on g4 __ */
                s4a[ijkw].j_smoothed[0] = s2[ijk].j_smoothed[0];
                s4a[ijkw].j_smoothed[1] = s2[ijk].j_smoothed[1];
                s4a[ijkw].j_smoothed[2] = s2[ijk].j_smoothed[2];
            }
        }
    }

    /* __ fill the buffer uw on g4 __ */
    GhostsSendRecv(ghosts, sx, s4a, GHOST_J_G4);


    /* __ calculate the smoothing : nested loops on the grid points of subdomain __ */
    for (i = 1; i < sx->n[0]+3; i++)
    {
        for (j = 1; j < sx->n[1]+3; j++)
        {
            for (k = 1; k < sx->n[2]+3; k++)
            {
                /* __ index on g4 __ */
                ijk = IDX(i, j, k, p0, p1, p2);

                /* __ index on g4 __ */
                ijka = IDX(i-1, j  , k  , p0, p1, p2);
                ijkb = IDX(i+1, j  , k  , p0, p1, p2);
                ijkc = IDX(i  , j-1, k  , p0, p1, p2);
                ijkd = IDX(i  , j+1, k  , p0, p1, p2);
                ijke = IDX(i  , j  , k-1, p0, p1, p2);
                ijkf = IDX(i  , j  , k+1, p0, p1, p2);

                /* __ index on g4 __ */
                ijkg = IDX(i-1, j-1, k  , p0, p1, p2);
                ijkh = IDX(i-1, j+1, k  , p0, p1, p2);
                ijki = IDX(i-1, j  , k-1, p0, p1, p2);
                ijkj = IDX(i-1, j  , k+1, p0, p1, p2);
                ijkk = IDX(i  , j-1, k-1, p0, p1, p2);
                ijkl = IDX(i  , j-1, k+1, p0, p1, p2);
                ijkm = IDX(i  , j+1, k-1, p0, p1, p2);
                ijkn = IDX(i  , j+1, k+1, p0, p1, p2);
                ijko = IDX(i+1, j  , k-1, p0, p1, p2);
                ijkp = IDX(i+1, j  , k+1, p0, p1, p2);
                ijkq = IDX(i+1, j-1, k  , p0, p1, p2);
                ijkr = IDX(i+1, j+1, k  , p0, p1, p2);

                /* __ index on g4 __ */
                ijks = IDX(i-1, j-1, k-1, p0, p1, p2);
                ijkt = IDX(i-1, j-1, k+1, p0, p1, p2);
                ijku = IDX(i-1, j+1, k-1, p0, p1, p2);
                ijkv = IDX(i-1, j+1, k+1, p0, p1, p2);
                ijkw = IDX(i+1, j-1, k-1, p0, p1, p2);
                ijkx = IDX(i+1, j-1, k+1, p0, p1, p2);
                ijky = IDX(i+1, j+1, k-1, p0, p1, p2);
                ijkz = IDX(i+1, j+1, k+1, p0, p1, p2);


                /* ---------------------------------------------------------- */
                /*                      SMOOTHING Current                    */
                /* ---------------------------------------------------------- */

				for (c = 0; c < 3; c++) {
					s4b[ijk].j_smoothed[c] = k2 * (s4a[ijk].j_smoothed[c])
					                                        + k3*(s4a[ijka].j_smoothed[c] + s4a[ijkb].j_smoothed[c]
					                                            + s4a[ijkc].j_smoothed[c] + s4a[ijkd].j_smoothed[c]
					                                            + s4a[ijke].j_smoothed[c] + s4a[ijkf].j_smoothed[c])

					                                        + k4*(s4a[ijkg].j_smoothed[c] + s4a[ijkh].j_smoothed[c]
					                                            + s4a[ijki].j_smoothed[c] + s4a[ijkj].j_smoothed[c]
					                                            + s4a[ijkk].j_smoothed[c] + s4a[ijkl].j_smoothed[c]
					                                            + s4a[ijkm].j_smoothed[c] + s4a[ijkn].j_smoothed[c]
					                                            + s4a[ijko].j_smoothed[c] + s4a[ijkp].j_smoothed[c]
					                                            + s4a[ijkq].j_smoothed[c] + s4a[ijkr].j_smoothed[c])

					                                        + k5*(s4a[ijks].j_smoothed[c] + s4a[ijkt].j_smoothed[c]
					                                            + s4a[ijku].j_smoothed[c] + s4a[ijkv].j_smoothed[c]
					                                            + s4a[ijkw].j_smoothed[c] + s4a[ijkx].j_smoothed[c]
					                                            + s4a[ijky].j_smoothed[c] + s4a[ijkz].j_smoothed[c]);
				} // end of Current smoothing


            }
        }
    }

    /* __ reduce sw in tw : nested loops on the grid points of subdomain __ */
    for (i = 0; i < n0; i++)
    {
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n2; k++)
            {
                /* __ index on g2 __ */
                ijk = IDX(i, j, k, n0, n1, n2);

                /* __ index on g4 __ */
                ijkw = IDX(i+1, j+1, k+1, p0, p1, p2);

                // unpack s4 into s2
                s2[ijk].j_smoothed[0] = s4b[ijkw].j_smoothed[0];
                s2[ijk].j_smoothed[1] = s4b[ijkw].j_smoothed[1];
                s2[ijk].j_smoothed[2] = s4b[ijkw].j_smoothed[2];
            }
        }
    }


    /* --------------------------------------------------------------------- */
    /*                       Applying boundary conditions                    */
    /* --------------------------------------------------------------------- */

    HeckleBCFieldApply(hbc, sx, s2, BC_VAR_NV);
    free(s4a);
    free(s4b);
}








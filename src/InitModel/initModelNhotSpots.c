

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"


#define PARA 0
#define PERP 1

#define MAXNBEAMS 16


typedef struct s_nhotspots
{
    int N;
    double fi[MAXNBEAMS];
    double psi[MAXNBEAMS];
    double b[MAXNBEAMS];
    double T[NS+1];
    double axis0[MAXNBEAMS];
    double axis1[MAXNBEAMS];
    double widthRatio[MAXNBEAMS];
    double thick[MAXNBEAMS];
    double center[MAXNBEAMS][3];
    double tRise;
    double bDrive[MAXNBEAMS];
} NhotSpots;




/* this is for private use only */
NhotSpots NhotSpotsParams;





/*---------------------------------------------------------------------------
    nhotspots_start()
  ---------------------------------------------------------------------------
    AIM : read the parameters for nhotspots model
 ---------------------------------------------------------------------------*/
void nhotspots_start(struct sti *si, struct stx *sx, char *dir)
{
    int g, h;
    char sfh[80];
    char junk[100];
    FILE *fp;

    (void)si;


    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "NhotSpots.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL) {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscanint(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.N), &(sx->irun));
    if (NhotSpotsParams.N > MAXNBEAMS) {
        printf("too many beams in this initialization : increase MAXNBEAMS\n");
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.fi[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.psi[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.b[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.T[0]), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.T[1]), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.axis0[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.axis1[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.widthRatio[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.thick[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.center[h][0]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.center[h][1]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.center[h][2]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.tRise), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NhotSpotsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NhotSpotsParams.bDrive[h]), &(sx->irun));
    }


/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : N hotSpots ______\n");

   for (h = 0; h < NhotSpotsParams.N; h++) {
       printf("... Beam #       :%12d\n", h);
       printf("fi  [0, 360]     :%12.4f\n", NhotSpotsParams.fi[h]);
       printf("psi [0,  90]     :%12.4f\n", NhotSpotsParams.psi[h]);
       printf("B field          :%12.4f\n", NhotSpotsParams.b[h]);
       printf("\n");

       for (g = 0; g < 3; g++) {
           printf("beam center  [%1d] :%12.4f\n", g, NhotSpotsParams.center[h][g]);
       }
       printf("axis X           :%12.4f\n", NhotSpotsParams.axis0[h]);
       printf("axis Y           :%12.4f\n", NhotSpotsParams.axis1[h]);
       printf("width ratio      :%12.4f\n", NhotSpotsParams.widthRatio[h]);
       printf("thickness (z)    :%12.4f\n", NhotSpotsParams.thick[h]);
       printf("\n");
   }

   for (g = 0; g < NS+1; g++) {
       printf("___ specie %1d ___\n", g);
       printf("... temperature  :%12.4f\n", NhotSpotsParams.T[g]);
       printf("\n");
   }

   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ polynom used for the interpolation __ */
double polynomNhS(double x)
{
double w, y;


w = -6.0*fabs(x)*fabs(x)*fabs(x)*fabs(x)*fabs(x)
    +15.0*x*x*x*x
    -10.0*fabs(x)*fabs(x)*fabs(x)
    +1;

y = (fabs(x) <= 1.0) ? w : 0.0;

return y;

}



/* __ rotate the coordinates with fi & psi angles __ */
void rotateCoordinatesNhS(double pos[3], int h, double* wf, double* wp, double* tx, double* ty, double *tz)
{
//    double xc, yc, zc;
double wx, wy, wz;


/* __ set fi & psi angles __ */
*wf = NhotSpotsParams.fi[h] *PI/180.0;
*wp = NhotSpotsParams.psi[h]*PI/180.0;

/* __ rotation with fi angle __ */
wx = (pos[0]-NhotSpotsParams.center[h][0])*cos(*wf)
    +(pos[1]-NhotSpotsParams.center[h][1])*sin(*wf);
wy =-(pos[0]-NhotSpotsParams.center[h][0])*sin(*wf)
    +(pos[1]-NhotSpotsParams.center[h][1])*cos(*wf);
wz =  pos[2]-NhotSpotsParams.center[h][2];

/* __ rotation with psi angle __ */
*tx =  wx;
*ty =  wy*cos(*wp)+wz*sin(*wp);
*tz = -wy*sin(*wp)+wz*cos(*wp);


//   /* __ coordinates relative to the ellipsis center __ */
//   xc = pos[0]-NhotSpotsParams.center[h][0];
//   yc = pos[1]-NhotSpotsParams.center[h][1];
//   zc = pos[2]-NhotSpotsParams.center[h][2];
//   
//   /* __ backward rotation with -psi angle */
//   wx =  xc;
//   wy =  yc*cos(*wp)+zc*sin(*wp);
//   wz = -yc*sin(*wp)+zc*cos(*wp);
//   
//   /* __ backward rotation with -fi angle */
//   *tx =  wx*cos(*wf)+wy*sin(*wf);
//   *ty = -wx*sin(*wf)+wy*cos(*wf);
//   *tz =  wz;

}



/*---------------------------------------------------------------------------
    nhotspotsDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for the given specie (ispe)
    at the given position (pos)
 ---------------------------------------------------------------------------*/
double nhotspotsDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe)
{
    (void)si;
    (void)sx;
    (void)pos;
    (void) ispe;


    return 1.0;

}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    nhotspotsMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void nhotspotsMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3])
{
    double wf, wp;
    double tx, ty, tz;
    double uX, uY, uT;
    double vX, vY, vT;
    double wa[MAXNBEAMS];
    double wb[MAXNBEAMS];
    double Z[MAXNBEAMS];
    double u[3][MAXNBEAMS];
    int g, h;
    (void)sx;
    (void)pos;


    /* __ loop on the shells __ */
    for (h = 0; h < NhotSpotsParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesNhS(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ ellipsis parameters __ */
        uX = tx/NhotSpotsParams.axis0[h];
        uY = ty/NhotSpotsParams.axis1[h];

        /* __ ellipsis arg __ */
        uT = sqrt(uX*uX+uY*uY);

        /* __ set arg of the polynom (in X,Y rotated coordinates) __ */
        wa[h] = (uT-1)/NhotSpotsParams.widthRatio[h];

        /* __ set arg of the polynom for Z extent __ */
        wb[h] = (si->n[2] > 1) ? tz/NhotSpotsParams.thick[h] : 0.0;

        /* __ vector along magnetic field in rotated coordinates __ */
        vX = +uY/NhotSpotsParams.axis1[h];
        vY = -uX/NhotSpotsParams.axis0[h];
        vT = sqrt(vX*vX+vY*vY);

        /* __ normalisation : angle between B & azimuthal direction (given by v) __ */
        Z[h] = (uT < EPS8) ? 0.0 : NhotSpotsParams.axis1[h]*vT/uT;

        /* __ unit vector along magnetic field lines in lab coordinates __ */
        u[1][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*sin(wf)+(vY/vT)*cos(wf)*cos(wp);
        u[0][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*cos(wf)-(vY/vT)*sin(wf)*cos(wp);
        u[2][h] = (vT < EPS8) ? 0.0 :                 +(vY/vT)        *sin(wp);
        }

    /* __ magnetic field __ */
    for (g = 0; g < 3; g++) {
        B[g] = 0.0;

        for (h = 0; h < NhotSpotsParams.N; h++) {
            B[g] += NhotSpotsParams.b[h]*polynomNhS(wa[h])*u[g][h]*Z[h]*polynomNhS(wb[h]);
        }
    }

}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    nhotspotsElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void nhotspotsElectric(struct sti *si, struct stx *sx,
                      double pos[3], double E[3])
{
    (void)si;
    (void)sx;
    (void)pos;

    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
    nhotspotsCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void nhotspotsCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3])
{
    (void)si;
    (void)sx;
    (void)pos;

    J[0] = 0.;
    J[1] = 0.;
    J[2] = 0.;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    nhotspotsTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void nhotspotsTemperature(struct sti *si, struct stx *sx, double pos[3], int ispe, double T[2])
{
    (void)si;
    (void)sx;
    (void)pos;



    /* __ electron temperatures __ */
    if (ispe == 0) {
        T[PARA] = NhotSpotsParams.T[0];
        T[PERP] = NhotSpotsParams.T[0];
    }

    /* __  proton temperatures __ */
    else if (ispe == 1) {
        T[PARA] = NhotSpotsParams.T[1];
        T[PERP] = NhotSpotsParams.T[1];
    }

    else {
        printf("what the fuck is this ispe value ?\n");
    }

}
/*===========================================================================*/







/*---------------------------------------------------------------------------
    nhotspotsCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void nhotspotsCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3])
{
    (void)si;
    (void)sx;
    (void)pos;
    (void)ispe;

    curdrift[0] = 0.;
    curdrift[1] = 0.;
    curdrift[2] = 0.;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
    nhotspotsDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void nhotspotsDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3])
{
    (void)si;
    (void)sx;
    (void)pos;
    (void) ispe;


    vdrift[0] = 0.0;
    vdrift[1] = 0.0;
    vdrift[2] = 0.0;

}
/*===========================================================================*/



void nhotspotsDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{
    double argt;
    double timeA;
    double pos[3];
    double uX, uY, uT;
    double vX, vY, vT;
    double wf, wp;
    double tx, ty, tz;
    double bDrive[3];
    double wa[MAXNBEAMS];
    double wb[MAXNBEAMS];
    double Z[MAXNBEAMS];
    double u[3][MAXNBEAMS];
    int nxg1, nyg1, nzg1;
    int ijk;
    int h, i, j, k;

    /* __ # of grid points on G1 __ */
    nxg1 = sx.n[0]+1;
    nyg1 = sx.n[1]+1;
    nzg1 = sx.n[2]+1;

    /* __ time dependent amplitude factor __ */
    argt = ((float)it*si.ts)/NhotSpotsParams.tRise-1.0;
    timeA = (argt < 0.0) ? polynomNhS(argt) : 1.0;

    /* __ loop on the N shells __ */
    for (h = 0; h < NhotSpotsParams.N; h++) {

        /* __ loop on each grid points of g1 to drive magnetic field __ */
        for (i = 0; i < nxg1; i++) {
            for (j = 0; j < nyg1; j++) {
                for (k = 0; k < nzg1; k++) {
                    /* __ coordinate of grid point __ */
                    pos[0] = (i + sx.i0[0]) * (si.dl[0]);
                    pos[1] = (j + sx.i0[1]) * (si.dl[1]);
                    pos[2] = (k + sx.i0[2]) * (si.dl[2]);

                    /* __ [tx, ty, tz] : rotated coordinates __ */
                    rotateCoordinatesNhS(pos, h, &wf, &wp, &tx, &ty, &tz);

                    /* __ ellipsis parameters __ */
                    uX = tx/NhotSpotsParams.axis0[h];
                    uY = ty/NhotSpotsParams.axis1[h];

                    /* __ ellipsis arg __ */
                    uT = sqrt(uX*uX+uY*uY);

                    /* __ set arg of the polynom (in X,Y rotated coordinates) __ */
                    wa[h] = (uT-1)/NhotSpotsParams.widthRatio[h];

                    /* __ set arg of the polynom for Z extent __ */
                    wb[h] = (si.n[2] > 1) ? tz/NhotSpotsParams.thick[h] : 0.0;

                    /* __ vector along magnetic field in rotated coordinates __ */
                    vX = +uY/NhotSpotsParams.axis0[h];
                    vY = -uX/NhotSpotsParams.axis1[h];
                    vT = sqrt(vX*vX+vY*vY);

                    /* __ normalisation : angle between B & azimuthal direction (given by v) __ */
                    Z[h] = (uT < EPS8) ? 0.0 : NhotSpotsParams.axis1[h]*vT/uT;

                    /* __ unit vector along magnetic field lines in lab coordinates __ */
                    u[1][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*sin(wf)+(vY/vT)*cos(wf)*cos(wp);
                    u[0][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*cos(wf)-(vY/vT)*sin(wf)*cos(wp);
                    u[2][h] = (vT < EPS8) ? 0.0 :                 +(vY/vT)        *sin(wp);

                    /* __ magnetic field __ */
                    bDrive[0] = NhotSpotsParams.bDrive[h]*polynomNhS(wa[h])*u[0][h]*Z[h]*polynomNhS(wb[h]);
                    bDrive[1] = NhotSpotsParams.bDrive[h]*polynomNhS(wa[h])*u[1][h]*Z[h]*polynomNhS(wb[h]);
                    bDrive[2] = NhotSpotsParams.bDrive[h]*polynomNhS(wa[h])*u[2][h]*Z[h]*polynomNhS(wb[h]);

                    /* __ grid point coordinate __ */
                    ijk = IDX(i, j, k, nxg1, nyg1, nzg1);

                    /* __ then add the driven magnetic field __ */
                    switch (ipc) {
                        case 0:
                            s1[ijk].c[0] += bDrive[0];
                            s1[ijk].c[1] += bDrive[1];
                            s1[ijk].c[2] += bDrive[2];
                            break;
                        case 1:
                            s1[ijk].b[0] += bDrive[0];
                            s1[ijk].b[1] += bDrive[1];
                            s1[ijk].b[2] += bDrive[2];
                            break;
                        /* __ no reason to get there __ */
                        default : IAMDEAD(sx.r);
                    }
                }
            }
        }
    }
}


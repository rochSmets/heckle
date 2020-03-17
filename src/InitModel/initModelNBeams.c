

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


typedef struct s_nbeams
{
    int N;
    double fi[MAXNBEAMS];
    double psi[MAXNBEAMS];
    double b[MAXNBEAMS];
    double n[NS+1][MAXNBEAMS];
    double V[NS+1][MAXNBEAMS];
    double T[NS+1][MAXNBEAMS];
    double axisX[MAXNBEAMS];
    double axisY[MAXNBEAMS];
    double widthRatio[MAXNBEAMS];
    double thick[MAXNBEAMS];
    double center[MAXNBEAMS][3];
} NBeams;




/* this is for private use only */
NBeams NBeamsParams;





/*---------------------------------------------------------------------------
    nbeams_start()
  ---------------------------------------------------------------------------
    AIM : read the parameters for nbeams model
 ---------------------------------------------------------------------------*/
void nbeams_start(struct sti *si, struct stx *sx, char *dir)
{
    int g, h;
    char sfh[80];
    char junk[100];
    FILE *fp;

    (void)si;


    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "Nbeams.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL) {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscanint(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.N), &(sx->irun));
    if (NBeamsParams.N > MAXNBEAMS) {
        printf("too many beams in this initialization : increase MAXNBEAMS\n");
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.fi[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.psi[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.b[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.n[1][h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.n[2][0]), &(sx->irun));
    for (h = 1; h < NBeamsParams.N; h++) {
        NBeamsParams.n[2][h] = NBeamsParams.n[2][0];
    }

    for (h = 0; h < NBeamsParams.N; h++) {
        NBeamsParams.V[0][h] = 0.0;
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.V[1][h]), &(sx->irun));
    }

    for (h = 0; h < NBeamsParams.N; h++) {
        NBeamsParams.V[2][h] = 0.0;
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.T[0][0]), &(sx->irun));
    for (h = 1; h < NBeamsParams.N; h++) {
        NBeamsParams.T[0][h] = NBeamsParams.T[0][0];
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.T[1][h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.T[2][0]), &(sx->irun));
    for (h = 1; h < NBeamsParams.N; h++) {
        NBeamsParams.T[2][h] = NBeamsParams.T[2][0];
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.axisX[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.axisY[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.widthRatio[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.thick[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.center[h][0]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.center[h][1]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NBeamsParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NBeamsParams.center[h][2]), &(sx->irun));
    }


/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : N beams _________\n");

   for (h = 0; h < NBeamsParams.N; h++) {
       printf("... Beam #       :%12d\n", h);
       printf("fi  [0, 360]     :%12.4f\n", NBeamsParams.fi[h]);
       printf("psi [0,  90]     :%12.4f\n", NBeamsParams.psi[h]);
       printf("B field          :%12.4f\n", NBeamsParams.b[h]);
       printf("\n");

       for (g = 0; g < 3; g++) {
           printf("beam center  [%1d] :%12.4f\n", g, NBeamsParams.center[h][g]);
       }
       printf("axis X           :%12.4f\n", NBeamsParams.axisX[h]);
       printf("axis Y           :%12.4f\n", NBeamsParams.axisY[h]);
       printf("width ratio      :%12.4f\n", NBeamsParams.widthRatio[h]);
       printf("thickness (z)    :%12.4f\n", NBeamsParams.thick[h]);
       printf("\n");

       for (g = 0; g < NS+1; g++) {
           printf("___ specie %1d ___\n", g);
           printf("... density      :%12.4f\n", NBeamsParams.n[g][h]);
           printf("... fluid vel.   :%12.4f\n", NBeamsParams.V[g][h]);
           printf("... temperature  :%12.4f\n", NBeamsParams.T[g][h]);
           printf("\n");
       }
   }

   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ polynom used for the interpolation __ */
double polynomNBeams(double x)
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
void rotateCoordinatesNB(double pos[3], int h, double* wf, double* wp, double* tx, double* ty, double *tz)
{
double wx, wy, wz;


/* __ set fi & psi angles __ */
*wf = NBeamsParams.fi[h] *PI/180.0;
*wp = NBeamsParams.psi[h]*PI/180.0;

/* __ rotation with fi angle __ */
wx = (pos[0]-NBeamsParams.center[h][0])*cos(*wf)
    +(pos[1]-NBeamsParams.center[h][1])*sin(*wf);
wy =-(pos[0]-NBeamsParams.center[h][0])*sin(*wf)
    +(pos[1]-NBeamsParams.center[h][1])*cos(*wf);
wz =  pos[2]-NBeamsParams.center[h][2];

/* __ rotation with psi angle __ */
*tx =  wx;
*ty =  wy*cos(*wp)+wz*sin(*wp);
*tz = -wy*sin(*wp)+wz*cos(*wp);

}



/*---------------------------------------------------------------------------
    nbeamsDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for the given specie (ispe)
    at the given position (pos)
 ---------------------------------------------------------------------------*/
double nbeamsDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe)
{
    double wf, wp;
    double tx, ty, tz;
    double uX, uY;
    double wa[MAXNBEAMS];
    double n_;
    int h;
    (void)si;
    (void)sx;
    (void)pos;


    /* __ loop on the N shells __ */
    for (h = 0; h < NBeamsParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesNB(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ ellipsis parameters __ */
        uX = tx/NBeamsParams.axisX[h];
        uY = ty/NBeamsParams.axisY[h];

        /* __ set arg of the polynom __ */
        wa[h] = sqrt(pow(uX, 2)+pow(uY, 2));
        }

    if (ispe == 0) {
        n_ = NBeamsParams.n[2][0];

        for (h = 0; h < NBeamsParams.N; h++) {
            n_ += NBeamsParams.n[1][h]*polynomNBeams(wa[h]);
        }

        return n_;
    }

    else if (ispe == 1) {
        n_ = 0.0;

        for (h = 0; h < NBeamsParams.N; h++) {
            n_ += NBeamsParams.n[1][h]*polynomNBeams(wa[h]);
        }

        return n_;
    }

    else if (ispe == 2) {
        return NBeamsParams.n[2][0];
    }

    else {
        printf("what the fuck is this ispe value ?\n");
        return 0;
    }

}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    nbeamsMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void nbeamsMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3])
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
    (void)si;
    (void)sx;
    (void)pos;


    /* __ loop on the 2 shells __ */
    for (h = 0; h < NBeamsParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesNB(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ ellipsis parameters __ */
        uX = tx/NBeamsParams.axisX[h];
        uY = ty/NBeamsParams.axisY[h];

        /* __ ellipsis arg __ */
        uT = sqrt(uX*uX+uY*uY);

        /* __ set arg of the polynom (in X,Y rotated coordinates) __ */
        wa[h] = (uT-1)/NBeamsParams.widthRatio[h];

        /* __ set arg of the polynom for Z extent __ */
        wb[h] = tz/NBeamsParams.thick[h];

        /* __ vector along magnetic field in rotated coordinates __ */
        vX = +uY/NBeamsParams.axisY[h];
        vY = -uX/NBeamsParams.axisX[h];
        vT = sqrt(vX*vX+vY*vY);

        /* __ normalisation : angle between B & azimuthal direction (given by v) __ */
        Z[h] = (uT < EPS8) ? 0.0 : NBeamsParams.axisY[h]*vT/uT;

        /* __ unit vector along magnetic field lines in lab coordinates __ */
        u[1][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*sin(wf)+(vY/vT)*cos(wf)*cos(wp);
        u[0][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*cos(wf)-(vY/vT)*sin(wf)*cos(wp);
        u[2][h] = (vT < EPS8) ? 0.0 :                 +(vY/vT)        *sin(wp);
        }

    /* __ magnetic field __ */
    for (g = 0; g < 3; g++) {
        B[g] = 0.0;

        for (h = 0; h < NBeamsParams.N; h++) {
            B[g] += NBeamsParams.b[h]*polynomNBeams(wa[h])*u[g][h]*Z[h]*polynomNBeams(wb[h]);
        }
    }

}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    nbeamsElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void nbeamsElectric(struct sti *si, struct stx *sx,
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
    nbeamsCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void nbeamsCurrent(struct sti *si, struct stx *sx,
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
    nbeamsTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void nbeamsTemperature(struct sti *si, struct stx *sx, double pos[3], int ispe, double T[2])
{
    double wf, wp;
    double tx, ty, tz;
    double uX, uY;
    double wa[MAXNBEAMS];
    double n_[MAXNBEAMS];
    double nTot;
    double P;
    int h;
    (void)si;
    (void)sx;
    (void)pos;


    /* __ loop on the 2 shells __ */
    for (h = 0; h < NBeamsParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesNB(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ ellipsis parameters __ */
        uX = tx/NBeamsParams.axisX[h];
        uY = ty/NBeamsParams.axisY[h];

        /* __ set arg of the polynom __ */
        wa[h] = sqrt(pow(uX, 2)+pow(uY, 2));
        }

    /* __ electron temperatures __ */
    if (ispe == 0) {
        T[PARA] = NBeamsParams.T[0][0];
        T[PERP] = NBeamsParams.T[0][0];
    }

    /* __  proton temperatures __ */
    else if (ispe == 1) {
        nTot = 0.0;
        for (h = 0; h < NBeamsParams.N; h++) {
            n_[h] = NBeamsParams.n[1][h]*polynomNBeams(wa[h]);
            nTot += n_[h];
        }

        if (nTot != 0.0) {

            P = 0.0;

            for (h = 0; h < NBeamsParams.N; h++) {
                P += n_[h]*NBeamsParams.T[1][h];
            }

            T[PARA] = P/nTot;
            T[PERP] = T[PARA];
        }

        else {
            T[PARA] = 0.0;
            T[PERP] = 0.0;
        }
    }

    /* __  alpha temperatures __ */
    else if (ispe == 2) {
        T[PARA] = NBeamsParams.T[2][0];
        T[PERP] = NBeamsParams.T[2][0];
    }

    else {
        printf("what the fuck is this ispe value ?\n");
    }

}
/*===========================================================================*/







/*---------------------------------------------------------------------------
    nbeamsCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void nbeamsCurDrift(struct sti *si, struct stx *sx, double pos[3],
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
    nbeamsDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void nbeamsDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3])
{
    double wf, wp;
    double tx, ty, tz;
    double uX, uY, uT;
    double vX, vY, vT;
    double wa[MAXNBEAMS];
    double u[3][MAXNBEAMS];
    int g, h;
    (void)si;
    (void)sx;
    (void)pos;


    /* __ loop on the 2 shells __ */
    for (h = 0; h < NBeamsParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesNB(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ ellipsis parameters __ */
        uX = tx/NBeamsParams.axisX[h];
        uY = ty/NBeamsParams.axisY[h];

        /* __ ellipsis arg __ */
        uT = sqrt(uX*uX+uY*uY);

        /* __ set arg of the polynom __ */
        wa[h] = 2*uT-1;

        /* __ vector along magnetic field in rotated coordinates __ */
        vX = +uX/NBeamsParams.axisX[h];
        vY = +uY/NBeamsParams.axisY[h];
        vT = sqrt(vX*vX+vY*vY);

        /* __ unit vector along magnetic field lines in lab coordinates __ */
        u[0][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*cos(wf)-(vY/vT)*sin(wf)*cos(wp);
        u[1][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*sin(wf)+(vY/vT)*cos(wf)*cos(wp);
        u[2][h] = (vT < EPS8) ? 0.0 :                 +(vY/vT)        *sin(wp);
        }

    /* __ electron velocity __ */
    if (ispe == 0) {
        vdrift[0] = 0.0;
        vdrift[1] = 0.0;
        vdrift[2] = 0.0;
    }

    /* __ proton velocity __ */
    else if (ispe == 1) {
        for (g = 0; g < 3; g++) {
            vdrift[g] = 0.0;

            for (h = 0; h < NBeamsParams.N; h++) {
                vdrift[g] += NBeamsParams.V[1][h]*polynomNBeams(wa[h])*u[g][h];
            }

        }
    }

    /* __ alpha velocity __ */
    else if (ispe == 2) {
        vdrift[0] = 0.0;
        vdrift[1] = 0.0;
        vdrift[2] = 0.0;
    }

    else {
        printf("what the fuck is this ispe value ?\n");
    }

}
/*===========================================================================*/



void nbeamsDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

}


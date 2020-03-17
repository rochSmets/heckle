

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


typedef struct s_lasers
{
    int N;
    double fi[MAXNBEAMS];
    double psi[MAXNBEAMS];
    double b[MAXNBEAMS];
    double n[NS+1][MAXNBEAMS];
    double V[NS+1][MAXNBEAMS];
    double T[NS+1][MAXNBEAMS];
    double center[MAXNBEAMS][2];
    double ellipsis[MAXNBEAMS];
    double bAxis0[MAXNBEAMS];
    double bWidthRatio[MAXNBEAMS];
    double tRise;
    double jDrive[MAXNBEAMS];
    double jAxis0[MAXNBEAMS];
    double bDrive[MAXNBEAMS];
    double cAxis0[MAXNBEAMS];
    double cWidthRatio[MAXNBEAMS];
} NLasers;



/* __ for private use only __ */
NLasers NLasersParams;



/* __ read the parameters for lasers model __ */
void nlasers_start(struct sti *si, struct stx *sx, char *dir)
{
    int g, h;
    char sfh[80];
    char junk[100];
    FILE *fp;

    (void)si;


    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "Nlasers.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL) {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscanint(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.N), &(sx->irun));
    if (NLasersParams.N > MAXNBEAMS) {
        printf("too many beams in this initialization : increase MAXNBEAMS\n");
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.fi[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.psi[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.b[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.n[1][h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.n[2][0]), &(sx->irun));
    for (h = 1; h < NLasersParams.N; h++) {
        NLasersParams.n[2][h] = NLasersParams.n[2][0];
    }

    for (h = 0; h < NLasersParams.N; h++) {
        NLasersParams.V[0][h] = 0.0;
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        NLasersParams.V[0][h] = 0.0;
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.V[1][h]), &(sx->irun));
        NLasersParams.V[2][h] = 0.0;
    }

    for (h = 0; h < NLasersParams.N; h++) {
        NLasersParams.V[2][h] = 0.0;
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.T[0][0]), &(sx->irun));
    for (h = 1; h < NLasersParams.N; h++) {
        NLasersParams.T[0][h] = NLasersParams.T[0][0];
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.T[1][h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.T[2][0]), &(sx->irun));
    for (h = 1; h < NLasersParams.N; h++) {
        NLasersParams.T[2][h] = NLasersParams.T[2][0];
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.center[h][0]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.center[h][1]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.ellipsis[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.bAxis0[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.bWidthRatio[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.tRise), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.jDrive[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.jAxis0[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.bDrive[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.cAxis0[h]), &(sx->irun));
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < NLasersParams.N; h++) {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(NLasersParams.cWidthRatio[h]), &(sx->irun));
    }

/* __ display the topology parameters __ */
if (sx->r == 0)
   {
   printf("________________ magnetic topology : N lasers ________\n");

   printf("\n");
   printf("tRise            :%12.4f\n", NLasersParams.tRise);
   printf("\n");

   for (h = 0; h < NLasersParams.N; h++) {
       printf("... Beam #       :%12d\n", h);
       printf("fi  [0, 360]     :%12.4f\n", NLasersParams.fi[h]);
       printf("psi [0,  90]     :%12.4f\n", NLasersParams.psi[h]);
       printf("B field          :%12.4f\n", NLasersParams.b[h]);
       printf("\n");

       for (g = 0; g < 2; g++) {
           printf("beam center  [%1d] :%12.4f\n", g, NLasersParams.center[h][g]);
       }
       printf("ellipticity      :%12.4f\n", NLasersParams.ellipsis[h]);
       printf("B axis 0         :%12.4f\n", NLasersParams.bAxis0[h]);
       printf("B width ratio    :%12.4f\n", NLasersParams.bWidthRatio[h]);
       printf("\n");

       printf("J drive (z)      :%12.4f\n", NLasersParams.jDrive[h]);
       printf("J axis 0         :%12.4f\n", NLasersParams.jAxis0[h]);
       printf("\n");

       printf("B drive          :%12.4f\n", NLasersParams.bDrive[h]);
       printf("B axis X 0       :%12.4f\n", NLasersParams.cAxis0[h]);
       printf("B width ratio    :%12.4f\n", NLasersParams.cWidthRatio[h]);
       printf("\n");

       for (g = 0; g < NS+1; g++) {
           printf("___ specie %1d ___\n", g);
           printf("... density      :%12.4f\n", NLasersParams.n[g][h]);
           printf("... fluid vel.   :%12.4f\n", NLasersParams.V[g][h]);
           printf("... temperature  :%12.4f\n", NLasersParams.T[g][h]);
           printf("\n");
       }
   }

   printf("______________________________________________________\n");
   printf("\n");
   }

}


/* __ polynom used for the interpolation __ */
double polynomLasers(double x)
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
void rotateCoordinatesL(double pos[3], int h, double* wf, double* wp, double* tx, double* ty, double *tz)
{
double wx, wy, wz;


/* __ set fi & psi angles __ */
*wf = NLasersParams.fi[h] *PI/180.0;
*wp = NLasersParams.psi[h]*PI/180.0;

/* __ rotation with fi angle __ */
wx = (pos[0]-NLasersParams.center[h][0])*cos(*wf)
    +(pos[1]-NLasersParams.center[h][1])*sin(*wf);
wy =-(pos[0]-NLasersParams.center[h][0])*sin(*wf)
    +(pos[1]-NLasersParams.center[h][1])*cos(*wf);
wz =  pos[2];

/* __ rotation with psi angle __ */
*tx =  wx;
*ty =  wy*cos(*wp)+wz*sin(*wp);
*tz = -wy*sin(*wp)+wz*cos(*wp);

}


/* __ density for the given specie (ispe) at the given position (pos) __ */
double nlasersDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe)
{
    double wf, wp;
    double tx, ty, tz;
    double axis0, axis1;
    double uX, uY;
    double wa[MAXNBEAMS];
    double n_;
    int h;
    (void)si;
    (void)sx;
    (void)pos;


    /* __ loop on the N shells __ */
    for (h = 0; h < NLasersParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesL(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ semi axis (major  minor) __ */
        axis0 = NLasersParams.bAxis0[h];
        axis1 = axis0*NLasersParams.ellipsis[h];

        /* __ ellipsis parameters __ */
        uX = tx/axis0;
        uY = ty/axis1;

        /* __ set arg of the polynom __ */
        wa[h] = sqrt(pow(uX, 2)+pow(uY, 2));
        }

    if (ispe == 0) {
        n_ = NLasersParams.n[2][0];

        for (h = 0; h < NLasersParams.N; h++) {
            n_ += NLasersParams.n[1][h]*polynomLasers(wa[h]);
        }

        return n_;
    }

    else if (ispe == 1) {
        n_ = 0.0;

        for (h = 0; h < NLasersParams.N; h++) {
            n_ += NLasersParams.n[1][h]*polynomLasers(wa[h]);
        }

        return n_;
    }

    else if (ispe == 2) {
        return NLasersParams.n[2][0];
    }

    else {
        printf("what the fuck is this ispe value ?\n");
        return 0;
    }

}


/* __ returns the magnetic field at the given position (pos) __ */
void nlasersMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3])
{
    double wf, wp;
    double tx, ty, tz;
    double axis0, axis1;
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


    /* __ loop on the N shells __ */
    for (h = 0; h < NLasersParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesL(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ semi axis (major  minor) __ */
        axis0 = NLasersParams.bAxis0[h];
        axis1 = axis0*NLasersParams.ellipsis[h];

        /* __ ellipsis parameters __ */
        uX = tx/axis0;
        uY = ty/axis1;

        /* __ ellipsis arg __ */
        uT = sqrt(uX*uX+uY*uY);

        /* __ set arg of the polynom (in X,Y rotated coordinates) __ */
        wa[h] = (uT-1)/NLasersParams.bWidthRatio[h];

        /* __ vector along magnetic field in rotated coordinates __ */
        vX = +uY/axis1;
        vY = -uX/axis0;
        vT = sqrt(vX*vX+vY*vY);

        /* __ normalisation : angle between B & azimuthal direction (given by v) __ */
        Z[h] = (uT < EPS8) ? 0.0 : axis1*vT/uT;

        /* __ unit vector along magnetic field lines in lab coordinates __ */
        u[1][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*sin(wf)+(vY/vT)*cos(wf)*cos(wp);
        u[0][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*cos(wf)-(vY/vT)*sin(wf)*cos(wp);
        u[2][h] = (vT < EPS8) ? 0.0 :                 +(vY/vT)        *sin(wp);
        }

    /* __ magnetic field __ */
    for (g = 0; g < 3; g++) {
        B[g] = 0.0;

        for (h = 0; h < NLasersParams.N; h++) {
            B[g] += NLasersParams.b[h]*polynomLasers(wa[h])*u[g][h]*Z[h];
        }
    }

}


/* __ electric field at the given position (pos) __ */
void nlasersElectric(struct sti *si, struct stx *sx,
                      double pos[3], double E[3])
{
    (void)si;
    (void)sx;
    (void)pos;

    E[0] = 0.;
    E[1] = 0.;
    E[2] = 0.;
}


/* __ current density at the given position (pos) __ */
void nlasersCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3])
{
    (void)si;
    (void)sx;
    (void)pos;

    J[0] = 0.0;
    J[1] = 0.0;
    J[2] = 0.0;
}


/* __ temperature for the given specie (ispe) at the given position (pos) __ */
void nlasersTemperature(struct sti *si, struct stx *sx, double pos[3], int ispe, double T[2])
{
    double wf, wp;
    double tx, ty, tz;
    double axis0, axis1;
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
    for (h = 0; h < NLasersParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesL(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ semi axis (major  minor) __ */
        axis0 = NLasersParams.bAxis0[h];
        axis1 = axis0*NLasersParams.ellipsis[h];

        /* __ ellipsis parameters __ */
        uX = tx/axis0;
        uY = ty/axis1;

        /* __ set arg of the polynom __ */
        wa[h] = sqrt(pow(uX, 2)+pow(uY, 2));
        }

    /* __ electron temperatures __ */
    if (ispe == 0) {
        T[PARA] = NLasersParams.T[0][0];
        T[PERP] = NLasersParams.T[0][0];
    }

    /* __  proton temperatures __ */
    else if (ispe == 1) {
        nTot = 0.0;
        for (h = 0; h < NLasersParams.N; h++) {
            n_[h] = NLasersParams.n[1][h]*polynomLasers(wa[h]);
            nTot += n_[h];
        }

        if (nTot != 0.0) {
            P = 0.0;
            for (h = 0; h < NLasersParams.N; h++) {
                P += n_[h]*NLasersParams.T[1][h];
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
        T[PARA] = NLasersParams.T[2][0];
        T[PERP] = NLasersParams.T[2][0];
    }

    else {
        printf("what the fuck is this ispe value ?\n");
    }

}


/* __ components of the drift vel. associated with the current __ */
void nlasersCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3])
{
    (void)si;
    (void)sx;
    (void)pos;
    (void)ispe;

    curdrift[0] = 0.0;
    curdrift[1] = 0.0;
    curdrift[2] = 0.0;
}


/* __ drift velocity at a given position (pos) __ */
void nlasersDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3])
{
    double wf, wp;
    double tx, ty, tz;
    double axis0, axis1;
    double uX, uY, uT;
    double vX, vY, vT;
    double wa[MAXNBEAMS];
    double u[3][MAXNBEAMS];
    int g, h;
    (void)si;
    (void)sx;
    (void)pos;


    /* __ loop on the 2 shells __ */
    for (h = 0; h < NLasersParams.N; h++)
        {
        /* __ [tx, ty, tz] : rotated coordinates __ */
        rotateCoordinatesL(pos, h, &wf, &wp, &tx, &ty, &tz);

        /* __ semi axis (major  minor) __ */
        axis0 = NLasersParams.bAxis0[h];
        axis1 = axis0*NLasersParams.ellipsis[h];

        /* __ ellipsis parameters __ */
        uX = tx/axis0;
        uY = ty/axis1;

        /* __ ellipsis arg __ */
        uT = sqrt(uX*uX+uY*uY);

        /* __ set arg of the polynom __ */
        wa[h] = 2*uT-1;

        /* __ vector along magnetic field in rotated coordinates __ */
        vX = +uX/axis0;
        vY = +uY/axis1;
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

            for (h = 0; h < NLasersParams.N; h++) {
                vdrift[g] += NLasersParams.V[1][h]*polynomLasers(wa[h])*u[g][h];
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


/* __ driver for the axial current and magnetic field __ */
void nlasersDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{
    double wf, wp;
    double tx, ty, tz;
    double axis0, axis1;
    double argt;
    double timeA;
    double pos[3];
    double uX, uY, uT;
    double vX, vY, vT;
    double driverJ;
    double jDrive[3];
    double bDrive[3];
    double wa[MAXNBEAMS];
    double wb[MAXNBEAMS];
    double Z[MAXNBEAMS];
    double u[3][MAXNBEAMS];
    int nxg1, nyg1, nzg1;
    int nxg2, nyg2, nzg2;
    int g, h, i, j, k;
    int ijk;


    /* __ # of grid points on G1 __ */
    nxg1 = sx.n[0]+1;
    nyg1 = sx.n[1]+1;
    nzg1 = sx.n[2]+1;

    /* __ # of grid points on G2 __ */
    nxg2 = sx.n[0]+2;
    nyg2 = sx.n[1]+2;
    nzg2 = sx.n[2]+2;

    /* __ time dependent amplitude factor __ */
    argt = ((float)it*si.ts)/NLasersParams.tRise-1.0;
    timeA = (argt < 0.0) ? polynomLasers(argt) : 1.0;

    /* __ loop on the N shells __ */
    for (h = 0; h < NLasersParams.N; h++)
    {
        /* __ loop on each grid points of g2 to drive current __ */
        for (i = 0; i < nxg2; i++)
        {
            for (j = 0; j < nyg2; j++)
            {
                for (k = 0; k < nzg2; k++)
                {
                    /* __ coordinate of grid point __ */
                    pos[0] = (i-0.5 + sx.i0[0]) * (si.dl[0]);
                    pos[1] = (j-0.5 + sx.i0[1]) * (si.dl[1]);
                    pos[2] = (k-0.5 + sx.i0[2]) * (si.dl[2]);

                    /* __ [tx, ty, tz] : rotated coordinates __ */
                    rotateCoordinatesL(pos, h, &wf, &wp, &tx, &ty, &tz);

                    /* __ semi axis (major  minor) __ */
                    axis0 = NLasersParams.jAxis0[h];
                    axis1 = axis0*NLasersParams.ellipsis[h];

                    /* __ ellipsis parameters __ */
                    uX = tx/axis0;
                    uY = ty/axis1;

                    /* __ ellipsis arg __ */
                    uT = sqrt(uX*uX+uY*uY);

                    /* __ driver amplitude depending on space and time __ */
                    driverJ = (uT < 1.0) ? -timeA : 0.0;

                    /* __ driven current __ */
                    jDrive[0] = - driverJ*NLasersParams.jDrive[h]*sin(wp)*sin(wf);
                    jDrive[1] =   driverJ*NLasersParams.jDrive[h]*sin(wp)*cos(wf);
                    jDrive[2] =   driverJ*NLasersParams.jDrive[h]*cos(wp);

                    /* __ grid point coordinate __ */
                    ijk = IDX(i, j, k, nxg2, nyg2, nzg2);

                    /* __ then add the driven current __ */
                    s2[ijk].j[0] += jDrive[0];
                    s2[ijk].j[1] += jDrive[1];
                    s2[ijk].j[2] += jDrive[2];
                }
            }
        }

        /* __ loop on each grid points of g1 to drive magnetic field __ */
        for (i = 0; i < nxg1; i++)
        {
            for (j = 0; j < nyg1; j++)
            {
                for (k = 0; k < nzg1; k++)
                {
                    /* __ coordinate of grid point __ */
                    pos[0] = (i + sx.i0[0]) * (si.dl[0]);
                    pos[1] = (j + sx.i0[1]) * (si.dl[1]);
                    pos[2] = (k + sx.i0[2]) * (si.dl[2]);

                    /* __ [tx, ty, tz] : rotated coordinates __ */
                    rotateCoordinatesL(pos, h, &wf, &wp, &tx, &ty, &tz);

                    /* __ semi axis (major  minor) __ */
                    axis0 = NLasersParams.cAxis0[h];
                    axis1 = axis0*NLasersParams.ellipsis[h];

                    /* __ ellipsis parameters __ */
                    uX = tx/axis0;
                    uY = ty/axis1;

                    /* __ ellipsis arg __ */
                    uT = sqrt(uX*uX+uY*uY);

                    /* __ set arg of the polynom (in X,Y rotated coordinates) __ */
                    wa[h] = (uT-1)/NLasersParams.cWidthRatio[h];

                    /* __ vector along magnetic field in rotated coordinates __ */
                    vX = +uY/axis0;
                    vY = -uX/axis1;
                    vT = sqrt(vX*vX+vY*vY);

                    /* __ normalisation : angle between B & azimuthal direction (given by v) __ */
                    Z[h] = (uT < EPS8) ? 0.0 : axis1*vT/uT;

                    /* __ unit vector along magnetic field lines in lab coordinates __ */
                    u[1][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*sin(wf)+(vY/vT)*cos(wf)*cos(wp);
                    u[0][h] = (vT < EPS8) ? 0.0 : +(vX/vT)*cos(wf)-(vY/vT)*sin(wf)*cos(wp);
                    u[2][h] = (vT < EPS8) ? 0.0 :                 +(vY/vT)        *sin(wp);

                    /* __ magnetic field __ */
                    bDrive[0] = NLasersParams.bDrive[h]*polynomLasers(wa[h])*u[0][h]*Z[h];
                    bDrive[1] = NLasersParams.bDrive[h]*polynomLasers(wa[h])*u[1][h]*Z[h];
                    bDrive[2] = NLasersParams.bDrive[h]*polynomLasers(wa[h])*u[2][h]*Z[h];

                    /* __ grid point coordinate __ */
                    ijk = IDX(i, j, k, nxg1, nyg1, nzg1);

                    /* __ then add the driven magnetic field __ */
                    switch (ipc)
                        {
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

    if (sx.r == 0 && ipc == 0) {
        printf("________________ driven system : %12.8f ________\n", timeA);
    }
}


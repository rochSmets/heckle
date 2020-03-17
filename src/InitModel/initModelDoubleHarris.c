

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"


#define TPARA 0
#define TPERP 1


typedef struct s_doubleharris
{
    double nb;
    double l;
    double teti;
    double perturb;
} Harris;




/* this is for private use only */
static Harris DoubleHarrisParams;





/*---------------------------------------------------------------------------
    doubleharrisStart()
  ---------------------------------------------------------------------------
    AIM : read the parameters for the Double Harris Sheet model
 ---------------------------------------------------------------------------*/
void doubleharris_start(struct sti *si, struct stx *sx, char *dir)
{
    char sfh[80];
    char junk[100];
    FILE *fp;



    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "doubleharris.txt");

    // unused but should be defined
    (void) si;

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL)
    {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(DoubleHarrisParams.nb), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(DoubleHarrisParams.l), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(DoubleHarrisParams.teti), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(DoubleHarrisParams.perturb), &(sx->irun));

    fclose(fp);

    /* __ display the topology parameters __ */
    if (sx->r == 0)
    {
       printf("___ magnetic topology : double harris ________________\n");
       printf("background density       : %8.4f\n", DoubleHarrisParams.nb);
       printf("sheet thickness          : %8.4f\n", DoubleHarrisParams.l);
       printf("electron/ion temeprature : %8.4f\n", DoubleHarrisParams.teti);
       printf("mag. perturb.            : %8.4f\n", DoubleHarrisParams.perturb);
       printf("______________________________________________________\n");
       printf("\n");
    }
}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    doubleharrisDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for all species at the given position (pos)
 ---------------------------------------------------------------------------*/
double doubleharrisDensity(struct sti *si, struct stx *sx, double pos[3], int ispe)
{
    double y1, y2;
    double nb, n, l;

    //unused but should be defined
    (void)sx;
    (void)ispe;

    nb = DoubleHarrisParams.nb;
    l  = DoubleHarrisParams.l;

    /* distances from the 2 current sheets */
    y1 = (pos[1] - 0.25*  si->l[1]) / l;
    y2 = (pos[1] - 0.75 * si->l[1]) / l;

    /* harris sheet density */
    n = nb + 1./(cosh(y1)*cosh(y1)) + 1./(cosh(y2)*cosh(y2));
    //n = 1;
    return n;
}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    doubleharrisMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void doubleharrisMagnetic(struct sti *si, struct stx *sx,
                                double pos[3], double B[3])
{
    double x0, y1, y2;
    const double w1 = 0.1, w2 = 2.0;
    double w3, w4, w5;
    double l;

    (void)sx;

    l = DoubleHarrisParams.l;

    /* distances from the 2 current sheets */
    y1 = (pos[1] - 0.25 * si->l[1]) / l;
    y2 = (pos[1] - 0.75 * si->l[1]) / l;

    /* distance from the middle of the box in x direction */
    x0 = (pos[0] - 0.5 * si->l[0]) / l;

    /* constants for the magnetic perturbation from seiji */
    w3 = exp(-(x0*x0 + y1*y1) / (w2*w2));
    w4 = exp(-(x0*x0 + y2*y2) / (w2*w2));
    w5 = 2.0*w1/w2;

    B[0] = tanh(y1) - tanh(y2) - 1.0 + (-w5*y1*w3) + (+w5*y2*w4);
    B[1] = ( (w5 * x0 * w3) + ( -w5 * x0 * w4));
    B[2] = 0.0;
}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    doubleharrisElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void doubleharrisElectric(struct sti *si, struct stx *sx,
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
    doubleharrisCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void doubleharrisCurrent(struct sti *si, struct stx *sx,
                               double pos[3], double J[3])
{
    double y1, y2;
    double l;

    (void)sx;

    l = DoubleHarrisParams.l;

    y1 = (pos[1] - 0.25*si->l[1]);
    y2 = (pos[1] - 0.75*si->l[1]);

    J[0] = 0.;
    J[1] = 0.;
    J[2] = - 1./(l * cosh(y1) * cosh(y1)) + 1./(l * cosh(y2) * cosh(y2));
}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    doubleharrisTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void doubleharrisTemperature(struct sti *si, struct stx *sx,
                                   double pos[3], int ispe, double T[2])
{
    double pmag;
    double teti;

    (void)si;
    (void)sx;
    (void)pos;

    teti = DoubleHarrisParams.teti;

    /* asymptotic magnetic field pressure */
    pmag = 0.5;

    if (ispe == 0) // electrons
    {
        T[TPARA] = pmag * teti / ((1 + teti));
        T[TPERP] = pmag * teti / ((1 + teti));
    }

    else //if (ispe == 1) // protons
    {
        T[TPARA] = pmag / ((1 + teti));
        T[TPERP] = pmag / ((1 + teti));
    }
}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    doubleharrisCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void doubleharrisCurDrift(struct sti *si, struct stx *sx, double pos[3],
                                int ispe, double curdrift[3])
{
    double n;
    int s;
    double T[2], J[3];
    double Ttot;

    Ttot = 0.;
    for (s=0; s < NS+1; s++)
    {
        doubleharrisTemperature(si, sx, pos, s, T);
        Ttot += T[TPARA] + 2 * T[TPERP];
    }

    doubleharrisCurrent(si, sx, pos, J);
    n = doubleharrisDensity(si, sx, pos, ispe);

    doubleharrisTemperature(si, sx, pos, ispe, T);

    curdrift[0] = 1./Ttot * ( T[TPARA] + 2 * T[TPERP]) * J[0]/(n);
    curdrift[1] = 1./Ttot * ( T[TPARA] + 2 * T[TPERP]) * J[1]/(n);
    curdrift[2] = 1./Ttot * ( T[TPARA] + 2 * T[TPERP]) * J[2]/(n);

    if (ispe == 0) /* electrons charge = -1 */
    {
        curdrift[0] *= -1;
        curdrift[1] *= -1;
        curdrift[2] *= -1;
    }
}
/*===========================================================================*/



/*---------------------------------------------------------------------------
    doubleharrisDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void doubleharrisDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3])
{
    (void)si;
    (void)sx;
    (void)pos;
    (void)ispe;

    vdrift[0] = 0.;
    vdrift[1] = 0.;
    vdrift[2] = 0.;
}
/*===========================================================================*/


void doubleharrisDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

}


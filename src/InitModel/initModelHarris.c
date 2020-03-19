
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
//#include "plasma.h"


#define TPARA 0
#define TPERP 1


typedef struct s_harris
{
    double nb;
    double l;
    double tetp;
    double teta;
    double perturb;
    double bg;
} Harris;




/* this is for private use only */
static Harris HarrisParams;





/*---------------------------------------------------------------------------
    harrisStart()
  ---------------------------------------------------------------------------
    AIM : read the parameters for the Double Harris Sheet model
 ---------------------------------------------------------------------------*/
void harris_start(struct sti *si, struct stx *sx, char *dir)
{
    char sfh[80];
    char junk[100];
    FILE *fp;


    // unused but should be defined
    (void) si;

    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "harris.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL)
    {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(HarrisParams.nb), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(HarrisParams.l), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(HarrisParams.tetp), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(HarrisParams.teta), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(HarrisParams.perturb), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(HarrisParams.bg), &(sx->irun));

    fclose(fp);

    /* __ display the topology parameters __ */
    if (sx->r == 0)
    {
       printf("________________ magnetic topology : harris ________\n");
       printf("background density       : %8.4f\n", HarrisParams.nb);
       printf("sheet thickness          : %8.4f\n", HarrisParams.l);
       printf("electron/proton temp.    : %8.4f\n", HarrisParams.tetp);
       printf("electron/alpha temp.     : %8.4f\n", HarrisParams.teta);
       printf("mag. perturb.            : %8.4f\n", HarrisParams.perturb);
       printf("guide field              : %8.4f\n", HarrisParams.bg);
       printf("______________________________________________________\n");
       printf("\n");
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
    harrisDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for all species at the given position (pos)
 ---------------------------------------------------------------------------*/
double harrisDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe)
{
    double y0;
    double n;

    //unused but should be defined
    (void)sx;
    (void)ispe;


    if (ispe == 1)
    {
        y0 = (pos[1] - 0.5*  si->l[1]) / HarrisParams.l;
        n = 1./(cosh(y0)*cosh(y0));
    }
    else if (ispe == 2)
    {
        n = HarrisParams.nb;
    }
    else
    {
        printf("what the fuck is this specie ?\n");
    }

    return n;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    harrisMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void harrisMagnetic(struct sti *si, struct stx *sx,
                    double pos[3], double B[3])
{
    double x0, y0;
    const double w1 = HarrisParams.perturb, w2 = 2.0;
    double w3, w5;
    double l;

    (void)sx;

    l = HarrisParams.l;

    /* center of the domain */
    y0 = (pos[1] - 0.5*  si->l[1]) / l;

    /* distance from the middle of the box in x direction */
    x0 = (pos[0] - 0.5 * si->l[0]) / l;

    /* constants for the magnetic perturbation from seiji */
    w3 = exp(-(x0*x0 + y0*y0) / (w2*w2));
    w5 = 2.0*w1/w2;

    B[0] = tanh(y0) + (-w5*y0*w3);
    B[1] = w5*x0*w3;
    B[2] = HarrisParams.bg;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    harrisElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void harrisElectric(struct sti *si, struct stx *sx,
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
    harrisCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void harrisCurrent(struct sti *si, struct stx *sx,
                     double pos[3], double J[3])
{
    double y0;
    double l;

    (void)sx;

    l = HarrisParams.l;

    y0 = (pos[1] - 0.50*si->l[1]);

    J[0] = 0.;
    J[1] = 0.;
    J[2] = - 1./(l * cosh(y0) * cosh(y0));
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    harrisTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void harrisTemperature(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double T[2])
{
    double pmag;
    double tetp, teta;

    (void)si;
    (void)sx;
    (void)pos;

    tetp = HarrisParams.tetp;
    teta = HarrisParams.teta;

    pmag = 0.5;

    if (ispe == 0)
    {
        T[TPARA] = pmag * tetp / ((1 + tetp));
        T[TPERP] = pmag * tetp / ((1 + tetp));
    }
    else if (ispe == 1)
    {
        T[TPARA] = pmag / ((1 + tetp));
        T[TPERP] = pmag / ((1 + tetp));
    }
    else if (ispe == 2)
    {
        T[TPARA] = pmag * tetp / ((1 + tetp) * teta);
        T[TPERP] = pmag * tetp / ((1 + tetp) * teta);
    }
    else
    {
        printf("what the fuck is this ispe value ?\n");
    }
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
    harrisCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void harrisCurDrift(struct sti *si, struct stx *sx, double pos[3],
                    int ispe, double curdrift[3])
{
    double n;
    int s;
    double T[2], J[3];
    double Ttot;

    Ttot = 0.;
    for (s=0; s < NS+1; s++)
    {
        harrisTemperature(si, sx, pos, s, T);
        Ttot += T[TPARA] + 2 * T[TPERP];
    }

    harrisCurrent(si, sx, pos, J);
    n = harrisDensity(si, sx, pos, ispe);

    harrisTemperature(si, sx, pos, ispe, T);

    curdrift[0] = 1./Ttot * ( T[TPARA] + 2 * T[TPERP]) * J[0]/(n);
    curdrift[1] = 1./Ttot * ( T[TPARA] + 2 * T[TPERP]) * J[1]/(n);
    curdrift[2] = 1./Ttot * ( T[TPARA] + 2 * T[TPERP]) * J[2]/(n);

    if (ispe == 0) /* electrons charge = -1 */
    {
        curdrift[0] *= -1;
        curdrift[1] *= -1;
        curdrift[2] *= -1;
    }
    else if (ispe == 2) /* uncomment below for alphas drifting as protons */
    {
        //curdrift[0] = 0.0;
        //curdrift[1] = 0.0;
        //curdrift[2] = 0.0;
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
    harrisDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void harrisDrift(struct sti *si, struct stx *sx,
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







void harrisDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

}

















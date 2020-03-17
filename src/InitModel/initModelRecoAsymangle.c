
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


typedef struct s_asymangle
{
    double n1, n2;
    double l;
    double b1, b2;
    double te;
    double beta1;
    double alpha1,alpha2;
    double perturb;
} AsymAngle;




/* this is for private use only */
static AsymAngle AsymAngleParams;





/*---------------------------------------------------------------------------
    AsymAngleStart()
  ---------------------------------------------------------------------------
    AIM : read the parameters for the AsymAngle model
 ---------------------------------------------------------------------------*/
void asymangle_start(struct sti *si, struct stx *sx, char *dir)
{
    char sfh[80];
    char junk[100];
    FILE *fp;


    // unused but should be defined
    (void) si;

    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "recoasymangle.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL)
    {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.n1), &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.n2), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.b1), &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.b2), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.l), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.te), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.beta1), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.alpha1), &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.alpha2), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(AsymAngleParams.perturb),
                                            &(sx->irun));

    fclose(fp);

    /* __ display the topology parameters __ */
    if (sx->r == 0)
    {
       printf("___magnetic topology : Asymmetric reconnection orientation _\n");
       printf("density 1                : %8.4f\n", AsymAngleParams.n1);
       printf("density 2                : %8.4f\n", AsymAngleParams.n2);
       printf("magnetic field 1         : %8.4f\n", AsymAngleParams.b1);
       printf("magnetic field 2         : %8.4f\n", AsymAngleParams.b2);
       printf("magnetic shear angle     : %8.4f\n", AsymAngleParams.alpha1);
       printf("magnetic shear angle     : %8.4f\n", AsymAngleParams.alpha2);
       printf("electron temperature     : %8.4f\n", AsymAngleParams.te);
       printf("beta 1                   : %8.4f\n", AsymAngleParams.beta1);
       printf("sheet thickness          : %8.4f\n", AsymAngleParams.l);
       printf("mag. perturb.            : %8.4f\n", AsymAngleParams.perturb);
       printf("______________________________________________________\n");
       printf("\n");
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
    asymangleDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for all species at the given position (pos)
 ---------------------------------------------------------------------------*/
double asymangleDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe)
{
    double y0;
    double n1, n2, n, l;

    //unused but should be defined
    (void)sx;
    (void)ispe;

    n1 = AsymAngleParams.n1;
    n2 = AsymAngleParams.n2;
    l  = AsymAngleParams.l;

    /* center of the domain */
    y0 = (pos[1] - 0.5*  si->l[1]) / l;

    /* harris sheet density */
    n = n1 + (n2-n1) * 0.5 * (1 + tanh(y0));
    //n = 1;
    return n;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
    asymangleMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleMagnetic(struct sti *si, struct stx *sx,
                    double pos[3], double B[3])
{
    double x0, y0;
    const double w1 = 0.1, w2 = 2.0;
    double w3, w5;
    double b1, b2, bmod;
    double alpha1, alpha2;
    double alphay;
    double l;

    (void)sx;

    l      = AsymAngleParams.l;
    b1     = AsymAngleParams.b1;
    b2     = AsymAngleParams.b2;
    alpha1 = AsymAngleParams.alpha1  * acos(-1)/180.;
    alpha2 = AsymAngleParams.alpha2  * acos(-1)/180.;


    /* center of the domain */
    y0 = (pos[1] - 0.5*  si->l[1]) / l;

    /* distance from the middle of the box in x direction */
    x0 = (pos[0] - 0.5 * si->l[0]) / l;


    // variation of the magnetic modulus
    bmod = b1 + 0.5*(b2-b1)*(1. + tanh(y0));

    // variation of the magnetic orientation (from alpha1 to alpha2)
    alphay = alpha1 + (alpha2-alpha1) * 0.5 * (1. + tanh(y0));


    /* constants for the magnetic perturbation from seiji */
    w3 = exp(-(x0*x0 + y0*y0) / (w2*w2));
    w5 = 2.0*w1/w2;


    B[0] = bmod * cos(alphay) + (-w5*y0*w3);
    B[1] =  ((w5 * x0 * w3));
    B[2] = bmod * sin(alphay);
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    asymangleElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleElectric(struct sti *si, struct stx *sx,
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
    asymangleCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleCurrent(struct sti *si, struct stx *sx,
                      double pos[3], double J[3])
{
    double y0;
    double l;
    double b1, b2, bmod, bmodp;
    double alpha1, alpha2;
    double alphay, alphayp;

    (void)sx;


    l      = AsymAngleParams.l;
    b1     = AsymAngleParams.b1;
    b2     = AsymAngleParams.b2;
    alpha1 = AsymAngleParams.alpha1 * acos(-1)/180.;
    alpha2 = AsymAngleParams.alpha2 * acos(-1)/180.;

    y0 = (pos[1] - 0.50*si->l[1]);

    // variation of the magnetic orientation and its y derivative
    alphay  = alpha1 + (alpha2 - alpha1) * 0.5 * (1. + tanh(y0/l));
    alphayp = (alpha2 - alpha1) /(2.*l)/(cosh(y0/l)*cosh(y0/l));

    // variation of the magnetic modulus and its y derivative
    bmod  = b1 + 0.5*(b2-b1)*(1. + tanh(y0/l));
    bmodp = (b2-b1)/(2.*l) / (cosh(y0/l)*cosh(y0/l));

    J[0] = bmodp * sin(alphay) +  bmod * cos(alphay) * alphayp;

    J[1] = 0.;

    J[2] = - bmodp* cos(alphay) + bmod * sin(alphay) * alphayp;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    asymangleTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void asymangleTemperature(struct sti *si, struct stx *sx,
                          double pos[3], int ispe, double T[2])
{
    double te;
    double bsq;
    double b1, beta1, n;
    double Ptot;

    (void)si;
    (void)sx;
    (void)pos;

    double b[3];

    // get the :
    // constant electron temperature 'te'
    // magnetic field amplitude
    // total density

    te    = AsymAngleParams.te;
    beta1 = AsymAngleParams.beta1;
    b1    = AsymAngleParams.b1;
    n     = asymangleDensity(si, sx, pos, ispe);
    asymangleMagnetic(si, sx, pos, b);
    bsq = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];


    // then calculate the total pressure, knowing beta on side 1
    Ptot = (1+beta1) * b1*b1/2.;

    // now the total temperature is T = 1/n(y) * (Ptot - B2/2)  - te;

    if (ispe == 0) // electrons
    {
        T[TPARA] = te;
        T[TPERP] = te;
    }

    else if (ispe == 1) // protons
    {
        T[TPARA] = 1./n * (Ptot - bsq/2.) - te;
        T[TPERP] = 1./n * (Ptot - bsq/2.) - te;

        if(T[TPARA] <= 0)
        {
            if (sx->r == 0)
            {
                printf("ERROR - Ion temperature negative\n");
            }
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
    }
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
    asymangleCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void asymangleCurDrift(struct sti *si, struct stx *sx, double pos[3],
                       int ispe, double curdrift[3])
{
    double n;
    int s;
    double T[2], J[3];
    double Ttot;

    Ttot = 0.;
    for (s=0; s < NS+1; s++)
    {
        asymangleTemperature(si, sx, pos, s, T);
        Ttot += T[TPARA] + 2 * T[TPERP];
    }

    asymangleCurrent(si, sx, pos, J);
    n = asymangleDensity(si, sx, pos, ispe);

    asymangleTemperature(si, sx, pos, ispe, T);

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
    asymangleDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void asymangleDrift(struct sti *si, struct stx *sx,
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






void asymangleDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

}















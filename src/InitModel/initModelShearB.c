
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


typedef struct s_shearB
{
    double l;
    double te;
    double ti;
    double perturb;
    double bg;
} ShearB;




/* __ this is for private use only __ */
static ShearB ShearBParams;




void shearB_start(struct sti *si, struct stx *sx, char *dir)
{
    char sfh[80];
    char junk[100];
    FILE *fp;


    /* __ unused but should be defined __ */
    (void) si;

    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "shearB.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL)
    {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(ShearBParams.l), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(ShearBParams.te), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(ShearBParams.ti), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(ShearBParams.perturb), &(sx->irun));

    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    fscandbl(__FILE__, __LINE__, fp, sx->r, &(ShearBParams.bg), &(sx->irun));

    fclose(fp);

    /* __ display the topology parameters __ */
    if (sx->r == 0)
    {
       printf("________________ magnetic topology : ShearB ________\n");
       printf("sheet thickness          : %8.4f\n", ShearBParams.l);
       printf("electron temeprature     : %8.4f\n", ShearBParams.te);
       printf("proton temeprature       : %8.4f\n", ShearBParams.ti);
       printf("mag. perturb.            : %8.4f\n", ShearBParams.perturb);
       printf("guide field              : %8.4f\n", ShearBParams.bg);
       printf("______________________________________________________\n");
       printf("\n");
    }
}



double shearBDensity(struct sti *si, struct stx *sx, double pos[3], int ispe)
{
    /* __ unused but should be defined __ */
    (void)si;
    (void)sx;
    (void)pos;
    (void)ispe;

    return 1.0;
}



void shearBMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3])
{
    double x0, y0;
    const double w1 = ShearBParams.perturb, w2 = 2.0;
    double w3, w5;
    double l;


    /* __ unused but should be defined __ */
    (void)sx;

    l = ShearBParams.l;

    /* __ center of the domain __ */
    y0 = (pos[1]-0.5*si->l[1])/l;

    /* __ distance from the middle of the box in x direction __ */
    x0 = (pos[0]-0.5*si->l[0])/l;

    /* __ constants for the magnetic perturbation from seiji __ */
    w3 = exp(-(x0*x0+y0*y0)/(w2*w2));
    w5 = 2.0*w1/w2;

    B[0] = tanh(y0)-w5*y0*w3;
    B[1] = w5*x0*w3;
    B[2] = ShearBParams.bg;
}



void shearBElectric(struct sti *si, struct stx *sx, double pos[3], double E[3])
{
    /* __ unused but should be defined __ */
    (void)si;
    (void)sx;
    (void)pos;


    E[0] = 0.0;
    E[1] = 0.0;
    E[2] = 0.0;
}



void shearBCurrent(struct sti *si, struct stx *sx, double pos[3], double J[3])
{
    double y0;
    double l;


    /* __ unused but should be defined __ */
    (void)sx;

    l = ShearBParams.l;

    y0 = pos[1]-0.5*si->l[1];

    J[0] = 0.;
    J[1] = 0.;
    J[2] = -1.0/(l*cosh(y0)*cosh(y0));
}



void shearBTemperature(struct sti *si, struct stx *sx, double pos[3], int ispe, double T[2])
{
    double pTot;
    double bg, te, ti;


    /* __ unused but should be defined __ */
    (void)si;
    (void)sx;
    (void)pos;

    bg = ShearBParams.bg;
    te = ShearBParams.te;
    ti = ShearBParams.ti;

    /* __ asymptotic total pressure __ */
    pTot = 0.5*(1.0+bg*bg)+te+ti;

    /* __ for the electrons __ */
    if (ispe == 0)
    {
        T[TPARA] = te;
        T[TPERP] = te;
    }

    /* __ for the protons __ */
    else
    {
        T[TPARA] = pTot-te-0.5*(1.0+bg*bg);
        T[TPERP] = pTot-te-0.5*(1.0+bg*bg);
    }
}



void shearBCurDrift(struct sti *si, struct stx *sx, double pos[3], int ispe, double curdrift[3])
{
    int s;
    double T[2], J[3];
    double Ttot;


    Ttot = 0.;
    for (s = 0; s < NS+1; s++)
    {
        shearBTemperature(si, sx, pos, s, T);
        Ttot += T[TPARA]+2*T[TPERP];
    }

    shearBCurrent(si, sx, pos, J);
    shearBTemperature(si, sx, pos, ispe, T);

    curdrift[0] = J[0]*(T[TPARA]+2*T[TPERP])/Ttot;
    curdrift[1] = J[1]*(T[TPARA]+2*T[TPERP])/Ttot;
    curdrift[2] = J[2]*(T[TPARA]+2*T[TPERP])/Ttot;

    /* __ electrons charge = -1 __ */
    if (ispe == 0)
    {
        curdrift[0] *= -1;
        curdrift[1] *= -1;
        curdrift[2] *= -1;
    }
}



void shearBDrift(struct sti *si, struct stx *sx, double pos[3], int ispe, double vdrift[3])
{
    (void)si;
    (void)sx;
    (void)pos;
    (void)ispe;


    vdrift[0] = 0.0;
    vdrift[1] = 0.0;
    vdrift[2] = 0.0;
}

void shearBDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

}


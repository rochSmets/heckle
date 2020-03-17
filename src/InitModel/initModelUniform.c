

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
//#include "plasma.h"


#define PARA 0
#define PERP 1


typedef struct s_uniform
{
    double b[3];
    double n[NS+1];
    double V[NS+1][3];
    double T[NS+1][2];
} Uniform;




/* this is for private use only */
Uniform UniformParams;





/*---------------------------------------------------------------------------
    uniformStart()
  ---------------------------------------------------------------------------
    AIM : read the parameters for Uniform model
 ---------------------------------------------------------------------------*/
void uniform_start(struct sti *si, struct stx *sx, char *dir)
{
    int g, h;
    char sfh[80];
    char junk[100];
    FILE *fp;


    (void)si;

    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "uniform.txt");

    /* __ read the topology parameters __ */
    fp = fopen(sfh, "r");
    if (fp == NULL)
    {
        printf("\n\nproblem in opening file %s\n", sfh);
        exit(-1);
    }

    /* __ dc magnetic field components __ */
    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (h = 0; h < 3; h++)
    {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(UniformParams.b[h]), &(sx->irun));
    }

    /* __ specie densities __ */
    fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
    for (UniformParams.n[0] = 0.0, g = 1; g < NS+1; g++)
    {
        fscandbl(__FILE__, __LINE__, fp, sx->r, &(UniformParams.n[g]), &(sx->irun));
        UniformParams.n[0] += UniformParams.n[g]*si->qs[g];
    }

    /* __ specie fluid velocities __ */
    for (h = 0; h < 3; h++)
    {
        fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
        for (g = 1; g < NS+1; g++)
        {
            fscandbl(__FILE__, __LINE__, fp, sx->r, &(UniformParams.V[g][h]), &(sx->irun));
        }
    }

    /* __ specie temperatures __ */
    for (h = 0; h < 2; h++)
    {
        fscanstr(__FILE__, __LINE__, fp, sx->r, junk, &(sx->irun));
        for (g = 0; g < NS+1; g++)
        {
            fscandbl(__FILE__, __LINE__, fp, sx->r, &(UniformParams.T[g][h]), &(sx->irun));
        }
    }


    fclose(fp);

    /* __ node 0 displays the topology parameters __ */
    if (sx->r == 0)
    {
       printf("________________ magnetic topology : uniform _________\n");
       for (h = 0; h < 3; h++) {
           printf("b field [%1d]      :%12.4f\n", h, UniformParams.b[h]);
       }
       printf("\n");

       for (g = 0; g < NS+1; g++) {
           printf("specie %1d\n", g);
           printf("... density      :%12.4f\n", UniformParams.n[g]);
           for (h = 0; h < 3; h++)
           {
               printf("... fluid V[%1d]   :%12.4f\n", h, UniformParams.V[g][h]);
           }
           printf("... parallel T   :%12.4f\n", UniformParams.T[g][0]);
           printf("... perp. T      :%12.4f\n", UniformParams.T[g][1]);
           printf("\n");
       }

//     printf("______________________________________________________\n");
//     printf("\n");
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
    uniformDensity()
  ---------------------------------------------------------------------------
    AIM : returns the density for all species at the given position (pos)
 ---------------------------------------------------------------------------*/
double uniformDensity(struct sti *si, struct stx *sx,
                     double pos[3], int ispe)
{
    (void)si;
    (void)sx;
    (void)pos;

    return UniformParams.n[ispe];
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    uniformMagnetic()
  ---------------------------------------------------------------------------
    AIM : returns the magnetic field at the given position (pos)
 ---------------------------------------------------------------------------*/
void uniformMagnetic(struct sti *si, struct stx *sx,
                      double pos[3], double B[3])
{
    (void)si;
    (void)sx;
    (void)pos;

    B[0] = UniformParams.b[0];
    B[1] = UniformParams.b[1];
    B[2] = UniformParams.b[2];
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
    uniformElectric()
  ---------------------------------------------------------------------------
    AIM : returns the electric field at the given position (pos)
 ---------------------------------------------------------------------------*/
void uniformElectric(struct sti *si, struct stx *sx,
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
    uniformCurrent()
  ---------------------------------------------------------------------------
    AIM : returns the current density at the given position (pos)
 ---------------------------------------------------------------------------*/
void uniformCurrent(struct sti *si, struct stx *sx,
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
    uniformTemperature()
  ---------------------------------------------------------------------------
    AIM : returns the temperature
 ---------------------------------------------------------------------------*/
void uniformTemperature(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double T[2])
{
    (void)si;
    (void)sx;
    (void)pos;

    T[PARA] = UniformParams.T[ispe][PARA];
    T[PERP] = UniformParams.T[ispe][PERP];

}
/*===========================================================================*/







/*---------------------------------------------------------------------------
    uniformCurDrift()
  ---------------------------------------------------------------------------
    AIM : returns the component of the drift vel. associated with the current
 ---------------------------------------------------------------------------*/
void uniformCurDrift(struct sti *si, struct stx *sx, double pos[3],
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
    uniformDrift()
  ---------------------------------------------------------------------------
    AIM : returns a drift velocity at a given position (pos)
 ---------------------------------------------------------------------------*/
void uniformDrift(struct sti *si, struct stx *sx,
                 double pos[3], int ispe, double vdrift[3])
{
    (void)si;
    (void)sx;
    (void)pos;
    (void)ispe;

    vdrift[0] = UniformParams.V[ispe][0];
    vdrift[1] = UniformParams.V[ispe][1];
    vdrift[2] = UniformParams.V[ispe][2];
}
/*===========================================================================*/







void uniformDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

}


/*
 * This is the module that plugs the model chosen by the user among all
 * available models.
 *
 * This is done by assigning model specific functions to generic model
 * function pointers
 *
 */



#include <initModel.h>
#include <structures.h>
#include <initModelHarris.h>
#include <initModelDoubleHarris.h>
#include <initModelUniform.h>
#include <initModelRecoAsymangle.h>
#include <initModelNBeams.h>
#include <initModelShearB.h>
#include <initModelNLasers.h>
#include <stdio.h>
#include <stdlib.h>
#include <Particle/particle.h>




/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                               PRIVATE                                 //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */

typedef struct s_model
{
    int model;

    void (*start)(struct sti *si, struct stx *sx, char *dir);

    void (*magnetic)(struct sti *si, struct stx *sx,
                       double pos[3], double B[3]);

    void (*electric)(struct sti *si, struct stx *sx,
                       double pos[3], double E[3]);

    void (*current)(struct sti *si, struct stx *sx,
                      double pos[3], double J[3]);

    void (*drift)(struct sti *si, struct stx *sx,
                    double pos[3], int ispe, double vdrift[3]);

    void (*curdrift)(struct sti *si, struct stx *sx,
                       double pos[3], int ispe, double curdrif[3]);

    double (*density)(struct sti *si, struct stx *sx, double pos[3], int speciesID);

    void (*temperature)(struct sti *si, struct stx *st,
                        double pos[3], int ispe, double T[2]);

    void (*drive)(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it);

} InitModel;



/* this structure is for private use only. */
InitModel MyModel;






/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                               PUBLIC                                  //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */






void initModelStart(struct sti *si, struct stx *sx, char *dir)
{

    switch(si->InitModelID)
    {

        /* ---------------------------------------------------------------- */
        /*                 FULLY PERIODIC HARRIS SHEET                      */
        /* ---------------------------------------------------------------- */
        case INITMODEL_DOUBLEHARRIS:

            MyModel.model       = INITMODEL_DOUBLEHARRIS;
            MyModel.start       = doubleharris_start;
            MyModel.magnetic    = doubleharrisMagnetic;
            MyModel.electric    = doubleharrisElectric;
            MyModel.current     = doubleharrisCurrent;
            MyModel.drift       = doubleharrisDrift;
            MyModel.curdrift    = doubleharrisCurDrift;
            MyModel.density     = doubleharrisDensity;
            MyModel.temperature = doubleharrisTemperature;
            MyModel.drive       = doubleharrisDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*                       HARRIS SHEET                               */
        /* ---------------------------------------------------------------- */
        case INITMODEL_HARRIS:

            MyModel.model       = INITMODEL_HARRIS;
            MyModel.start       = harris_start;
            MyModel.magnetic    = harrisMagnetic;
            MyModel.electric    = harrisElectric;
            MyModel.current     = harrisCurrent;
            MyModel.drift       = harrisDrift;
            MyModel.curdrift    = harrisCurDrift;
            MyModel.density     = harrisDensity;
            MyModel.temperature = harrisTemperature;
            MyModel.drive       = harrisDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*                       UNIFORM DENSITY                            */
        /* ---------------------------------------------------------------- */
        case INITMODEL_UNIFORM:
            MyModel.model       = INITMODEL_UNIFORM;
            MyModel.start       = uniform_start;
            MyModel.magnetic    = uniformMagnetic;
            MyModel.electric    = uniformElectric;
            MyModel.current     = uniformCurrent;
            MyModel.drift       = uniformDrift;
            MyModel.curdrift    = uniformCurDrift;
            MyModel.density     = uniformDensity;
            MyModel.temperature = uniformTemperature;
            MyModel.drive       = uniformDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               ASYMMETRIC RECONNECTION ANGLE                      */
        /* ---------------------------------------------------------------- */
        case INITMODEL_ASYMANGLE:
            MyModel.model       = INITMODEL_ASYMANGLE;
            MyModel.start       = asymangle_start;
            MyModel.magnetic    = asymangleMagnetic;
            MyModel.electric    = asymangleElectric;
            MyModel.current     = asymangleCurrent;
            MyModel.drift       = asymangleDrift;
            MyModel.curdrift    = asymangleCurDrift;
            MyModel.density     = asymangleDensity;
            MyModel.temperature = asymangleTemperature;
            MyModel.drive       = asymangleDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               N BEAMS                                            */
        /* ---------------------------------------------------------------- */
        case INITMODEL_NBEAMS:
            MyModel.model       = INITMODEL_NBEAMS;
            MyModel.start       = nbeams_start;
            MyModel.magnetic    = nbeamsMagnetic;
            MyModel.electric    = nbeamsElectric;
            MyModel.current     = nbeamsCurrent;
            MyModel.drift       = nbeamsDrift;
            MyModel.curdrift    = nbeamsCurDrift;
            MyModel.density     = nbeamsDensity;
            MyModel.temperature = nbeamsTemperature;
            MyModel.drive       = nbeamsDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               SHEAR B                                            */
        /* ---------------------------------------------------------------- */
        case INITMODEL_SHEARB:
            MyModel.model       = INITMODEL_SHEARB;
            MyModel.start       = shearB_start;
            MyModel.magnetic    = shearBMagnetic;
            MyModel.electric    = shearBElectric;
            MyModel.current     = shearBCurrent;
            MyModel.drift       = shearBDrift;
            MyModel.curdrift    = shearBCurDrift;
            MyModel.density     = shearBDensity;
            MyModel.temperature = shearBTemperature;
            MyModel.drive       = shearBDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               LASERS                                             */
        /* ---------------------------------------------------------------- */
        case INITMODEL_NLASERS:
            MyModel.model       = INITMODEL_NLASERS;
            MyModel.start       = nlasers_start;
            MyModel.magnetic    = nlasersMagnetic;
            MyModel.electric    = nlasersElectric;
            MyModel.current     = nlasersCurrent;
            MyModel.drift       = nlasersDrift;
            MyModel.curdrift    = nlasersCurDrift;
            MyModel.density     = nlasersDensity;
            MyModel.temperature = nlasersTemperature;
            MyModel.drive       = nlasersDrive;
        break;

#if 0
        /* ---------------------------------------------------------------- */
        /*               EXAMPLE OF NON_DEFINED MODEL                       */
        /* ---------------------------------------------------------------- */
        case INITMODEL_TODEFINE:
            MyModel.model       = INITMODEL_TODEFINE;
            MyModel.start       = todefine_start;
            MyModel.magnetic    = todefineMagnetic;
            MyModel.electric    = todefineElectric;
            MyModel.current     = todefineCurrent;
            MyModel.drift       = todefineDrift;
            MyModel.curdrift    = todefineCurDrift;
            MyModel.density     = todefineDensity;
            MyModel.temperature = todefineTemperature;
        break;

#endif

        default:
        break;
    }

    MyModel.start(si, sx, dir);
}
/*===========================================================================*/






void initModelMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3])
{
    MyModel.magnetic(si, sx, pos, B);
}

/*===========================================================================*/





void initModelElectric(struct sti *si, struct stx *sx, double pos[3], double E[3])
{
    MyModel.electric(si, sx, pos, E);
}

/*===========================================================================*/






void initModelCurrent(struct sti *si, struct stx *sx, double pos[3], double J[3])
{
    MyModel.current(si, sx, pos, J);
}


/*===========================================================================*/




double initModelDensity(struct sti *si, struct stx *sx, double pos[3], int speciesID)
{
    return MyModel.density(si, sx, pos, speciesID);
}



/*===========================================================================*/


void initModelTemperature(struct sti *si, struct stx *sx,
                      double pos[3], int ispe, double T[2])
{
    MyModel.temperature(si, sx, pos, ispe, T);
}

/*===========================================================================*/





void initModelCurDrift(struct sti *si, struct stx *sx,
                   double pos[3], int ispe, double curdrift[3])
{
    MyModel.curdrift(si, sx, pos, ispe, curdrift);
}



/*===========================================================================*/


void initModelDrift(struct sti *si,
                struct stx *sx,
                double pos[3],
                int ispe,
                double vdrift[3])
{
    MyModel.drift(si, sx, pos, ispe, vdrift);
}


void initModelDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

    MyModel.drive(si, sx, s0, s1, s2, sp, ipc, it);

}


void normal(struct sti *si, struct stx *sx)
{
    double pos[3], ns[NS+1];
    double dlxyz;
    double tnrp;        // total number of real particles
    int i, j, k, s;


    /* __ volume of a cell __ */
    dlxyz = si->dl[0]*si->dl[1]*si->dl[2];

    // fake value... just to do not writes bizarres values in species.h5
    si->ws[0] = 0.0;

    /* __ loop on the part. __ */
    for (s = 1; s < NS+1; s++)
    {
        /* set the total number of real particles to 0 first*/
        tnrp = 0.0;

        /* then integrate it over the whole domain */
        for (i = 0; i < si->n[0]; i++)
        {
            for (j = 0; j < si->n[1]; j++)
            {
               for (k = 0; k < si->n[2]; k++)
               {
                    /* __ center of the cell __ */
                    pos[0] = (i +0.5) * si->dl[0];
                    pos[1] = (j +0.5) * si->dl[1];
                    pos[2] = (k +0.5) * si->dl[2];

                    /* __ local density for part. s __ */
                    ns[s] = initModelDensity(si, sx, pos, s);

                    /* __ integrated local "mass" of part. s __ */
                    tnrp += ns[s] * dlxyz;
                }
            }
        }

        /* __ weight of part. "s" __ */
        si->ws[s] = (si->ns[s] == 0) ? EPS6 : tnrp / (si->ns[s] * dlxyz);
    }
    //printf("si->ws = %f \t %f \t %f\n", si->ws[0], si->ws[1], si->ws[2]);

}

// .............................................................................
//
// the functions polynom & rotateCoordinates are needed for both NBeams & Lasers
//
// .............................................................................








#if 0
double model_j2vfac(Species *spe, double rho)
{
    return model.model_j2vfac(spe, rho);
}


/*===========================================================================*/


int model_species_kind(int ispe)
{
    return model.model_species_kind(ispe);
}



/*===========================================================================*/



int model_species_number(void)
{
    return model.model_species_number();
}


/*===========================================================================*/


int model_species_npart(int ispe)
{
    return model.model_species_npart(ispe);
}
#endif



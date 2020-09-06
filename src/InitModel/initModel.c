
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
#include <initModelNhotSpots.h>
#include <stdio.h>
#include <stdlib.h>
#include <Particle/particle.h>




/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                               PRIVATE                                 //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */

typedef struct init_model
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
InitModel MyInitModel;






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

            MyInitModel.model       = INITMODEL_DOUBLEHARRIS;
            MyInitModel.start       = doubleharris_start;
            MyInitModel.magnetic    = doubleharrisMagnetic;
            MyInitModel.electric    = doubleharrisElectric;
            MyInitModel.current     = doubleharrisCurrent;
            MyInitModel.drift       = doubleharrisDrift;
            MyInitModel.curdrift    = doubleharrisCurDrift;
            MyInitModel.density     = doubleharrisDensity;
            MyInitModel.temperature = doubleharrisTemperature;
            MyInitModel.drive       = doubleharrisDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*                       HARRIS SHEET                               */
        /* ---------------------------------------------------------------- */
        case INITMODEL_HARRIS:

            MyInitModel.model       = INITMODEL_HARRIS;
            MyInitModel.start       = harris_start;
            MyInitModel.magnetic    = harrisMagnetic;
            MyInitModel.electric    = harrisElectric;
            MyInitModel.current     = harrisCurrent;
            MyInitModel.drift       = harrisDrift;
            MyInitModel.curdrift    = harrisCurDrift;
            MyInitModel.density     = harrisDensity;
            MyInitModel.temperature = harrisTemperature;
            MyInitModel.drive       = harrisDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*                       UNIFORM DENSITY                            */
        /* ---------------------------------------------------------------- */
        case INITMODEL_UNIFORM:
            MyInitModel.model       = INITMODEL_UNIFORM;
            MyInitModel.start       = uniform_start;
            MyInitModel.magnetic    = uniformMagnetic;
            MyInitModel.electric    = uniformElectric;
            MyInitModel.current     = uniformCurrent;
            MyInitModel.drift       = uniformDrift;
            MyInitModel.curdrift    = uniformCurDrift;
            MyInitModel.density     = uniformDensity;
            MyInitModel.temperature = uniformTemperature;
            MyInitModel.drive       = uniformDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               ASYMMETRIC RECONNECTION ANGLE                      */
        /* ---------------------------------------------------------------- */
        case INITMODEL_ASYMANGLE:
            MyInitModel.model       = INITMODEL_ASYMANGLE;
            MyInitModel.start       = asymangle_start;
            MyInitModel.magnetic    = asymangleMagnetic;
            MyInitModel.electric    = asymangleElectric;
            MyInitModel.current     = asymangleCurrent;
            MyInitModel.drift       = asymangleDrift;
            MyInitModel.curdrift    = asymangleCurDrift;
            MyInitModel.density     = asymangleDensity;
            MyInitModel.temperature = asymangleTemperature;
            MyInitModel.drive       = asymangleDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               N BEAMS                                            */
        /* ---------------------------------------------------------------- */
        case INITMODEL_NBEAMS:
            MyInitModel.model       = INITMODEL_NBEAMS;
            MyInitModel.start       = nbeams_start;
            MyInitModel.magnetic    = nbeamsMagnetic;
            MyInitModel.electric    = nbeamsElectric;
            MyInitModel.current     = nbeamsCurrent;
            MyInitModel.drift       = nbeamsDrift;
            MyInitModel.curdrift    = nbeamsCurDrift;
            MyInitModel.density     = nbeamsDensity;
            MyInitModel.temperature = nbeamsTemperature;
            MyInitModel.drive       = nbeamsDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               SHEAR B                                            */
        /* ---------------------------------------------------------------- */
        case INITMODEL_SHEARB:
            MyInitModel.model       = INITMODEL_SHEARB;
            MyInitModel.start       = shearB_start;
            MyInitModel.magnetic    = shearBMagnetic;
            MyInitModel.electric    = shearBElectric;
            MyInitModel.current     = shearBCurrent;
            MyInitModel.drift       = shearBDrift;
            MyInitModel.curdrift    = shearBCurDrift;
            MyInitModel.density     = shearBDensity;
            MyInitModel.temperature = shearBTemperature;
            MyInitModel.drive       = shearBDrive;
        break;



        /* ---------------------------------------------------------------- */
        /*               N LASERS                                           */
        /* ---------------------------------------------------------------- */
        case INITMODEL_NLASERS:
            MyInitModel.model       = INITMODEL_NLASERS;
            MyInitModel.start       = nlasers_start;
            MyInitModel.magnetic    = nlasersMagnetic;
            MyInitModel.electric    = nlasersElectric;
            MyInitModel.current     = nlasersCurrent;
            MyInitModel.drift       = nlasersDrift;
            MyInitModel.curdrift    = nlasersCurDrift;
            MyInitModel.density     = nlasersDensity;
            MyInitModel.temperature = nlasersTemperature;
            MyInitModel.drive       = nlasersDrive;
        break;



       /* ---------------------------------------------------------------- */
       /*               N HOTSPOTS                                         */
       /* ---------------------------------------------------------------- */
       case INITMODEL_NHOTSPOTS:
           MyInitModel.model       = INITMODEL_NHOTSPOTS;
           MyInitModel.start       = nhotspots_start;
           MyInitModel.magnetic    = nhotspotsMagnetic;
           MyInitModel.electric    = nhotspotsElectric;
           MyInitModel.current     = nhotspotsCurrent;
           MyInitModel.drift       = nhotspotsDrift;
           MyInitModel.curdrift    = nhotspotsCurDrift;
           MyInitModel.density     = nhotspotsDensity;
           MyInitModel.temperature = nhotspotsTemperature;
           MyInitModel.drive       = nhotspotsDrive;
       break;



#if 0
        /* ---------------------------------------------------------------- */
        /*               EXAMPLE OF NON_DEFINED MODEL                       */
        /* ---------------------------------------------------------------- */
        case INITMODEL_TODEFINE:
            MyInitModel.model       = INITMODEL_TODEFINE;
            MyInitModel.start       = todefine_start;
            MyInitModel.magnetic    = todefineMagnetic;
            MyInitModel.electric    = todefineElectric;
            MyInitModel.current     = todefineCurrent;
            MyInitModel.drift       = todefineDrift;
            MyInitModel.curdrift    = todefineCurDrift;
            MyInitModel.density     = todefineDensity;
            MyInitModel.temperature = todefineTemperature;
        break;

#endif

        default:
        break;
    }

    MyInitModel.start(si, sx, dir);
}
/*===========================================================================*/






void initModelMagnetic(struct sti *si, struct stx *sx, double pos[3], double B[3])
{
    MyInitModel.magnetic(si, sx, pos, B);
}

/*===========================================================================*/





void initModelElectric(struct sti *si, struct stx *sx, double pos[3], double E[3])
{
    MyInitModel.electric(si, sx, pos, E);
}

/*===========================================================================*/






void initModelCurrent(struct sti *si, struct stx *sx, double pos[3], double J[3])
{
    MyInitModel.current(si, sx, pos, J);
}


/*===========================================================================*/




double initModelDensity(struct sti *si, struct stx *sx, double pos[3], int speciesID)
{
    return MyInitModel.density(si, sx, pos, speciesID);
}



/*===========================================================================*/


void initModelTemperature(struct sti *si, struct stx *sx,
                      double pos[3], int ispe, double T[2])
{
    MyInitModel.temperature(si, sx, pos, ispe, T);
}

/*===========================================================================*/





void initModelCurDrift(struct sti *si, struct stx *sx,
                   double pos[3], int ispe, double curdrift[3])
{
    MyInitModel.curdrift(si, sx, pos, ispe, curdrift);
}



/*===========================================================================*/


void initModelDrift(struct sti *si,
                struct stx *sx,
                double pos[3],
                int ispe,
                double vdrift[3])
{
    MyInitModel.drift(si, sx, pos, ispe, vdrift);
}


void initModelDrive(struct sti si, struct stx sx, Grid0 *s0, ST1 *s1, ST2 *s2, Particle *sp[NS+1], int ipc, int it)
{

    MyInitModel.drive(si, sx, s0, s1, s2, sp, ipc, it);

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



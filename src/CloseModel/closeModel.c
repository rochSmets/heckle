
/*
 * This is the module that plugs the model chosen by the user among all
 * available models.
 *
 * This is done by assigning model specific functions to generic model
 * function pointers
 *
 */



#include "structures.h"
#include "closeModel.h"
#include "closeModelIsotherm.h"
#include "closeModelPolytrop.h"
#include "closeModelFullPressureImplicit.h"
#include "closeModelFullPressureSubCycling.h"
#include "stdio.h"
#include "stdlib.h"





/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                               PRIVATE                                 //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */

typedef struct s_model
{
    int model;

    void (*pressure)(int it, const STI * const si, const STX * const sx,
            struct st1 *s1, struct st2 *s2, HeckleBC *hbc, Ghosts *ghosts, int ipc);

} CloseModel;



/* this structure is for private use only. */
CloseModel MyModel;






/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                               PUBLIC                                  //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */






void closeModelStart(struct sti *si, struct stx *sx, char *dir)
{

    switch(si->CloseModelID)
    {

        /* ---------------------------------------------------------------- */
        /*                 ISOTHERM CLOSURE FOR ELECTRONS                   */
        /* ---------------------------------------------------------------- */
        case CLOSE_ISOTHERM:

            MyModel.model       = CLOSE_ISOTHERM;
            MyModel.pressure    = isothermPressure;
        break;

#ifdef FULLP
        /* ---------------------------------------------------------------- */
        /*                 FULLP CLOSURE FOR ELECTRONS (IMPLICIT)           */
        /* ---------------------------------------------------------------- */
        case CLOSE_FULLIMPL:

            MyModel.model       = CLOSE_FULLIMPL;
            MyModel.pressure    = implicitPressure;
        break;

        /* ---------------------------------------------------------------- */
        /*                 FULLP CLOSURE FOR ELECTRONS (SUBCYCLING)         */
        /* ---------------------------------------------------------------- */
        case CLOSE_FULLSUB:

            MyModel.model       = CLOSE_FULLSUB;
            MyModel.pressure    = subCycledPressure;
        break;
#endif
        case CLOSE_POLYTROP:
            MyModel.model       = CLOSE_POLYTROP;
            MyModel.pressure    = polytropPressure;
        break;
#if 0
        /* ---------------------------------------------------------------- */
        /*               EXAMPLE OF NON_DEFINED MODEL                       */
        /* ---------------------------------------------------------------- */
        case CLOSE_TODEFINE:
            MyModel.model       = CLOSE_TODEFINE;
            MyModel.pressure    = todefinePressure;
        break;

#endif

        default:
        break;
    }

}
/*===========================================================================*/






void closeModelPressure(int it,
                        const STI * const si,
                        const STX * const sx,
                        struct st1 *s1,
                        struct st2 *s2,
                        HeckleBC *hbc,
                        Ghosts *ghosts,
                        int ipc)
{
    MyModel.pressure(it, si, sx, s1, s2, hbc, ghosts, ipc);
}


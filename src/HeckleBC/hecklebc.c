
#include <stdlib.h>
#include <stdio.h>

#include <structures.h>
#include <particlebc.h>

#include <perfcond.h>
#include <periodic.h>


#include <bc_constants.h>

//#include <bc_open.h>


/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */








typedef struct s_hecklebc
{

    // boundary condition ID
    int bc[3];

    /* Particle boundary condition handle */
    PartBC   *partbc;


    /* those pointers will point to the appropriate kind of
       boundary condition. There is one pointer for each
       type of quantity.

       The functions for a specifc type of BC are defined
       in its associated module, ex. periodic.c or perfcond.c
    */


    // TODO rename osvs, nv etc. in more 'english names'
    void (*fix_osvs[3])  (const STX* const sx, ST2 *s2);
    void (*fix_nv[3])    (const STX* const sx, ST2 *s2);
    void (*fix_e[3])     (const STX* const sx, ST2 *s2);
    void (*fix_f[3])     (const STX* const sx, ST2 *s2);
    void (*fix_j[3])     (const STX* const sx, ST2 *s2);
    void (*fix_P[3])     (const STX* const sx, ST2 *s2);
#ifdef FULLP
    void (*fix_P_winske[3])     (const STX* const sx, ST2 *s2);
    void (*fix_Driver_winske[3])     (const STX* const sx, ST2 *s2);
#endif

} HeckleBC;








/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */






/*---------------------------------------------------------------------------
  HeckleBCInit()
  ---------------------------------------------------------------------------
  AIM : Initialize the Heckle Boundary Conditions module
 ---------------------------------------------------------------------------*/
HeckleBC* HeckleBCInit(const STI* const si)
{
    HeckleBC* self;
    int c;

    self = malloc(sizeof *self);


    for (c=0; c <3; c++)
    {
        // internal copy of the boundary conditions
        self->bc[c] = si->bc[c];
    }


    /* ------------------------------------------------------------- */
    /*                         X BOUNDARIES                          */
    /* ------------------------------------------------------------- */


    switch(self->bc[0])
    {
        case BC_TYPE_PERIODIC:
            self->fix_osvs[0]  = BC_periodic_os_vs_x;
            self->fix_nv[0]    = BC_periodic_nv_x;
            self->fix_e[0]     = BC_periodic_e_x;
            self->fix_f[0]     = BC_periodic_f_x;
            self->fix_j[0]     = BC_periodic_j_x;
            self->fix_P[0]     = BC_periodic_P_x;
#ifdef FULLP
            self->fix_P_winske[0]     = BC_periodic_P_winske_x;
            self->fix_Driver_winske[0]     = BC_periodic_Driver_winske_x;
#endif
        break;


        case BC_TYPE_PERFCOND:
            self->fix_osvs[0]  = BC_perfcond_os_vs_x;
            self->fix_nv[0]    = BC_perfcond_nv_x;
            self->fix_e[0]     = BC_perfcond_e_x;
            self->fix_f[0]     = BC_perfcond_f_x;
            self->fix_j[0]     = BC_perfcond_j_x;
            self->fix_P[0]     = BC_perfcond_P_x;
#ifdef FULLP
            self->fix_P_winske[0]     = BC_perfcond_P_winske_x;
            self->fix_Driver_winske[0]     = BC_perfcond_Driver_winske_x;
#endif
        break;
    } // end switch bc 0


    /* ------------------------------------------------------------- */
    /*                         Y BOUNDARIES                          */
    /* ------------------------------------------------------------- */

    switch(self->bc[1])
    {
        case BC_TYPE_PERIODIC:
            self->fix_osvs[1]  = BC_periodic_os_vs_y;
            self->fix_nv[1]    = BC_periodic_nv_y;
            self->fix_e[1]     = BC_periodic_e_y;
            self->fix_f[1]     = BC_periodic_f_y;
            self->fix_j[1]     = BC_periodic_j_y;
            self->fix_P[1]     = BC_periodic_P_y;
#ifdef FULLP
            self->fix_P_winske[1]     = BC_periodic_P_winske_y;
            self->fix_Driver_winske[1]     = BC_periodic_Driver_winske_y;
#endif
        break;


        case BC_TYPE_PERFCOND:
            self->fix_osvs[1]  = BC_perfcond_os_vs_y;
            self->fix_nv[1]    = BC_perfcond_nv_y;
            self->fix_e[1]     = BC_perfcond_e_y;
            self->fix_f[1]     = BC_perfcond_f_y;
            self->fix_j[1]     = BC_perfcond_j_y;
            self->fix_P[1]     = BC_perfcond_P_y;
#ifdef FULLP
            self->fix_P_winske[1]     = BC_perfcond_P_winske_y;
            self->fix_Driver_winske[1]     = BC_perfcond_Driver_winske_y;
#endif
        break;
    } // end switch bc 1


    /* ------------------------------------------------------------- */
    /*                         Z BOUNDARIES                          */
    /* ------------------------------------------------------------- */


    switch(self->bc[2])
    {
        case BC_TYPE_PERIODIC:
            self->fix_osvs[2]  = BC_periodic_os_vs_z;
            self->fix_nv[2]    = BC_periodic_nv_z;
            self->fix_e[2]     = BC_periodic_e_z;
            self->fix_f[2]     = BC_periodic_f_z;
            self->fix_j[2]     = BC_periodic_j_z;
            self->fix_P[2]     = BC_periodic_P_z;
#ifdef FULLP
            self->fix_P_winske[2]     = BC_periodic_P_winske_z;
            self->fix_Driver_winske[2]     = BC_periodic_Driver_winske_z;
#endif
        break;


        case BC_TYPE_PERFCOND:
            self->fix_osvs[2]  = BC_perfcond_os_vs_z;
            self->fix_nv[2]    = BC_perfcond_nv_z;
            self->fix_e[2]     = BC_perfcond_e_z;
            self->fix_f[2]     = BC_perfcond_f_z;
            self->fix_j[2]     = BC_perfcond_j_z;
            self->fix_P[2]     = BC_perfcond_P_z;
#ifdef FULLP
            self->fix_P_winske[2]     = BC_perfcond_P_winske_z;
            self->fix_Driver_winske[2]     = BC_perfcond_Driver_winske_z;
#endif
        break;
    } // end switch bc 2



    // initialize moments and particle boundary condition modules
    // with appropriate boundary condition codes
    self->partbc = PartBCInit(self->bc[0], self->bc[1], self->bc[2], 50000);


    return self;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  HeckleBCDelete()
  ---------------------------------------------------------------------------
  AIM : delete a heckleBC object
 ---------------------------------------------------------------------------*/
void HeckleBCDelete(HeckleBC *self)
{
    if (self)
    {
        PartBCDelete(self->partbc);
        free(self);
    }
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  HeckleBCFieldApply()
  ---------------------------------------------------------------------------
  AIM : apply boundary condition for the selected quantity
 ---------------------------------------------------------------------------*/
void HeckleBCFieldApply(HeckleBC *self, const STX* const sx, ST2 *s2, int qtyID)
{
    int c;

    switch (qtyID)
    {
        case BC_VAR_OSVS:
            for (c=0; c < 3; c++)
            {
                self->fix_osvs[c](sx, s2);
            }
        break;


        case BC_VAR_NV:
            for (c=0; c < 3; c++)
            {
                self->fix_nv[c](sx, s2);
            }
        break;



        case BC_VAR_E:
            for(c=0; c < 3; c++)
            {
                self->fix_e[c](sx, s2);
            }
        break;


        case BC_VAR_F:
            for(c=0; c < 3; c++)
            {
                self->fix_f[c](sx, s2);
            }
        break;


        case BC_VAR_J:
            for(c=0; c < 3; c++)
            {
                self->fix_j[c](sx, s2);
            }
        break;


        case BC_VAR_P:
            for(c=0; c < 3; c++)
            {
                self->fix_P[c](sx, s2);
            }
        break;
#ifdef FULLP
        case BC_VAR_P_FULLP:
            for(c=0; c < 3; c++)
            {
                self->fix_P_winske[c](sx, s2);
            }
        break;

        case BC_VAR_DRIVER_FULLP:
             for(c=0; c < 3; c++)
             {
                self->fix_Driver_winske[c](sx, s2);
             }
        break;
#endif
    }
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  HeckleBCGetPartBC()
  ---------------------------------------------------------------------------
  AIM : returns a handle on the particle boundary condition
 ---------------------------------------------------------------------------*/
PartBC* HeckleBCGetPartBC(HeckleBC *self)
{
    return self->partbc;
}
/*===========================================================================*/

















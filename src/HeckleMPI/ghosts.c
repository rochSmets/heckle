

#include <stdlib.h>


#include<mpi.h>



#include <ghostg2mom.h>
#include <ghostg2fields.h>
#include <ghostg2p.h>
#include <ghostg2pwinske.h>
#include <ghostg4mom.h>
#include <ghostg4current.h>
#include <ghostg4winske.h>
#include <structures.h>





/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */



#define GHOST_NS_VS          0
#define GHOST_NV             1
#define GHOST_P              2
#define GHOST_E              3
#define GHOST_F              4
#define GHOST_J              5

#define GHOST_P_FULLP       6
#define GHOST_DRIVER_FULLP  7
#define GHOST_P_FULLP_G4    8
#define GHOST_J_G4           9


typedef struct ghosts_s
{
    GhostG2Mom     *ghost_ns_vs;
    GhostG2Fields  *ghost_e;
    GhostG2Fields  *ghost_f;
    GhostG2Fields  *ghost_j;

    // TODO to implement
    GhostG4Mom     *ghost_nv;
    GhostG4Current *ghost_current;
    GhostG2P       *ghost_Pi;
    #ifdef FULLP
    GhostG4Winske       *ghost_g4winske;
    GhostG2Winske       *ghost_P_Winske;
    GhostG2Winske       *ghost_Driver_Winske;
    #endif
} Ghosts;






/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */







/*---------------------------------------------------------------------------
  GhostsInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
Ghosts* GhostsInit(const STI * const si, const STX * const sx)
{

    Ghosts *self;

    self = malloc(sizeof *self);


    self->ghost_ns_vs = GhostG2MomInit(sx, si->n);

    self->ghost_e     = GhostG2FieldsInit(sx, GHOSTFIELD_E);
    self->ghost_f     = GhostG2FieldsInit(sx, GHOSTFIELD_F);
    self->ghost_j     = GhostG2FieldsInit(sx, GHOSTFIELD_J);


    // TODO implement those function
    self->ghost_nv    = GhostG4MomInit(sx);
    self->ghost_current    = GhostG4CurrentInit(sx);
    self->ghost_Pi    = GhostG2PInit(sx);
    #ifdef FULLP
    self->ghost_g4winske         = GhostG4WinskeInit(sx);
    self->ghost_P_Winske         = GhostG2WinskeInit(sx, GHOSTFIELD_P);
    self->ghost_Driver_Winske    = GhostG2WinskeInit(sx, GHOSTFIELD_DRIVER);
    #endif
    return self;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  GhostsDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostsDelete(Ghosts* self)
{
    if (self)
    {
        GhostG2MomDelete(self->ghost_ns_vs);

        GhostG2FieldsDelete(self->ghost_e);
        GhostG2FieldsDelete(self->ghost_f);
        GhostG2FieldsDelete(self->ghost_j);

        GhostG2PDelete(self->ghost_Pi);
        #ifdef FULLP
        GhostG4WinskeDelete(self->ghost_g4winske);
        GhostG2WinskeDelete(self->ghost_P_Winske);
        GhostG2WinskeDelete(self->ghost_Driver_Winske);
        #endif
        GhostG4MomDelete(self->ghost_nv);
        GhostG4CurrentDelete(self->ghost_current);

    }
    free(self);
}
/*===========================================================================*/












/*---------------------------------------------------------------------------
  GhostSendRecv()
  ---------------------------------------------------------------------------
  AIM : Communicate ghost points for the specified quantity
 ---------------------------------------------------------------------------*/
void GhostsSendRecv(Ghosts* self, const STX* const sx, void *data, int qtyID)
{
    ST2 *s2;
    ST4 *s4;

    s2 = (ST2*)data;
    s4 = (ST4*)data;


    switch( qtyID)
    {
        case GHOST_NS_VS:
            GhostG2MomSendRecv(self->ghost_ns_vs, sx, s2);
        break;

        case GHOST_NV:
            GhostG4MomSendRecv(self->ghost_nv, sx, s4);
        break;

        case GHOST_J_G4:
            GhostG4CurrentSendRecv(self->ghost_current, sx, s4);
        break;

        case GHOST_P:
            GhostG2PSendRecv(self->ghost_Pi, sx, s2);
        break;

        #ifdef FULLP
        case GHOST_P_FULLP_G4:
            GhostG4WinskeSendRecv(self->ghost_g4winske, sx, s4);
        break;

        case GHOST_P_FULLP:
            GhostG2WinskeSendRecv(self->ghost_P_Winske, sx, s2);
        break;

        case GHOST_DRIVER_FULLP:
            GhostG2WinskeSendRecv(self->ghost_Driver_Winske, sx, s2);
        break;
        #endif

        case GHOST_E:
            GhostG2FieldsSendRecv(self->ghost_e, sx, s2);
        break;

        case GHOST_F:
            GhostG2FieldsSendRecv(self->ghost_f, sx, s2);
        break;

        case GHOST_J:
            GhostG2FieldsSendRecv(self->ghost_j, sx, s2);
        break;

        default:
            printf("Error - Invalid ghost type %s %d\n", __FILE__, __LINE__);
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
}












































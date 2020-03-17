

#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>

#include <mpi.h>


static MPI_Datatype mpighost;



/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */





/* this structure holds the indices of ghost points
   i,j,k represent the dimensions while 0 and 1 mean
   the first and last points on the border to be transfered


    ex. for n and V we want to transfer the column 0 and 1
    to the neighbor on the left, i0=0 and i1 = 1.

       we we'd have liked to send column 1 i0=i1=1 */
typedef struct ghostindices_s
{
    int i0, i1;
    int j0, j1;
    int k0, k1;
}GhostIndices;






/* this is the structure holding the data for the ghost moments
   each ghost node has the density n and the 3 components of the
   velocity V, all this is defined for all ion species.
*/
typedef struct ghostpoint_s
{
    double n;
    double vx,vy,vz;
} GhostPoint;






typedef struct ghostg4mom_s
{
    // send and receive buffers for moments n and V
    GhostPoint *sendbuff[27];
    GhostPoint *recvbuff[27];

    // number of ghost points to be exchanged
    int counter[27];

    // indices of ghost points
    GhostIndices isend[27];
    GhostIndices irecv[27];

} GhostG4Mom;









/*---------------------------------------------------------------------------
  findGhostIndices()
  ---------------------------------------------------------------------------
  AIM : define ghost indices for the density and the bulk velocity. Those are
  the first and second points on each border
 ---------------------------------------------------------------------------*/
static void findGhostIndices(GhostG4Mom *self, const STX* const sx)
{

    int a, b, c;
    int t;
    int i,j,k;

    /* First define the indices of ghost nodes for each neighbor */

    for (a = -1; a <= 1; a++)
    {
        for (b = -1; b <= 1; b++)
        {
            for (c = -1; c <= 1; c++)
            {
                t = (1+c)+3*((1+b)+3*(1+a));

                switch (a)
                {
                    case -1 :
                        self->isend[t].i0 = 3;
                        self->isend[t].i1 = 3;
                        self->irecv[t].i0 = sx->n[0]+3;
                        self->irecv[t].i1 = sx->n[0]+3;
                    break;


                    case 0 :
                        self->isend[t].i0 = 1;
                        self->isend[t].i1 = sx->n[0]+2;
                        self->irecv[t].i0 = 1;
                        self->irecv[t].i1 = sx->n[0]+2;
                    break;


                    case 1 :
                        self->isend[t].i0 = sx->n[0];
                        self->isend[t].i1 = sx->n[0];
                        self->irecv[t].i0 = 0;
                        self->irecv[t].i1 = 0;
                    break;
                } // end switch a






                switch (b)
                {
                    case -1 :
                        self->isend[t].j0 = 3;
                        self->isend[t].j1 = 3;
                        self->irecv[t].j0 = sx->n[1]+3;
                        self->irecv[t].j1 = sx->n[1]+3;
                    break;


                    case 0 :
                        self->isend[t].j0 = 1;
                        self->isend[t].j1 = sx->n[1]+2;
                        self->irecv[t].j0 = 1;
                        self->irecv[t].j1 = sx->n[1]+2;
                    break;


                    case 1 :
                        self->isend[t].j0 = sx->n[1];
                        self->isend[t].j1 = sx->n[1];
                        self->irecv[t].j0 = 0;
                        self->irecv[t].j1 = 0;
                    break;
                } // end switch a





                switch (c)
                {
                    case -1 :
                        self->isend[t].k0 = 3;
                        self->isend[t].k1 = 3;
                        self->irecv[t].k0 = sx->n[2]+3;
                        self->irecv[t].k1 = sx->n[2]+3;
                    break;


                    case 0 :
                        self->isend[t].k0 = 1;
                        self->isend[t].k1 = sx->n[2]+2;
                        self->irecv[t].k0 = 1;
                        self->irecv[t].k1 = sx->n[2]+2;
                    break;


                    case 1 :
                        self->isend[t].k0 = sx->n[2];
                        self->isend[t].k1 = sx->n[2];
                        self->irecv[t].k0 = 0;
                        self->irecv[t].k1 = 0;
                    break;
                } // end switch a

            } // end loop on c
        } // end loop on b
    } // end loop on a



    // now count how many elements we need to exhange with each neighbor

    // first initialize counter to 0
    for (t = 0; t < 27; t++)
    {
        self->counter[t] = 0;

        for (i = self->isend[t].i0; i <= self->isend[t].i1; i++)
        {
            for (j = self->isend[t].j0; j <= self->isend[t].j1; j++)
            {
                for (k = self->isend[t].k0; k <= self->isend[t].k1; k++)
                {
                    self->counter[t]++;
                }
            }
        }
    }
}
/*===========================================================================*/








/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */







/*---------------------------------------------------------------------------
  GhostG4MomInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
GhostG4Mom* GhostG4MomInit(const STX * const sx)
{
    int t;
    GhostG4Mom* self = malloc(sizeof *self);


    // determine the indices of ghost nodes
    findGhostIndices(self, sx);


    // at this point, we know ghost node indices and how many of them we have
    // so we can allocate memory for buffers for incoming/outgoing data

    for (t=0; t < 27; t++)
    {
        self->sendbuff[t] = malloc(self->counter[t] * sizeof * self->sendbuff[t]);
        self->recvbuff[t] = malloc(self->counter[t] * sizeof * self->recvbuff[t]);
    }



    // now define the MPI Datatype
    int count;
    int blocklength[4];
    MPI_Aint displacements[4];
    MPI_Datatype types[4];
    MPI_Aint startaddress, tmpaddress;
    GhostPoint gptmp;

    count = 4;

    blocklength[0] = 1; //  n
    blocklength[1] = 1; //  vx
    blocklength[2] = 1; //  vy
    blocklength[3] = 1; //  vz

    MPI_Get_address(&gptmp, &startaddress);
    MPI_Get_address(&gptmp.n, &tmpaddress);
    displacements[0] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vx, &tmpaddress);
    displacements[1] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vy, &tmpaddress);
    displacements[2] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vz, &tmpaddress);
    displacements[3] = tmpaddress - startaddress;

    types[0] = MPI_DOUBLE;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;
    types[3] = MPI_DOUBLE;

    MPI_Type_create_struct(count, blocklength ,displacements, types, &mpighost);
    MPI_Type_commit(&mpighost);


    return self;
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  GhostG4MomDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG4MomDelete(GhostG4Mom* self)
{
    int t;

    if (self)
    {
        for (t=0; t < 27; t++)
        {
            free(self->sendbuff[t]);
            free(self->recvbuff[t]);
        }

        free(self);
    }
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  GhostG4MomSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG4MomSendRecv(GhostG4Mom* self, const STX * const sx, ST4 * s4)
{
    int t;
    int i,j,k,ijk;
    int is0, is1, js0, js1, ks0, ks1;
    int ir0, ir1, jr0, jr1, kr0, kr1;
    int cpt = 0;
    MPI_Status st;
    int n0, n1, n2;

    n0 = sx->n[0]+4;
    n1 = sx->n[1]+4;
    n2 = sx->n[2]+4;


    for (t=0; t < 27; t++)
    {
        cpt = 0;


        is0 = self->isend[t].i0;
        is1 = self->isend[t].i1;
        js0 = self->isend[t].j0;
        js1 = self->isend[t].j1;
        ks0 = self->isend[t].k0;
        ks1 = self->isend[t].k1;

        // fill the buffer
        for (i = is0; i <= is1; i++)
        {
            for (j = js0; j <= js1; j++)
            {
                for (k = ks0; k <= ks1; k++)
                {
                    ijk = IDX(i,j,k,n0,n1,n2);

                    self->sendbuff[t][cpt].n  = s4[ijk].n;
                    self->sendbuff[t][cpt].vx = s4[ijk].v[0];
                    self->sendbuff[t][cpt].vy = s4[ijk].v[1];
                    self->sendbuff[t][cpt].vz = s4[ijk].v[2];
                    cpt++; // one more ghost point stored
                }// end k loop
            } // end j loop
        } // end i loop
    } // end neighbor loop


    // ok now all ghost nodes have been stored in the outgoing buffers
    // proceed to sending to neighbors


    for (t=0; t < 27; t++)
    {
        if (t != 13) // I'm not sending to myself
        {
            MPI_Sendrecv(self->sendbuff[t],         // send buffer
                         self->counter[t],          // # of ghost nodes to send
                         mpighost,                  // MPI datatype
                         sx->nt[t],                 // destination process,
                         t,                         // message send tag,
                         self->recvbuff[t],         // receiving buffer
                         self->counter[t],          // # of ghost nodes to recv
                         mpighost,                  // recv datatype
                         sx->nf[t],                 // reception process,
                         t,                         // message recv tag
                         MPI_COMM_WORLD,
                         &st);
        }
    } //end neighbor loop


    // ok now all the received data is stored in recvbuff
    // we need to unpack it in

    // I've not received anything from myself, so avoid unpacking t=13 recvbuffs
    // and also do not unpack data if receiving neighbor is MPI_PROC_NULL
    // because that is the border of the simulation domain and those points
    // will be fixed by boundary conditions
    for (t=0; t < 27; t++)
    {
        if (t != 13 && sx->nf[t] != MPI_PROC_NULL)
        {
            cpt = 0;

            ir0 = self->irecv[t].i0;
            ir1 = self->irecv[t].i1;
            jr0 = self->irecv[t].j0;
            jr1 = self->irecv[t].j1;
            kr0 = self->irecv[t].k0;
            kr1 = self->irecv[t].k1;

            for (i = ir0; i <= ir1; i++)
            {
                for (j = jr0; j <= jr1; j++)
                {
                    for (k = kr0; k <= kr1; k++)
                    {
                        ijk = IDX(i,j,k,n0,n1,n2);

                        s4[ijk].n    = self->recvbuff[t][cpt].n;
                        s4[ijk].v[0] = self->recvbuff[t][cpt].vx;
                        s4[ijk].v[1] = self->recvbuff[t][cpt].vy;
                        s4[ijk].v[2] = self->recvbuff[t][cpt].vz;

                        cpt++; // oe more ghost node unpacked
                    } // end k loop
                } // end j loop
            } // end i loop
        } // end of not 13 and not proc_null
    } // end neighbor loop

}
/*===========================================================================*/

























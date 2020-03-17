

#include <stdlib.h>
#include <stdio.h>


#include <structures.h>
#include <defines.h>

#include <mpi.h>


static MPI_Datatype mpighost;



// fonction de reset des packs (entre E et J par ex.)
// fonction a appeler dans init et qui switch entre prendre 2 points sur chaque cote (moments)
// et prendre que le point ext a recevoir et interieur a envoyer (E,J).
// voir ce qui est fait pour n, v, E et J lors du unpack et de l'application des BCs



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
    double n[NS];
    double vx[NS],vy[NS],vz[NS];
} GhostPoint;






typedef struct ghostg2mom_s
{
    // global number of cells
    int ng[3];

    // send and receive buffers for moments n and V
    GhostPoint *sendbuff[27];
    GhostPoint *recvbuff[27];

    // number of ghost points to be exchanged
    int counter[27];

    // indices of ghost points
    GhostIndices isend[27];
    GhostIndices irecv[27];

} GhostG2Mom;









/*---------------------------------------------------------------------------
  findGhostIndices()
  ---------------------------------------------------------------------------
  AIM : define ghost indices for the density and the bulk velocity. Those are
  the first and second points on each border
 ---------------------------------------------------------------------------*/
static void findGhostIndices(GhostG2Mom *self, const STX* const sx)
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
                        self->isend[t].i0 = 0;
                        self->isend[t].i1 = 1;
                        self->irecv[t].i0 = sx->n[0];
                        self->irecv[t].i1 = sx->n[0]+1;
                    break;


                    case 0 :
                        self->isend[t].i0 = 0;
                        self->isend[t].i1 = sx->n[0]+1;
                        self->irecv[t].i0 = 0;
                        self->irecv[t].i1 = sx->n[0]+1;
                    break;


                    case 1 :
                        self->isend[t].i0 = sx->n[0];
                        self->isend[t].i1 = sx->n[0]+1;
                        self->irecv[t].i0 = 0;
                        self->irecv[t].i1 = 1;
                    break;
                } // end switch a






                switch (b)
                {
                    case -1 :
                        self->isend[t].j0 = 0;
                        self->isend[t].j1 = 1;
                        self->irecv[t].j0 = sx->n[1];
                        self->irecv[t].j1 = sx->n[1]+1;
                    break;


                    case 0 :
                        self->isend[t].j0 = 0;
                        self->isend[t].j1 = sx->n[1]+1;
                        self->irecv[t].j0 = 0;
                        self->irecv[t].j1 = sx->n[1]+1;
                    break;


                    case 1 :
                        self->isend[t].j0 = sx->n[1];
                        self->isend[t].j1 = sx->n[1]+1;
                        self->irecv[t].j0 = 0;
                        self->irecv[t].j1 = 1;
                    break;
                } // end switch a





                switch (c)
                {
                    case -1 :
                        self->isend[t].k0 = 0;
                        self->isend[t].k1 = 1;
                        self->irecv[t].k0 = sx->n[2];
                        self->irecv[t].k1 = sx->n[2]+1;
                    break;


                    case 0 :
                        self->isend[t].k0 = 0;
                        self->isend[t].k1 = sx->n[2]+1;
                        self->irecv[t].k0 = 0;
                        self->irecv[t].k1 = sx->n[2]+1;
                    break;


                    case 1 :
                        self->isend[t].k0 = sx->n[2];
                        self->isend[t].k1 = sx->n[2]+1;
                        self->irecv[t].k0 = 0;
                        self->irecv[t].k1 = 1;
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
  GhostG2MomInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications (n, V)
 ---------------------------------------------------------------------------*/
GhostG2Mom* GhostG2MomInit(const STX * const sx, const int n[3])
{
    int t;
    GhostG2Mom* self = malloc(sizeof *self);

    self->ng[0] = n[0];
    self->ng[1] = n[1];
    self->ng[2] = n[2];


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

    blocklength[0] = NS; //  n
    blocklength[1] = NS; //  vx
    blocklength[2] = NS; //  vy
    blocklength[3] = NS; //  vz

    MPI_Get_address(&gptmp, &startaddress);
    MPI_Get_address(&gptmp.n[0], &tmpaddress);
    displacements[0] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vx[0], &tmpaddress);
    displacements[1] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vy[0], &tmpaddress);
    displacements[2] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vz[0], &tmpaddress);
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
  GhostG2MomDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG2MomDelete(GhostG2Mom* self)
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
  GhostG2SendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG2MomSendRecv(GhostG2Mom* self, const STX * const sx, ST2 * s2)
{
    int t;
    int i,j,k,ijk;
    int is0, is1, js0, js1, ks0, ks1;
    int ir0, ir1, jr0, jr1, kr0, kr1;
    int ispe;
    int cpt = 0;
    MPI_Status st;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;

#if 1
    for (t=0; t < 27; t++)
    {
        cpt = 0;


        is0 = self->isend[t].i0;
        is1 = self->isend[t].i1;
        js0 = self->isend[t].j0;
        js1 = self->isend[t].j1;
        ks0 = self->isend[t].k0;
        ks1 = self->isend[t].k1;
/*
        ir0 = self->irecv[t].i0;
        ir1 = self->irecv[t].i1;
        jr0 = self->irecv[t].j0;
        jr1 = self->irecv[t].j1;
        kr0 = self->irecv[t].k0;
        kr1 = self->irecv[t].k1;
*/


        // fill the buffer
        for (i = is0; i <= is1; i++)
        {
            for (j = js0; j <= js1; j++)
            {
                for (k = ks0; k <= ks1; k++)
                {
                    ijk = IDX(i,j,k,n0,n1,n2);

                    // loop over ion species
                    for (ispe = 1; ispe < NS+1; ispe++)
                    {
                        //printf("cpt = %d\n",cpt);
                        self->sendbuff[t][cpt].n[ispe-1]  = s2[ijk].os[ispe];
                        self->sendbuff[t][cpt].vx[ispe-1] = s2[ijk].vs[ispe][0];
                        self->sendbuff[t][cpt].vy[ispe-1] = s2[ijk].vs[ispe][1];
                        self->sendbuff[t][cpt].vz[ispe-1] = s2[ijk].vs[ispe][2];
                    }// end specie loop
                    cpt++; // one more ghost point stored
                }// end k loop
            } // end j loop
        } // end i loop

    } // end neighbor loop
#endif
    //MPI_Abort(MPI_COMM_WORLD,-1);
        // ok now all ghost nodes have been stored in the outgoing buffers
        // proceed to sending to neighbors

    for (t=0; t < 27; t++)
    {
#if 1
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
#endif
    } //end neighbor loop


    // ok now all the received data is stored in recvbuff
    // we need to unpack it in
#if 1
    // I've not received anything from myself, so avoid unpacking t=13 recvbuffs
    // and also do not unpack data if receiving neighbor is MPI_PROC_NULL
    // because that is the border of the simulation domain and those points
    // will be fixed by boundary conditions
    for (t=0; t < 27; t++)
    {
        if (t != 13 && sx->nf[t] != MPI_PROC_NULL)
        {
            cpt = 0;
/*
            is0 = self->isend[t].i0;
            is1 = self->isend[t].i1;
            js0 = self->isend[t].j0;
            js1 = self->isend[t].j1;
            ks0 = self->isend[t].k0;
            ks1 = self->isend[t].k1;
*/
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

                        for (ispe=1; ispe < NS+1; ispe++)
                        {
                            s2[ijk].os[ispe]    += self->recvbuff[t][cpt].n[ispe-1];
                            s2[ijk].vs[ispe][0] += self->recvbuff[t][cpt].vx[ispe-1];
                            s2[ijk].vs[ispe][1] += self->recvbuff[t][cpt].vy[ispe-1];
                            s2[ijk].vs[ispe][2] += self->recvbuff[t][cpt].vz[ispe-1];
                        } // end species loop
                        cpt++; // oe more ghost node unpacked
                    } // end k loop
                } // end j loop
            } // end i loop
        } // end of not 13 and not proc_null
    } // end neighbor loop
#endif

#if 1
    /* --------------------------------------------------------------------- */
    /*              SPECIAL TREATMENT FOR 1D AND 2D                          */
    /* --------------------------------------------------------------------- */


    /* __ need special attention if # of cell = 1 in x direction __ */
    if (self->ng[0] == 1)
    {
        /* __ in the x direction : appropriate nested loops __ */
        for (j = 0; j < n1; j++)
        {
            for (k = 0; k < n2; k++)
            {
                for (ispe=1; ispe < NS+1; ispe++)
                {
                    /* __ index on g2 __ */
                    ijk0 = IDX(0, j, k, n0, n1, n2);
                    ijk  = IDX(1, j, k, n0, n1, n2);
                    ijk1 = IDX(2, j, k, n0, n1, n2);

                    /* __ update the edges of the subdomain in x direction __ */
                    s2[ijk0].os[ispe]    = s2[ijk].os[ispe];
                    s2[ijk0].vs[ispe][0] = s2[ijk].vs[ispe][0];
                    s2[ijk0].vs[ispe][1] = s2[ijk].vs[ispe][1];
                    s2[ijk0].vs[ispe][2] = s2[ijk].vs[ispe][2];

                    s2[ijk1].os[ispe]    = s2[ijk].os[ispe];
                    s2[ijk1].vs[ispe][0] = s2[ijk].vs[ispe][0];
                    s2[ijk1].vs[ispe][1] = s2[ijk].vs[ispe][1];
                    s2[ijk1].vs[ispe][2] = s2[ijk].vs[ispe][2];

                }
            }
        }
    }// end if ng==1


    /* __ need special attention if # of cell = 1 in x direction __ */
    if (self->ng[1] == 1)
    {
        /* __ in the x direction : appropriate nested loops __ */
        for (i = 0; i < n0; i++)
        {
            for (k = 0; k < n2; k++)
            {
                for (ispe=1; ispe < NS+1; ispe++)
                {
                    /* __ index on g2 __ */
                    ijk0 = IDX(i, 0, k, n0, n1, n2);
                    ijk  = IDX(i, 1, k, n0, n1, n2);
                    ijk1 = IDX(i, 2, k, n0, n1, n2);

                    /* __ update the edges of the subdomain in x direction __ */
                    s2[ijk0].os[ispe]    = s2[ijk].os[ispe];
                    s2[ijk0].vs[ispe][0] = s2[ijk].vs[ispe][0];
                    s2[ijk0].vs[ispe][1] = s2[ijk].vs[ispe][1];
                    s2[ijk0].vs[ispe][2] = s2[ijk].vs[ispe][2];

                    s2[ijk1].os[ispe]    = s2[ijk].os[ispe];
                    s2[ijk1].vs[ispe][0] = s2[ijk].vs[ispe][0];
                    s2[ijk1].vs[ispe][1] = s2[ijk].vs[ispe][1];
                    s2[ijk1].vs[ispe][2] = s2[ijk].vs[ispe][2];

                }
            }
        }
    }// end if ng==1



    /* __ need special attention if # of cell = 1 in x direction __ */
    if (self->ng[2] == 1)
    {
        /* __ in the x direction : appropriate nested loops __ */
        for (i = 0; i < n0; i++)
        {
            for (j = 0; j < n1; j++)
            {
                for (ispe=1; ispe < NS+1; ispe++)
                {
                    /* __ index on g2 __ */
                    ijk0 = IDX(i, j, 0, n0, n1, n2);
                    ijk  = IDX(i, j, 1, n0, n1, n2);
                    ijk1 = IDX(i, j, 2, n0, n1, n2);

                    /* __ update the edges of the subdomain in x direction __ */
                    s2[ijk0].os[ispe]    = s2[ijk].os[ispe];
                    s2[ijk0].vs[ispe][0] = s2[ijk].vs[ispe][0];
                    s2[ijk0].vs[ispe][1] = s2[ijk].vs[ispe][1];
                    s2[ijk0].vs[ispe][2] = s2[ijk].vs[ispe][2];

                    s2[ijk1].os[ispe]    = s2[ijk].os[ispe];
                    s2[ijk1].vs[ispe][0] = s2[ijk].vs[ispe][0];
                    s2[ijk1].vs[ispe][1] = s2[ijk].vs[ispe][1];
                    s2[ijk1].vs[ispe][2] = s2[ijk].vs[ispe][2];
                }
            }
        }
    }// end if ng==1

#endif

}
/*===========================================================================*/
























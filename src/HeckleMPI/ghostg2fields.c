
#include <stdio.h>
#include <stdlib.h>


//#include<ghostg2fields.h>
#include<structures.h>
#include<mpi.h>



#define GHOSTFIELD_E 0
#define GHOSTFIELD_F 1
#define GHOSTFIELD_J 2



static MPI_Datatype mpighost;




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




/* this is the structure holding the data for the ghost point */
typedef struct ghostpoint_s
{
    double vx,vy,vz;
} GhostPoint;






typedef struct ghostg2fields_s
{
    // send and receive buffers for moments n and V
    GhostPoint *sendbuff[27];
    GhostPoint *recvbuff[27];

    // number of ghost points to be exchanged
    int counter[27];

    // indices of ghost points
    GhostIndices isend[27];
    GhostIndices irecv[27];


    void (*GetV)(const ST2 * const s2,
              const STX * const sx,
              int i, int j, int k,
              double *vx, double *vy, double *vz);

    void (*SetV)(double vx, double vy, double vz,
              int i, int j, int k,
              const STX* sx,
              ST2 *s2);


} GhostG2Fields;







/*---------------------------------------------------------------------------
  findGhostIndicesFields()
  ---------------------------------------------------------------------------
  AIM : define ghost indices for the fields (electric, current)
 ---------------------------------------------------------------------------*/
void findGhostIndices(GhostG2Fields *self, const STX* const sx)
{

    int a, b, c;
    int t;
    int i, j, k;

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
                        self->isend[t].i0 = 1;
                        self->isend[t].i1 = 1;
                        self->irecv[t].i0 = sx->n[0]+1;
                        self->irecv[t].i1 = sx->n[0]+1;
                    break;


                    case 0 :
                        self->isend[t].i0 = 1;
                        self->isend[t].i1 = sx->n[0];
                        self->irecv[t].i0 = 1;
                        self->irecv[t].i1 = sx->n[0];
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
                        self->isend[t].j0 = 1;
                        self->isend[t].j1 = 1;
                        self->irecv[t].j0 = sx->n[1]+1;
                        self->irecv[t].j1 = sx->n[1]+1;
                    break;


                    case 0 :
                        self->isend[t].j0 = 1;
                        self->isend[t].j1 = sx->n[1];
                        self->irecv[t].j0 = 1;
                        self->irecv[t].j1 = sx->n[1];
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
                        self->isend[t].k0 = 1;
                        self->isend[t].k1 = 1;
                        self->irecv[t].k0 = sx->n[2]+1;
                        self->irecv[t].k1 = sx->n[2]+1;
                    break;


                    case 0 :
                        self->isend[t].k0 = 1;
                        self->isend[t].k1 = sx->n[2];
                        self->irecv[t].k0 = 1;
                        self->irecv[t].k1 = sx->n[2];
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


    for (t=0; t < 27; t++)
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





#define ex(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[0]
#define ey(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[1]
#define ez(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[2]


#define fx(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].f[0]
#define fy(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].f[1]
#define fz(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].f[2]


#define jx(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].j[0]
#define jy(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].j[1]
#define jz(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].j[2]



void GetE(const ST2 * const s2,
          const STX * const sx,
          int i, int j, int k,
          double *vx, double *vy, double *vz)
{
   int nyg2, nzg2;

   nyg2 = sx->n[1]+2;
   nzg2 = sx->n[2]+2;

    *vx = ex(i,j,k);
    *vy = ey(i,j,k);
    *vz = ez(i,j,k);
}



void GetF(const ST2 * const s2,
          const STX * const sx,
          int i, int j, int k,
          double *vx, double *vy, double *vz)
{
    int nyg2, nzg2;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;

     *vx = fx(i,j,k);
     *vy = fy(i,j,k);
     *vz = fz(i,j,k);
}




void GetJ(const ST2 * const s2,
          const STX * const sx,
          int i, int j, int k,
          double *vx, double *vy, double *vz)
{
    int nyg2, nzg2;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;

     *vx = jx(i,j,k);
     *vy = jy(i,j,k);
     *vz = jz(i,j,k);
}




void SetE(double vx, double vy, double vz,
          int i, int j, int k,
          const STX* sx,
          ST2 *s2)
{
    int nyg2, nzg2;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;


    ex(i,j,k) = vx;
    ey(i,j,k) = vy;
    ez(i,j,k) = vz;
}





void SetF(double vx, double vy, double vz,
          int i, int j, int k,
          const STX* sx,
          ST2 *s2)
{
    int nyg2, nzg2;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;

    fx(i,j,k) = vx;
    fy(i,j,k) = vy;
    fz(i,j,k) = vz;
}




void SetJ(double vx, double vy, double vz,
          int i, int j, int k,
          const STX* sx,
          ST2 *s2)
{
    int nyg2, nzg2;

    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;

    jx(i,j,k) = vx;
    jy(i,j,k) = vy;
    jz(i,j,k) = vz;
}






/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */







/*---------------------------------------------------------------------------
  GhostG2FieldsInit()
  ---------------------------------------------------------------------------
  AIM :  Initialize the module for ghost communications
 ---------------------------------------------------------------------------*/
GhostG2Fields* GhostG2FieldsInit(const STX * const sx, int fieldID)
{

    int t;
    GhostG2Fields *self = malloc(sizeof *self);


    // determine the indices of ghost nodes
    findGhostIndices(self, sx);


    // at this point, we know ghost node indices and how many of them we have
    // so we can allocate memory for buffers for incoming/outgoing data


    for (t=0; t < 27; t++)
    {
        self->sendbuff[t] = malloc(self->counter[t] * sizeof * self->sendbuff[t]);
        self->recvbuff[t] = malloc(self->counter[t] * sizeof * self->recvbuff[t]);
    }


    switch (fieldID)
    {

        case GHOSTFIELD_E:
            self->GetV = GetE;
            self->SetV = SetE;
        break;

        case GHOSTFIELD_F:
            self->GetV = GetF;
            self->SetV = SetF;
        break;

        case GHOSTFIELD_J:
            self->GetV = GetJ;
            self->SetV = SetJ;
        break;


        default:
            printf("error - ghostfield ID %d invalid\n", fieldID);
            MPI_Abort(MPI_COMM_WORLD, -1);
            exit(-1);
    }




    // now define the MPI Datatype
    int count;
    int blocklength[3];
    MPI_Aint displacements[3];
    MPI_Datatype types[3];
    MPI_Aint startaddress, tmpaddress;
    GhostPoint gptmp;

    count = 3;

    blocklength[0] = 1; //  x
    blocklength[1] = 1; //  y
    blocklength[2] = 1; //  z

    MPI_Get_address(&gptmp, &startaddress);
    MPI_Get_address(&gptmp.vx, &tmpaddress);
    displacements[0] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vy, &tmpaddress);
    displacements[1] = tmpaddress - startaddress;

    MPI_Get_address(&gptmp.vz, &tmpaddress);
    displacements[2] = tmpaddress - startaddress;


    types[0] = MPI_DOUBLE;
    types[1] = MPI_DOUBLE;
    types[2] = MPI_DOUBLE;

    MPI_Type_create_struct(count, blocklength ,displacements, types, &mpighost);
    MPI_Type_commit(&mpighost);


    return self;
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  GhostG2FieldsDelete()
  ---------------------------------------------------------------------------
  AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG2FieldsDelete(GhostG2Fields *self)
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
  GhostG2FieldsSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG2FieldsSendRecv(GhostG2Fields *self, const STX * const sx, ST2 * s2)
{
    int t;
    int i,j,k;
    int is0, is1, js0, js1, ks0, ks1;
    int ir0, ir1, jr0, jr1, kr0, kr1;
    int cpt = 0;
    MPI_Status st;
   // int n0, n1, n2;

  //  n0 = sx->n[0]+2;
    //n1 = sx->n[1]+2;
   // n2 = sx->n[2]+2;


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

                    self->GetV(s2, sx, i,j,k,
                               &self->sendbuff[t][cpt].vx,
                               &self->sendbuff[t][cpt].vy,
                               &self->sendbuff[t][cpt].vz);

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
                        self->SetV(self->recvbuff[t][cpt].vx,
                                   self->recvbuff[t][cpt].vy,
                                   self->recvbuff[t][cpt].vz,
                                   i,j,k,
                                   sx, s2);

                        cpt++; // oe more ghost node unpacked
                    } // end k loop
                } // end j loop
            } // end i loop
        } // end of not 13 and not proc_null
    } // end neighbor loop

}
/*===========================================================================*/







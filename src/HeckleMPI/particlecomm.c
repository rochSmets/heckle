
#include<particlecomm.h>
#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>




/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */


// this is a global variable for this module
// it is the mpi datatype of a particle
static MPI_Datatype mpiparticle;









/*---------------------------------------------------------------------------
  PpackInit()
  ---------------------------------------------------------------------------
  AIM : initialize a Ppack object and allocate memory for particle arrays
  of a given 'size'.
  PARAMETERS : 'size' is the size of the particle buffers, not the actual
                number of those particles packed into te buffer (npart)
 ---------------------------------------------------------------------------*/
Ppack *PpackInit(int size)
{
    Ppack *self;

    self = malloc(sizeof *self);

    if (self == NULL)
    {
        printf("Error (Ppack) - memory allocation impossible -"
               " %s in %s line %d\n", __func__, __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }

    // size of the particle buffer
    // and number of particles in it set to 0
    self->size  = size;
    self->npart = 0;

    self->parr = malloc(size * sizeof *self->parr);

    return self;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PpackDelete()
  ---------------------------------------------------------------------------
  AIM : delete a Ppack object
 ---------------------------------------------------------------------------*/
void PpackDelete(Ppack *self)
{

    if (self)
    {
        free(self->parr);
        free(self);
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PpackRealloc()
  ---------------------------------------------------------------------------
  AIM : NOT IMPLEMENTED - this function reallocates memory for particle buffers
  of Ppack objects in case they are too small
 ---------------------------------------------------------------------------*/
void PpackRealloc(Ppack *self, int partn)
{
    // reallocation of the pack if needed
    (void)self; // to turn off Not Used warnings
    (void)partn; // to turn off Not Used Warnings
}
/*===========================================================================*/











/*---------------------------------------------------------------------------
  PpackIt()
  ---------------------------------------------------------------------------
  AIM : Put a given particle in the buffer of the Ppack object 'self'
 ---------------------------------------------------------------------------*/
void PpackIt(Ppack *self, int ispe, Particle *p)
{

#if 0 // not implemented - buffer reallocation in case outgoing buffers are too small
    if (self->npart == self->size)
    {
        PpackRealloc(self, -1);
    }
#endif

    self->parr[self->npart].r[0] = p->r[0];
    self->parr[self->npart].r[1] = p->r[1];
    self->parr[self->npart].r[2] = p->r[2];

    self->parr[self->npart].s[0] = p->s[0];
    self->parr[self->npart].s[1] = p->s[1];
    self->parr[self->npart].s[2] = p->s[2];

    self->parr[self->npart].v[0] = p->v[0];
    self->parr[self->npart].v[1] = p->v[1];
    self->parr[self->npart].v[2] = p->v[2];

    self->parr[self->npart].w[0] = p->w[0];
    self->parr[self->npart].w[1] = p->w[1];
    self->parr[self->npart].w[2] = p->w[2];

    self->parr[self->npart].b[0] = p->b[0];
    self->parr[self->npart].b[1] = p->b[1];
    self->parr[self->npart].b[2] = p->b[2];

    self->parr[self->npart].ispe = ispe;
    self->parr[self->npart].i    = p->i;

    self->npart++;
}
/*===========================================================================*/












/*---------------------------------------------------------------------------
  PpackGetPart()
  ---------------------------------------------------------------------------
  AIM : Returns the particle at index 'iparticle' of the Ppack 'self'
  in the particle 'p'. Also returns it species type 'ispe'. The returned
  values are copies of what is in the Ppack, this is not a 'pop' routine.
 ---------------------------------------------------------------------------*/
void PpackGetPart(Ppack *self, int iparticle, Particle *p, int *ispe)
{
        p->r[0] = self->parr[iparticle].r[0];
        p->r[1] = self->parr[iparticle].r[1];
        p->r[2] = self->parr[iparticle].r[2];

        p->s[0] = self->parr[iparticle].s[0];
        p->s[1] = self->parr[iparticle].s[1];
        p->s[2] = self->parr[iparticle].s[2];

        p->v[0] = self->parr[iparticle].v[0];
        p->v[1] = self->parr[iparticle].v[1];
        p->v[2] = self->parr[iparticle].v[2];

        p->w[0] = self->parr[iparticle].w[0];
        p->w[1] = self->parr[iparticle].w[1];
        p->w[2] = self->parr[iparticle].w[2];

        p->b[0] = self->parr[self->npart].b[0];
        p->b[1] = self->parr[self->npart].b[1];
        p->b[2] = self->parr[self->npart].b[2];

        p->i    = self->parr[iparticle].i;

        *ispe   = self->parr[iparticle].ispe;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  Ppack_debug()
  ---------------------------------------------------------------------------
  AIM : this function aims at displaying informations about Ppack objects
 ---------------------------------------------------------------------------*/
void Ppack_debug(Ppack *self)
{
    if (self == NULL)
    {
        printf("ERROR - invalid Ppack object pointer (NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    else
    {
        printf("Ppack object holds %d particles and can store %d\n",
               self->npart, self->size);
        fflush(stdout);
    }
}
/*===========================================================================*/





/*---------------------------------------------------------------------------
  PartComms_debug()
  ---------------------------------------------------------------------------
  AIM : this function displays some stuff to check objects PartComms
 ---------------------------------------------------------------------------*/
void PartComms_debug(PartComms *self)
{
    int ipack;

    if (self == NULL)
    {
        printf("ERROR - invalid PartComms object pointer (NULL)\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }


    printf("Incoming buffers\n");

    for (ipack = 0; ipack < 27; ipack++)
    {
        Ppack_debug(self->packsin[ipack]);
    }

    printf("Outgoing buffers\n");

    for (ipack = 0; ipack < 27; ipack++)
    {
        Ppack_debug(self->packsout[ipack]);
    }
}
/*===========================================================================*/





void leavingP(PartComms *self, const STX * const sx, int p)
{
    int ipack;
    int ipart;

    if (sx->r == p)
    {
        for (ipack=0; ipack < 27; ipack++)
        {
            printf("Outgoing Pack %d (neighbor %d) : %d particles \n", ipack,
                   sx->nt[ipack], self->packsout[ipack]->npart);

            for (ipart=0; ipart < self->packsout[ipack]->npart; ipart++)
            {
                printf("Particle %04d : index %07d, positions : %f %f %f"
                       " | %f %f %f"
                       " | species : %d\n",
                       ipart, self->packsout[ipack]->parr[ipart].i,
                       self->packsout[ipack]->parr[ipart].r[0],
                       self->packsout[ipack]->parr[ipart].r[1],
                       self->packsout[ipack]->parr[ipart].r[2],
                       self->packsout[ipack]->parr[ipart].s[0],
                       self->packsout[ipack]->parr[ipart].s[1],
                       self->packsout[ipack]->parr[ipart].s[2],
                       self->packsout[ipack]->parr[ipart].ispe);
            }
        }
    }
}



void arrivingP(PartComms *self, const STX * const sx, int p)
{
    int ipack;
    int ipart;

    if (sx->r == p)
    {
        for (ipack=0; ipack < 27; ipack++)
        {
            printf("Incoming Pack %d (neighbor %d) : %d particles \n", ipack,
                   sx->nf[ipack], self->packsin[ipack]->npart);

            for (ipart=0; ipart < self->packsin[ipack]->npart; ipart++)
            {
                printf("Particle %04d : index %07d, positions : %8.4f %8.4f %8.4f"
                       " | %f %f %f"
                       " | species : %d\n",
                       ipart, self->packsin[ipack]->parr[ipart].i,
                       self->packsin[ipack]->parr[ipart].r[0],
                       self->packsin[ipack]->parr[ipart].r[1],
                       self->packsin[ipack]->parr[ipart].r[2],
                       self->packsin[ipack]->parr[ipart].s[0],
                       self->packsin[ipack]->parr[ipart].s[1],
                       self->packsin[ipack]->parr[ipart].s[2],
                       self->packsin[ipack]->parr[ipart].ispe);
            }
        }
    }
}







/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PUBLIC FUNCTIONS                             //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */









/*---------------------------------------------------------------------------
  PartCommsInit()
  ---------------------------------------------------------------------------
  AIM : initialize the particle communication module.
        is called once at the begining of the simulation
        returns the communication module handle
 ---------------------------------------------------------------------------*/
PartComms* PartCommsInit(void)
{
   int ipack;
   PartComms *self = malloc(sizeof *self);
   int ppacksize;

   // variable for the MPI derived datatype
   int count;
   int blocklength[7];
   MPI_Aint displacements[7];
   MPI_Datatype types[7];
   MPI_Aint startaddress, tmpaddress;
   ParticleMPI pmpi;

   // this is a hard-coded estimation that is totally arbitrary and hopefully
   // way larger than necessary. A more accurate estimate would be helpful
   // although we do not want to have something that involves collective
   // communications.
   ppacksize = 50000;


   // allocate memory for particle packs, incoming and outgoing
   for (ipack=0; ipack < 27; ipack++)
   {
       self->packsin[ipack]  = PpackInit(ppacksize);
       self->packsout[ipack] = PpackInit(ppacksize);
   }


   // now create the MPI derived datatype

   count = 7;

   blocklength[0] = 3;          // r[0..2]
   blocklength[1] = 3;          // s[0..2]
   blocklength[2] = 3;          // v[0..2]
   blocklength[3] = 3;          // w[0..2]
   blocklength[4] = 3;          // b[0..2]
   blocklength[5] = 1;          // i
   blocklength[6] = 1;          // ispe

   MPI_Get_address(&pmpi, &startaddress);
   MPI_Get_address(&pmpi.r[0], &tmpaddress);
   displacements[0] = tmpaddress - startaddress;

   MPI_Get_address(&pmpi.s[0], &tmpaddress);
   displacements[1] = tmpaddress - startaddress;

   MPI_Get_address(&pmpi.v[0], &tmpaddress);
   displacements[2] = tmpaddress - startaddress;

   MPI_Get_address(&pmpi.w[0], &tmpaddress);
   displacements[3] = tmpaddress - startaddress;

   MPI_Get_address(&pmpi.b[0], &tmpaddress);
   displacements[4] = tmpaddress - startaddress;

   MPI_Get_address(&pmpi.i, &tmpaddress);
   displacements[5] = tmpaddress - startaddress;

   MPI_Get_address(&pmpi.ispe, &tmpaddress);
   displacements[6] = tmpaddress - startaddress;


   types[0] = MPI_DOUBLE;
   types[1] = MPI_DOUBLE;
   types[2] = MPI_DOUBLE;
   types[3] = MPI_DOUBLE;
   types[4] = MPI_INT;
   types[5] = MPI_INT;
   types[6] = MPI_INT;

    MPI_Type_create_struct(count,
                           blocklength,
                           displacements,
                           types,
                           &mpiparticle);

    MPI_Type_commit(&mpiparticle);

   return self;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  PartCommsDelete()
  ---------------------------------------------------------------------------
  AIM : delete the communication module
 ---------------------------------------------------------------------------*/
void PartCommsDelete(PartComms *self)
{
    int ipack;

    if (self)
    {
        for (ipack = 0; ipack < 27; ipack++)
        {
            if (ipack != 13)
            {
                PpackDelete(self->packsin[ipack]);
                PpackDelete(self->packsout[ipack]);
            }
        } // end ipack loop
    } // end if self
}
/*===========================================================================*/











/*---------------------------------------------------------------------------
  PartCommsPack()
  ---------------------------------------------------------------------------
  AIM : pack the particle 'p' of species 'ispe' into the pack associated
        with particles leaving for the partition 'ineighb'
 ---------------------------------------------------------------------------*/
void PartCommsPack(PartComms *self, Particle *p, int ispe, int ineighb)
{
    Ppack *ppack;

    // convenience pointer
    ppack = self->packsout[ineighb];

#if 0
    if (ppack->npart  == ppack->size)
    {
        // reallocation
    }
#endif

    // pack the particle 'p' of species 'ispe' into 'ppack'
    PpackIt(ppack, ispe, p);
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PartCommsResetPacks()
  ---------------------------------------------------------------------------
  AIM : This function resets the number of particles in incoming/outgoing
  buffers. It should be always used after send/recv operations.
 ---------------------------------------------------------------------------*/
void PartCommsResetPacks(PartComms *self)
{
    int ipack;

    if (self)
    {
        for (ipack=0; ipack < 27; ipack++)
        {
            if (ipack != 27)
            {
                self->packsin[ipack]->npart  = 0;
                self->packsout[ipack]->npart = 0;
            }
        }
    }
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  PartCommsSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive particles from and in outgoing and incoming internal
  buffers. Received particles are placed in the particle array of the
  correct species.
 ---------------------------------------------------------------------------*/
void PartCommsSendRecv(PartComms *self, STI *si, STX *sx, Particle *sp[NS+1])
{
    int ineighb;
    MPI_Status st;
    struct stp p;
    int ispe;
    int ip;


    //leavingP(self, sx, 3);



    // loop over the 26 neighbors (27 cells in the array)
   for (ineighb=0; ineighb < 27; ineighb++)
   {
       // we exclude 13, which is the current processor
       if (ineighb != 13)
       {
           /* send and receive particles from neighbor */

           //   !!! WARNING !!!
           //   WE NOW ASSUME THERE IS ENOUGH ROOM IN PPACK INCOMING BUFFERS
           //   FOR ALL INCOMING PARTICLES
           //   !!! WARNING !!!

           MPI_Sendrecv(self->packsout[ineighb]->parr,      // Outgoing buffer
                        self->packsout[ineighb]->npart,     // # of particles sent
                        mpiparticle,                        // MPI derived type
                        sx->nt[ineighb],                    // destination
                        ineighb,                            // tag
                        self->packsin[ineighb]->parr,       // Incoming buffer
                        self->packsin[ineighb]->size,       // inc. buf. size
                        mpiparticle,                        // inc. datatype
                        sx->nf[ineighb],                    // source neighb.
                        ineighb,                            // tag
                        MPI_COMM_WORLD,                     // communicator
                        &st);                               // satus

           // now count how many particles we have received
           MPI_Get_count(&st, mpiparticle, &self->packsin[ineighb]->npart);


           // ok now unpack all particles into the particle array
           for (int iparticle = 0; iparticle < self->packsin[ineighb]->npart; iparticle++)
           {
                // get a particle from the buffer
                // and also its species
                PpackGetPart(self->packsin[ineighb],
                             iparticle,
                             &p,
                             &ispe);


                //   !!! WARNING !!!
                //   WE NOW ASSUME THERE IS ENOUGH ROOM IN PARTICLE ARRAY
                //   FOR ALL INCOMING PARTICLES
                //   !!! WARNING !!!

                // now we know its species
                // we can have the index of the particle array of that species
                // where the particle will be put
                ip = sx->ns[ispe];
                if (ip >= si->nm) {
                    fprintf(stderr, "you are loading to many particles ! try to increase the si.nm value\n");
                    exit(-1);
                }

                // put the particle (SHOULD CHECK SIZE HERE)
                sp[ispe][ip] = p;

                // there is now one more particle in that species's part. array
                sx->ns[ispe]++;

           }// end loop over incoming particle buffer
       }// end if not 13
   } // end loop over neighbors

  // arrivingP(self, sx, 3);

   // reset the number of particles in incom./outgo. particle buffers.
   PartCommsResetPacks(self);

}
/*===========================================================================*/













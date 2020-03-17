#ifndef PPACK_H
#define PPACK_H

#include<structures.h>
#include<particle.h>


//  this is a structure for particles
// it is different from stp only because it lacks ijk.
// if ijk is not needed anymore, this structure should be struct stp
typedef struct  s_particlempi
{
    double r[3];
    double s[3];
    double v[3];
    double w[3];
    int b[3];
    int i;
    int ispe;
} ParticleMPI;




// this is a particle pack
// it contains an of ParticleMPI items
// of size 'size' and of 'npart' elements
// (size >= npart)
typedef struct s_ppack
{
    ParticleMPI *parr;
    int size;
    int npart;
} Ppack;




typedef struct s_partcomms
{
    Ppack *packsin[27];
    Ppack *packsout[27];
} PartComms;




/*---------------------------------------------------------------------------
  PartCommsInit()
  ---------------------------------------------------------------------------
  AIM : initialize the particle communication module.
        is called once at the begining of the simulation
        returns the communication module handle
 ---------------------------------------------------------------------------*/
PartComms* PartCommsInit(void);
/*===========================================================================*/







/*---------------------------------------------------------------------------
  PartCommsDelete()
  ---------------------------------------------------------------------------
  AIM : delete the communication module
 ---------------------------------------------------------------------------*/
void PartCommsDelete(PartComms *self);
/*===========================================================================*/






/*---------------------------------------------------------------------------
  PartCommsPack()
  ---------------------------------------------------------------------------
  AIM : pack the particle 'p' of species 'ispe' into the pack associated
        with particles leaving for the partition 'ineighb'
 ---------------------------------------------------------------------------*/
void PartCommsPack(PartComms *self, Particle *p, int ispe, int ineighb);
/*===========================================================================*/







/*---------------------------------------------------------------------------
  PartCommsSendRecv()
  ---------------------------------------------------------------------------
  AIM : Send and receive particles from and in outgoing and incoming internal
  buffers. Received particles are placed in the particle array of the
  correct species.
 ---------------------------------------------------------------------------*/
void PartCommsSendRecv(PartComms *self, STI *si, STX *sx, Particle *sp[NS+1]);
/*===========================================================================*/





#endif // PPACK_H


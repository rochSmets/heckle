#ifndef PARTICLEBC_H
#define PARTICLEBC_H

#include<structures.h>
#include<particle.h>
#include<particlecomm.h>


#define ISOUT 1
#define ISIN  0

#define DELETED 2
#define NOTDELETED 3



typedef struct s_partBC
{
    /* 'indices' contains the index of all the particles that are
        detected by the module as 'leaving particles'. Leaving means
        they have crossed the domain boundaries (processor). Some of them
        will actually have crossed the simulation boundaries and boundary
        conditions will be applied to them, others will just move to an
        adjacent processor.
    */
    int *indices;

    /* 'species' holds the species index of all leaving particles */
    int *species;

    /* 'size' is the length of the arrays indices and species, while 'npart'
        is the actual number of items they contains. npart <= size */
    int size;
    int npart;

    /* this is the MPI Particle Communication object */
    PartComms *comms;

    // boundary condition routines

    // NOTE : in some cases those routines will have to delete particles from
    // the particle array (e.g. open boundaries). In this case, 'indices'
    // and 'npart' should be changed accordingly because routines called
    // after will work on these particles.

    int (*bc_x)(const STI* const si,
                const STX* const sx,
                Particle *sp[NS+1],
                int ispecies,
                int ipart,
                int ipc);

    int (*bc_y)(const STI* const si,
                const STX* const sx,
                Particle *sp[NS+1],
                int ispecies,
                int ipart,
                int ipc);

    int (*bc_z)(const STI* const si,
                const STX* const sx,
                Particle *sp[NS+1],
                int ispecies,
                int ipart,
                int ipc);


    /* private */

    int (*A) (const STI* const si, const STX* const sx, double x0, double x1);
    int (*B) (const STI* const si, const STX* const sx, double y0, double y1);
    int (*C) (const STI* const si, const STX* const sx, double z0, double z1);

} PartBC;






/*---------------------------------------------------------------------------
  PartBCInit()
  ---------------------------------------------------------------------------
  AIM : Initializes the object PartBC
 ---------------------------------------------------------------------------*/
PartBC* PartBCInit(int bcx,
                   int bcy, int bcz,
                   int size);






/*---------------------------------------------------------------------------
  PartBCDelete()
  ---------------------------------------------------------------------------
  AIM : delete a PartBC object
 ---------------------------------------------------------------------------*/
void PartBCDelete(PartBC *self);






/*---------------------------------------------------------------------------
  PartBCisOut()
  ---------------------------------------------------------------------------
  AIM : returns 1 if the particle is outside the domain and 0 if not
 ---------------------------------------------------------------------------*/
int PartBCisOut(PartBC *self, const STI * const si,
                const STX* const sx, double pos[3]);






/*---------------------------------------------------------------------------
  PartBCStoreIndex()
  ---------------------------------------------------------------------------
  AIM : stores the index of a particle detected outside the domain
 ---------------------------------------------------------------------------*/
void PartBCStore(PartBC *self, int ip, int ispe);





/*---------------------------------------------------------------------------
  PartBCApply()
  ---------------------------------------------------------------------------
  AIM :
 ---------------------------------------------------------------------------*/
void PartBCApply(PartBC *self, STI * si, STX *sx, Particle *sp[NS+1], int ipc);





/*---------------------------------------------------------------------------
  PartBCReset()
  ---------------------------------------------------------------------------
  AIM : Reset the buffer of leaving particles
 ---------------------------------------------------------------------------*/
void PartBCReset(PartBC *self);





#endif // PARTICLEBC_H












/* --------------------------------------------------------------------- */
/*  MODULE     : PartBC                                                  */
/*  AIM        : Deal with Particle Boundary conditions                  */
/* --------------------------------------------------------------------- */



#include <stdio.h>
#include<stdlib.h>
#include<math.h>



#include<particlecomm.h>
#include<structures.h>
#include <particle.h>
#include<particlebc.h>
#include <hecklebc.h>
#include <bc_constants.h>
#include<stdlib.h>







/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */










/*---------------------------------------------------------------------------
  PartBCisOK()
  ---------------------------------------------------------------------------
  AIM : debug function used to check the object PartBC
 ---------------------------------------------------------------------------*/
int PartBCisOK(const STX * const sx, PartBC *self)
{
    int errorcode = 1;
    int iproc;

    if (self == NULL)
    {
        printf("process %d > ERROR - PartBC * == NULL in %s file %s at %d\n",
               sx->r, __func__, __FILE__, __LINE__);
        errorcode = 0;
    }


    else
    {
        for (iproc = 0; iproc < sx->s; iproc ++)
        {
            if (sx->r == iproc)
            {
                // now cheks indices and species arrays
                printf("process %d > indices address : %p\n", sx->r, self->indices);
                printf("process %d > species address : %p\n", sx->r, self->species);
                printf("process %d > size            : %d\n", sx->r, self->size);
                printf("process %d > npart           : %d\n", sx->r, self->npart);
            }
        }
    }

    return errorcode;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PartBCDeleteParticle()
  ---------------------------------------------------------------------------
  AIM : delete a particle from the particle array
 ---------------------------------------------------------------------------*/
void PartBCDeleteParticle(PartBC *self, Particle *sp[NS+1], STX* sx, int ipart)
{
    int id;
    int index, species;



    index   = self->indices[ipart];
    species = self->species[ipart];



    // swap with the last particle with the last one of the SAME SPECIES
    sp[species][index] = sp[species][sx->ns[species]-1];



    // now deal with indices, npart and species


    // for all particles stored in our buffer
    for (id = ipart+1; id < self->npart; id++)
    {
        // if the one we've just swaped, i.e. the last of the species 'species'
        // is in our buffer
        if (self->indices[id] == sx->ns[species]-1
                && self->species[id] == species)
        {
            //if (sx->r == 3)
            //    printf("oh oh deleted(%d) : indice %d points to %d now moved to %d\n",
            //           ipart, id, sx->ns[species]-1,index);
            // we need to re-direct 'indice' to its new position 'index'
            self->indices[id] = index;
        }
    }

    // there is now one less particle in the array if that species
    sx->ns[species]--;
}
/*===========================================================================*/











/*---------------------------------------------------------------------------
  PartBCPeriodicX()
  ---------------------------------------------------------------------------
  AIM : Periodic boundary condition along the direction X
 ---------------------------------------------------------------------------*/
int PartBCPeriodicX(const STI * const si,
                     const STX * const sx,
                     Particle *sp[NS+1],
                     int ispecies,
                     int ipart,
                     int ipc)
{
    struct stp *pcur;

    // unused but must be defined
    (void)sx;
    (void)ipc;

    pcur = &sp[ispecies][ipart];


    if (pcur->r[0] < 0.0)
    {
        pcur->r[0] += si->l[0];
        pcur->s[0] += si->l[0];
        pcur->b[0] += -1;
    }
    if (pcur->r[0] >= si->l[0])
    {
        pcur->r[0] -= si->l[0];
        pcur->s[0] -= si->l[0];
        pcur->b[0] += 1;
    }

    return NOTDELETED;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  PartBCPeriodicY()
  ---------------------------------------------------------------------------
  AIM : Periodic boundary condition along the direction Y
 ---------------------------------------------------------------------------*/
int PartBCPeriodicY(const STI * const si,
                     const STX * const sx,
                     Particle *sp[NS+1],
                     int ispecies,
                     int ipart,
                     int ipc)
{
    struct stp *pcur;

    // unused but must be defined
    (void)sx;
    (void)ipc;

    pcur = &sp[ispecies][ipart];

    if (pcur->r[1] < 0.0)
    {
        pcur->r[1] += si->l[1];
        pcur->s[1] += si->l[1];
        pcur->b[1] += -1;
    }
    if (pcur->r[1] >= si->l[1])
    {
        pcur->r[1] -= si->l[1];
        pcur->s[1] -= si->l[1];
        pcur->b[1] += 1;
    }

    return NOTDELETED;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  PartBCPeriodicZ()
  ---------------------------------------------------------------------------
  AIM : Perdiodic boundary condition along the direciton Z
 ---------------------------------------------------------------------------*/
int  PartBCPeriodicZ(const STI * const si,
                     const STX * const sx,
                     Particle *sp[NS+1],
                     int ispecies,
                     int ipart,
                     int ipc)
{
    struct stp *pcur;

    // unused but must be defined
    (void)sx;
    (void)ipc;

    pcur = &sp[ispecies][ipart];


    if (pcur->r[2] < 0.0)
    {
        pcur->r[2] += si->l[2];
        pcur->s[2] += si->l[2];
        pcur->b[2] += -1;
    }
    if (pcur->r[2] >= si->l[2])
    {
        pcur->r[2] -= si->l[2];
        pcur->s[2] -= si->l[2];
        pcur->b[2] += 1;
    }

    return NOTDELETED;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  APeriodic()
  ---------------------------------------------------------------------------
  AIM : calculates the index 'a' for periodic boundaries
  for periodic boundaries, the index uses the 'old' position 'x0' since new
  position will be at the other end of the box and the index will not belong to
  [-1,1]
 ---------------------------------------------------------------------------*/
int APeriodic(const STI *const si, const STX * const sx, double x0, double x1)
{
    double xw;
    int a;

    // must be defined, but is not used
    (void) x1;

    xw = (x0 - sx->i0[0]*si->dl[0] ) /sx->l[0];
    a = (int) floor(xw);

    return a;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  BPeriodic()
  ---------------------------------------------------------------------------
  AIM : calculates the index 'b' for periodic boundaries
  for periodic boundaries, the index uses the 'old' position 'x0' since new
  position will be at the other end of the box and the index will not belong to
  [-1,1]
 ---------------------------------------------------------------------------*/
int BPeriodic(const STI *const si, const STX * const sx, double y0, double y1)
{
    double yw;
    int b;

    // must be defined but is not used
    (void)y1;

    yw = (y0 - sx->i0[1]*si->dl[1] ) /sx->l[1];
    b = (int) floor(yw);

    return b;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  CPeriodic()
  ---------------------------------------------------------------------------
  AIM : calculates the index 'c' for periodic boundaries
  for periodic boundaries, the index uses the 'old' position 'x0' since new
  position will be at the other end of the box and the index will not belong to
  [-1,1]
 ---------------------------------------------------------------------------*/
int CPeriodic(const STI *const si, const STX * const sx, double z0, double z1)
{
    double zw;
    int c;

    // must be defined but is not used
    (void)z1;

    zw = (z0 - sx->i0[2]*si->dl[2] ) /sx->l[2];
    c = (int) floor(zw);

    return c;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PartBCReflectionX()
  ---------------------------------------------------------------------------
  AIM : Reflection of the particle reaching a X boundary
 ---------------------------------------------------------------------------*/
int PartBCReflectionX(const STI * const si,
                       const STX * const sx,
                       Particle *sp[NS+1],
                       int ispecies,
                       int ipart,
                       int ipc)
{
    struct stp *pcur;
    // unused but must be defined
    (void)sx;
    (void)ipc;

    pcur = &sp[ispecies][ipart];

    /* contrary to periodic boundary conditions
       reflection should be applied also on the particle at corrector step */

    switch(ipc)
    {

    case 0: // predictor
        if (pcur->r[0] < 0.0)
        {
            pcur->r[0] = -pcur->r[0];
            pcur->v[0] *= -1;
        }
        if (pcur->r[0] >= si->l[0])
        {
            pcur->r[0] = -pcur->r[0] + 2 * si->l[0];
            pcur->v[0] *= -1;
        }
        break;


    case 1: // corrector
        if (pcur->s[0] < 0.0)
        {
            pcur->s[0] = -pcur->s[0];
            pcur->w[0] *= -1;
        }
        if (pcur->s[0] >= si->l[0])
        {
            pcur->s[0] = -pcur->s[0] + 2 * si->l[0];
            pcur->w[0] *= -1;
        }
        break;


    }

    return NOTDELETED;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  PartBCReflectionY()
  ---------------------------------------------------------------------------
  AIM : Reflection of a particle reaching a Y boundary
 ---------------------------------------------------------------------------*/
int PartBCReflectionY(const STI * const si,
                       const STX * const sx,
                       Particle *sp[NS+1],
                       int ispecies,
                       int ipart,
                       int ipc)
{
    struct stp *pcur;
    // unused but must be defined
    (void)sx;
    (void)ipc;

    pcur = &sp[ispecies][ipart];

    /* contrary to periodic boundary conditions
       reflection should be applied also on the particle at corrector step */

    switch(ipc)
    {

    case 0: // predictor
        if (pcur->r[1] < 0.0)
        {
            pcur->r[1] = -pcur->r[1];
            pcur->v[1] *= -1;
        }
        if (pcur->r[1] >= si->l[1])
        {
            pcur->r[1] = -pcur->r[1] + 2 * si->l[1];
            pcur->v[1] *= -1;
        }
        break;


    case 1: // corrector
        if (pcur->s[1] < 0.0)
        {
            pcur->s[1] = -pcur->s[1];
            pcur->w[1] *= -1;
        }
        if (pcur->s[1] >= si->l[1])
        {
            pcur->s[1] = -pcur->s[1] + 2 * si->l[1];
            pcur->w[1] *= -1;
        }
        break;


    }

    return NOTDELETED;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PartBCReflectionZ()
  ---------------------------------------------------------------------------
  AIM : Reflection of a particle reaching a Z boundary
 ---------------------------------------------------------------------------*/
int PartBCReflectionZ(const STI * const si,
                       const STX * const sx,
                       Particle *sp[NS+1],
                       int ispecies,
                       int ipart,
                       int ipc)
{
    struct stp *pcur;
    // unused but must be defined
    (void)sx;
    (void)ipc;

    pcur = &sp[ispecies][ipart];

    /* contrary to periodic boundary conditions
       reflection should be applied also on the particle at corrector step */

    switch(ipc)
    {

    case 0: // predictor
        if (pcur->r[2] < 0.0)
        {
            pcur->r[2] = -pcur->r[2];
            pcur->v[2] *= -1;
        }
        if (pcur->r[2] >= si->l[2])
        {
            pcur->r[2] = -pcur->r[2] + 2 * si->l[2];
            pcur->v[2] *= -1;
        }
        break;


    case 1: // corrector
        if (pcur->s[2] < 0.0)
        {
            pcur->s[2] = -pcur->s[2];
            pcur->w[2] *= -1;
        }
        if (pcur->s[2] >= si->l[2])
        {
            pcur->s[2] = -pcur->s[2] + 2 * si->l[2];
            pcur->w[2] *= -1;
        }
        break;
    }

    return NOTDELETED;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  ANonPer()
  ---------------------------------------------------------------------------
  AIM : calculates the index 'a' for non-periodic boundaries
  for non-periodic boundaries, the index uses the 'new' position 'x1' since old
  position will be outside of the box
 ---------------------------------------------------------------------------*/
int ANonPer(const STI *const si, const STX * const sx, double x0, double x1)
{
    double xw;
    int a;

    // must be defined, but is not used
    (void) x0;

    xw = (x1 - sx->i0[0]*si->dl[0] ) /sx->l[0];
    a = (int) floor(xw);

    return a;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  BNonPer()
  ---------------------------------------------------------------------------
  AIM : calculates the index 'b' for non-periodic boundaries
  for non-periodic boundaries, the index uses the 'new' position 'y1' since old
  position will be outside of the box
 ---------------------------------------------------------------------------*/
int BNonPer(const STI *const si, const STX * const sx, double y0, double y1)
{
    double yw;
    int b;

    // must be defined but is not used
    (void)y0;

    yw = (y1 - sx->i0[1]*si->dl[1] ) /sx->l[1];
    b = (int) floor(yw);

    return b;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  CNonPer()
  ---------------------------------------------------------------------------
  AIM : calculates the index 'c' for non-periodic boundaries
  for non-periodic boundaries, the index uses the 'new' position 'z1' since old
  position will be outside of the box
 ---------------------------------------------------------------------------*/
int CNonPer(const STI *const si, const STX * const sx, double z0, double z1)
{
    double zw;
    int c;

    // must be defined but is not used
    (void)z0;

    zw = (z1 - sx->i0[2]*si->dl[2] ) /sx->l[2];
    c = (int) floor(zw);

    return c;
}
/*===========================================================================*/








/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */







/*---------------------------------------------------------------------------
  PartBCInit()
  ---------------------------------------------------------------------------
  AIM : Initializes the object PartBC
 ---------------------------------------------------------------------------*/
PartBC* PartBCInit(int bcx, int bcy, int bcz, int size)
{
    PartBC *self = malloc(sizeof *self);

    // initializes the arrays storing information for leaving particles

    self->indices = malloc(size * sizeof *self->indices);
    if (self->indices == NULL)
    {
        printf("Error - cannot allocate memory in %s of file %s at line %d\n",
               __func__, __FILE__, __LINE__);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }




    self->species = malloc(size * sizeof *self->species);
    if (self->species == NULL)
    {
        printf("Error - cannot allocate memory in %s of file %s at line %d\n",
               __func__, __FILE__, __LINE__);

        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }

    self->size = size;
    self->npart = 0;


    // initializes the MPI communication module
    self->comms = PartCommsInit();


    // now setup boundary conditions

    switch(bcx)
    {
    case BC_TYPE_PERIODIC:
        self->bc_x = PartBCPeriodicX;
        self->A    = APeriodic;
        break;

    case BC_TYPE_PERFCOND:
        self->bc_x = PartBCReflectionX;
        self->A    = ANonPer;
        break;

    default:
        printf("Error - Boundary condition not implemented : %s in %s at %d\n",
               __func__, __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }



    switch(bcy)
    {
    case BC_TYPE_PERIODIC:
        self->bc_y = PartBCPeriodicY;
        self->B    = BPeriodic;
        break;

    case BC_TYPE_PERFCOND:
        self->bc_y = PartBCReflectionY;
        self->B    = BNonPer;
        break;

    default:
        printf("Error - Boundary condition not implemented : %s in %s at %d\n",
               __func__, __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }


    switch(bcz)
    {
    case BC_TYPE_PERIODIC:
        self->bc_z = PartBCPeriodicZ;
        self->C    = CPeriodic;
        break;

    case BC_TYPE_PERFCOND:
        self->bc_z = PartBCReflectionZ;
        self->C    = CNonPer;
        break;

    default:
        printf("Error - Boundary condition not implemented : %s in %s at %d\n",
               __func__, __FILE__, __LINE__);
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }

//    PartBCisOK(si, sx, self);
    return self;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PartBCDelete()
  ---------------------------------------------------------------------------
  AIM : delete a PartBC object
 ---------------------------------------------------------------------------*/
void PartBCDelete(PartBC *self)
{
    if(self)
    {
        if (self->indices) free(self->indices);
        if (self->species) free(self->species);
        if (self->comms) PartCommsDelete(self->comms);

        free(self);
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  PartBCisOut()
  ---------------------------------------------------------------------------
  AIM : returns 1 if the particle is outside the domain and 0 if not
 ---------------------------------------------------------------------------*/
int PartBCisOut(PartBC *self, const STI * const si,
                const STX* const sx, double pos[3])
{
    double xw, yw, zw;
    int a, b, c;
    int t;

    // unused but must be defined
    (void)self;

    xw = (pos[0] - sx->i0[0] * si->dl[0])/sx->l[0];
    yw = (pos[1] - sx->i0[1] * si->dl[1])/sx->l[1];
    zw = (pos[2] - sx->i0[2] * si->dl[2])/sx->l[2];

    a = (int)floor(xw);
    b = (int)floor(yw);
    c = (int)floor(zw);

    t = (1+c)+3*((1+b)+3*(1+a));

    return (t == 13) ? ISIN : ISOUT;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  PartBCStoreIndex()
  ---------------------------------------------------------------------------
  AIM : stores the index of a particle detected outside the domain
 ---------------------------------------------------------------------------*/
void PartBCStore(PartBC *self, int ip, int ispe)
{
#if 0
    if (self->npart == self->size)
    {
        // realloc memory
    }
#endif



    self->indices[self->npart] = ip;
    self->species[self->npart] = ispe;
    self->npart++;
}
/*===========================================================================*/





/*---------------------------------------------------------------------------
  PartBCReset()
  ---------------------------------------------------------------------------
  AIM : Reset the buffer of leaving particles
 ---------------------------------------------------------------------------*/
void PartBCReset(PartBC *self)
{
    self->npart = 0;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  PartBCApply()
  ---------------------------------------------------------------------------
  AIM :
 ---------------------------------------------------------------------------*/
void PartBCApply(PartBC *self, STI * si, STX *sx, Particle *sp[NS+1], int ipc)
{

    int ipart;
    int ispe, ip;
    int t;
    int deleted;
    int a, b, c;
    double x0, y0, z0, x1, y1, z1;


    // for all particles detected as outside, first apply
    // boundary condition
    for (ipart = 0; ipart < self->npart; ipart++)
    {
        ispe = self->species[ipart];
        ip   = self->indices[ipart];



        // save the position before applying the boundary condition

        x0 = sp[ispe][ip].r[0];
        y0 = sp[ispe][ip].r[1];
        z0 = sp[ispe][ip].r[2];


        // then apply boundary conditions.


        deleted = self->bc_x(si, sx, sp, self->species[ipart],
                             self->indices[ipart], ipc);

        if (deleted == DELETED)
            continue; //part. has been deleted, nothing else to do, proceed to next

        deleted = self->bc_y(si, sx, sp, self->species[ipart],
                             self->indices[ipart], ipc);

        if (deleted == DELETED)
            continue;//part. has been deleted, nothing else to do, proceed to next

        deleted = self->bc_z(si, sx, sp, self->species[ipart],
                             self->indices[ipart], ipc);

        if (deleted == DELETED)
            continue;//part. has been deleted, nothing else to do, proceed to next




        x1 = sp[ispe][ip].r[0];
        y1 = sp[ispe][ip].r[1];
        z1 = sp[ispe][ip].r[2];



        a = self->A(si, sx, x0, x1);
        b = self->B(si, sx, y0, y1);
        c = self->C(si, sx, z0, z1);

        t = (1+c)+3*((1+b)+3*(1+a));


        if (t != 13) // if the particle predictor position is still out
        {
            // pack the particle

            PartCommsPack(self->comms,
                          &sp[self->species[ipart]][self->indices[ipart]],
                          self->species[ipart],
                          t);


            // delete the particle from the buffer.
            //'indices' might have changed since push() because of bx*Apply()
            // but everything should be ok.
            PartBCDeleteParticle(self, sp, sx, ipart);

        }
    }

   PartCommsSendRecv(self->comms, si, sx, sp);

}
/*===========================================================================*/





























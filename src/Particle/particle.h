#ifndef PARTICLE_H
#define PARTICLE_H

#include <structures.h>


/* _____ struct for the particles position and velocity _____ */
typedef struct stp
{
    int i;             /* __ index of particle __ */
    double r[3];       /* __ part. position __ */
    double s[3];       /* __ proj. part. position (correction) __ */
    double v[3];       /* __ part. velocity __ */
    double w[3];       /* __ proj. part. velocity (correction) __ */
    int b[3];          /* __ # of crossing of the domain __ */
    int ijk;           /* __ index of the grid point in g2 __ */
} Particle;




typedef struct s_particle_dbg ParticleDBG;






/*---------------------------------------------------------------------------
  ParticleDBG_NanInf()
  ---------------------------------------------------------------------------
  AIM : checks if a particle has NaNs in its attributes.
 ---------------------------------------------------------------------------*/
int ParticleDBG_NanInf(const Particle * const sp, ParticleDBG *pdbg);





/*---------------------------------------------------------------------------
  ParticleDBG_IsInside()
  ---------------------------------------------------------------------------
  AIM : this is a debug routine that tests whether all particles are within
  the simulation domain.
 ---------------------------------------------------------------------------*/
int ParticleDBG_IsInside(const STI * const si, const Particle * const sp);







/*---------------------------------------------------------------------------
  pdbg_init()
  ---------------------------------------------------------------------------
  AIM : initialize the particle debug object
 ---------------------------------------------------------------------------*/
ParticleDBG* pdbg_init(void);





/*---------------------------------------------------------------------------
  pdbg_init()
  ---------------------------------------------------------------------------
  AIM : deletes the particle debug object
 ---------------------------------------------------------------------------*/
void pdbg_delete(ParticleDBG* pdbg);






#endif // PARTICLE_H




#include <math.h>
#include <particle.h>
#include <structures.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>



struct s_particle_dbg
{
    int nan_r;
    int nan_s;
    int nan_v;
    int nan_w;

    int inf_r;
    int inf_s;
    int inf_v;
    int inf_w;

};






/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           DEBUG FUNCTIONS                             //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */





/*---------------------------------------------------------------------------
  ParticleDBG_NanInf()
  ---------------------------------------------------------------------------
  AIM : checks if a particle has NaNs in its attributes.
 ---------------------------------------------------------------------------*/
int ParticleDBG_NanInf(const Particle * const sp, ParticleDBG *pdbg)
{
    int ok = 1;
    int c;


    for (c=0; c < 3; c++)
    {
        if (isnan(sp->r[c]) == 1)
        {
            ok = 0;
            pdbg->nan_r ++;
        }

        if (isnan(sp->s[c]) == 1)
        {
            ok  = 0;
            pdbg->nan_s ++;
        }

        if(isnan(sp->v[c]) == 1)
        {
            ok = 0;
            pdbg->nan_v ++;
        }

        if(isnan(sp->w[c]) == 1)
        {
            ok = 0;
            pdbg->nan_w ++;
        }

        if (isinf(sp->r[c]) == 1)
        {
            ok = 0;
            pdbg->inf_r ++;
        }

        if (isinf(sp->s[c]) == 1)
        {
            ok = 0;
            pdbg->inf_s ++;
        }

        if (isinf(sp->v[c]) == 1)
        {
            ok = 0;
            pdbg->inf_v ++;
        }

        if (isinf(sp->w[c]) == 1)
        {
            ok = 0;
            pdbg->inf_w ++;
        }
    }

    return ok;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  ParticleDBG_IsInside()
  ---------------------------------------------------------------------------
  AIM : this is a debug routine that tests whether all particles are within
  the simulation domain.
 ---------------------------------------------------------------------------*/
int ParticleDBG_IsInside(const STI * const si, const Particle * const sp)
{
    int ok = 1;

    /**/
    if (sp->r[0] < 0)
    {
        ok = 0;
    }
    if (sp->r[0] >= si->l[0])
    {
        ok = 0;
    }
    if (sp->r[1] <  0)
    {
        ok = 0;
    }
    if (sp->r[1] >=  si->l[1])
    {
        ok = 0;
    }
    return ok;
}
/*===========================================================================*/




/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */



/*---------------------------------------------------------------------------
  pdbg_init()
  ---------------------------------------------------------------------------
  AIM : initialize the particle debug object
 ---------------------------------------------------------------------------*/
ParticleDBG* pdbg_init(void)
{
    ParticleDBG* pdbg;

    pdbg = malloc(sizeof *pdbg);

    if (pdbg != NULL)
    {
        pdbg->inf_r = 0;
        pdbg->inf_s = 0;
        pdbg->inf_v = 0;
        pdbg->inf_w = 0;

        pdbg->nan_r = 0;
        pdbg->nan_s = 0;
        pdbg->nan_v = 0;
        pdbg->nan_w = 0;
    }
    else
    {
        fprintf(stderr, "ERROR - Invalid ParticleDBG pointer\n");
        exit(-1);
    }

    return pdbg;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  pdbg_init()
  ---------------------------------------------------------------------------
  AIM : deletes the particle debug object
 ---------------------------------------------------------------------------*/
void pdbg_delete(ParticleDBG* pdbg)
{
    if (pdbg)
    {
        free(pdbg);
    }
    else
    {
        fprintf(stderr,"ERROR - Invalid ParticleDBG pointer : %s at %d\n",
                __FILE__,__LINE__);
    }
}
/*===========================================================================*/













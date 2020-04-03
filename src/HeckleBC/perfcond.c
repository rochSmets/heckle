
#ifndef PERFCOND_H
#define PERFCOND_H



#include <structures.h>
#include <defines.h>
#include <hecklebc.h>


/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */





/*---------------------------------------------------------------------------
  BC_periodic_ns_vs_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_os_vs_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    double tmp;
    int ispe;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);


                for (ispe=1; ispe <= NS; ispe++)
                {
                    tmp = s2[ijk0].os[ispe];
                    s2[ijk0].os[ispe] += s2[ijk1].os[ispe];
                    s2[ijk1].os[ispe] += tmp;

                    for (c=0; c < 3 ;c++)
                    {
                        tmp = s2[ijk0].vs[ispe][c];
                        s2[ijk0].vs[ispe][c] += s2[ijk1].vs[ispe][c];
                        s2[ijk1].vs[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);


                for (ispe=1; ispe <= NS; ispe++)
                {
                    tmp = s2[ijk0].os[ispe];
                    s2[ijk0].os[ispe] += s2[ijk1].os[ispe];
                    s2[ijk1].os[ispe] += tmp;

                    for (c=0; c < 3 ;c++)
                    {
                        tmp = s2[ijk0].vs[ispe][c];
                        s2[ijk0].vs[ispe][c] += s2[ijk1].vs[ispe][c];
                        s2[ijk1].vs[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  BC_perfcond_os_vs_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_os_vs_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    double tmp;
    int ispe;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);


                for (ispe=1; ispe <= NS; ispe++)
                {
                    tmp = s2[ijk0].os[ispe];
                    s2[ijk0].os[ispe] += s2[ijk1].os[ispe];
                    s2[ijk1].os[ispe] += tmp;

                    for (c=0; c < 3 ;c++)
                    {
                        tmp = s2[ijk0].vs[ispe][c];
                        s2[ijk0].vs[ispe][c] += s2[ijk1].vs[ispe][c];
                        s2[ijk1].vs[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);


                for (ispe=1; ispe <= NS; ispe++)
                {
                    tmp = s2[ijk0].os[ispe];
                    s2[ijk0].os[ispe] += s2[ijk1].os[ispe];
                    s2[ijk1].os[ispe] += tmp;

                    for (c=0; c < 3 ;c++)
                    {
                        tmp = s2[ijk0].vs[ispe][c];
                        s2[ijk0].vs[ispe][c] += s2[ijk1].vs[ispe][c];
                        s2[ijk1].vs[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  BC_perfcond_os_vs_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_os_vs_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;
    double tmp;
    int ispe;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);


                for (ispe=1; ispe <= NS; ispe++)
                {
                    tmp = s2[ijk0].os[ispe];
                    s2[ijk0].os[ispe] += s2[ijk1].os[ispe];
                    s2[ijk1].os[ispe] += tmp;

                    for (c=0; c < 3 ;c++)
                    {
                        tmp = s2[ijk0].vs[ispe][c];
                        s2[ijk0].vs[ispe][c] += s2[ijk1].vs[ispe][c];
                        s2[ijk1].vs[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);


                for (ispe=1; ispe <= NS; ispe++)
                {
                    tmp = s2[ijk0].os[ispe];
                    s2[ijk0].os[ispe] += s2[ijk1].os[ispe];
                    s2[ijk1].os[ispe] += tmp;

                    for (c=0; c < 3 ;c++)
                    {
                        tmp = s2[ijk0].vs[ispe][c];
                        s2[ijk0].vs[ispe][c] += s2[ijk1].vs[ispe][c];
                        s2[ijk1].vs[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  BC_perfcond_nv_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_nv_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;
#if 1

    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

                s2[ijk0].ns[0] =  s2[ijk1].ns[0];

                s2[ijk0].vi[0] = -s2[ijk1].vi[0];
                s2[ijk0].vi[1] =  s2[ijk1].vi[1];
                s2[ijk0].vi[2] =  s2[ijk1].vi[2];

                s2[ijk0].j_smoothed[0] = -s2[ijk1].j_smoothed[0];
                s2[ijk0].j_smoothed[1] =  s2[ijk1].j_smoothed[1];
                s2[ijk0].j_smoothed[2] =  s2[ijk1].j_smoothed[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

                s2[ijk0].ns[0] =  s2[ijk1].ns[0];

                s2[ijk0].vi[0] = -s2[ijk1].vi[0];
                s2[ijk0].vi[1] =  s2[ijk1].vi[1];
                s2[ijk0].vi[2] =  s2[ijk1].vi[2];

                s2[ijk0].j_smoothed[0] = -s2[ijk1].j_smoothed[0];
                s2[ijk0].j_smoothed[1] =  s2[ijk1].j_smoothed[1];
                s2[ijk0].j_smoothed[2] =  s2[ijk1].j_smoothed[2];
            } //end k loop
        } // end j loop
    }
#endif
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  BC_perfcond_nv_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_nv_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;
#if 1

    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

                s2[ijk0].ns[0] =  s2[ijk1].ns[0];

                s2[ijk0].vi[0] =  s2[ijk1].vi[0];
                s2[ijk0].vi[1] = -s2[ijk1].vi[1];
                s2[ijk0].vi[2] =  s2[ijk1].vi[2];

                s2[ijk0].j_smoothed[0] =  s2[ijk1].j_smoothed[0];
                s2[ijk0].j_smoothed[1] = -s2[ijk1].j_smoothed[1];
                s2[ijk0].j_smoothed[2] =  s2[ijk1].j_smoothed[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

                s2[ijk0].ns[0] =  s2[ijk1].ns[0];

                s2[ijk0].vi[0] =  s2[ijk1].vi[0];
                s2[ijk0].vi[1] = -s2[ijk1].vi[1];
                s2[ijk0].vi[2] =  s2[ijk1].vi[2];

                s2[ijk0].j_smoothed[0] =  s2[ijk1].j_smoothed[0];
                s2[ijk0].j_smoothed[1] = -s2[ijk1].j_smoothed[1];
                s2[ijk0].j_smoothed[2] =  s2[ijk1].j_smoothed[2];
            } //end k loop
        } // end j loop
    }
#endif
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  BC_perfcond_nv_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_nv_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);


                s2[ijk0].ns[0] =  s2[ijk1].ns[0];

                s2[ijk0].vi[0] =  s2[ijk1].vi[0];
                s2[ijk0].vi[1] =  s2[ijk1].vi[1];
                s2[ijk0].vi[2] = -s2[ijk1].vi[2];

                s2[ijk0].j_smoothed[0] =  s2[ijk1].j_smoothed[0];
                s2[ijk0].j_smoothed[1] =  s2[ijk1].j_smoothed[1];
                s2[ijk0].j_smoothed[2] = -s2[ijk1].j_smoothed[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

                s2[ijk0].ns[0] =   s2[ijk1].ns[0];

                s2[ijk0].vi[0] =   s2[ijk1].vi[0];
                s2[ijk0].vi[1] =   s2[ijk1].vi[1];
                s2[ijk0].vi[2] = - s2[ijk1].vi[2];

                s2[ijk0].j_smoothed[0] =   s2[ijk1].j_smoothed[0];
                s2[ijk0].j_smoothed[1] =   s2[ijk1].j_smoothed[1];
                s2[ijk0].j_smoothed[2] = - s2[ijk1].j_smoothed[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  BC_perfcond_e_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_e_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;


    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

                s2[ijk0].e[0] =  s2[ijk1].e[0];
                s2[ijk0].e[1] = -s2[ijk1].e[1];
                s2[ijk0].e[2] = -s2[ijk1].e[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

                s2[ijk0].e[0] = s2[ijk1].e[0];
                s2[ijk0].e[1] = -s2[ijk1].e[1];
                s2[ijk0].e[2] = -s2[ijk1].e[2];
            } //end k loop
        } // end j loop
    }

}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  BC_perfcond_e_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_e_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

                s2[ijk0].e[0] = -s2[ijk1].e[0];
                s2[ijk0].e[1] =  s2[ijk1].e[1];
                s2[ijk0].e[2] = -s2[ijk1].e[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

                s2[ijk0].e[0] = -s2[ijk1].e[0];
                s2[ijk0].e[1] =  s2[ijk1].e[1];
                s2[ijk0].e[2] = -s2[ijk1].e[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  BC_perfcond_e_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_e_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);

                s2[ijk0].e[0] = -s2[ijk1].e[0];
                s2[ijk0].e[1] = -s2[ijk1].e[1];
                s2[ijk0].e[2] =  s2[ijk1].e[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

                s2[ijk0].e[0] = -s2[ijk1].e[0];
                s2[ijk0].e[1] = -s2[ijk1].e[1];
                s2[ijk0].e[2] =  s2[ijk1].e[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/











/*---------------------------------------------------------------------------
  BC_perfcond_f_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_f_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

                s2[ijk0].f[0] =  s2[ijk1].f[0];
                s2[ijk0].f[1] = -s2[ijk1].f[1];
                s2[ijk0].f[2] = -s2[ijk1].f[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

                s2[ijk0].f[0] =  s2[ijk1].f[0];
                s2[ijk0].f[1] = -s2[ijk1].f[1];
                s2[ijk0].f[2] = -s2[ijk1].f[2];
            } //end k loop
        } // end j loop
    }

}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_f_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_f_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

                s2[ijk0].f[0] = -s2[ijk1].f[0];
                s2[ijk0].f[1] =  s2[ijk1].f[1];
                s2[ijk0].f[2] = -s2[ijk1].f[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

                s2[ijk0].f[0] = -s2[ijk1].f[0];
                s2[ijk0].f[1] =  s2[ijk1].f[1];
                s2[ijk0].f[2] = -s2[ijk1].f[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_f_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_f_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);

                s2[ijk0].f[0] = -s2[ijk1].f[0];
                s2[ijk0].f[1] = -s2[ijk1].f[1];
                s2[ijk0].f[2] =  s2[ijk1].f[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

                s2[ijk0].f[0] = -s2[ijk1].f[0];
                s2[ijk0].f[1] = -s2[ijk1].f[1];
                s2[ijk0].f[2] =  s2[ijk1].f[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/












/*---------------------------------------------------------------------------
  BC_perfcond_j_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_j_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

                s2[ijk0].j[0] = -s2[ijk1].j[0];
                s2[ijk0].j[1] =  s2[ijk1].j[1];
                s2[ijk0].j[2] =  s2[ijk1].j[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

                s2[ijk0].j[0] = -s2[ijk1].j[0];
                s2[ijk0].j[1] =  s2[ijk1].j[1];
                s2[ijk0].j[2] =  s2[ijk1].j[2];
            } //end k loop
        } // end j loop
    }

}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_j_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_j_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

                s2[ijk0].j[0] =  s2[ijk1].j[0];
                s2[ijk0].j[1] = -s2[ijk1].j[1];
                s2[ijk0].j[2] =  s2[ijk1].j[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

                s2[ijk0].j[0] =  s2[ijk1].j[0];
                s2[ijk0].j[1] = -s2[ijk1].j[1];
                s2[ijk0].j[2] =  s2[ijk1].j[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_j_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_j_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);

                s2[ijk0].j[0] =  s2[ijk1].j[0];
                s2[ijk0].j[1] =  s2[ijk1].j[1];
                s2[ijk0].j[2] = -s2[ijk1].j[2];
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

                s2[ijk0].j[0] =  s2[ijk1].j[0];
                s2[ijk0].j[1] =  s2[ijk1].j[1];
                s2[ijk0].j[2] = -s2[ijk1].j[2];
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  BC_perfcond_P_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_P_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    double tmp;
    int ispe;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

                for (ispe=1; ispe <= NS; ispe++)
                {
                    for (c=0; c < 6; c++)
                    {
                        tmp = s2[ijk0].ps[ispe][c];
                        s2[ijk0].ps[ispe][c] += s2[ijk1].ps[ispe][c];
                        s2[ijk1].ps[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

                for (ispe=1; ispe <= NS; ispe++)
                {
                    for (c=0; c < 6 ;c++)
                    {
                        tmp = s2[ijk0].ps[ispe][c];
                        s2[ijk0].ps[ispe][c] += s2[ijk1].ps[ispe][c];
                        s2[ijk1].ps[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  BC_perfcond_P_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_P_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    double tmp;
    int ispe;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

                for (ispe=1; ispe <= NS; ispe++)
                {
                    for (c=0; c < 6 ;c++)
                    {
                        tmp = s2[ijk0].ps[ispe][c];
                        s2[ijk0].ps[ispe][c] += s2[ijk1].ps[ispe][c];
                        s2[ijk1].ps[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

                for (ispe=1; ispe <= NS; ispe++)
                {
                    for (c=0; c < 6; c++)
                    {
                        tmp = s2[ijk0].ps[ispe][c];
                        s2[ijk0].ps[ispe][c] += s2[ijk1].ps[ispe][c];
                        s2[ijk1].ps[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_P_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_P_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;
    double tmp;
    int ispe;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);

                for (ispe=1; ispe <= NS; ispe++)
                {
                    for (c=0; c < 6 ;c++)
                    {
                        tmp = s2[ijk0].ps[ispe][c];
                        s2[ijk0].ps[ispe][c] += s2[ijk1].ps[ispe][c];
                        s2[ijk1].ps[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

                for (ispe=1; ispe <= NS; ispe++)
                {
                    for (c=0; c < 6 ;c++)
                    {
                        tmp = s2[ijk0].ps[ispe][c];
                        s2[ijk0].ps[ispe][c] += s2[ijk1].ps[ispe][c];
                        s2[ijk1].ps[ispe][c] += tmp;
                    }
                } // end species loop
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/

/*===========================================================================*/

#ifdef FULLP
/*---------------------------------------------------------------------------
  BC_perfcond_P_winske_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_P_winske_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

				for (c = 0; c < 6; c++) {
					s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
				}
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

				for (c = 0; c < 6; c++) {
					s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
				}
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  BC_perfcond_P_winske_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_P_winske_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

				for (c = 0; c < 6; c++) {
					s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
				}
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

				for (c = 0; c < 6; c++) {
					s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
				}
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_P_winske_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_P_winske_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);

				for (c = 0; c < 6; c++) {
					s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
				}
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

				for (c = 0; c < 6; c++) {
					s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
				}
            } //end k loop
        } // end j loop
    }
}


/*---------------------------------------------------------------------------
  BC_perfcond_Driver_winske_x()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the x direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_Driver_winske_x(const STX* const sx, ST2 *s2)
{
    int j, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                           LEFT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the left neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_LEFT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(0,j,k, n0, n1, n2);
                ijk1 = IDX(1,j,k, n0, n1, n2);

                for (c = 0; c < 6; c++) {
                    s2[ijk0].dFull[c] = s2[ijk1].dFull[c];
                }
            } //end k loop
        } // end j loop
    }



    /* ------------------------------------------------------------- */
    /*                          RIGHT NEIGHBOR                       */
    /* ------------------------------------------------------------- */

    // first check whether the right neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_RIGHT] == MPI_PROC_NULL)
    {
        for (j = 0;  j < n1; j++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(sx->n[0]+1,j,k, n0, n1, n2);
                ijk1 = IDX(sx->n[0]  ,j,k, n0, n1, n2);

                for (c = 0; c < 6; c++) {
                    s2[ijk0].dFull[c] = s2[ijk1].dFull[c];
                }
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  BC_perfcond_Driver_winske_y()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the y direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_Driver_winske_y(const STX* const sx, ST2 *s2)
{
    int i, k;
    int n0, n1, n2;
    int ijk0, ijk1;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BOTTOM NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the bottom neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BOTTOM] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,0,k, n0, n1, n2);
                ijk1 = IDX(i,1,k, n0, n1, n2);

                for (c = 0; c < 6; c++) {
                    s2[ijk0].dFull[c] = s2[ijk1].dFull[c];
                }
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                           TOP NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the top neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_TOP] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (k=0; k < n2; k++)
            {
                ijk0 = IDX(i,sx->n[1]+1,k, n0, n1, n2);
                ijk1 = IDX(i,sx->n[1]  ,k, n0, n1, n2);

                for (c = 0; c < 6; c++) {
                    s2[ijk0].dFull[c] = s2[ijk1].dFull[c];
                }
            } //end k loop
        } // end j loop
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  BC_perfcond_Driver_winske_z()
  ---------------------------------------------------------------------------
  AIM : perfcond boundary condition along the z direction
 ---------------------------------------------------------------------------*/
void BC_perfcond_Driver_winske_z(const STX* const sx, ST2 *s2)
{
    int i, j;
    int n0, n1, n2;
    int ijk0, ijk1;
    int c;

    n0 = sx->n[0]+2;
    n1 = sx->n[1]+2;
    n2 = sx->n[2]+2;


    /* ------------------------------------------------------------- */
    /*                        BACK NEIGHBOR                          */
    /* ------------------------------------------------------------- */

    // first check whether the back neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_BACK] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,0, n0, n1, n2);
                ijk1 = IDX(i,j,1, n0, n1, n2);

                for (c = 0; c < 6; c++) {
                    s2[ijk0].dFull[c] = s2[ijk1].dFull[c];
                }
            } //end k loop
        } // end j loop
    }

    /* ------------------------------------------------------------- */
    /*                         FRONT NEIGHBOR                        */
    /* ------------------------------------------------------------- */

    // first check whether the front neighbor is the left side of the simulation
    // domain, only in that case BC are applied.
    if (sx->nt[NEIGHBOR_FRONT] == MPI_PROC_NULL)
    {
        for (i = 0;  i < n0; i++)
        {
            for (j=0; j < n1; j++)
            {
                ijk0 = IDX(i,j,sx->n[2]+1, n0, n1, n2);
                ijk1 = IDX(i,j,sx->n[2]  , n0, n1, n2);

                for (c = 0; c < 6; c++) {
                    s2[ijk0].ps[0][c] = s2[ijk1].ps[0][c];
                    s2[ijk0].dFull[c] = s2[ijk1].dFull[c];
                    s2[ijk0].dHalf[c] = s2[ijk1].dHalf[c];
                }
            } //end k loop
        } // end j loop
    }
}
#endif
/*===========================================================================*/

#endif // PERFCOND_H












#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"


/* __ update ghost cells from neighbors for density & velocity ______________ */
void ghost(struct sti si, struct stx sx, double *tw,
           double i0x, double i1x, double i0y, double i1y, double i0z, double i1z)
{
double *ws[27], *wr[27];
double *sw;
int a, b, c, i, j, k, t;
int ijk;
int ijk0, ijk1;
int n0, n1, n2;
int is0, is1, js0, js1, ks0, ks1;
int ir0, ir1, jr0, jr1, kr0, kr1;
int wc[27], qc;
int nn2;
MPI_Status st;


/* __ # of grid points __ */
n0 = sx.n[0]+2;
n1 = sx.n[1]+2;
n2 = sx.n[2]+2;

/* __ # of grid points on g2 __ */
nn2 = n0*n1*n2;


/* __ memory allocation for send & receive buffers : loop on all the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ not more than the g2 size in subdomain __ */
    ws[t] = (double *)malloc(nn2*sizeof(double));
    wr[t] = (double *)malloc(nn2*sizeof(double));
    }

/* __ memory allocation for dump on the edges of the domain __ */
sw = (double *)malloc(nn2*sizeof(double));


/* --------------------------------------------------------------------- */
/*                        FILL SEND BUFFERS                              */
/* --------------------------------------------------------------------- */



/* __ pack the data : nested loops on the neighbors __ */
for (a = -1; a <= +1; a++)
    {
    for (b = -1; b <= +1; b++)
        {
        for (c = -1; c <= +1; c++)
            {
            /* __ set the tag values __ */
            t = (1+c)+3*((1+b)+3*(1+a));

            switch (a)
                   {
                   case -1 : /* __ left : indexes of the data to load __ */
                             is0 = 0;
                             is1 = 1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             is0 = 0;
                             is1 = sx.n[0]+1;
                             break;

                   case +1 : /* __ right : indexes of the data to load __ */
                             is0 = sx.n[0];
                             is1 = sx.n[0]+1;
                             break;

                   /* __ no reason to get there __ */
                   default : IAMDEAD(sx.r);
                   }

            switch (b)
                   {
                   case -1 : /* __ bottom : indexes of the data to load __ */
                             js0 = 0;
                             js1 = 1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             js0 = 0;
                             js1 = sx.n[1]+1;
                             break;

                   case +1 : /* __ top : indexes of the data to load __ */
                             js0 = sx.n[1];
                             js1 = sx.n[1]+1;
                             break;

                   /* __ no reason to get there __ */
                   default : IAMDEAD(sx.r);
                   }

            switch (c)
                   {
                   case -1 : /* __ back : indexes of the data to load __ */
                             ks0 = 0;
                             ks1 = 1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             ks0 = 0;
                             ks1 = sx.n[2]+1;
                             break;

                   case +1 : /* __ front : indexes of the data to load __ */
                             ks0 = sx.n[2];
                             ks1 = sx.n[2]+1;
                             break;

                   /* __ no reason to get there __ */
                   default : IAMDEAD(sx.r);
               }

            /* __ set the counter for the packed data for send/receive __ */
            wc[t] = 0;

            /* __ nested loops on the appropriate indexes __ */
            for (i = is0; i <= is1; i++)
                {
                for (j = js0; j <= js1; j++)
                    {
                    for (k = ks0; k <= ks1; k++)
                        {
                        /* __ index on g2 __ */
                        ijk = IDX(i, j, k, n0, n1, n2);

                        /* __ fill the buffer to send __ */
                        ws[t][wc[t]++] = tw[ijk];
                        }
                    }
                }
            }
        }
    }


/* --------------------------------------------------------------------- */
/*                        NOW SEND/RECEIVE DATA                          */
/* --------------------------------------------------------------------- */


/* __ send/receive the data __ */
for (t = 0; t < 27; t++)
    {
    if (t != 13)
       {
       MPI_Sendrecv(ws[t], wc[t], MPI_DOUBLE, sx.nt[t], t,
                    wr[t], wc[t], MPI_DOUBLE, sx.nf[t], t, MPI_COMM_WORLD, &st);
       }
    }



/* --------------------------------------------------------------------- */
/*                        NOW UPDATE GHOST NODES                         */
/* --------------------------------------------------------------------- */




/* __ update tw : nested loops on the neighbors __ */
for (a = -1; a <= +1; a++)
    {
    for (b = -1; b <= +1; b++)
        {
        for (c = -1; c <= +1; c++)
            {
            /* __ set the tag value __ */
            t = (1+c)+3*((1+b)+3*(1+a));

            switch (a)
                   {
                   case -1 : /* __ left : indexes of the data to dump __ */
                             ir0 = sx.n[0];
                             ir1 = sx.n[0]+1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to dump __ */
                             ir0 = 0;
                             ir1 = sx.n[0]+1;
                             break;

                   case +1 : /* __ right : indexes of the data to dump __ */
                             ir0 = 0;
                             ir1 = 1;
                             break;

                   /* __ no reason to get there __ */
                   default : IAMDEAD(sx.r);
                   }

            switch (b)
                   {
                   case -1 : /* __ bottom : indexes of the data to dump __ */
                             jr0 = sx.n[1];
                             jr1 = sx.n[1]+1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to dump __ */
                             jr0 = 0;
                             jr1 = sx.n[1]+1;
                             break;

                   case +1 : /* __ top : indexes of the data to dump __ */
                             jr0 = 0;
                             jr1 = 1;
                             break;

                   /* __ no reason to get there __ */
                   default : IAMDEAD(sx.r);
                   }

            switch (c)
                   {
                   case -1 : /* __ back : indexes of the data to dump __ */
                             kr0 = sx.n[2];
                             kr1 = sx.n[2]+1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to dump __ */
                             kr0 = 0;
                             kr1 = sx.n[2]+1;
                             break;

                   case +1 : /* __ top : indexes of the data to dump __ */
                             kr0 = 0;
                             kr1 = 1;
                             break;

                   /* __ no reason to get there __ */
                   default : IAMDEAD(sx.r);
                   }

            if (t != 13 && sx.nf[t] != MPI_PROC_NULL)
               {
               /* __ set the counter __ */
               qc = 0;

               /* __ loop on the appropriate indexes __ */
               for (i = ir0; i <= ir1; i++)
                   {
                   for (j = jr0; j <= jr1; j++)
                       {
                       for (k = kr0; k <= kr1; k++)
                           {
                           /* __ index on g2 __ */
                           ijk = IDX(i, j, k, sx.n[0]+2, sx.n[1]+2, sx.n[2]+2);

                           /* __ fill tw from the buffer __ */
                           tw[ijk] += wr[t][qc++];
                           }
                       }
                   }
               }
            }
        }
    }




/* --------------------------------------------------------------------- */
/*              SPECIAL TREATMENT FOR 1D AND 2D AND 3D CASES             */
/* --------------------------------------------------------------------- */


/* __ need special attention if # of cell = 1 in x direction __ */
if (si.n[0] == 1)
   {
   /* __ in the x direction : appropriate nested loops __ */
   for (j = 0; j < n1; j++)
       {
       for (k = 0; k < n2; k++)
           {
           /* __ index on g2 __ */
           ijk0 = IDX(0, j, k, n0, n1, n2);
           ijk  = IDX(1, j, k, n0, n1, n2);
           ijk1 = IDX(2, j, k, n0, n1, n2);

           /* __ update the edges of the subdomain in x direction __ */
           tw[ijk0] = tw[ijk];
           tw[ijk1] = tw[ijk];
           }
       }
   }

/* __ need special attention if # of cell = 1 in y direction __ */
if (si.n[1] == 1)
   {
   /* __ in the y direction : appropriate nested loops __ */
   for (i = 0; i < n0; i++)
       {
       for (k = 0; k < n2; k++)
           {
           /* __ index on g2 __ */
           ijk0 = IDX(i, 0, k, n0, n1, n2);
           ijk  = IDX(i, 1, k, n0, n1, n2);
           ijk1 = IDX(i, 2, k, n0, n1, n2);

           /* __ update the edges of the subdomain in y direction __ */
           tw[ijk0] = tw[ijk];
           tw[ijk1] = tw[ijk];
           }
       }
   }

/* __ need special attention if # of cell = 1 in z direction __ */
if (si.n[2] == 1)
   {
   /* __ in the z direction : appropriate nested loops __ */
   for (i = 0; i < n0; i++)
       {
       for (j = 0; j < n1; j++)
           {
           /* __ index on g2 __ */
           ijk0 = IDX(i, j, 0, n0, n1, n2);
           ijk  = IDX(i, j, 1, n0, n1, n2);
           ijk1 = IDX(i, j, 2, n0, n1, n2);

           /* __ update the edges of the subdomain in z direction __ */
           tw[ijk0] = tw[ijk];
           tw[ijk1] = tw[ijk];
           }
       }
   }





/* --------------------------------------------------------------------- */
/*                    NON PERIODIC BOUNDARY CONDITION                    */
/* --------------------------------------------------------------------- */

/* __ set the border of the domain if bounded : initial dump value __ */
for (ijk = 0; ijk < nn2; ijk++)
    {
    sw[ijk] = 0.0;
    }

/* __ set the border of the domain if bounded : loop on the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ boundary condition if needed __ */
    if (sx.nt[t] == MPI_PROC_NULL)
       {
       switch (t)
              {
              case  4 : /* __ left : loop on the appropriate indexes __ */
                        for (j = 0; j < n1; j++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(0, j, k, n0, n1, n2);
                                ijk1 = IDX(1, j, k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                sw[ijk0] += tw[ijk1];
                                sw[ijk1] += tw[ijk0];
                                }
                            }
                        break ;

              case 22 : /* __ right : loop on the appropriate indexes __ */
                        for (j = 0; j < n1; j++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(sx.n[0]+1, j, k, n0, n1, n2);
                                ijk1 = IDX(sx.n[0]  , j, k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                sw[ijk0] += tw[ijk1];
                                sw[ijk1] += tw[ijk0];
                                }
                            }
                        break ;

              case 10 : /* __ bottom : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, 0, k, n0, n1, n2);
                                ijk1 = IDX(i, 1, k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                sw[ijk0] += tw[ijk1];
                                sw[ijk1] += tw[ijk0];
                                }
                            }
                        break ;

              case 16 : /* __ top : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, sx.n[1]+1, k, n0, n1, n2);
                                ijk1 = IDX(i, sx.n[1]  , k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                sw[ijk0] += tw[ijk1];
                                sw[ijk1] += tw[ijk0];
                                }
                            }
                        break ;

              case 12 : /* __ back : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (j = 0; j < n1; j++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, j, 0, n0, n1, n2);
                                ijk1 = IDX(i, j, 1, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                sw[ijk0] += tw[ijk1];
                                sw[ijk1] += tw[ijk0];
                                }
                            }
                        break ;

              case 14 : /* __ front : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (j = 0; j < n1; j++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, j, sx.n[2]+1, n0, n1, n2);
                                ijk1 = IDX(i, j, sx.n[2]  , n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                sw[ijk0] += tw[ijk1];
                                sw[ijk1] += tw[ijk0];
                                }
                            }
                        break ;
              }
       }
    }

/* __ add the dumped values at the edges of the domain __ */
for (ijk = 0; ijk < nn2; ijk++)
    {
    tw[ijk] += sw[ijk];
    }



/* --------------------------------------------------------------------- */
/*                TAKE CARE OF BOUNDARY CONDITIONS COEF                  */
/* --------------------------------------------------------------------- */



/* __ update ghost cells @ the edges of the domain : loop on the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ boundary condition if needed __ */
    if (sx.nt[t] == MPI_PROC_NULL)
       {
       switch (t)
              {
              case  4 : /* __ left : loop on the appropriate indexes __ */
                        for (j = 0; j < n1; j++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(0, j, k, n0, n1, n2);
                                ijk1 = IDX(1, j, k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                tw[ijk0]  = tw[ijk1]*i0x;
                                }
                            }
                        break ;

              case 22 : /* __ right : loop on the appropriate indexes __ */
                        for (j = 0; j < sx.n[1]+2; j++)
                            {
                            for (k = 0; k < sx.n[2]+2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(sx.n[0]+1, j, k, n0, n1, n2);
                                ijk1 = IDX(sx.n[0]  , j, k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                tw[ijk0]  = tw[ijk1]*i1x;
                                }
                            }
                        break ;

              case 10 : /* __ bottom : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, 0, k, n0, n1, n2);
                                ijk1 = IDX(i, 1, k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                tw[ijk0]  = tw[ijk1]*i0y;
                                }
                            }
                        break ;

              case 16 : /* __ top : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, sx.n[1]+1, k, n0, n1, n2);
                                ijk1 = IDX(i, sx.n[1]  , k, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                tw[ijk0]  = tw[ijk1]*i1y;
                                }
                            }
                        break ;

              case 12 : /* __ back : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (j = 0; j < n1; j++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, j, 0, n0, n1, n2);
                                ijk1 = IDX(i, j, 1, n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                tw[ijk0]  = tw[ijk1]*i0z;
                                }
                            }
                        break ;

              case 14 : /* __ front : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (j = 0; j < n1; j++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i, j, sx.n[2]+1, n0, n1, n2);
                                ijk1 = IDX(i, j, sx.n[2]  , n0, n1, n2);

                                /* __ set the ghost cell value __ */
                                tw[ijk0]  = tw[ijk1]*i1z;
                                }
                            }
                        break ;
              }
       }
    }

/* __ clean-up the pointers __ */
for (t = 0; t < 27; t++)
    {
    free(ws[t]);
    free(wr[t]);
    }

free(sw);


}


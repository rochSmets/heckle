
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"

#ifdef OUTDATED

/* _____ fill the boundaries of g2 fields from neighbors ____________________ */
void fillg2(struct sti si, struct stx sx, double *tw,
            double i0x, double i1x, double i0y, double i1y, double i0z, double i1z)
{
double *ws[27], *wr[27];
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

/* __ loop on all the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ memory allocation for send & receive __ */
    ws[t] = (double *)malloc(nn2*sizeof(double));
    wr[t] = (double *)malloc(nn2*sizeof(double));
    }

/* __ pack and send/receive the data : nested loops on all the neighbors __ */
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
                   case -1 : /* __ left : indexes of the data to load __ */
                             is0 = 1;
                             is1 = 1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             is0 = 1;
                             is1 = sx.n[0];
                             break;

                   case +1 : /* __ right : indexes of the data to load __ */
                             is0 = sx.n[0];
                             is1 = sx.n[0];
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (b)
                   {
                   case -1 : /* __ bottom : indexes of the data to load __ */
                             js0 = 1;
                             js1 = 1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             js0 = 1;
                             js1 = sx.n[1];
                             break;

                   case +1 : /* __ top : indexes of the data to load __ */
                             js0 = sx.n[1];
                             js1 = sx.n[1];
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (c)
                   {
                   case -1 : /* __ back : indexes of the data to load __ */
                             ks0 = 1;
                             ks1 = 1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             ks0 = 1;
                             ks1 = sx.n[2];
                             break;

                   case +1 : /* __ front : indexes of the data to load __ */
                             ks0 = sx.n[2];
                             ks1 = sx.n[2];
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            /* __ counter on packed ws for send/receive __ */
            wc[t] = 0;

            /* __ nested loops on the appropriate indexes __ */
            for (i = is0; i <= is1; i++)
                {
                for (j = js0; j <= js1; j++)
                    {
                    for (k = ks0; k <= ks1; k++)
                        {
                        /* __ set the index on g2 __ */
                        ijk = IDX(i, j, k, n0, n1, n2);

                        /* __ fill the buffer __ */
                        ws[t][wc[t]++] = tw[ijk];
                        }
                    }
                }
            }
        }
    }

/* __ send/receive the data __ */
for (t = 0; t < 27; t++)
    {
    if (t != 13)
       {
       MPI_Sendrecv(ws[t], wc[t], MPI_DOUBLE, sx.nt[t], t,
                    wr[t], wc[t], MPI_DOUBLE, sx.nf[t], t, MPI_COMM_WORLD, &st);
       }
    }

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
                             ir0 = sx.n[0]+1;
                             ir1 = sx.n[0]+1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to dump __ */
                             ir0 = 1;
                             ir1 = sx.n[0];
                             break;

                   case +1 : /* __ right : indexes of the data to dump __ */
                            ir0 = 0;
                            ir1 = 0;
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (b)
                   {
                   case -1 : /* __ bottom : indexes of the data to dump __ */
                             jr0 = sx.n[1]+1;
                             jr1 = sx.n[1]+1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to dump __ */
                             jr0 = 1;
                             jr1 = sx.n[1];
                             break;

                   case +1 : /* __ top : indexes of the data to dump __ */
                             jr0 = 0;
                             jr1 = 0;
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (c)
                   {
                   case -1 : /* __ back : indexes of the data to dump __ */
                             kr0 = sx.n[2]+1;
                             kr1 = sx.n[2]+1;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to dump __ */
                             kr0 = 1;
                             kr1 = sx.n[2];
                             break;

                   case +1 : /* __ front : indexes of the data to dump __ */
                             kr0 = 0;
                             kr1 = 0;
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            /* __ unpack from send/receive __ */
            if (t != 13 && sx.nf[t] != MPI_PROC_NULL)
               {
               /* __ set the counter __ */
               qc = 0;

               /* __ nested loop on the appropriate indexes __ */
               for (i = ir0; i <= ir1; i++)
                   {
                   for (j = jr0; j <= jr1; j++)
                       {
                       for (k = kr0; k <= kr1; k++)
                           {
                           /* __ set index on g2 __ */
                           ijk = IDX(i  , j  , k  , n0, n1, n2);

                           /* __ fill tw from the buffer __ */
                           tw[ijk] = wr[t][qc++];
                           }
                       }
                   }
               }
            }
        }
    }

/*___ update the ghost cells with dirichlet-neumann conditions : loop on all the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ boundary condition if needed : left __ */
    if (sx.nt[t] == MPI_PROC_NULL)
       {
       switch (t)
              {
              case  4 : /* __ left boundary condition : loop on the appropriate indexes __ */
                        for (j = 0; j < n1; j++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(0  , j  , k  , n0, n1, n2);
                                ijk1 = IDX(1  , j  , k  , n0, n1, n2);

                                /* __ set the outer cell value __ */
                                tw[ijk0] = tw[ijk1]*i0x;
                                }
                            }
                        break;

              case 22 : /* __ right boundary condition : loop on the appropriate indexes __ */
                        for (j = 0; j < n1; j++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(sx.n[0]+1, j  , k  , n0, n1, n2);
                                ijk1 = IDX(sx.n[0]  , j  , k  , n0, n1, n2);

                                /* __ set the outer cell value __ */
                                tw[ijk0] = tw[ijk1]*i1x;
                                }
                            }
                        break;

              case 10 : /* __ bottom boundary condition : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i  , 0  , k  , n0, n1, n2);
                                ijk1 = IDX(i  , 1  , k  , n0, n1, n2);

                                /* __ set the outer cell value __ */
                                tw[ijk0] = tw[ijk1]*i0y;
                                }
                            }
                        break;

              case 16 : /* __ top boundary condition : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (k = 0; k < n2; k++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i  , sx.n[1]+1 , k  , n0, n1, n2);
                                ijk1 = IDX(i  , sx.n[1]   , k  , n0, n1, n2);

                                /* __ set the outer cell value __ */
                                tw[ijk0] = tw[ijk1]*i1y;
                                }
                            }
                        break;

              case 12 : /* __ back boundary condition : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (j = 0; j < n1; j++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i  , j  , 0  , n0, n1, n2);
                                ijk1 = IDX(i  , j  , 1  , n0, n1, n2);

                                /* __ set the outer cell value __ */
                                tw[ijk0] = tw[ijk1]*i0z;
                                }
                            }
                        break;

              case 14 : /* __ front boundary condition : loop on the appropriate indexes __ */
                        for (i = 0; i < n0; i++)
                            {
                            for (j = 0; j < n1; j++)
                                {
                                /* __ outbound & inbound indexes __ */
                                ijk0 = IDX(i  , j  , sx.n[2]+1 , n0, n1, n2);
                                ijk1 = IDX(i  , j  , sx.n[2]   , n0, n1, n2);

                                /* __ set the outer cell value __ */
                                tw[ijk0] = tw[ijk1]*i1z;
                                }
                            }
                        break;

              default : break ;
              }
       }
    }

/* __ loop on the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ clean-up the pointers __ */
    free(ws[t]);
    free(wr[t]);
    }

}
#endif // OUTDATED




/* _____ fill the boundaries of g4 fields from neighbors ____________________ */
void fillg4(struct sti si, struct stx sx, double *tw,
            double i0x, double i1x, double i0y, double i1y, double i0z, double i1z)
{
double *ws[27], *wr[27];
int a, b, c, i, j, k, t;
int ijk;
int ijk0, ijk1;
int p0, p1, p2;
int is0, is1, js0, js1, ks0, ks1;
int ir0, ir1, jr0, jr1, kr0, kr1;
int wc[27], qc;
int nn4;
MPI_Status st;


/* __ # of grid points __ */
p0 = sx.n[0]+4;
p1 = sx.n[1]+4;
p2 = sx.n[2]+4;


/* __ # of grid points on g2 __ */
nn4 = p0*p1*p2;

/* __ loop on all the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ memory allocation for send & receive __ */
    ws[t] = (double *)malloc(nn4*sizeof(double));
    wr[t] = (double *)malloc(nn4*sizeof(double));
    }

/* __ pack and send/receive the data : nested loops on the neighbors __ */
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
                   case -1 : /* __ left : indexes of the data to send __ */
                             is0 = 3;
                             is1 = 3;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to send __ */
                             is0 = 1;
                             is1 = sx.n[0]+2;
                             break;

                   case +1 : /* __ right : indexes of the data to send __ */
                             is0 = sx.n[0];
                             is1 = sx.n[0];
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (b)
                   {
                   case -1 : /* __ bottom : indexes of the data to send __ */
                             js0 = 3;
                             js1 = 3;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to send __ */
                             js0 = 1;
                             js1 = sx.n[1]+2;
                             break;

                   case +1 : /* __ top : indexes of the data to send __ */
                             js0 = sx.n[1];
                             js1 = sx.n[1];
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (c)
                   {
                   case -1 : /* __ back : indexes of the data to send __ */
                             ks0 = 3;
                             ks1 = 3;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to send __ */
                             ks0 = 1;
                             ks1 = sx.n[2]+2;
                             break;

                   case +1 : /* __ front : indexes of the data to send __ */
                             ks0 = sx.n[2];
                             ks1 = sx.n[2];
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            /* __ set the counter to pack ws for send/receive __ */
            wc[t] = 0;

            /* __ nested loops on the appropriate indexes __ */
            for (i = is0; i <= is1; i++)
                {
                for (j = js0; j <= js1; j++)
                    {
                    for (k = ks0; k <= ks1; k++)
                        {
                        /* __ index on g4 __ */
                        ijk = IDX(i, j, k, p0, p1, p2);

                        /* __ fill the buffer to send __ */
                        ws[t][wc[t]++] = tw[ijk];
                        }
                    }
                }

            /* __ send/receive the data __ */
            if (t != 13)
               {
               MPI_Sendrecv(ws[t], wc[t], MPI_DOUBLE, sx.nt[t], t,
                            wr[t], wc[t], MPI_DOUBLE, sx.nf[t], t,
                            MPI_COMM_WORLD, &st);
               }
            }
        }
    }

/* __ update uw : nested loops on the neighbors __ */
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
                   case -1 : /* __ left : indexes of the data to load __ */
                             ir0 = sx.n[0]+3;
                             ir1 = sx.n[0]+3;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             ir0 = 1;
                             ir1 = sx.n[0]+2;
                             break;

                   case +1 : /* __ right : indexes of the data to load __ */
                             ir0 = 0;
                             ir1 = 0;
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (b)
                   {
                   case -1 : /* __ bottom : indexes of the data to load __ */
                             jr0 = sx.n[1]+3;
                             jr1 = sx.n[1]+3;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             jr0 = 1;
                             jr1 = sx.n[1]+2;
                             break;

                   case +1 : /* __ top : indexes of the data to load __ */
                             jr0 = 0;
                             jr1 = 0;
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            switch (c)
                   {
                   case -1 : /* __ back : indexes of the data to load __ */
                             kr0 = sx.n[2]+3;
                             kr1 = sx.n[2]+3;
                             break;

                   case  0 : /* __ inside the subdomain : indexes of the data to load __ */
                             kr0 = 1;
                             kr1 = sx.n[2]+2;
                             break;

                   case +1 : /* __ front : indexes of the data to load __ */
                             kr0 = 0;
                             kr1 = 0;
                             break;

                   /* __ no reason to get there __ */
                   default : printf("no reason to get there !\n");
                   }

            /* __ set the counter to unpack from send/receive __ */
            qc = 0;

            /* __ if not on the subdomain boundary __ */
            if (t != 13 && sx.nf[t] != MPI_PROC_NULL)
               {
               /* __ loop on the appropriate indexes __ */
               for (i = ir0; i <= ir1; i++)
                   {
                   for (j = jr0; j <= jr1; j++)
                       {
                       for (k = kr0; k <= kr1; k++)
                           {
                           /* __ index on g4 __ */
                           ijk = IDX(i, j, k, p0, p1, p2);

                           /* __ fill uw from the buffer __ */
                           tw[ijk] = wr[t][qc++];
                           }
                       }
                   }
               }
            }
        }
    }

/* __ duplicate uw if on the border of the domain : loop on the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ boundary condition if needed __ */
    if (sx.nt[t] == MPI_PROC_NULL)
       {
       switch (t)
              {
              case  4 : /* __ left : loop on the dividing plan __ */
                        for (j = 0; j < p1; j++)
                            {
                            for (k = 0; k < p2; k++)
                                {
                                /* __ set outbounds & inbound indexes __ */
                                ijk0 = IDX(0, j, k, p0, p1, p2);
                                ijk1 = IDX(1, j, k, p0, p1, p2);

                                /* __ dump the outer domain values __ */
                                tw[ijk0] = tw[ijk1]*i0x;
                                }
                            }
                        break;

              case 22 : /* __ right : loop on the dividing plan __ */
                        for (j = 0; j < p1; j++)
                            {
                            for (k = 0; k < p2; k++)
                                {
                                /* __ set outbounds & inbound indexes __ */
                                ijk0 = IDX(sx.n[0]+3, j, k, p0, p1, p2);
                                ijk1 = IDX(sx.n[0]+2, j, k, p0, p1, p2);

                                /* __ dump the outer domain values __ */
                                tw[ijk0] = tw[ijk1]*i1x;
                                }
                            }
                        break;

              case 10 : /* __ bottom : loop on the dividing plan __ */
                        for (i = 0; i < p0; i++)
                            {
                            for (k = 0; k < p2; k++)
                                {
                                /* __ set outbounds & inbound indexes __ */
                                ijk0 = IDX(i, 0, k, p0, p1, p2);
                                ijk1 = IDX(i, 1, k, p0, p1, p2);

                                /* __ dump the outer domain values __ */
                                tw[ijk0] = tw[ijk1]*i0y;
                                }
                            }
                        break;

              case 16 : /* __ top : loop on the dividing plan __ */
                        for (i = 0; i < p0; i++)
                            {
                            for (k = 0; k < p2; k++)
                                {
                                /* __ set outbounds & inbound indexes __ */
                                ijk0 = IDX(i, sx.n[1]+3, k, p0, p1, p2);
                                ijk1 = IDX(i, sx.n[1]+2, k, p0, p1, p2);

                                /* __ dump the outer domain values __ */
                                tw[ijk0] = tw[ijk1]*i1y;
                                }
                            }
                        break;

               case 12 :/* __ back : loop on the dividing plan __ */
                        for (i = 0; i < p0; i++)
                            {
                            for (j = 0; j < p1; j++)
                                {
                                /* __ set outbounds & inbound indexes __ */
                                ijk0 = IDX(i, j, 0, p0, p1, p2);
                                ijk1 = IDX(i, j, 1, p0, p1, p2);

                                /* __ dump the outer domain values __ */
                                tw[ijk0] = tw[ijk1]*i0z;
                                }
                            }
                        break;

              case 14 : /* __ front : loop on the dividing plan __ */
                        for (i = 0; i < p0; i++)
                            {
                            for (j = 0; j < p1; j++)
                                {
                                /* __ set outbounds & inbound indexes __ */
                                ijk0 = IDX(i, j, sx.n[2]+3, p0, p1, p2);
                                ijk1 = IDX(i, j, sx.n[2]+2, p0, p1, p2);

                                /* __ dump the outer domain values __ */
                                tw[ijk0] = tw[ijk1]*i1z;
                                }
                            }
                        break;

              default : break;

              }
       }
    }

/* __ loop on the neighbors __ */
for (t = 0; t < 27; t++)
    {
    /* __ clean-up the pointers __ */
    free(ws[t]);
    free(wr[t]);
    }
}


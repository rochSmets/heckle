
/*

  This code tests the module 'GhostG2Mom'.
  It fills a grid g2 with 1 even on ghost nodes.
  Then it checks whether there is :

  - 2 on faces
  - 4 on edges
  - 8 on corners

*/


#define PASSED 0
#define FAILED 1


#include <structures.h>
#include <ghostg4mom.h>
#include <stdlib.h>
#include <stdio.h>
#include <subdomains.h>





/*---------------------------------------------------------------------------
  initsti()
  ---------------------------------------------------------------------------
  AIM : initialize a structure STI. The most important parameters here for
  this test are the number of grid points.
 ---------------------------------------------------------------------------*/
void initsti(STI *si, int nx, int ny, int nz)
{


    // these are the important params.
    si->n[0]   = nx;
    si->n[1]   = ny;
    si->n[2]   = nz;


    // a faire varier
    si->bc[0]  = 1.;
    si->bc[1]  = 1.;
    si->bc[2]  = 1.;
    //////


    // these are initialized but not used...
    si->l[0]   = 10.;
    si->l[1]   = 10.;
    si->l[2]   = 10.;
    si->rst    = 0;
    si->ts = 1e-3;
    si->nt = 100000;
    si->tf = 100;
    si->tp = 100;
    si->tt = 100;
    si->nm = 10000000;
    si->feed = 0;
    si->coll = 0;
    si->me   = 0.02;
    si->kapp = 0.;
    si->rsty = 1e-4;
    si->tr = 10000;
    si->ns[0] = 0;
    si->ns[1] = 1000000;
    si->ns[2] = 0;
    si->ms[0] = 0;
    si->ms[1] = 1;
    si->qs[0] = -1;
    si->qs[1] =  1;
    si->ws[0] = 0.00;
    si->ws[1] = 0.01;
    si->dl[0] = si->l[0] / si->n[0];
    si->dl[1] = si->l[1] / si->n[1];
    si->dl[2] = si->l[2] / si->n[2];
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  test_ghostg4mominit()
  ---------------------------------------------------------------------------
  AIM : this routine tests the initialization of the ghost module.
 ---------------------------------------------------------------------------*/
int test_ghostg4mominit(void)
{
    int res = FAILED;
    STI si;
    STX sx;
    GhostG4Mom *gf;

    MPI_Comm_rank(MPI_COMM_WORLD, &sx.r);
    MPI_Comm_size(MPI_COMM_WORLD, &sx.s);

    // initialize the STI parameter structure
    initsti(&si,100,100,100);

    // perform domain decomposition
    subdomains(si, &sx);


    // now initialize ghost points


    gf = GhostG4MomInit(&sx);
    if (gf == NULL)
    {
        return FAILED;
    }
    GhostG4MomDelete(gf);


    res = PASSED;

    return res;
}
/*===========================================================================*/





// mon nico... j'ai change les defines pour qu'ils utilisent nxg2... et eviter les warnings

#define dens(iw, jw, kw) s4[(kw)+(nzg4)*((jw)+(nyg4)*((iw)+(nxg4)*(0)))].n

#define vx(iw, jw, kw) s4[(kw)+(nzg4)*((jw)+(nyg4)*((iw)+(nxg4)*(0)))].v[0]
#define vy(iw, jw, kw) s4[(kw)+(nzg4)*((jw)+(nyg4)*((iw)+(nxg4)*(0)))].v[1]
#define vz(iw, jw, kw) s4[(kw)+(nzg4)*((jw)+(nyg4)*((iw)+(nxg4)*(0)))].v[2]



//#define dens(i,j,k) s4[k + (j)*nzg4 + (i)*(nzg4*nyg4)].n
//
//#define vx(i,j,k)   s4[k + (j)*nzg4 + (i)*(nzg4*nyg4)].v[0]
//#define vy(i,j,k)   s4[k + (j)*nzg4 + (i)*(nzg4*nyg4)].v[1]
//#define vz(i,j,k)   s4[k + (j)*nzg4 + (i)*(nzg4*nyg4)].v[2]




/*---------------------------------------------------------------------------
  fills4()
  ---------------------------------------------------------------------------
  AIM : fills the grid g4.
 ---------------------------------------------------------------------------*/
void fills4(STX *sx, ST4 *s4)
{
    int i,j,k;
    int nxg4, nyg4, nzg4;
    nxg4 = sx->n[0]+4;
    nyg4 = sx->n[1]+4;
    nzg4 = sx->n[2]+4;


    // fill the whole thing with -12
    for (i=0; i < sx->n[0]+4; i++)
    {
        for (j=0; j < sx->n[1]+4; j++)
        {
            for (k=0; k < sx->n[2]+4; k++)
            {
                dens(i,j,k) = -12;
                vx(i,j,k)   = -12;
                vy(i,j,k)   = -12;
                vz(i,j,k)   = -12;
            }
        }
    }

    // fill inside points with rank+1
    for (i=1; i < sx->n[0]+3; i++)
    {
        for (j=1; j < sx->n[1]+3; j++)
        {
            for (k=1; k < sx->n[2]+3; k++)
            {
                dens(i,j,k) = sx->r+1;
                vx(i,j,k)   = sx->r+1;
                vy(i,j,k)   = sx->r+1;
                vz(i,j,k)   = sx->r+1;
            }
        }
    }
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  checkghosts()
  ---------------------------------------------------------------------------
  AIM : routine called after SendRecv() that checks if the ghost nodes have the
  value 'rank+1' with 'rank' the rank of the neighbor proc that sent the data.
 ---------------------------------------------------------------------------*/
int checkghosts(const STX * const sx,
                const ST4 * const s4)
{
    int a,b,c;
    int i,j,k,i0,i1,j0,j1,k0,k1;
    int t;
    int nxg4, nyg4, nzg4;
    nxg4 = sx->n[0]+4;
    nyg4 = sx->n[1]+4;
    nzg4 = sx->n[2]+4;
    int ret = FAILED;


    // a, b and c are the direction in which data is sent
    // e.g. (-1,0,0) means processor has sent to LEFT and therefore
    // we receive in nx+1
    for (a=-1; a <=1; a++)
    {
        for (b = -1; b <= 1; b++)
        {
            for (c=-1; c <= 1; c++)
            {
                t = (1+c)+3*((1+b)+3*(1+a));


                switch (a)
                {
                    case -1:
                        i0 = sx->n[0]+3;
                        i1 = sx->n[0]+3;
                    break;

                    case 0:
                        i0 = 1;
                        i1 = sx->n[0]+2;
                    break;

                    case 1:
                        i0 = 0;
                        i1 = 0;
                    break;
                }




                switch(b)
                {
                    case -1:
                        j0 = sx->n[1]+3;
                        j1 = sx->n[1]+3;
                    break;

                    case 0:
                        j0 = 1;
                        j1 = sx->n[1]+2;
                    break;

                    case 1:
                        j0 = 0;
                        j1 = 0;
                    break;
                }




                switch(c)
                {
                    case -1:
                        k0 = sx->n[2]+3;
                        k1 = sx->n[2]+3;
                    break;

                    case 0:
                        k0 = 1;
                        k1 = sx->n[2];
                    break;

                    case 1:
                        k0 = 0;
                        k1 = 0;
                    break;
                }


#if 1
                printf("MERDE\n");
                for (i=i0; i <= i1; i++)
                {
                    for (j=j0; j<= j1; j++)
                    {
                        for (k=k0; k <= k1; k++)
                        {
                              if (dens(i,j,k) != sx->nf[t]+1)
                              {
                                  printf("Error N %f != %d %d %d %d | %d %d %d\n",
                                         dens(i,j,k),sx->nf[t]+1,i,j,k, a, b, c);
                                  return ret;
                              }


                              if (vx(i,j,k) != sx->nf[t]+1)
                              {
                                  printf("Error Vx %f != %d %d %d %d | %d %d %d\n",
                                         vx(i,j,k),sx->nf[t]+1,i,j,k, a, b, c);
                                  return ret;
                              }

                              if (vy(i,j,k) != sx->nf[t]+1)
                              {
                                  printf("Error Vy %f != %d %d %d %d | %d %d %d\n",
                                         vy(i,j,k),sx->nf[t]+1,i,j,k, a, b, c);
                                  return ret;
                              }

                              if (vz(i,j,k) != sx->nf[t]+1)
                              {
                                  printf("Error Vz %f != %d %d %d %d | %d %d %d\n",
                                         vz(i,j,k),sx->nf[t]+1,i,j,k, a, b, c);
                                  return ret;
                              }
                        }// end k
                    }// end j
                }//end i
#endif

            }// end c
        } // end b
    } //end a



    return PASSED;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  test_ghostg4momsendrecv()
  ---------------------------------------------------------------------------
  AIM : routine that creates a GhostG4Mom object, fills the g4 grid,
  exchange ghosts points and calls checkghosts.
 ---------------------------------------------------------------------------*/
int test_ghostg4momsendrecv(const STI* const si)
{
    int res = FAILED;
    STX sx;
    GhostG4Mom *gf;
    ST4 *s4=NULL;
    int nxg4, nyg4, nzg4;

    MPI_Comm_rank(MPI_COMM_WORLD, &sx.r);
    MPI_Comm_size(MPI_COMM_WORLD, &sx.s);

    subdomains(*si, &sx);

    if (sx.r == 0)
    {
        printf("X direction : %d\n", sx.d[0]);
        printf("Y direction : %d\n", sx.d[1]);
        printf("Z direction : %d\n", sx.d[2]);
    }

    nxg4 = sx.n[0]+4;
    nyg4 = sx.n[1]+4;
    nzg4 = sx.n[2]+4;


    s4 = malloc(nxg4*nyg4*nzg4 * sizeof*s4);
    if (s4 == NULL)
        printf("S4 NULL\n");


    gf = GhostG4MomInit(&sx);

    fills4(&sx, s4);


    GhostG4MomSendRecv(gf, &sx, s4);


    res = checkghosts(&sx, s4);

    free(s4);

    return res;
}
/*===========================================================================*/






#define NDOMAINS 10




int min(int results[][2], int size)
{
    int m;
    m = results[0][0];

    for (int id=0; id < NDOMAINS; id++)
    {
        for (int i=0; i < size; i++)
        {
            m = (results[id][i] < m) ? results[id][i] : m;
        }
    }

    return m;
}










int main(int argc, char **argv)
{
    int results[NDOMAINS][2];
    int i;
    int ncpus, rank, r;
    int ret;
    STI si;
    int domains[NDOMAINS][3] = {{100,100,100},
                                {50,100,10},
                                {20,10,10},
                                {10,10,1000},
                                {2,10,100},
                                {100,2,2},
                                {20,4,4},
                                {20,20,10},
                                {44,10,33},
                                {10,10,10}};

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (ncpus > 10)
    {
        if (rank == 0)
        {
            printf("Too many cpus, max allowed # is 2\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }


    for (int id=0; id < NDOMAINS; id++)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        initsti(&si, domains[id][0], domains[id][1], domains[id][2]);
        results[id][0] = test_ghostg4mominit();
        results[id][1] = test_ghostg4momsendrecv(&si);


        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) printf("Domain (%04d,%04d,%04d), Checking results...\n",
                              domains[id][0], domains[id][1], domains[id][2]);
        MPI_Barrier(MPI_COMM_WORLD);
        for (r =0; r < ncpus; r++)
        {
            for (i=0; i < 2; i++)
            {
                if (r == rank)
                    printf("proc %d > test %d : %s\n",rank,i,
                           (results[id][i]==FAILED) ? "FAILED" : "PASSED");
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) printf("----------\n \n");

    }


    ret = min(results, 4);

    MPI_Finalize();

    return ret;
}









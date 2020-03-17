
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
#include <ghostg2mom.h>
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
    si->bc[0]  = 0.;
    si->bc[1]  = 0.;
    si->bc[2]  = 0.;
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
  test_ghostg2mominit()
  ---------------------------------------------------------------------------
  AIM : this routine tests the initialization of the ghost module.
 ---------------------------------------------------------------------------*/
int test_ghostg2mominit(void)
{
    int res = FAILED;
    STI si;
    STX sx;
    GhostG2Mom *gf;

    MPI_Comm_rank(MPI_COMM_WORLD, &sx.r);
    MPI_Comm_size(MPI_COMM_WORLD, &sx.s);

    // initialize the STI parameter structure
    initsti(&si,100,100,100);

    // perform domain decomposition
    subdomains(si, &sx);


    // now initialize ghost points


    gf = GhostG2MomInit(&sx, si.n);
    if (gf == NULL)
    {
        return FAILED;
    }
    GhostG2MomDelete(gf);


    res = PASSED;

    return res;
}
/*===========================================================================*/





// mon nico... j'ai change les defines pour qu'ils utilisent nxg2... et eviter les warnings

#define os(iw, jw, kw, isp) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].os[isp]

#define vx(iw, jw, kw, isp) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].vs[isp][0]
#define vy(iw, jw, kw, isp) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].vs[isp][1]
#define vz(iw, jw, kw, isp) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].vs[isp][2]


//#define os(i,j,k,ispe) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].os[ispe]
//
//#define vx(i,j,k,ispe) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].vs[ispe][0]
//#define vy(i,j,k,ispe) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].vs[ispe][1]
//#define vz(i,j,k,ispe) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].vs[ispe][2]




/*---------------------------------------------------------------------------
  fills2()
  ---------------------------------------------------------------------------
  AIM : fills the grid g2.
 ---------------------------------------------------------------------------*/
void fills2(STX *sx, ST2 *s2)
{
    int i,j,k;
    int nxg2, nyg2, nzg2;
    nxg2 = sx->n[0]+2;
    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;


    for (i=0; i < sx->n[0]+2; i++)
    {
        for (j=0; j < sx->n[1]+2; j++)
        {
            for (k=0; k < sx->n[2]+2; k++)
            {
                os(i,j,k,1) = 1;
                vx(i,j,k,1) = 1;
                vy(i,j,k,1) = 1;
                vz(i,j,k,1) = 1;
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
                const ST2 * const s2)
{
    int i,j,k;
    int nxg2, nyg2, nzg2;
    nxg2 = sx->n[0]+2;
    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;
    int ret = FAILED;
    int ii[4], jj[4], kk[4];

    // first check faces, where values should be 2
    // then check edges, where values should be 4
    // then corners where the values should be 8


    // ok start with faces, there are 4 'i = cst faces'
    ii[0] = 0;
    ii[1] = 1;
    ii[2] = sx->n[0];
    ii[3] = sx->n[0]+1;
    for (int iii=0; iii < 4; iii++)
    {
        i = ii[iii];
        for (j=2; j < sx->n[1]; j++)
        {
            for (k=2; k < sx->n[2]; k++)
            {
                if (os(i,j,k,1) != 2)
                {
                    printf("Face OS Error : %d %f %d %d %d\n",
                           sx->r,os(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vx(i,j,k,1) != 2)
                {
                    printf("Face VX Error : %d %f %d %d %d\n",
                           sx->r,vx(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vy(i,j,k,1) != 2)
                {
                    printf("Face Error : %d %f %d %d %d\n",
                           sx->r,vy(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vz(i,j,k,1) != 2)
                {
                    printf("Face Error : %d %f %d %d %d\n",
                           sx->r,vz(i,j,k,1), i,j,k);
                    return ret;
                }
            }
        }
    }



    // now the 4 'J=cst' faces
    jj[0] = 0;
    jj[1] = 1;
    jj[2] = sx->n[1];
    jj[3] = sx->n[1]+1;
    for (int jjj=0; jjj < 4; jjj++)
    {
        j = jj[jjj];
        for (i=2; i < sx->n[0]; i++)
        {
            for (k=2; k < sx->n[2]; k++)
            {
                if (os(i,j,k,1) != 2)
                {
                    printf("Face OS Error : %d %f %d %d %d\n",
                           sx->r,os(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vx(i,j,k,1) != 2)
                {
                    printf("Face VX Error : %d %f %d %d %d\n",
                           sx->r,vx(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vy(i,j,k,1) != 2)
                {
                    printf("Face Error : %d %f %d %d %d\n",
                           sx->r,vy(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vz(i,j,k,1) != 2)
                {
                    printf("Face Error : %d %f %d %d %d\n",
                           sx->r,vz(i,j,k,1), i,j,k);
                    return ret;
                }
            }
        }
    }





    // now the 4 'k=cst' faces
    kk[0] = 0;
    kk[1] = 1;
    kk[2] = sx->n[2];
    kk[3] = sx->n[2]+1;
    for (int kkk=0; kkk < 4; kkk++)
    {
        k = kk[kkk];
        for (i=2; i < sx->n[0]; i++)
        {
            for (j=2; j < sx->n[1]; j++)
            {
                if (os(i,j,k,1) != 2)
                {
                    printf("Face OS Error : %d %f %d %d %d\n",
                           sx->r,os(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vx(i,j,k,1) != 2)
                {
                    printf("Face VX Error : %d %f %d %d %d\n",
                           sx->r,vx(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vy(i,j,k,1) != 2)
                {
                    printf("Face Error : %d %f %d %d %d\n",
                           sx->r,vy(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vz(i,j,k,1) != 2)
                {
                    printf("Face Error : %d %f %d %d %d\n",
                           sx->r,vz(i,j,k,1), i,j,k);
                    return ret;
                }
            }
        }
    }



    // ok now edges

    // edges along the J direction
    ii[0] = 0;
    ii[1] = 1;
    ii[2] = sx->n[0];
    ii[3] = sx->n[0]+1;
    kk[0] = 0;
    kk[1] = 1;
    kk[2] = sx->n[2];
    kk[3] = sx->n[2]+1;

    for (j = 2; j < sx->n[1]; j++)
    {
        for (int iii=0; iii < 4; iii++)
        {
            i = ii[iii];
            for (int kkk=0; kkk < 4; kkk++)
            {
                k = kk[kkk];
                if (os(i,j,k,1) != 4)
                {
                    printf("Edge OS Error : %d %f %d %d %d\n",
                           sx->r,os(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vx(i,j,k,1) != 4)
                {
                    printf("Edge VX Error : %d %f %d %d %d\n",
                           sx->r,vx(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vy(i,j,k,1) != 4)
                {
                    printf("Edge Error : %d %f %d %d %d\n",
                           sx->r,vy(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vz(i,j,k,1) != 4)
                {
                    printf("Edge Error : %d %f %d %d %d\n",
                           sx->r,vz(i,j,k,1), i,j,k);
                    return ret;
                }
            }
        }
    }




    //i edges
    jj[0] = 0;
    jj[1] = 1;
    jj[2] = sx->n[1];
    jj[3] = sx->n[1]+1;
    kk[0] = 0;
    kk[1] = 1;
    kk[2] = sx->n[2];
    kk[3] = sx->n[2]+1;

    for (i = 2; i < sx->n[0]; i++)
    {
        for (int jjj=0; jjj < 4; jjj++)
        {
            j = jj[jjj];
            for (int kkk=0; kkk < 4; kkk++)
            {
                k = kk[kkk];
                if (os(i,j,k,1) != 4)
                {
                    printf("Edge OS Error : %d %f %d %d %d\n",
                           sx->r,os(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vx(i,j,k,1) != 4)
                {
                    printf("Edge VX Error : %d %f %d %d %d\n",
                           sx->r,vx(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vy(i,j,k,1) != 4)
                {
                    printf("Edge Error : %d %f %d %d %d\n",
                           sx->r,vy(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vz(i,j,k,1) != 4)
                {
                    printf("Edge Error : %d %f %d %d %d\n",
                           sx->r,vz(i,j,k,1), i,j,k);
                    return ret;
                }
            }
        }
    }





    // corners, there are 8 of them
    ii[0] = 0;
    ii[1] = 1;
    ii[2] = sx->n[0];
    ii[3] = sx->n[0]+1;

    jj[0] = 0;
    jj[1] = 1;
    jj[2] = sx->n[1];
    jj[3] = sx->n[1]+1;

    kk[0] = 0;
    kk[1] = 1;
    kk[2] = sx->n[2];
    kk[3] = sx->n[2]+1;


    for (int iii = 0; iii < 4; iii++)
    {
        for (int jjj=0; jjj < 4; jjj++)
        {
            for (int kkk=0; kkk < 4; kkk++)
            {
                if (os(i,j,k,1) != 8)
                {
                    printf("Corner OS Error : %d %f %d %d %d\n",
                           sx->r,os(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vx(i,j,k,1) != 8)
                {
                    printf("Corner VX Error : %d %f %d %d %d\n",
                           sx->r,vx(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vy(i,j,k,1) != 8)
                {
                    printf("Corner Error : %d %f %d %d %d\n",
                           sx->r,vy(i,j,k,1), i,j,k);
                    return ret;
                }
                if (vz(i,j,k,1) != 8)
                {
                    printf("Corner Error : %d %f %d %d %d\n",
                           sx->r,vz(i,j,k,1), i,j,k);
                    return ret;
                }
            }
        }
    }



    return PASSED;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  test_ghostg2momsendrecv()
  ---------------------------------------------------------------------------
  AIM : routine that creates a GhostG2Mom object, fills the g2 grid,
  exchange ghosts points and calls checkghosts.
 ---------------------------------------------------------------------------*/
int test_ghostg2momsendrecv(const STI* const si)
{
    int res = FAILED;
    STX sx;
    GhostG2Mom *gf;
    ST2 *s2=NULL;
    int nxg2, nyg2, nzg2;

    MPI_Comm_rank(MPI_COMM_WORLD, &sx.r);
    MPI_Comm_size(MPI_COMM_WORLD, &sx.s);

    subdomains(*si, &sx);

    if (sx.r == 0)
    {
        printf("X direction : %d\n", sx.d[0]);
        printf("Y direction : %d\n", sx.d[1]);
        printf("Z direction : %d\n", sx.d[2]);
    }

    nxg2 = sx.n[0]+2;
    nyg2 = sx.n[1]+2;
    nzg2 = sx.n[2]+2;


    s2 = malloc(nxg2*nyg2*nzg2 * sizeof*s2);
    if (s2 == NULL)
        printf("S2 NULL\n");


    gf = GhostG2MomInit(&sx, si->n);

    fills2(&sx, s2);


    GhostG2MomSendRecv(gf, &sx, s2);


    res = checkghosts(&sx, s2);

    free(s2);

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
        results[id][0] = test_ghostg2mominit();
        results[id][1] = test_ghostg2momsendrecv(&si);


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








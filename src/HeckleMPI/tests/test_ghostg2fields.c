
/*

  This code tests the module 'GhostG2Fields'.
  It fills a grid g2 with some value inside that depends on the processor
  and -12 on ghost nodes. Then it calls sendrecv and checks wether ghost
  nodes have the right value.


  IMPORTANT : the code has just been tested with fully periodic domains
  need to add a IF proc_null for non periodic BCs.


  This code uses 'subdomain' to decompose the domain.


  returns PASSED (0) if ok, FAILED (1) if there is an error.

*/


#define PASSED 0
#define FAILED 1


#include <structures.h>
#include <ghostg2fields.h>
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
  test_ghostg2fieldsinit()
  ---------------------------------------------------------------------------
  AIM : this routine tests the initialization of the ghost module.
 ---------------------------------------------------------------------------*/
int test_ghostg2fieldsinit(void)
{
    int res = FAILED;
    STI si;
    STX sx;
    GhostG2Fields *gf;

    MPI_Comm_rank(MPI_COMM_WORLD, &sx.r);
    MPI_Comm_size(MPI_COMM_WORLD, &sx.s);

    // initialize the STI parameter structure
    initsti(&si,100,100,100);

    // perform domain decomposition
    subdomains(si, &sx);


    // now initialize ghost points


    gf = GhostG2FieldsInit(&sx, GHOSTFIELD_E);
    if (gf == NULL)
    {
        return FAILED;
    }
    GhostG2FieldsDelete(gf);



    gf = GhostG2FieldsInit(&sx, GHOSTFIELD_F);
    if (gf == NULL)
    {
        return FAILED;
    }
    GhostG2FieldsDelete(gf);



    gf = GhostG2FieldsInit(&sx, GHOSTFIELD_J);
    if (gf == NULL)
    {
        return FAILED;
    }


    res = PASSED;



    return res;
}
/*===========================================================================*/




// mon nico... j'ai change les defines pour qu'ils utilisent nxg2... et eviter les warnings

#define ex(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].e[0]
#define ey(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].e[1]
#define ez(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].e[2]

#define fx(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].f[0]
#define fy(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].f[1]
#define fz(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].f[2]

#define jx(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].j[0]
#define jy(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].j[1]
#define jz(iw, jw, kw) s2[(kw)+(nzg2)*((jw)+(nyg2)*((iw)+(nxg2)*(0)))].j[2]


//#define ex(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[0]
//#define ey(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[1]
//#define ez(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].e[2]
//
//#define fx(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].f[0]
//#define fy(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].f[1]
//#define fz(i,j,k) s2[k + (j)*nzg2 + (i)*(nzg2*nyg2)].f[2]
//
//#define jx(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].j[0]
//#define jy(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].j[1]
//#define jz(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].j[2]




/*---------------------------------------------------------------------------
  fills2()
  ---------------------------------------------------------------------------
  AIM : fills the grid g2. Points inside the domain are filled with 'rank+1'
  while ghost points are filled with -12
 ---------------------------------------------------------------------------*/
void fills2(STX *sx, ST2 *s2, int qty)
{
    int i,j,k;
    int nxg2, nyg2, nzg2;
    nxg2 = sx->n[0]+2;
    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;


    for (i=1; i < sx->n[0]+1; i++)
    {
        for (j=1; j < sx->n[1]+1; j++)
        {
            for (k=1; k < sx->n[2]+1; k++)
            {
                switch (qty)
                {
                    case GHOSTFIELD_E:
                        ex(i,j,k) = sx->r+1;
                        ey(i,j,k) = sx->r+1;
                        ez(i,j,k) = sx->r+1;
                    break;

                    case GHOSTFIELD_F:
                        fx(i,j,k) = sx->r+1;
                        fy(i,j,k) = sx->r+1;
                        fz(i,j,k) = sx->r+1;
                    break;

                    case GHOSTFIELD_J:
                        jx(i,j,k) = sx->r+1;
                        jy(i,j,k) = sx->r+1;
                        jz(i,j,k) = sx->r+1;
                    break;
                }
            }
        }
    }


    // now puts bad value on ghost points

    for (i=0; i < nxg2; i++)
    {
        for (j=0; j < nyg2; j++)
        {

            switch(qty)
            {
                case GHOSTFIELD_E:
                    ex(i,j,0) = -12;
                    ez(i,j,0) = -12;
                    ey(i,j,0) = -12;
                    ex(i,j,sx->n[2]+1) = -12;
                    ez(i,j,sx->n[2]+1) = -12;
                    ey(i,j,sx->n[2]+1) = -12;
                break;

                case GHOSTFIELD_F:
                    fx(i,j,0) = -12;
                    fz(i,j,0) = -12;
                    fy(i,j,0) = -12;
                    fx(i,j,sx->n[2]+1) = -12;
                    fz(i,j,sx->n[2]+1) = -12;
                    fy(i,j,sx->n[2]+1) = -12;
                break;


                case GHOSTFIELD_J:
                    jx(i,j,0) = -12;
                    jz(i,j,0) = -12;
                    jy(i,j,0) = -12;
                    jx(i,j,sx->n[2]+1) = -12;
                    jz(i,j,sx->n[2]+1) = -12;
                    jy(i,j,sx->n[2]+1) = -12;
                break;
            }// end switch qty
        }
    }


    for (i=0; i < nxg2; i++)
    {
        for (k=0; k < nzg2; k++)
        {
            switch (qty)
            {
                case GHOSTFIELD_E:
                    ex(i,0,k) = -12;
                    ez(i,0,k) = -12;
                    ey(i,0,k) = -12;
                    ex(i,sx->n[1]+1,k) = -12;
                    ez(i,sx->n[1]+1,k) = -12;
                    ey(i,sx->n[1]+1,k) = -12;
                break;

                case GHOSTFIELD_F:
                    fx(i,0,k) = -12;
                    fz(i,0,k) = -12;
                    fy(i,0,k) = -12;
                    fx(i,sx->n[1]+1,k) = -12;
                    fz(i,sx->n[1]+1,k) = -12;
                    fy(i,sx->n[1]+1,k) = -12;
                break;

                case GHOSTFIELD_J:
                    jx(i,0,k) = -12;
                    jz(i,0,k) = -12;
                    jy(i,0,k) = -12;
                    jx(i,sx->n[1]+1,k) = -12;
                    jz(i,sx->n[1]+1,k) = -12;
                    jy(i,sx->n[1]+1,k) = -12;
                break;
            }
        }
    }


    for (j=0; j < nyg2; j++)
    {
        for (k=0; k < nzg2; k++)
        {
            switch(qty)
            {
                case GHOSTFIELD_E:
                    ex(0,j,k) = -12;
                    ez(0,j,k) = -12;
                    ey(0,j,k) = -12;
                    ex(sx->n[0]+1,j,k) = -12;
                    ez(sx->n[0]+1,j,k) = -12;
                    ey(sx->n[0]+1,j,k) = -12;
                break;

                case GHOSTFIELD_F:
                    fx(0,j,k) = -12;
                    fz(0,j,k) = -12;
                    fy(0,j,k) = -12;
                    fx(sx->n[0]+1,j,k) = -12;
                    fz(sx->n[0]+1,j,k) = -12;
                    fy(sx->n[0]+1,j,k) = -12;
                break;

                case GHOSTFIELD_J:
                    jx(0,j,k) = -12;
                    jz(0,j,k) = -12;
                    jy(0,j,k) = -12;
                    jx(sx->n[0]+1,j,k) = -12;
                    jz(sx->n[0]+1,j,k) = -12;
                    jy(sx->n[0]+1,j,k) = -12;
                break;
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
                const ST2 * const s2,
                int qty)
{
    int a,b,c;
    int i,j,k,i0,i1,j0,j1,k0,k1;
    int t;
    int nxg2, nyg2, nzg2;
    nxg2 = sx->n[0]+2;
    nyg2 = sx->n[1]+2;
    nzg2 = sx->n[2]+2;
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
                        i0 = sx->n[0]+1;
                        i1 = sx->n[0]+1;
                    break;

                    case 0:
                        i0 = 1;
                        i1 = sx->n[0];
                    break;

                    case 1:
                        i0 = 0;
                        i1 = 0;
                    break;
                }




                switch(b)
                {
                    case -1:
                        j0 = sx->n[1]+1;
                        j1 = sx->n[1]+1;
                    break;

                    case 0:
                        j0 = 1;
                        j1 = sx->n[1];
                    break;

                    case 1:
                        j0 = 0;
                        j1 = 0;
                    break;
                }




                switch(c)
                {
                    case -1:
                        k0 = sx->n[2]+1;
                        k1 = sx->n[2]+1;
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



                for (i=i0; i <= i1; i++)
                {
                    for (j=j0; j<= j1; j++)
                    {
                        for (k=k0; k <= k1; k++)
                        {
                            switch (qty)
                            {
                                case GHOSTFIELD_E:
                                    if (ex(i,j,k) != sx->nf[t]+1)
                                    {
                                        if (sx->r == 0)
                                            printf("Error Vx %d %d %d %d"
                                                   " %d %d %f %d %d %d\n",
                                                   sx->nf[t],t,sx->r,i,j,k,
                                              ex(i,j,k), a, b, c);
                                        return ret;
                                    }
                                    if (ey(i,j,k) != sx->nf[t]+1)
                                    {
                                        if (sx->r == 0)
                                        printf("Error Vy %d %d %d %d"
                                               " %d %d %f %d %d %d\n",sx->nf[t],t,sx->r,i,j,k,
                                              ex(i,j,k), a, b, c);
                                        return ret;
                                    }
                                    if (ez(i,j,k) != sx->nf[t]+1)
                                    {
                                        if (sx->r == 0)
                                        printf("Error Vz %d %d %d %d"
                                               "% d %d %f %d %d %d\n",sx->nf[t],t,sx->r,i,j,k,
                                              ex(i,j,k), a, b, c);
                                        return ret;
                                    }
                                break;


                            case GHOSTFIELD_F:
                                if (fx(i,j,k) != sx->nf[t]+1)
                                {
                                    if (sx->r == 0)
                                        printf("Error Vx %d %d %d %d"
                                               " %d %d %f %d %d %d\n",
                                               sx->nf[t],t,sx->r,i,j,k,
                                          fx(i,j,k), a, b, c);
                                    return ret;
                                }
                                if (fy(i,j,k) != sx->nf[t]+1)
                                {
                                    if (sx->r == 0)
                                    printf("Error Vy %d %d %d %d"
                                           " %d %d %f %d %d %d\n",sx->nf[t],t,sx->r,i,j,k,
                                          fx(i,j,k), a, b, c);
                                    return ret;
                                }
                                if (fz(i,j,k) != sx->nf[t]+1)
                                {
                                    if (sx->r == 0)
                                    printf("Error Vz %d %d %d %d"
                                           "% d %d %f %d %d %d\n",sx->nf[t],t,sx->r,i,j,k,
                                          fx(i,j,k), a, b, c);

                                    return ret;
                                }
                            break;




                            case GHOSTFIELD_J:
                                if (jx(i,j,k) != sx->nf[t]+1)
                                {
                                    if (sx->r == 0)
                                        printf("Error Vx %d %d %d %d"
                                               " %d %d %f %d %d %d\n",
                                               sx->nf[t],t,sx->r,i,j,k,
                                          jx(i,j,k), a, b, c);
                                    return ret;
                                }
                                if (jy(i,j,k) != sx->nf[t]+1)
                                {
                                    if (sx->r == 0)
                                    printf("Error Vy %d %d %d %d"
                                           " %d %d %f %d %d %d\n",sx->nf[t],t,sx->r,i,j,k,
                                          fx(i,j,k), a, b, c);
                                    return ret;
                                }
                                if (jz(i,j,k) != sx->nf[t]+1)
                                {
                                    if (sx->r == 0)
                                    printf("Error Vz %d %d %d %d"
                                           "% d %d %f %d %d %d\n",sx->nf[t],t,sx->r,i,j,k,
                                          jz(i,j,k), a, b, c);

                                    return ret;
                                }
                            break;

                            default:
                                return FAILED;


                            }
                        } // end k
                    } //  end j
                } // end i

            }// end c
        } // end b
    } //end a


    return PASSED;
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  test_ghostg2fieldssendrecv()
  ---------------------------------------------------------------------------
  AIM : routine that creates a GhostG2Fields object, fills the g2 grid,
  exchange ghosts points and calls checkghosts.
 ---------------------------------------------------------------------------*/
int test_ghostg2fieldssendrecv(const STI* const si, int qty)
{
    int res = FAILED;
    STX sx;
    GhostG2Fields *gf;
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


    gf = GhostG2FieldsInit(&sx, qty);

    fills2(&sx, s2, qty);


    GhostG2FieldsSendRecv(gf, &sx, s2);


    res = checkghosts(&sx, s2, qty);

    free(s2);

    return res;
}
/*===========================================================================*/




#define NDOMAINS 10




int min(int results[][4], int size)
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
    int results[NDOMAINS][4];
    int i;
    int ncpus, rank, r;
    int ret;
    STI si;
    int domains[NDOMAINS][3] = {{100,100,100},
                                {50,100,10},
                                {20,10,10},
                                {10,10,1000},
                                {2,10,100},
                                {100,2,1},
                                {20,2,1},
                                {20,20,10},
                                {44,10,33},
                                {203,1,1}};

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
        results[id][0] = test_ghostg2fieldsinit();
        results[id][1] = test_ghostg2fieldssendrecv(&si, GHOSTFIELD_E);
        results[id][2] = test_ghostg2fieldssendrecv(&si, GHOSTFIELD_F);
        results[id][3] = test_ghostg2fieldssendrecv(&si, GHOSTFIELD_J);


        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) printf("Domain (%04d,%04d,%04d), Checking results...\n",
                              domains[id][0], domains[id][1], domains[id][2]);
        MPI_Barrier(MPI_COMM_WORLD);
        for (r =0; r < ncpus; r++)
        {
            for (i=0; i < 4; i++)
            {
                if (r == rank)
                    printf("proc %d > test %d : %s\n",rank,i,
                           (results[id][i] == FAILED) ? "FAILED" : "PASSED");
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) printf("----------\n \n");

    }


    ret = min(results, 4);

    MPI_Finalize();

    return ret;
}








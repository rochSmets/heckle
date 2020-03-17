

/* this file is the unit test for perfcond.c
   it must check :

    - BC_perfcond_os_vs_x()
    - BC_perfcond_os_vs_y()
    - BC_perfcond_os_vs_z()


    - BC_perfcond_nv_x()
    - BC_perfcond_nv_y()
    - BC_perfcond_nv_z()

    - BC_perfcond_e_x()
    - BC_perfcond_e_y()
    - BC_perfcond_e_z()

    - BC_perfcond_f_x()
    - BC_perfcond_f_y()
    - BC_perfcond_f_z()

    - BC_perfcond_j_x()
    - BC_perfcond_j_y()
    - BC_perfcond_j_z()

    - BC_perfcond_P_x()
    - BC_perfcond_P_y()
    - BC_perfcond_P_z()


    the strategy is to intialize a domain STI
    and subdomains STX, then fill a g2 grid with some values except on
    the domain boundaries.
    Call the BC function and check afterwards whether the appropriate BC
    has been applied.

*/




#include <stdio.h>
#include <stdlib.h>


#include <structures.h>
#include <perfcond.h>

#define PASSED 0
#define FAILED 1



#define jx(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].j[0]
#define jy(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].j[1]
#define jz(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].j[2]


#define ex(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].e[0]
#define ey(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].e[1]
#define ez(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].e[2]

#define fx(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].f[0]
#define fy(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].f[1]
#define fz(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].f[2]

#define dens(ii,jj,kk) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].n
#define dens(ii,jj,kk, ispe) s2[kk + (jj)*nzg2 + (ii)*(nzg2*nyg2)].os[ispe]






int check_os(const STX* const sx, const ST2* const s2)
{
    int i,j,k;
    int nxg2, nyg2, nzg2;
    int ret = FAILED;

    nxg2 = sx->n[0];
    nyg2 = sx->n[1];
    nzg2 = sx->n[2];



}






int main(int argc, char *argv[])
{

    return 0;
}

















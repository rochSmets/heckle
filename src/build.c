
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>

#include "defines.h"
#include "structures.h"
#include "memory.h"
#include "read.h"


/* __ main __________________________________________________________________ */
int main(int ac, char **av)
{
struct sti si;
struct stf sf;
struct stq sq;
struct sto so;
struct stm *sm;
fftw_complex *bx, *BX, *by, *BY, *bz, *BZ;
fftw_complex *ax, *AX, *ay, *AY, *az, *AZ;
fftw_plan pxf, pyf, pzf;
fftw_plan pxb, pyb, pzb;
float kx, ky, k2;
float lx, ly;
int n1, nw[2], nn;
int nf, np, no;
int f, i, j, k, l, m, o, p, r, s;
int ijk, ijkw;
int ijk0;
int m0, m1, m2;
char nfile[11];
FILE *fp;


/* __ read the "heckle.txt" file to set the si structure __ */
readheckle(&si, av[1], 0);

/* __ # of grid points __ */
n1 = (si.n[0]+1)*(si.n[1]+1)*(si.n[2]+1);

/* __ # of grid points __ */
m0 = si.n[0]+1;
m1 = si.n[1]+1;
m2 = si.n[2]+1;

/* __ # of grid points in fourier space __ */
nw[0] = si.n[0];
nw[1] = si.n[1];

/* __ # of grid points __ */
nn = nw[0]*nw[1];

/* __ # of field, part. & orbits files __ */
nf = si.nt/si.tf+1;
np = si.nt/si.tp+1;
no = si.nt/si.tt+1;

/* __ read the "horbi.txt" file to set the so structure __ */
readhorbi(&so);

/* __ print some comments __ */
printf("\n");
printf("________________ heckle code _________________________\n");
printf("\n");
printf("________________ build the hf, hp & ho files _________\n");
printf("\n");
printf("# of cpu       : %8i\n", si.ncpu);
printf("# of f. files  : %8i\n", nf);
printf("# of p. files  : %8i\n", np);
printf("# of maps      : %8i\n", so.wo);
printf("# of t. dumps  : %8i\n", no);
printf("\n");
printf("________________ id. of all maps _____________________\n");
printf("\n");
for (m = 0; m < so.wo; m++)
    {
    printf("map # %8i : %8i    %8i\n", m, so.so[m], so.io[m]);
    }
printf("\n");

/* __ memory allocation for sf structure __ */
memorystf(&sf, si);

/* __ memory allocation __ */
bx = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
BX = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
by = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
BY = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
bz = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
BZ = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
ax = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
AX = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
ay = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
AY = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
az = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));
AZ = (fftw_complex *)fftw_malloc(nn*sizeof(fftw_complex));

/* __ create 3d plan __ */
pxf = fftw_plan_dft_2d(nw[0], nw[1], bx, BX, FFTW_FORWARD, FFTW_ESTIMATE);
pyf = fftw_plan_dft_2d(nw[0], nw[1], by, BY, FFTW_FORWARD, FFTW_ESTIMATE);
pzf = fftw_plan_dft_2d(nw[0], nw[1], bz, BZ, FFTW_FORWARD, FFTW_ESTIMATE);
pxb = fftw_plan_dft_2d(nw[0], nw[1], AX, ax, FFTW_BACKWARD, FFTW_ESTIMATE);
pyb = fftw_plan_dft_2d(nw[0], nw[1], AY, ay, FFTW_BACKWARD, FFTW_ESTIMATE);
pzb = fftw_plan_dft_2d(nw[0], nw[1], AZ, az, FFTW_BACKWARD, FFTW_ESTIMATE);

/* __ print some comments __ */
printf("________________ build the hf___.dat files ___________\n");
printf("\n");

/* __ read all the field files : loop on the files __ */
for (f = 0; f < nf; f++)
    {
    /* __ loop on the nodes __ */
    for (r = 0; r < si.ncpu; r++)
        {
        /* __ read the fragmented field files __ */
        readfield(si, &sf, r, f);
        }

    /* __ build the potential vector __ */
    for (i = 0; i < si.n[0]; i++)
        {
        for (j = 0; j < si.n[1]; j++)
            {
            /* __ set the index on large domain __ */
            ijk = IDX(i, j, 0, m0, m1, m2);

            /* __ set the index on large replicated domain __ */
            ijk0 = IDX(i , j , 0 , nw[0], nw[1], 1);

            bx[ijk0][0] = sf.b[0][ijk];
            bx[ijk0][1] = 0.0;

            by[ijk0][0] = sf.b[1][ijk];
            by[ijk0][1] = 0.0;

            bz[ijk0][0] = sf.b[2][ijk];
            bz[ijk0][1] = 0.0;
            }
        }

    /* __ 2d fft __ */
    fftw_execute(pxf);
    fftw_execute(pyf);
    fftw_execute(pzf);

    /* __ set the field for fft __ */
    for (i = 0; i < nw[0]; i++)
        {
        for (j = 0; j < nw[1]; j++)
            {
            /* __ set lengths of the box __ */
            lx = si.l[0];
            ly = si.l[1];

            /* __ set k values __ */
            kx = (i < nw[0]/2+1) ?         i*2.0*PI/lx :
                                   (i-nw[0])*2.0*PI/lx;
            ky = (j < nw[1]/2+1) ?         j*2.0*PI/ly :
                                   (j-nw[1])*2.0*PI/ly;

            /* __  set the indexes __ */
            ijkw = IDX(i, j, 0, nw[0], nw[1], 1);

            /* __ remove nyquist  __ */
            if ((si.n[0]%2 == 0 && i == si.n[0]) ||
                (si.n[1]%2 == 0 && j == si.n[1]))
               {
               BX[ijkw][0] = 0.0, BY[ijkw][0] = 0.0, BZ[ijkw][0] = 0.0;
               BX[ijkw][1] = 0.0, BY[ijkw][1] = 0.0, BZ[ijkw][1] = 0.0;
               }

            /* __ set k2 __ */
            k2 = kx*kx+ky*ky;

            /* __ consider the k=0 mode __ */
            AX[ijkw][0] = (k2 == 0) ? 0.0 :
                                      (-ky*BZ[ijkw][1])/k2;
            AX[ijkw][1] = (k2 == 0) ? 0.0 :
                                      (+ky*BZ[ijkw][0])/k2;

            AY[ijkw][0] = (k2 == 0) ? 0.0 :
                                      (+kx*BZ[ijkw][1])/k2;
            AY[ijkw][1] = (k2 == 0) ? 0.0 :
                                      (-kx*BZ[ijkw][0])/k2;

            AZ[ijkw][0] = (k2 == 0) ? 0.0 :
                                      (-kx*BY[ijkw][1]+ky*BX[ijkw][1])/k2;
            AZ[ijkw][1] = (k2 == 0) ? 0.0 :
                                      (+kx*BY[ijkw][0]-ky*BX[ijkw][0])/k2;
            }
        }

    /* __ 2d fft __ */
    fftw_execute(pxb);
    fftw_execute(pyb);
    fftw_execute(pzb);

    /* __ set the field from fft __ */
    for (i = 0; i < si.n[0]+1; i++)
        {
        for (j = 0; j < si.n[1]+1; j++)
            {
            for (k = 0; k < si.n[2]+1; k++)
                {
                /* __  set the indexes __ */
                ijk = IDX(i, j, k, m0, m1, m2);
                ijkw = IDX(i%nw[0], j%nw[1], 0, nw[0], nw[1], 1);

                sf.a[0][ijk] = ax[ijkw][0]/nn;
                sf.a[1][ijk] = ay[ijkw][0]/nn;
                sf.a[2][ijk] = az[ijkw][0]/nn;
                }
            }
        }

    /* __ make the file name __ */
    sprintf(nfile, "hf%04i.dat", f);

    /* __ open the file __ */
    fp = fopen(nfile, "wb");
    if (fp == NULL) printf("problem in opening file %s\n", nfile);

    /* __ write the fields in the hf file __ */
    for (l = 0; l < 3; l++)
        {
        fwrite(sf.b[l], sizeof(float), n1, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sf.a[l], sizeof(float), n1, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sf.e[l], sizeof(float), n1, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sf.j[l], sizeof(float), n1, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sf.vi[l], sizeof(float), n1, fp);
        }

    for (s = 0; s < NS+1; s++)
        {
        fwrite(sf.n[s], sizeof(float), n1, fp);
        }

    for (s = 0; s < NS+1; s++)
        {
        for (l = 0; l < 3; l++)
            {
            fwrite(sf.v[s][l], sizeof(float), n1, fp);
            }
        }

    for (s = 0; s < NS+1; s++)
        {
        for (l = 0; l < 6; l++)
            {
            fwrite(sf.p[s][l], sizeof(float), n1, fp);
            }
        }

    /* __ print some comments __ */
    printf("\r%9s      :    %3d %%", nfile, 100*f/(nf-1));
    fflush(stdout);

    /* __ close the file __ */
    fclose(fp);
    }

/* __ free the sf structure __ */
freestf(&sf);

/* __ free the vectors for fft __ */
fftw_free(bx);
fftw_free(BX);
fftw_free(by);
fftw_free(BY);
fftw_free(bz);
fftw_free(BZ);
fftw_free(ax);
fftw_free(AX);
fftw_free(ay);
fftw_free(AY);
fftw_free(az);
fftw_free(AZ);

/* __ destroy 2d plan __ */
fftw_destroy_plan(pxf);
fftw_destroy_plan(pyf);
fftw_destroy_plan(pzf);
fftw_destroy_plan(pxb);
fftw_destroy_plan(pyb);
fftw_destroy_plan(pzb);

/* __ print some comments __ */
printf("\n\n");

/* __ memory allocation for sq structure __ */
memorystq(&sq, si);

/* __ print some comments __ */
printf("________________ build the hp___.dat files ___________\n");
printf("\n");

/* __ read all the particles files : loop on the files __ */
for (p = 0; p < np; p++)
    {
    /* __ loop on the nodes __ */
    for (r = 0; r < si.ncpu; r++)
        {
        /* __ read the fragmented particles files __ */
        readpart(&si, &sq, r, p);
        }

    /* __ make the file name __ */
    sprintf(nfile, "hp%04i.dat", p);

    /* __ open the file __ */
    fp = fopen(nfile, "wb");
    if (fp == NULL) printf("problem in opening file %s\n", nfile);

    /* __ write the part. in the hp file : loop on the specie __ */
    for (s = 1; s < NS+1; s++)
        {
        /* __ write the part. in the hp files __ */
        fwrite(&si.ns[s], sizeof(int), 1, fp);

        if (si.ns[s] != 0) fwrite(sq.r[s][0], sizeof(float), si.ns[s], fp);
        if (si.ns[s] != 0) fwrite(sq.r[s][1], sizeof(float), si.ns[s], fp);
        if (si.ns[s] != 0) fwrite(sq.r[s][2], sizeof(float), si.ns[s], fp);

        if (si.ns[s] != 0) fwrite(sq.v[s][0], sizeof(float), si.ns[s], fp);
        if (si.ns[s] != 0) fwrite(sq.v[s][1], sizeof(float), si.ns[s], fp);
        if (si.ns[s] != 0) fwrite(sq.v[s][2], sizeof(float), si.ns[s], fp);
        }

    /* __ print some comments __ */
    printf("\r%9s      :    %3d %%", nfile, 100*p/(np-1));
    fflush(stdout);

    /* __ close the file __ */
    fclose(fp);
    }

/* __ print some comments __ */
printf("\n\n");

/* __ free the sq structure __ */
freestq(&sq);

/* __ print some comments __ */
printf("________________ build the ho___.dat files ___________\n");
printf("\n");

/* __ memory allocation for sm structure __ */
sm = (struct stm *) malloc(so.wo*sizeof(struct stm));

/* __ memory allocation for each map : loop on the maps __ */
for (o = 0; o < so.wo; o++)
    {
    /* __ memory allocation for each map __ */
    memorystm(sm+o, no);
    }

/* __ read all the orbits files : loop on the files __ */
for (o = 0; o < no; o++)
    {
    /* __ loop on the nodes __ */
    for (r = 0; r < si.ncpu; r++)
        {
        /* __ read the fragmented orbits files __ */
        readorbit(so, sm, r, o);
        }
    }

/* __ read all the orbits files : loop on the files __ */
for (m = 0; m < so.wo; m++)
    {
    /* __ make the file name __ */
    sprintf(nfile, "ho%04i.dat", m);

    /* __ open the file __ */
    fp = fopen(nfile, "wb");
    if (fp == NULL) printf("problem in opening file %s\n", nfile);

    /* __ write the maps in the ho files __ */
    for (l = 0; l < 3; l++)
        {
        fwrite(sm[m].r[l], sizeof(float), no, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sm[m].w[l], sizeof(float), no, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sm[m].b[l], sizeof(float), no, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sm[m].e[l], sizeof(float), no, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sm[m].j[l], sizeof(float), no, fp);
        }

    for (l = 0; l < 3; l++)
        {
        fwrite(sm[m].vi[l], sizeof(float), no, fp);
        }

    for (s = 0; s < NS+1; s++)
        {
        fwrite(sm[m].n[s], sizeof(float), no, fp);
        }

    for (s = 0; s < NS+1; s++)
        {
        for (l = 0; l < 3; l++)
            {
            fwrite(sm[m].v[s][l], sizeof(float), no, fp);
            }
        }

    for (s = 0; s < NS+1; s++)
        {
        for (l = 0; l < 6; l++)
            {
            fwrite(sm[m].p[s][l], sizeof(float), no, fp);
            }
        }

    /* __ print some comments __ */
    printf("\r%9s      :    %3d %%", nfile, 100*(m+1)/(so.wo));
    fflush(stdout);

    /* __ close the file __ */
    fclose(fp);
    }

/* __ print some comments __ */
printf("\n\n");

/* __ free the sm structure : loop on the maps __ */
for (o = 0; o < so.wo; o++)
    {
    /* __ free the sm structure __ */
    freestm(sm+o);
    }

/* __ free the sm structure __ */
free(sm);

/* __ free the so structure __ */
free(so.so);
free(so.io);
free(so.s);
free(so.m);

return 0;

}


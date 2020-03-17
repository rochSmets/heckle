
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "structures.h"
#include "memory.h"
#include "read.h"


/* __ main __________________________________________________________________ */
int main(int ac, char **av)
{
struct sti si;
struct stf sf;
int ijk;
int i, j;
char sfd[80];
char sfr[80];
char sff[80];
FILE *fp;


/* __ which run to select __ */
strcpy(sfd, "/media/roch/TOSHIBA/sherpa/shells/run/03b/");

/* __ which run to select __ */
strcpy(sfr, "hf0400.dat");

/* __ which file to load __ */
strcpy(sff, sfd);

printf("%s\n", sff);

/* __ read the "heckle.txt" file to set the si structure __ */
readheckle(&si, sff, 0);

/* __ which file to load __ */
strcpy(sff, sfd);
strcat(sff, sfr);

printf("%s\n", sff);

/* __ memory allocation for sf structure __ */
memorystf(&sf, si);

/* __ read the full particle file __ */
readfullfield(si, &sf, sff);

/* __ open the file __ */
fp = fopen("file.txt", "w");
if (fp == NULL) printf("problem in opening file file.txt\n");

fprintf(fp, "        x          y        bx        by        bz        ex        ey        ez\n");

/* __ build the potential vector __ */
for (i = 0; i < si.n[0]+1; i++)
    {
    for (j = 0; j < si.n[1]+1; j++)
        {
        /* __ set the index on large domain __ */
        ijk = IDX(i, j, 0, si.n[0]+1, si.n[1]+1, si.n[2]+1);

        fprintf(fp, "%9.4f  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", i*si.dl[0], j*si.dl[1], sf.b[0][ijk], sf.b[1][ijk], sf.b[2][ijk], sf.e[0][ijk], sf.e[1][ijk], sf.e[2][ijk]);
        }
    }

/* __ close the file __ */
fclose(fp);

/* __ free the sf structure __ */
freestf(&sf);

return 0;
}


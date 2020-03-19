
//#include <stdio.h>
//#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/utsname.h>
#include <inttypes.h>

#include "structures.h"
#include "collision.h"
#include "read.h"
#include "closeModel.h"


/* __ read the initial parameters to set sti structure ______________________ */
struct sti setup(struct stx sx, char *dir)
{
struct sti si;
//char hostname[80];
int g;
struct utsname buffer;


/* __ read the "heckle.txt" file __ */
readheckle(&si, dir, sx.s);


   errno = 0;
   if (uname(&buffer) != 0) {
      perror("uname");
      exit(EXIT_FAILURE);
   }



/* __ display parameters of the run __ */
if (sx.r == 0)
   {
   printf("\n");
   printf("________________ heckle code _________________________\n");
   printf("cpu #            :%12d\n", sx.s);

   if (si.rst == 0) {
       printf("restart          :     scratch\n");
   }

   else {
       printf("restart          :%12.4lf\n", si.time4rst);
   }

   printf("\n");
   printf("system name      : %s\n", buffer.sysname);
   printf("node name        : %s\n", buffer.nodename);
   printf("release          : %s\n", buffer.release);
   printf("version          : %s\n", buffer.version);
   printf("machine          : %s\n", buffer.machine);

   #ifdef _GNU_SOURCE
   printf("domain name      : %s\n", buffer.domainname);
   #endif

   printf("compilation date :%12s\n", __DATE__);
   printf("compilation time :%12s\n", __TIME__);
   printf("\n");
   printf("________________ simulation settings _________________\n");
   printf("# of cells (x)   :%12d\n", si.n[0]);
   printf("# of cells (y)   :%12d\n", si.n[1]);
   printf("# of cells (z)   :%12d\n", si.n[2]);
   printf("\n");
   printf("box size (x)     :%12.6lf\n", si.l[0]);
   printf("box size (y)     :%12.6lf\n", si.l[1]);
   printf("box size (z)     :%12.6lf\n", si.l[2]);
   printf("\n");
   printf("time step        :%12.6lf\n", si.ts);
   printf("# time steps     :%12d\n", si.nt);
   printf("total time       :%12.6lf\n", si.nt*si.ts);
   printf("\n");
   printf("resistivity      :%12.6lf\n", si.rsty);
   printf("hyperviscosity   :%12.6lf\n", si.hyvi);
   printf("\n");
   printf("electron mass    :%12.6lf\n", si.me);
   printf("\n");
   printf("thermal cond.    :%12.6lf\n", si.kapp);
   printf("\n");
   printf("max # part/node  :%12ld\n", si.nm);
   printf("min N value (Ohm):%12.6lf\n", si.nmin);
   printf("\n");

   switch (si.coll) {
       case 0 : /* __ no collisions __ */
           printf("collisions       :          no\n");
           break;

       case 1 : /* __ ion-ion collisions __ */
           printf("collisions       :     ion-ion\n");
           break;

       /* __ no reason to get there __ */
       default : printf("no reason to get there !\n");
   }

   switch (si.feed) {
       case 0 : /* __ const. # of particles __ */
           printf("feeded run       :          no\n");
           break;

       /* __ no reason to get there __ */
       default : printf("no reason to get there !\n");
   }

   switch (si.drive) {
       case 0 : /* __ no driver __ */
           printf("driven system    :          no\n");
           break;

       case 1 : /* __ driven system __ */
           printf("driven system    :         yes\n");
           break;

       /* __ no reason to get there __ */
       default : printf("no reason to get there !\n");
   }

   printf("\n");
   #ifdef NS
   g = NS;
   printf("number of specie :%12d\n", g);
   #else
   printf("number of specie : not defined !!!\n");
   exit(EXIT_FAILURE);
   #endif
   printf("\n");
   printf("________________ particles setting ___________________\n");
   for (g = 1; g < NS+1; g++) {
       printf("specie %1d\n", g);
       printf("... # of part    :%12zu\n" , si.ns[g]);
       printf("... mass         :%12.6lf\n", si.ms[g]);
       printf("... charge       :%12.6lf\n", si.qs[g]);
       printf("\n");
   }
   printf("\n");
   printf("________________ closure equation for electrons _____\n");
   if (si.CloseModelID == CLOSE_ISOTHERM) {
           printf("close model      : isothermal\n");
   }

   if (si.CloseModelID == CLOSE_POLYTROP) {
           printf("close model      : polytrop\n");
   }

   if (si.CloseModelID == CLOSE_FULLIMPL) {
       printf("close model      : full P (implicit)\n");
   }

   if (si.CloseModelID == CLOSE_FULLSUB) {
       printf("close model      : full P (subcycling)\n");
   }

   printf("\n");
   printf("________________ boundary conditions _________________\n");
   if (si.bc[0] == 0) printf("b-c (x)          :    periodic\n");
   if (si.bc[0] == 1) printf("b-c (x)          : perf. cond.\n");
   if (si.bc[0] == 2) printf("b-c (x)          :        open\n");
   if (si.bc[1] == 0) printf("b-c (y)          :    periodic\n");
   if (si.bc[1] == 1) printf("b-c (y)          : perf. cond.\n");
   if (si.bc[1] == 2) printf("b-c (y)          :        open\n");
   if (si.bc[2] == 0) printf("b-c (z)          :    periodic\n");
   if (si.bc[2] == 1) printf("b-c (z)          : perf. cond.\n");
   if (si.bc[2] == 2) printf("b-c (z)          :        open\n");
   printf("\n");
   printf("________________ files dumping _______________________\n");
   printf("dump field       :%12d\n", si.tf);
   printf("dump particle    :%12d\n", si.tp);
   printf("dump time        :%12d\n", si.tt);
   printf("dump restart     :%12d\n", si.tr);
   printf("\n");
   }

return si;

}



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "structures.h"
#include <particle.h>
#include "defines.h"
#include <stdint.h>


/* __ calculation of maximum value of a double vector _______________________ */
double maxdbl(double *w, int s)
{
double mw;
int i;


mw = w[0];

for (i = 1; i < s; i++)
    {
    if (w[i] > mw)
       {
       mw = w[i];
       }
    }

return mw;

}


/* __ calculation of minimum value of a double vector _______________________ */
double mindbl(double *w, int s)
{
double mw;
int i;


mw = w[0];

for (i = 1; i < s; i++)
    {
    if (w[i] < mw)
       {
       mw = w[i];
       }
    }

return mw;

}


/* __ calculation of minimum value of a int vector___________________________ */
int64_t minint64(int64_t *w, int s)
{
int64_t mw;
int i;


mw = w[1];

for (i = 1; i < s; i++)
    {
    if (w[i] < mw)
       {
       mw = w[i];
       }
    }

return mw;

}


/* __ calculation of maximum value of a int vector __________________________ */
int64_t maxint64(int64_t *w, int s)
{
int64_t mw;
int i;


mw = w[0];

for (i = 1; i < s; i++)
    {
    if (w[i] > mw)
       {
       mw = w[i];
       }
    }

return mw;

}


/* __ calculation of minimum value of a int vector___________________________ */
int minint(int *w, int s)
{
int mw;
int i;


mw = w[0];

for (i = 1; i < s; i++)
    {
    if (w[i] < mw)
       {
       mw = w[i];
       }
    }

return mw;

}


/* __ base of orthonormal vectors ___________________________________________ */
void ortho(double bw[3], double aw[3][3])
{
double rw, sw;


/* __ modulus of b __ */
rw = sqrt(bw[0]*bw[0]+bw[1]*bw[1]+bw[2]*bw[2]);

/* __ if modulus null : arbitrary base __ */
if (rw < EPS8)
   {
   /* __ unit vector along b __ */
   aw[0][0] = 1.0;
   aw[0][1] = 0.0;
   aw[0][2] = 0.0;

   /* __ unit vector perp to b __ */
   aw[1][0] = 0.0;
   aw[1][1] = 1.0;
   aw[1][2] = 0.0;

   /* __ last unit vector for direct triedr __ */
   aw[2][0] = 0.0;
   aw[2][1] = 0.0;
   aw[2][2] = 1.0;
   }

/* __ general case : need bz != 0 & bx != -by __ */
else
   {
   /* __ normalization for the second vector __ */
   //sw = sqrt(2.0*(bw[0]*bw[0]+bw[1]*bw[1]+bw[2]*bw[2]
   //              +bw[0]*bw[1]-bw[0]*bw[2]+bw[1]*bw[2]));


   /* __ unit vector along b __ */
   aw[0][0] = bw[0]/rw;
   aw[0][1] = bw[1]/rw;
   aw[0][2] = bw[2]/rw;


	/* __ unit vector perp to b __ */
	if (aw[0][0] == aw[0][1] && aw[0][1] == aw[0][2]) {
		aw[1][0] = aw[0][1];
		aw[1][1] = -aw[0][0];
		aw[1][2] = 0;
	} else {
		aw[1][0] = (aw[0][2] - aw[0][1]);   //(bw[2]+bw[1])/sw;
		aw[1][1] = (aw[0][0] - aw[0][2]);   //(bw[2]-bw[0])/sw;
		aw[1][2] = (aw[0][1] - aw[0][0]);   //-(bw[0]+bw[1])/sw;
	}

   sw = sqrt(aw[1][0]*aw[1][0] + aw[1][1]*aw[1][1] + aw[1][2]*aw[1][2]);
   aw[1][0] /= sw;
   aw[1][1] /= sw;
   aw[1][2] /= sw;

   /* __ last unit vector for direct triedr __ */
   aw[2][0] = aw[0][1]*aw[1][2]-aw[0][2]*aw[1][1];
   aw[2][1] = aw[0][2]*aw[1][0]-aw[0][0]*aw[1][2];
   aw[2][2] = aw[0][0]*aw[1][1]-aw[0][1]*aw[1][0];
   }
#if 0
    printf("|aw0| = %f\n", aw[0][0]*aw[0][0]+aw[0][1]*aw[0][1]+aw[0][2]*aw[0][2]);
    printf("|aw1| = %f\n", aw[1][0]*aw[1][0]+aw[1][1]*aw[1][1]+aw[1][2]*aw[1][2]);
    printf("|aw2| = %f\n", aw[2][0]*aw[2][0]+aw[2][1]*aw[2][1]+aw[2][2]*aw[2][2]);
#endif
}


/* __ fscanf function of a string ___________________________________________ */
void fscanstr(char *file, int line, FILE *fp, int r, char *buff, int *ir)
{
int z;


z = fscanf(fp, "%s", buff);

if (z != 1)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fscanf function of a int ______________________________________________ */
void fscanint(char *file, int line, FILE *fp, int r, int *buff, int *ir)
{
int z;


z = fscanf(fp, "%d", buff);

if (z != 1)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
printf("z = %d buff = %d\n", z, *buff);
   *ir = 1;
   }

}



/* __ fscanf function of a int ______________________________________________ */
void fscanint64(char *file, int line, FILE *fp, int r, int64_t *buff, int *ir)
{
int z;


//z = fscanf(fp, "%" PRId64, buff);
z = fscanf(fp, "%zu", buff);

if (z != 1)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   // printf("z = %d buff = %lld\n", z, *buff);
   *ir = 1;
   }

}




/* __ fscanf function of a double ___________________________________________ */
void fscandbl(char *file, int line, FILE *fp, int r, double *buff, int *ir)
{
int z;


z = fscanf(fp, "%lf", buff);

if (z != 1)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fread function for an int pointer _____________________________________ */
void freadint(char *file, int line, FILE *fp, int r, int length, int *buff, int *ir)
{
int z;


z = fread(buff, sizeof(int), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}




/* __ fread function for a float pointer ____________________________________ */
void freadflt(char *file, int line, FILE *fp, int r, int length, float *buff, int *ir)
{
int z;


z = fread(buff, sizeof(float), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fread function for a double pointer ___________________________________ */
void freaddbl(char *file, int line, FILE *fp, int r, int length, double *buff, int *ir)
{
int z;


z = fread(buff, sizeof(double), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fread function for a stx function _____________________________________ */
void freadstx(char *file, int line, FILE *fp, int r, int length, struct stx *buff, int *ir)
{
int z;


z = fread(buff, sizeof(struct stx), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fread function for a st1 function _____________________________________ */
void freadst1(char *file, int line, FILE *fp, int r, int length, struct st1 *buff, int *ir)
{
int z;


z = fread(buff, sizeof(struct st1), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fread function for a st2 function _____________________________________ */
void freadst2(char *file, int line, FILE *fp, int r, int length, struct st2 *buff, int *ir)
{
int z;


z = fread(buff, sizeof(struct st2), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


/* __ fread function for a stp function _____________________________________ */
void freadstp(char *file, int line, FILE *fp, int r, int length, struct stp *buff, int *ir)
{
int z;


z = fread(buff, sizeof(struct stp), length, fp);

if (z != length)
   {
   printf("file %s @ line %d on node %d\n", file, line, r);
   *ir = 1;
   }

}


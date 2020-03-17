
#ifndef MISC
#define MISC


#include <structures.h>
#include <particle.h>
#include <stdio.h>
#include <stdint.h>


///* _____ calculation of index on g0 grid ____________________________________ */
//int idxg0(int, int, int, struct stx);
//
///* _____ calculation of index on g1 grid ____________________________________ */
//int idxg1(int, int, int, struct stx);
//
///* _____ calculation of index on g2 grid ____________________________________ */
//int idxg2(int, int, int, struct stx);
//
///* _____ calculation of index on g4 grid ____________________________________ */
//int idxg4(int, int, int, struct stx);

/* _____ calculation of maximum value of a double vector ____________________ */
double maxdbl(double *, int);

/* _____ calculation of minimum value of a double vector_____________________ */
double mindbl(double *, int);

/* _____ calculation of minimum value of a int vector________________________ */
int minint64(int64_t *, int);

/* _____ calculation of maximum value of a int vector________________________ */
int maxint64(int64_t *, int);

/* _____ calculation of minimum value of a int vector________________________ */
int minint(int *, int);

/* _____ base of orthonormal vectors ________________________________________ */
void ortho(double [3], double [3][3]);

/* _____ fscanf function of a string ________________________________________ */
void fscanstr(char *, int, FILE *, int, char *, int *);

/* _____ fscanf function of an int __________________________________________ */
void fscanint(char *, int, FILE *, int, int *, int *);

/* _____ fscanf function of an int64_t ______________________________________ */
void fscanint64(char *file, int line, FILE *fp, int r, int64_t *buff, int *ir);

/* _____ fscanf function of a double ________________________________________ */
void fscandbl(char *, int, FILE *, int, double *, int *);

/* __ fread function for an int pointer _____________________________________ */
void freadint(char *, int, FILE *, int, int, int *, int *);

/* __ fread function for a float pointer ____________________________________ */
void freadflt(char *, int, FILE *, int, int, float *, int *);

/* __ fread function for a double pointer ___________________________________ */
void freaddbl(char *, int, FILE *, int, int, double *, int *);

/* __ fread function for a stx structure ____________________________________ */
void freadstx(char *, int, FILE *, int, int, struct stx *, int *);

/* __ fread function for a st1 structure ____________________________________ */
void freadst1(char *, int, FILE *, int, int, struct st1 *, int *);

/* __ fread function for a st2 structure ____________________________________ */
void freadst2(char *, int, FILE *, int, int, struct st2 *, int *);

/* __ fread function for a stp structure ____________________________________ */
void freadstp(char *, int, FILE *, int, int, struct stp *, int *);

#endif


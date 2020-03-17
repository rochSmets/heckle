
#ifndef FELAFEL
#define FELAFEL


/* _____ struct for the plasma topology : uniform _____ */
struct stt {
           double b[3];    /* __ magnetic field __ */
           double n[3];    /* __ density __ */
           double v[3];    /* __ specie drift velocity __ */
           double l;       /* __ X-position of velocity step __ */
           double t;       /* __ half-thickness of velocity step __ */
           double k;       /* __ total pressure __ */
           double e;       /* __ beta electronsi / beta total __ */
           double p;       /* __ beta protons / beta total __ */
           double i;       /* __ proton anisotropy factor __ */
           double j;       /* __ alpha anisotropy factor __ */
           };


#ifndef PYP
/* _____ rotational of magnetic field _______________________________________ */
void current(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ density ____________________________________________________________ */
void density(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ drift velocity _____________________________________________________ */
void drift(struct sti, struct stx, struct stt *, double [3], double [3][3], MPI_Comm);

/* _____ electric field _____________________________________________________ */
void electric(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ temperature and thermal velocity ___________________________________ */
void kinetic(struct sti, struct stx, struct stt *, double [3], double [3], double [3], MPI_Comm);

/* _____ magnetic field _____________________________________________________ */
void magnetic(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ calcul of the normalization coefficients ___________________________ */
void normal(struct sti *, struct stx, struct stt *, MPI_Comm);

/* _____ read the initial tangential discontiniuty configuration ____________ */
void topo(struct sti, struct stx, struct stt *, MPI_Comm);
#endif

#endif


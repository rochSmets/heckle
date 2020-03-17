
#ifndef YAOURT
#define YAOURT


/* _____ struct for the plasma topology : 1 beam _____ */
struct stt {
           double b0[2];   /* __ magnetic field __ */
           double fi;      /* __ tilt of the targets __ */
           double psi;     /* __ tilt of the targets __ */
           double nb[2];   /* __ proton density __ */
           double n0;      /* __ alpha density __ */
           double v0[2];   /* __ proton drift velocity __ */
           double k[2];    /* __ total pressure __ */
           double p[2];    /* __ beta protons / beta total __ */
           double e[2];    /* __ beta electrons / beta total __ */
           double lb[2];   /* __ initial radius of the bubble __ */
           double lr[2];   /* __ initial half-width of the ribbon __ */
           };


#ifndef PYP
/* _____ rotational of magnetic field _______________________________________ */
void current(struct sti, struct stx sx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ density ____________________________________________________________ */
void density(struct sti, struct stx sx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ drift velocity _____________________________________________________ */
void drift(struct sti, struct stx sx, struct stt *, double [3], double [3][3], MPI_Comm);

/* _____ electric field _____________________________________________________ */
void electric(struct sti, struct stx sx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ temperature and thermal velocity ___________________________________ */
void kinetic(struct sti, struct stx sx, struct stt *, double [3], double [3], double [3], MPI_Comm);

/* _____ magnetic field _____________________________________________________ */
void magnetic(struct sti, struct stx sx, struct stt *, double [3], double [3], MPI_Comm);

/* _____ calcul of the normalization coefficients ___________________________ */
void normal(struct sti *, struct stx sx, struct stt *, MPI_Comm);

/* __ polynom used for the interpolation __ */
double polynom(double);

/* _____ read the initial tangential discontiniuty configuration ____________ */
void topo(struct sti, struct stx, struct stt *, MPI_Comm);
#endif

#endif


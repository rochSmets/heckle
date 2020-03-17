
#ifndef GEM
#define GEM


/* _____ struct for the plasma topology : uniform _____ */
struct stt {
           double b;       /* __ magnetic field __ */
           double n;       /* __ proton density __ */
           double m;       /* __ alpha density __ */
           double l;       /* __ sheet thickness __ */
           double a;       /* __ temp electrons / temp protons __ */
           double p;       /* __ psi magnetic perturbation __ */
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
void topo(struct sti, struct stx *, struct stt *, MPI_Comm);
#endif

#endif



#ifndef SHELLS
#define SHELLS


/* _____ struct for the plasma topology : 1 beam _____ */
struct stt {
           double fi[2];   /* __ tilt of the targets __ */
           double psi[2];  /* __ tilt of the targets __ */
           double b0[2];   /* __ magnetic field __ */
           double nb[2];   /* __ proton density __ */
           double n0;      /* __ alpha density __ */
           double v0[2];   /* __ proton drift velocity __ */
           double te;      /* __ electrons temperature __ */
           double tp[2];   /* __ protons temperature __ */
           double ta;      /* __ alpha temperature __ */
           double ls[2];   /* __ initial radius of the shells __ */
           double lw[2];   /* __ initial half-width of the shells __ */
           double lz[2];   /* __ z extension of the shells __ */
           double xs[2];   /* __ x location of the shells __ */
           double ys[2];   /* __ y location of the shells __ */
           double zs[2];   /* __ z location of the shells __ */
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
void topo(struct sti, struct stx *, struct stt *, MPI_Comm);
#endif

#endif


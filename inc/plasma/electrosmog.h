
#ifndef ELECTROSMOG
#define ELECTROSMOG


/* __ struct for the plasma topology : munster _____ */
struct stt {
           double b[2];       /* __ asymptotic magnetic field __ */
           double m[2];       /* __ proton density __ */
           double l[2];       /* __ sheet thickness __ */
           double te;         /* __ electrons temerature __ */
           double th1;        /* __ proton temperature out of antenna __*/
           double d[2];       /* __ discontinuity positions __ */
           double a[2];       /* __ magnetic field angle with ? __*/
           double eb[3];      /* __ energy for b : x & y & z __ */
           double slb[2];     /* __ slope of b spectrum in x and y direction __  */
           int ix[2];         /* __ min & max (positive) modes x direction __ */
           int iy[2];         /* __ min & max (positive) modes y direction__ */
           int td;            /* __ # of time step for driving __ */
           int sd[2];         /* __ # of grid cell for driving __ */
           double c[2];       /* __ antenne center position x & y __ */
           double *bx;        /* __ magnetic energy density in x dir. __ */
           double *by;        /* __ magnetic energy density in y dir. __ */
           double *bz;        /* __ magnetic energy density in z dir. __ */
           double *px;        /* __ phase of the mode for x magnetic field __ */
           double *py;        /* __ phase of the mode for y magnetic field __ */
           double *pz;        /* __ phase of the mode for z magnetic field __ */
          };


#ifndef PYP
/* __ rotational of magnetic field __________________________________________ */
void current(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* __ density _______________________________________________________________ */
void density(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* __ drift velocity ________________________________________________________ */
void drift(struct sti, struct stx, struct stt *, double [3], double [3][3], MPI_Comm);

/* __ drive the magnetic field ______________________________________________ */
void drive(struct sti, struct stx, struct stt *, struct st1 *, struct st2 *, MPI_Comm);

/* __ electric field ________________________________________________________ */
void electric(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* __ temperature and thermal velocity ______________________________________ */
void kinetic(struct sti, struct stx, struct stt *, double [3], double [3], double [3], MPI_Comm);

/* __ magnetic field ________________________________________________________ */
void magnetic(struct sti, struct stx, struct stt *, double [3], double [3], MPI_Comm);

/* __ calcul of the normalization coefficients ______________________________ */
void normal(struct sti *, struct stx, struct stt *, MPI_Comm);

/* __ read the initial tangential discontiniuty configuration _______________ */
void topo(struct sti, struct stx *, struct stt *, MPI_Comm);
#endif

#endif


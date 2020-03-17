
#ifndef MUNSTER
#define MUNSTER


/* __ struct for the plasma topology : munster _____ */
struct stt {
           double b[3];    /* __ magnetic field __ */
           double ne;      /* __ electron density __ */
           double np;      /* __ proton density __ */
           double na;      /* __ alpha density __ */
           double te;      /* __ electron temperature __ */
           double tp;      /* __ proton temperature __ */
           double ta;      /* __ alpha temperature __ */
           double eb[2];   /* __ energy for b : x & y __ */
           double slb;     /* __ slope of b spectrum __ */
           double ee[3];   /* __ energy for e : x, y & z __ */
           double sle;     /* __ slope of e spectrum __ */
           int m[2];       /* __ min & max (positive) modes __ */
           int td;         /* __ # of time step for driving __ */
           int sd;         /* __ # of grid cell for driving __ */
           double *bx;     /* __ magnetic energy density in x dir. __ */
           double *by;     /* __ magnetic energy density in y dir. __ */
           double *ex;     /* __ electric energy density in x dir. __ */
           double *ey;     /* __ electric energy density in y dir. __ */
           double *ez;     /* __ electric energy density in z dir. __ */
           double *px;     /* __ phase of the mode for x magnetic field __ */
           double *py;     /* __ phase of the mode for y magnetic field __ */
           double *fx;     /* __ phase of the mode for x electric field __ */
           double *fy;     /* __ phase of the mode for y electric field __ */
           double *fz;     /* __ phase of the mode for z electric field __ */
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


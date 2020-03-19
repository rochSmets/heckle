
#ifndef STRUCTURES
#define STRUCTURES

#include <mpi.h>
#define MAXNUMPART                                                             \
  1500 // affreux !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include <stdint.h>

/* _____ struct for the simulation parameters _____ */
typedef struct sti {
  int n[3];           /* __ # of cells __ */
  double l[3];        /* __ box size __ */
  int bc[3];          /* __ boundary condition __ */
  double ts;          /* __ time step __ */
  int nt;             /* __ total number of time step __ */
  int tf;             /* __ # of time steps to dump a field file __ */
  int tp;             /* __ # of time steps to dump a part. file __ */
  int tt;             /* __ # of time steps to dump a time file __ */
  int64_t nm;         /* __ maximum number of particles __ */
  double nmin;        /* __ min value of the density for Ohm's law __ */
  int feed;           /* __ feed the simulation __ */
  int coll;           /* __ consider collisions __ */
  int drive;          /* __ if the simulation needs to be driven __ */
  double me;          /* __ electron mass __ */
  double kapp;        /* __ thermal conductivity __ */
  double rsty;        /* __ resistivity __ */
  double hyvi;        /* __ hyperviscosity __ */
  int rst;            /* __ use of a restart __ */
  int tr;             /* __ # of time step to dump a restart __ */
  double time4rst;    /* __ time for restart __ */
  int64_t ns[NS + 1]; /* __ # of part. __ */
  double ms[NS + 1];  /* __ mass of part. __ */
  double qs[NS + 1];  /* __ charge of part. __ */
  double ws[NS + 1];  /* __ weight of part. __ */
  double dl[3];       /* __ grid step __ */
  int ncpu;           /* __ # of cpu __ */
  int mpidom[3];      /* __ # of mpi domains in each direction __ */
  int InitModelID;    /* __ initial condition/model __ */
  int CloseModelID;   /* __ closure condition/model __ */
} STI;

/* _____ struct for the simulation parameters _____ */
typedef struct stc {
  double nuIonIon;   /* __ ion-ion collision frequency __ */
  double coulombLog; /* __ coulombian logarithm __ */
} Collision;

/* _____ struct for the fields on g1 grid _____ */
typedef struct st0 {
  int ppc;            /* __ # of part in the cell whatever specie __ */
  int npart[NS + 1];  /* __ # of part in the cell (collision) __ */
  int *ipart[NS + 1]; /* __ index of part in the cell (collision) __ */
} Grid0;

/* _____ struct for the fields on g1 grid _____ */
typedef struct st1 {
  double b[3]; /* __ magnetic field @ predictor step __ */
  double c[3]; /* __ magnetic field @ corrector step __ */
#ifdef FULLP
  double bKept[3]; /* __  value of the magnetic field  at n- 1/2__ */
#endif
} ST1;

/* _____ struct for the fields and moments on g2 grid _____ */
typedef struct st2 {
  double e[3];          /* __ electric field (predictor) __ */
  double f[3];          /* __ electric field (corrector) __ */
  double r;             /* __ resistivity __ */
  double s;             /* __ "entropy" defined as p.pow(n,gamma) __ */
  double te[6];         /* __ electron temperature tensor __ */
  double vi[3];         /* __ ion fluid velocity __ */
  double j[3];          /* __ current density __ */
  double j_smoothed[3]; /* __ current density __ */
  double ms[NS + 1];    /* __ specie charge density (even # half dt) __ */
  double ns[NS + 1];    /* __ specie charge density (odd # half dt) __ */
  double os[NS + 1];    /* __ specie charge density (even # half dt) __ */
  double vs[NS + 1][3]; /* __ specie fluid velocity __ */
  double ps[NS + 1][6]; /* __ specie full pressure tensor __ */
#ifdef FULLP
  double pKept[6]; /* __ keep pressure [n-1/2] for ipc = 0 __ */
  double dHalf[6]; /* __ driver @ half time step __ */
  double dFull[6]; /* __ driver @ full time step __ */
#endif
} ST2;

typedef struct st4 {
  double n;
  double v[3];
  double j_smoothed[3];
  double ps[6];
} ST4;

/* _____ struct for all fields on g1 _____ */
struct stf {
  float *b[3];         /* __ magnetic field __ */
  float *e[3];         /* __ electric field __ */
  float *a[3];         /* __ potential vector __ */
  float *j[3];         /* __ current density __ */
  float *vi[3];        /* __ ion fluid velocity __ */
  float *n[NS + 1];    /* __ specie density __ */
  float *v[NS + 1][3]; /* __ specie fluid velocity __ */
  float *p[NS + 1][6]; /* __ full specie pressure tensor __ */
};

/* _____ struct for the particles position and velocity _____ */
struct stq {
  float *r[NS + 1][3]; /* __ particle position __ */
  float *v[NS + 1][3]; /* __ particle velocity __ */
};

/* _____ struct for the part. orbits _____ */
struct sto {
  int wo;  /* __ # of orbits in the whole domain __ */
  int *so; /* __ s specie of part. (size wo) __ */
  int *io; /* __ i index of part. (size wo) __ */
  int no;  /* __ # of orbits in the sub-domain __ */
  int *s;  /* __ s specie of part. (size no) __ */
  int *m;  /* __ m index of part. (size no) __ */
};

/* _____ struct for the part. maps (orbits to follow) _____ */
struct stm {
  float *r[3];         /* __ orbit position __ */
  float *w[3];         /* __ orbit velocity __ */
  float *b[3];         /* __ magnetic field along the orbit __ */
  float *e[3];         /* __ electric field along the orbit __ */
  float *j[3];         /* __ current density field along... __ */
  float *vi[3];        /* __ ion fluid velocity field along... __ */
  float *n[NS + 1];    /* __ density field along... __ */
  float *v[NS + 1][3]; /* __ specie fluid velocity field along... __ */
  float *p[NS + 1][6]; /* __ full specie pressure tensor along... __ */
};

/* _____ struct for the time dump _____ */
struct std {
  double db;         /* __ divergence b __ */
  double dw;         /* __ pseudo divergence b __ */
  int pc;            /* __ min # of parts per cells __ */
  int64_t no;        /* __ min # of parts in subdomain __ */
  int64_t nm;        /* __ max # of parts in subdomain __ */
  double ma;         /* __ magnetic energy __ */
  double e0;         /* __ initial total energy __ */
  double fx;         /* __ electric fluctuation in X direction __ */
  double fy;         /* __ electric fluctuation in Y direction __ */
  double fz;         /* __ electric fluctuation in Z direction __ */
  double pb[NS + 1]; /* __ bulk energy __ */
  double ta[NS + 1]; /* __ parallel thermal energy __ */
  double te[NS + 1]; /* __ perpendicular thermal energy __ */
  double lii;        /* __ min of coulomb log(ion-ion) __ */
  double niimin;     /* __ min of nu(ion-ion) __ */
};

/* _____ struct for mpi _____ */
typedef struct stx {
  int s;              /* __ total # of nodes __ */
  int r;              /* __ index of nodes __ */
  int d[3];           /* __ # of subdomains in 3 directions __ */
  int64_t ns[NS + 1]; /* __ # of part __ */
  int i0[3];          /* __ index min (g1) in each directions __ */
  int i1[3];          /* __ index max (g1) in each directions __ */
  int n[3];           /* __ # of cells of subdomain __ */
  double l[3];        /* __ size of subdomain __ */
  int nt[27];         /* __ index of neighbors (send to) __ */
  int nf[27];         /* __ index of neighbors (get from) __ */
  int64_t no[NS + 1]; /* __ minimum # of part in subdomain __ */
  int64_t nm[NS + 1]; /* __ maximum # of part in subdomain __ */
  int irun;           /* __ index null if the code has to run __ */
  MPI_Comm com;       /* __ mpi communicator __ */
} STX;

#endif

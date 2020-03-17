
#ifndef WRITE
#define WRITE



typedef struct s_HeckleIOFields HeckleIOFields;
typedef struct s_HeckleIOSpecies HeckleIOSpecies;
typedef struct s_HeckleIORestart HeckleIORestart;




/*---------------------------------------------------------------------------
  HeckleIOInitFields()
  ---------------------------------------------------------------------------
  AIM : Initialize the heckleIOFields module
 ---------------------------------------------------------------------------*/
HeckleIOFields *HeckleIOInitFields(struct sti *si, struct stx *sx);


/*---------------------------------------------------------------------------
  HeckleIOInitSpecies()
  ---------------------------------------------------------------------------
  AIM : Initialize the heckleIOSpecies module
 ---------------------------------------------------------------------------*/
HeckleIOSpecies *HeckleIOInitSpecies(struct sti *si, struct stx *sx);


/*---------------------------------------------------------------------------   
  HeckleIOInitRestarts()                                                        
  ---------------------------------------------------------------------------   
  AIM : Initialize the HeckleIORestarts module                                  
 ---------------------------------------------------------------------------*/  
HeckleIORestart *HeckleIOInitRestarts(struct sti *si, struct stx *sx);



/*---------------------------------------------------------------------------
  HeckleIODeleteFields()
  ---------------------------------------------------------------------------
  AIM : Delete the heckleIOFields module handle
 ---------------------------------------------------------------------------*/
void HeckleIODeleteFields(HeckleIOFields *iof);


/*---------------------------------------------------------------------------
  HeckleIODeleteSpecies()
  ---------------------------------------------------------------------------
  AIM : Delete the heckleIOSpecies module handle
 ---------------------------------------------------------------------------*/
void HeckleIODeleteSpecies(HeckleIOSpecies *ios);



/*---------------------------------------------------------------------------   
  HeckleIODeleteRestarts()                                                      
  ---------------------------------------------------------------------------   
  AIM :                                                                         
 ---------------------------------------------------------------------------*/  
void HeckleIODeleteRestarts(HeckleIORestart *hior);


/*---------------------------------------------------------------------------
  writeFields()
  ---------------------------------------------------------------------------
  AIM : this routine writes the fields (electromagnetic and fluid moments)
  in a HDF5 file.
 ---------------------------------------------------------------------------*/

void writeFields(HeckleIOFields *hiof,  /* IO module handle  */
                 struct stx *sx,        /* MPI parameters    */
                 struct st1 *s1,        /* g1 grid           */
                 struct st2 *s2,        /* g2 grid           */
                 double time);          /* current time      */


/*---------------------------------------------------------------------------
  writeSpecies()
  ---------------------------------------------------------------------------
  AIM : this routine writes the species (position, velocity & id)
  in a HDF5 file.
 ---------------------------------------------------------------------------*/
void writeSpecies(HeckleIOSpecies *hios, /* IO module handle  */
                  struct sti *si,         /* run parameters    */
                  struct stx *sx,         /* MPI parameters    */
                  struct stp *sp[NS+1],   /* part positions, velocities & id */
                  double time);           /* current time      */


/*---------------------------------------------------------------------------   
  writeRestarts()                                                               
  ---------------------------------------------------------------------------   
  AIM : this routine writes the restarts in a HDF5 file.                        
 ---------------------------------------------------------------------------*/  
void writeRestarts(HeckleIORestart *hior,                                       
                 struct sti *si,                                                
                 struct stx *sx,                                                
                 struct st1 *s1,                                                
                 struct st2 *s2,                                                
                 struct stp *sp[NS+1],                                          
                 struct std sd,
                 double time);



/* _____ write the time dump ________________________________________________ */
void writedump(struct sti,
               struct stx,
               struct std,
               int);


/* _____ write the xplosed orbits ___________________________________________ */
void writeorbit(struct sti,
                struct stx *,
                struct st1 *,
                struct st2 *,
                struct stp *[NS+1],
                struct sto *,
                int);

/* _____ write the xplosed particles ________________________________________ */
void writepart(struct sti,
               struct stx,
               struct stp *[NS+1],
               int);

/* _____ write the xplosed restart __________________________________________ */
void writerestart(struct sti,
                  struct stx,
                  struct st1 *,
                  struct st2 *,
                  struct stp *[NS+1],
                  struct std,
                  int);

#endif



#ifndef READ
#define READ


/* __ read the full field file ______________________________________________ */
void readfullfield(struct sti, struct stf *, char *);

/* __ read the full particle file ___________________________________________ */
void readfullpart(struct sti *, struct stq *, char *);

/* _____ read the heckle.txt file ___________________________________________ */
void readheckle(struct sti *, char *, int);

/* _____ read the horbi.txt file ____________________________________________ */
void readhorbi(struct sto *);

/* _____ read the xplosed field files _______________________________________ */
void readfield(struct sti, struct stf *, int, int);

/* _____ read the xplosed particles files ___________________________________ */
void readpart(struct sti *, struct stq *, int, int);

/* __ read the xplosed orbits files _________________________________________ */
void readorbit(struct sto, struct stm *, int, int);

#endif


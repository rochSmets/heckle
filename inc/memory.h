
#ifndef MEMORY
#define MEMORY


/* _____ free the stf fields ________________________________________________ */
void freestf(struct stf *);

/* _____ free the stq particles _____________________________________________ */
void freestq(struct stq *);

/* __ free a stm structure  _________________________________________________ */
void freestm(struct stm *);

/* _____ memory allocation for the stf fields _______________________________ */
void memorystf(struct stf *, struct sti);

/* _____ memory allocation for the stq particles ____________________________ */
void memorystq(struct stq *, struct sti);

/* _____ memory allocation for a stm structure ______________________________ */
void memorystm(struct stm *, int);

#endif


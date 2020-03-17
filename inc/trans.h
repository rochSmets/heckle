
#ifndef TRANS
#define TRANS


/* ____ update moments from particles dynamics ______________________________ */
void trans(struct sti *, struct stx *, struct stt *, struct st1 *, struct st2 *, struct stp *[NS+1], int, int, MPI_Comm);

#endif


#ifndef DEAD
#define DEAD

#include <structures.h>
#include <particle.h>


/* __ report problem if part is not in the subdomain ________________________ */
void deadpart(struct sti,
              struct stx *,
              struct st0 *,
              struct st1 *,
              struct st2 *,
              struct stp *[NS+1],
              int,
              int,
              int,
              char *,
              int);

#endif


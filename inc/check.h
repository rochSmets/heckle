
#ifndef CHECK
#define CHECK


#include <structures.h>
#include <particle.h>

/* _____ check conservation of div b, p(z) and total energy _________________ */
void check(struct sti,
           struct stx,
           struct st0 *,
           struct st1 *,
           struct st2 *,
           struct stp *[NS+1],
           struct std *,
           int,
           clock_t [2],
           time_t [2]);

/* _____ calcul of divergence b _____________________________________________ */
void db(struct sti,
        struct stx,
        struct st0 *,
        struct st1 *,
        double [2],
        int *);

/* _____ calculate the level of electric fluctuation ________________________ */
void ef(struct stx,
        struct st2 *,
        double [3],
        double [3]);

/* _____ calcul of total energy _____________________________________________ */
void et(struct sti,
        struct stx,
        struct st1 *,
        struct st2 *,
        struct stp *[NS+1],
        double[9],
        double wle[3]);

#endif


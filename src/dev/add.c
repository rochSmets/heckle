
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "structures.h"
#include "plasma.h"
#include "defines.h"
#include "misc.h"


/* __ add part. when the density is too low _________________________________ */
void add(struct sti *si, struct stx *sx, struct stt *st, struct st2 *s2, struct stp *sp[NS+1], MPI_Comm com)
{
double sa, sn;
double ww;
double nb[NS+1];
double aw[3], ew[3];
double qw[3];
double bw[3];
double vw[3][3];
double rw, w0, w1, w2;
double dw[3][3];
int m0, m1, wa;
int nw;
int ijk;
int d, i, j, k, l, m, s;
int *na[NS+1], ta[NS+1];
MPI_Comm co;


/* __ duplicate the communicator __ */
MPI_Comm_dup(com, &co);

/* __ loop on the species __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ left boundary __ */
    i = 0;

    sa = 0.0;
    sn = 0.0;

    /* __ left boundary : loop on the j index __ */
    for (j = 0; j < sx->n[1]; j++)
        {
        /* __ left boundary : loop on the k index __ */
        for (k = 0; k < sx->n[2]; k++)
            {
            /* __ index on the g2 grid __ */
            ijk = idx(i+1, j+1, k+1, sx->n[0]+2, sx->n[1]+2, sx->n[2]+2);

            /* __ coordinates at the middle of the cell __ */
            qw[0] = (i+0.5+sx->i0[0])*(si->dl[0]);
            qw[1] = (j+0.5+sx->i0[1])*(si->dl[1]);
            qw[2] = (k+0.5+sx->i0[2])*(si->dl[2]);

            /* __ asymptotic density on the boundary __ */
            density(*si, *sx, st, qw, nb, co);

            /* __ sum for analytical & numerical density __ */
            sa += nb[s];
            sn += s2[ijk].ns[s];
            }
        }

    /* __ min & max m values of added part. __ */
    m0 = sx->ns[s];
    m1 = sx->ns[s]-1;

    /* __ # of part. to add if no left neighbor __ */
    ww = (sa-sn)/si->ws[s];
    nw = (ww > 0.0 && sx->i0[0] == 0) ? (int)round(ww) : 0;

    /* __ loop on the # of part to add __ */
    for (d = 0; d < nw; d++)
        {
        /* __ coordinates of grid point __ */
        qw[0] = (RNM)*(si->dl[0]);
        qw[1] = (RNM)*(si->l[1]);
        qw[2] = (RNM)*(si->l[2]);

        /* __ index of the added part. __ */
        m = ++m1;

        /* __ set position __ */
        sp[s][m].r[0] = (qw[0] == sx->n[0]*si->dl[0]) ?
                         sx->i1[0]*si->dl[0]-EPS6 :
                         sx->i0[0]*si->dl[0]+qw[0];
        sp[s][m].r[1] = (qw[1] == sx->n[1]*si->dl[1]) ?
                         sx->i1[1]*si->dl[1]-EPS6 :
                         sx->i0[1]*si->dl[1]+qw[1];
        sp[s][m].r[2] = (qw[2] == sx->n[2]*si->dl[2]) ?
                         sx->i1[2]*si->dl[2]-EPS6 :
                         sx->i0[2]*si->dl[2]+qw[2];

        /* __ set proj. position __ */
        sp[s][m].s[0] = sp[s][m].r[0];
        sp[s][m].s[1] = sp[s][m].r[1];
        sp[s][m].s[2] = sp[s][m].r[2];

        /* __ set correct value of qw in the large domain __ */
        qw[0] += sx->i0[0]*si->dl[0];
        qw[1] += sx->i0[1]*si->dl[1];
        qw[2] += sx->i0[2]*si->dl[2];

        /* __ set # of crossing of the boundary domain __ */
        sp[s][m].b[0] = 0;
        sp[s][m].b[1] = 0;
        sp[s][m].b[2] = 0;

        /* __ drift velocity __ */
        drift(*si, *sx, st, qw, dw, co);

        /* __ parallel & perp. temperature __ */
        kinetic(*si, *sx, st, qw, aw, ew, co);

        /* __ parallel velocity (box & muller algo.) __ */
        ww = RNM;
        ww = (ww != 1.0) ? ww : 1.0-EPS6;
        rw = -log(ww)*aw[s];
        ww = RNM;
        w0 = sqrt(2.0*rw/si->ms[s])*cos(ww*2.0*PI);

        /* __ perpendicular velocity __ */
        ww = RNM;
        ww = (ww != 1.0) ? ww : 1.0-EPS6;
        rw = -log(ww)*ew[s];
        ww = RNM;
        w1 = sqrt(2.0*rw/si->ms[s])*cos(ww*2.0*PI);
        w2 = sqrt(2.0*rw/si->ms[s])*sin(ww*2.0*PI);

        /* __ base of orthonormal vectors __ */
        ortho(bw, vw);

        /* __ loop in the 3 directions __ */
        for (l = 0; l < 3; l++)
            {
            /* __ set velocity __ */
            sp[s][m].v[l] = vw[0][l]*w0+vw[1][l]*w1+vw[2][l]*w2+dw[s][l];

            /* __ set proj. velocity __ */
            sp[s][m].w[l] = sp[s][m].v[l];
            }
    }

    /* __ update the # of part. __ */
    sx->ns[s] = (m0 > m1) ? m0 : m1;

    /* __ # of added part. __ */
    wa = m1-m0+1;

    /* __ buffer needed for the total # of part. in nodes < r __ */
    na[s] = (int *)malloc(sx->s*sizeof(int));

    /* __ set the total # of part. in nodes < r __ */
    MPI_Allgather(&wa, 1, MPI_INT, na[s], 1, MPI_INT, co);

    /* __ set counter for the total # of added part. in nodes < r __ */
    ta[s] = 0;

    /* __ total # of added part. in nodes < r __ */
    for (d = 0; d < sx->r; d++)
        {
        ta[s] += na[s][d];
        }

    /* __ set i index of part. : loop on added part.__ */
    for (m = m0; m < m1; m++)
        {
        /* __ set i index of part. __ */
        sp[s][m].i = si->ns[s]+ta[s]+m-m0;
        }

    /* __ total # of part. __ */
    MPI_Allreduce(&(sx->ns[s]), &(si->ns[s]), 1, MPI_INT, MPI_SUM, co);

    /* __ reduce max # of part. __ */
    MPI_Allreduce(&(sx->ns[s]), &(sx->nm[s]), 1, MPI_INT, MPI_MAX, co);

    /* __ clean-up the pointers __ */
    free(na[s]);
    }

/* __ clean up the communicator __ */
MPI_Comm_free(&co);

}


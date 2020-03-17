
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "structures.h"
#include "defines.h"
#include "iamdead.h"


/* __ manage part. when reaching the walls __________________________________ */
void wall(struct sti *si,
          struct stx *sx,
          struct st1 *s1,
          struct st2 *s2,
          struct stp *sp[NS+1],
          int ipc,
          int s,
          int it)
{
    double xw, yw, zw;
    int a, b, c, l, m, n, t;
    int uc[27], vc[27], wc[27];
    int ps;
    char *ws[27], *wr[27];
    MPI_Status st;
    int packsize_int, packsize_dbl;



    MPI_Pack_size(4 ,  MPI_INT,    MPI_COMM_WORLD, &packsize_int);
    MPI_Pack_size(12,  MPI_DOUBLE, MPI_COMM_WORLD, &packsize_dbl);


    /* __ maximum size of packet __ */
    //ps = (sx->n[0]*2 + sx->n[1]*2 + sx->n[2]*2)*2000;
    //ps = sx->nm[s]*sizeof(struct stp);
    //ps = sx->nm[s]*(packsize_dbl + packsize_int);


    /* this assumes that at most 200 per cell of the perimeter of the MPI domain
   may be transfered. This is an overallocation that may be too important
   and something smarter should be done here */
    //ps = (sx->n[0]*2 + sx->n[1]*2+sx->n[2]*2)*2000*(packsize_dbl + packsize_int);


    /* --------------------------------------------------------------------- */
    /*                  REFLECTION BOUNDARY CONDITION                        */
    /* --------------------------------------------------------------------- */

    /* __ predictor __ */
    if (ipc == 0)
    {
        /* __ loop on all the part. __ */
        for (m = 0; m < sx->ns[s]; m++)
        {
            /* __ if non-periodic in (x) on the boundary of the domain __ */
            if (si->bc[0] == 0)
            {
                /* __ left boundary __ */
                if (sp[s][m].r[0] < 0.0)
                {
                    /* __ particle reflexion __ */
                    sp[s][m].r[0] = -sp[s][m].r[0];
                    sp[s][m].v[0] = -sp[s][m].v[0];
                }
                /* __ right boundary __ */
                if (sp[s][m].r[0] >= si->l[0])
                {
                    /* __ particle reflexion __ */
                    sp[s][m].r[0] = -sp[s][m].r[0]+2.0*si->l[0];
                    sp[s][m].v[0] = -sp[s][m].v[0];
                }
            } /* End if x non periodic */


            /* __ if non-periodic in (y) on the boundary of the domain __ */
            if (si->bc[1] == 0)
            {
                /* __ bottom boundary __ */
                if (sp[s][m].r[1] < 0.0)
                {
                    /* __ particle reflexion __ */
                    sp[s][m].r[1] = -sp[s][m].r[1];
                    sp[s][m].v[1] = -sp[s][m].v[1];
                }
                /* __ top boundary __ */
                if (sp[s][m].r[1] >= si->l[1])
                {
                    /* __ particle reflexion __ */
                    sp[s][m].r[1] = -sp[s][m].r[1]+2.0*si->l[1];
                    sp[s][m].v[1] = -sp[s][m].v[1];
                }
            }

            /* __ if non-periodic in (z) on the boundary of the domain __ */
            if (si->bc[2] == 0)
            {
                /* __ back boundary __ */
                if (sp[s][m].r[2] < 0.0)
                {
                    /* __ particle reflexion __ */
                    sp[s][m].r[2] = -sp[s][m].r[2];
                    sp[s][m].v[2] = -sp[s][m].v[2];
                }
                /* __ front boundary __ */
                if (sp[s][m].r[2] >= si->l[2])
                {
                    /* __ particle reflexion __ */
                    sp[s][m].r[2] = -sp[s][m].r[2]+2.0*si->l[2];
                    sp[s][m].v[2] = -sp[s][m].v[2];
                }
            }
        }
    } /* End predictor */



    /* __ corrector __ */
    if (ipc == 1)
    {
        /* __ loop on all the part. __ */
        for (m = 0; m < sx->ns[s]; m++)
        {
            /* __ if non-periodic in (x) on the boundary of the domain __ */
            if (si->bc[0] == 0)
            {
                /* __ left boundary __ */
                if (sp[s][m].s[0] < 0.0)
                {
                    /* __ particle reflexion __ */
                    sp[s][m].s[0] = -sp[s][m].s[0];
                    sp[s][m].w[0] = -sp[s][m].w[0];
                }
                /* __ right boundary __ */
                if (sp[s][m].s[0] >= si->l[0])
                {
                    /* __ particle reflexion __ */
                    sp[s][m].s[0] = -sp[s][m].s[0]+2.0*si->l[0];
                    sp[s][m].w[0] = -sp[s][m].w[0];
                }
            }

            /* __ if non-periodic in (y) on the boundary of the domain __ */
            if (si->bc[1] == 0)
            {
                /* __ bottom boundary __ */
                if (sp[s][m].s[1] < 0.0)
                {
                    /* __ particle reflexion __ */
                    sp[s][m].s[1] = -sp[s][m].s[1];
                    sp[s][m].w[1] = -sp[s][m].w[1];
                }
                /* __ top boundary __ */
                if (sp[s][m].s[1] >= si->l[1])
                {
                    /* __ particle reflexion __ */
                    sp[s][m].s[1] = -sp[s][m].s[1]+2.0*si->l[1];
                    sp[s][m].w[1] = -sp[s][m].w[1];
                }
            }

            /* __ if non-periodic in (z) on the boundary of the domain __ */
            if (si->bc[2] == 0)
            {
                /* __ back boundary __ */
                if (sp[s][m].s[2] < 0.0)
                {
                    /* __ particle reflexion __ */
                    sp[s][m].s[2] = -sp[s][m].s[2];
                    sp[s][m].w[2] = -sp[s][m].w[2];
                }
                /* __ front boundary __ */
                if (sp[s][m].s[2] >= si->l[2])
                {
                    /* __ particle reflexion __ */
                    sp[s][m].s[2] = -sp[s][m].s[2]+2.0*si->l[2];
                    sp[s][m].w[2] = -sp[s][m].w[2];
                }
            }
        }
    }


    /* --------------------------------------------------------------------- */
    /*                       PARTICLES OUTSIDE THE DOMAIN                    */
    /* --------------------------------------------------------------------- */

    /* first count the particles leaving and get the max over all processors
       so that we can allocate buffers to that number */
    int cpt_out_loc = 0;
    int cpt_out_glo = 0;
    for (m = 0; m < sx->ns[s]; m++)
    {
        /* __ part "position" needed to send to the neighbor __ */
        xw = (sp[s][m].r[0]-sx->i0[0]*si->dl[0])/sx->l[0];
        yw = (sp[s][m].r[1]-sx->i0[1]*si->dl[1])/sx->l[1];
        zw = (sp[s][m].r[2]-sx->i0[2]*si->dl[2])/sx->l[2];

        /* __ set the index of the neighbor to send the part. __ */
        a = (int)floor(xw);
        b = (int)floor(yw);
        c = (int)floor(zw);


        /* __ set the tag value __ */
        t = (1+c)+3*((1+b)+3*(1+a));


        /* __ condition for the part. to be out of the subdomain __ */
        if (t != 13)
        {
            cpt_out_loc++;
        }

    }

    MPI_Allreduce(&cpt_out_loc, &cpt_out_glo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    ps =cpt_out_glo*(packsize_dbl + packsize_int);


    /* __ loop on all the neighbors __ */
    for (t = 0; t < 27; t++)
    {
        /* __ memory allocation for send & receive packet __ */
        ws[t] = (char *)malloc(ps);
        wr[t] = (char *)malloc(ps);


#if BUG
        if (ws[t] == NULL)
        {
            printf("BORDEL wt NULL %d %ld %d %ld\n", ps,
                   sx->nm[s]*sizeof (struct stp), sx->nm[s], sizeof(struct stp));
        }

        if (wr[t] == NULL)
        {
            printf("BORDEL wr[%d] NULL %d\n", t, ps);
        }
#endif

        /* __ set the size of the packet __ */
        uc[t] = 0;
        vc[t] = 0;
    }


    /* __ prepare the packet if part. is out of subdomain : loop on all the part. __ */
    for (m = 0; m < sx->ns[s]; m++)
    {

        /* __ part "position" needed to send to the neighbor __ */
        xw = (sp[s][m].r[0]-sx->i0[0]*si->dl[0])/sx->l[0];
        yw = (sp[s][m].r[1]-sx->i0[1]*si->dl[1])/sx->l[1];
        zw = (sp[s][m].r[2]-sx->i0[2]*si->dl[2])/sx->l[2];

        /* __ set the index of the neighbor to send the part. __ */
        a = (int)floor(xw);
        b = (int)floor(yw);
        c = (int)floor(zw);


#ifdef BUG
        if ((a < -1 || a > +1) && it == 1)
        {
            deadpart(*si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
        if ((b < -1 || b > +1) && it == 1)
        {
            deadpart(*si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
        if ((c < -1 || c > +1) && it == 1)
        {
            deadpart(*si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
        }
#endif

        /* __ set the tag value __ */
        t = (1+c)+3*((1+b)+3*(1+a));


        /* __ condition for the part. to be out of the subdomain __ */
        if (t != 13)
        {
#ifdef BUG
            if (m >= sx->nm[s])
            {
                deadpart(*si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
            }
#endif

            /* __ set the counter __ */
            vc[t]++;

            /* __ need to shift the part. if periodic in (x) direction __ */
            if (si->bc[0] != 0)
            {
                if (sp[s][m].r[0] < 0.0)
                {
                    sp[s][m].r[0] += si->l[0];
                    sp[s][m].s[0] += si->l[0];
                    sp[s][m].b[0] -= 1;
                }
                if (sp[s][m].r[0] >= si->l[0])
                {
                    sp[s][m].r[0] -= si->l[0];
                    sp[s][m].s[0] -= si->l[0];
                    sp[s][m].b[0] += 1;
                }
            }

            /* __ need to shift the part. if periodic in (y) direction __ */
            if (si->bc[1] != 0)
            {
                if (sp[s][m].r[1] < 0.0)
                {
                    sp[s][m].r[1] += si->l[1];
                    sp[s][m].s[1] += si->l[1];
                    sp[s][m].b[1] -= 1;
                }
                if (sp[s][m].r[1] >= si->l[1])
                {
                    sp[s][m].r[1] -= si->l[1];
                    sp[s][m].s[1] -= si->l[1];
                    sp[s][m].b[1] += 1;
                }
            }

            /* __ need to shift the part. if periodic in (z) direction __ */
            if (si->bc[2] != 0)
            {
                if (sp[s][m].r[2] < 0.0)
                {
                    sp[s][m].r[2] += si->l[2];
                    sp[s][m].s[2] += si->l[2];
                    sp[s][m].b[2] -= 1;
                }
                if (sp[s][m].r[2] >= si->l[2])
                {
                    sp[s][m].r[2] -= si->l[2];
                    sp[s][m].s[2] -= si->l[2];
                    sp[s][m].b[2] += 1;
                }
            }

            /* __ pack the data __ */
            MPI_Pack(&sp[s][m].i,  1, MPI_INT   , ws[t], ps, &uc[t], MPI_COMM_WORLD);
            MPI_Pack( sp[s][m].r,  3, MPI_DOUBLE, ws[t], ps, &uc[t], MPI_COMM_WORLD);
            MPI_Pack( sp[s][m].s,  3, MPI_DOUBLE, ws[t], ps, &uc[t], MPI_COMM_WORLD);
            MPI_Pack( sp[s][m].v,  3, MPI_DOUBLE, ws[t], ps, &uc[t], MPI_COMM_WORLD);
            MPI_Pack( sp[s][m].w,  3, MPI_DOUBLE, ws[t], ps, &uc[t], MPI_COMM_WORLD);
            MPI_Pack( sp[s][m].b,  3, MPI_INT   , ws[t], ps, &uc[t], MPI_COMM_WORLD);

            /* __ decrease the counter __ */
            n = --sx->ns[s];

            /* __ put the last part. in the hole __ */
            sp[s][m].i = sp[s][n].i;

            /* __ loop on the 3 directions __ */
            for (l = 0; l < 3; l++)
            {
                /* __ put the last part. in the hole __ */
                sp[s][m].r[l] = sp[s][n].r[l];
                sp[s][m].s[l] = sp[s][n].s[l];
                sp[s][m].v[l] = sp[s][n].v[l];
                sp[s][m].w[l] = sp[s][n].w[l];
                sp[s][m].b[l] = sp[s][n].b[l];
            }

            /* __ test swapped (last in stack) part. __ */
            m--;
        }
    }

    /* __ send/receive the part. : loop on all the neighbors __ */
    for (t = 0; t < 27; t++)
    {
        /* __ except myself __ */
        if (t != 13)
        {
            /* __ send/receive the # of part. __ */
            MPI_Sendrecv(&vc[t], 1, MPI_INT, sx->nt[t], t,
                         &wc[t], 1, MPI_INT, sx->nf[t], t, MPI_COMM_WORLD, &st);

            /* __ send/receive the part. __ */
            MPI_Sendrecv(ws[t], uc[t], MPI_PACKED, sx->nt[t], t,
                         wr[t], ps   , MPI_PACKED, sx->nf[t], t, MPI_COMM_WORLD, &st);

            /* __ set the size of the packet __ */
            uc[t] = 0;

            /* __ index of the first added part __ */
            m = sx->ns[s];

            /* __ loop to add the part. __ */
            while (wc[t] > 0 && sx->nf[t] != MPI_PROC_NULL)
            {
                /* __ unpack data associated to the part. __ */
                MPI_Unpack(wr[t], ps, &uc[t], &sp[s][m].i, 1, MPI_INT   , MPI_COMM_WORLD);
                MPI_Unpack(wr[t], ps, &uc[t],  sp[s][m].r, 3, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Unpack(wr[t], ps, &uc[t],  sp[s][m].s, 3, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Unpack(wr[t], ps, &uc[t],  sp[s][m].v, 3, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Unpack(wr[t], ps, &uc[t],  sp[s][m].w, 3, MPI_DOUBLE, MPI_COMM_WORLD);
                MPI_Unpack(wr[t], ps, &uc[t],  sp[s][m].b, 3, MPI_INT   , MPI_COMM_WORLD);

                /* __ update the # of part. __ */
                m = ++sx->ns[s];

                /* __ update the counter of # of part. to unpack __ */
                wc[t]--;
            }
        }
    }

    /* __ loop on all the neighbors __ */
    for (t = 0; t < 27; t++)
    {
        /* __ clean up the pointers __ */
        free(ws[t]);
        free(wr[t]);
    }

}


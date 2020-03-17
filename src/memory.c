
#include <stdlib.h>
#include <stdio.h>

#include "structures.h"


/* __ memory allocation for a stf structure _________________________________ */
void memorystf(struct stf *sf, struct sti si)
{
int nn1;
int l, s;


/* __ # of cells on g1 __ */
nn1 = (si.n[0]+1)*(si.n[1]+1)*(si.n[2]+1);

/* __ memory allocation for fields depending on specie __ */
for (s = 0; s < NS+1; s++)
    {
    /* __ specie density __ */
    sf->n[s] = (float *) malloc(nn1*sizeof(float));

    /* __ specie velocity __ */
    for (l = 0; l < 3; l++)
        {
        sf->v[s][l] = (float *) malloc(nn1*sizeof(float));
        }

    /* __ specie pressure __ */
    for (l = 0; l < 6; l++)
        {
        sf->p[s][l] = (float *) malloc(nn1*sizeof(float));
        }
    }

/* __ memory allocation for electromagnetic & fluid fields __ */
for (l = 0; l < 3; l++)
    {
    sf->b[l] = (float *) malloc(nn1*sizeof(float));
    sf->e[l] = (float *) malloc(nn1*sizeof(float));
    sf->a[l] = (float *) malloc(nn1*sizeof(float));
    sf->j[l] = (float *) malloc(nn1*sizeof(float));
    sf->vi[l] = (float *) malloc(nn1*sizeof(float));
    }

}


/* _____ memory allocation for a stq structure ______________________________ */
void memorystq(struct stq *sq, struct sti si)
{
int l, s;


/* __ memory allocation 1d fields of sq structure : loop on the first direction __ */
for (s = 1; s < NS+1; s++)
    {
    for (l = 0; l < 3; l++)
        {
        /* __ for the 1d fields __ */
        sq->r[s][l] = (float *) malloc(si.nm*sizeof(float));
        sq->v[s][l] = (float *) malloc(si.nm*sizeof(float));
        }
    }

}


/* _____ memory allocation for a stm structure ______________________________ */
void memorystm(struct stm *sm, int no)
{
int l, s;


/* __ memory allocation for field depending on specie __ */
for (s = 0; s < NS+1; s++)
    {
    /* __ specie density __ */
    sm->n[s] = (float *) malloc(no*sizeof(float));

    /* __ specie velocity __ */
    for (l = 0; l < 3; l++)
        {
        sm->v[s][l] = (float *) malloc(no*sizeof(float));
        }

    /* __ specie pressure __ */
    for (l = 0; l < 6; l++)
        {
        sm->p[s][l] = (float *) malloc(no*sizeof(float));
        }
    }

/* __ memory allocation for electromagnetic & fluid fields __ */
for (l = 0; l < 3; l++)
    {
    /* __ for the 1d fields __ */
    sm->r[l] = (float *) malloc(no*sizeof(float));
    sm->w[l] = (float *) malloc(no*sizeof(float));
    sm->b[l] = (float *) malloc(no*sizeof(float));
    sm->e[l] = (float *) malloc(no*sizeof(float));
    sm->j[l] = (float *) malloc(no*sizeof(float));
    sm->vi[l] = (float *) malloc(no*sizeof(float));
    }

}


/* __ free a stf structure __________________________________________________ */
void freestf(struct stf *sf)
{
int l, s;


/* __ free the fields depending on specie __ */
for (s = 0; s < NS+1; s++)
    {
    /* __ specie density __ */
    free(sf->n[s]);

    /* __ specie velocity __ */
    for (l = 0; l < 3; l++)
        {
        free(sf->v[s][l]);
        }

    /* __ specie pressure __ */
    for (l = 0; l < 6; l++)
        {
        free(sf->p[s][l]);
        }
    }

/* __ free the electromagnetic & fluid fields __ */
for (l = 0; l < 3; l++)
    {
    free(sf->b[l]);
    free(sf->e[l]);
    free(sf->a[l]);
    free(sf->j[l]);
    free(sf->vi[l]);
    }

}


/* __ free a stq structure  _________________________________________________ */
void freestq(struct stq *sq)
{
int l, s;


/* __ free 1d fields of sq structure : loop on the first direction __ */
for (s = 1; s < NS+1; s++)
    {
    for (l = 0; l < 3; l++)
        {
        /* __ for the 1d fields __ */
        free(sq->r[s][l]);
        free(sq->v[s][l]);
        }
    }

}


/* __ free a stm structure  _________________________________________________ */
void freestm(struct stm *sm)
{
int l, s;


/* __ free 1d fields of sm structure : loop on the first direction __ */
for (l = 0; l < 3; l++)
    {
    /* __ for the 1d fields __ */
    free(sm->r[l]);
    free(sm->w[l]);
    free(sm->b[l]);
    free(sm->e[l]);
    free(sm->j[l]);
    free(sm->vi[l]);
    }

/* __ free the moments depending on specie __ */
for (s = 0; s < NS+1; s++)
    {
    /* __ specie density __ */
    free(sm->n[s]);

    /* __ specie velocity __ */
    for (l = 0; l < 3; l++)
        {
        free(sm->v[s][l]);
        }

    /* __ specie pressure __ */
    for (l = 0; l < 6; l++)
        {
        free(sm->p[s][l]);
        }
    }

}


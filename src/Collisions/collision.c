
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "structures.h"
#include "defines.h"
#include "collision.h"
#include "collisionIonIon.h"
#include "ions.h"
#include "ghosts.h"
#include "hecklebc.h"

#define SWAP(type, a, b) do { type t=(a); (a)=(b); (b)=t; } while (0)



void collision1(Collision sc,
                STI si,
                STX sx,
                struct std *sd,
                Grid0 *s0,
                ST1 *s1,
                ST2 *s2,
                Particle *sp[NS+1],
                Ghosts *ghosts,
                HeckleBC *hbc,
                int ipc)
{


    ionIonCollision(sc,
                    si,
                    sx,
                    sd,
                    s0,
                    s1,
                    s2,
                    sp,
                    ghosts,
                    hbc,
                    ipc);

}



/* __ read the initial parameters to set Collision structure ________________ */
Collision setCollision(struct sti si,
                       struct stx sx,
                       struct std *sd,
                       char *dir)
{
    Collision sc;


    if (si.coll != 0) {
        /* __ read the "collision.txt" file __ */
        readCollision(&sc, dir);

        if (sx.r == 0 && si.coll != 0) {
            printf("\n");
            printf("________________ parameters for ion-ion collisions ___\n");
            printf("ion-ion coll. nu :%12.6lf\n", sc.nuIonIon);
            printf("coulomb log      :%12.6lf\n", sc.coulombLog);
            printf("\n");
        }
    }

//  /* __ set initial values for collision diagnostics __ */
//  if (si.coll == 1 || si.coll == 3 || si.coll ==  5 || si.coll == 7) {
//      sd->lii = sc.coulombLog;
//      sd->niimin = sc.nuIonIon;
//  }

    switch(si.coll) {
        case 1 : // ion-ion collisions
            collision = collision1;
        break;

        case 2 : // ion-electron collisions
            printf("not implemented yet !\n");
            exit(EXIT_FAILURE);
        break;

        case 4 : // ion-neutral collisions
            printf("not implemented yet !\n");
            exit(EXIT_FAILURE);
        break;

        case 3 : // ion-ion & ion-electron collisions
            printf("not implemented yet !\n");
            exit(EXIT_FAILURE);
        break;

        case 5 : // ion-ion & ion-neutral collisions
            printf("not implemented yet !\n");
            exit(EXIT_FAILURE);
        break;

        case 6 : // ion-electron & ion-neutral collisions
            printf("not implemented yet !\n");
            exit(EXIT_FAILURE);
        break;

        case 7 : // ion-ion, ion-electron & ion-neutral collisions
            printf("not implemented yet !\n");
            exit(EXIT_FAILURE);
        break;
    }

    return sc;

}



/* __ read the collision.txt file __ */
void readCollision(Collision *sc,
                   char *dir)
{
    FILE *fp;
    int ir;
    char sfh[80];
    char junk[80];

    /* __ build the file name __ */
    strcpy(sfh, dir);
    strcat(sfh, "collision.txt");

    /* __ open the collision.txt file __ */
    fp = fopen(sfh, "r+");
    if (fp == NULL) printf("\n\nproblem in opening file %s\n", sfh);

    /* collision frequency */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(sc->nuIonIon), &ir);

    /* coulomb logarithm */
    fscanstr(__FILE__, __LINE__, fp, 0, junk, &ir);
    fscandbl(__FILE__, __LINE__, fp, 0, &(sc->coulombLog), &ir);

    /* __ close the file __ */
    fclose(fp);

    return;
}


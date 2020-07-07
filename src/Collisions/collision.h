
#ifndef __COLLISION_H__
#define __COLLISION_H__


#include<structures.h>
#include<particle.h>
#include<particlebc.h>
#include <ghosts.h>
#include <hecklebc.h>



/*---------------------------------------------------------------------------
  collision()
  ---------------------------------------------------------------------------
  AIM : wrapper for all kind of collisions
 ---------------------------------------------------------------------------*/
void collisionModel (Collision sc,
                     STI si,
                     STX sx,
                     struct std *sd,
                     Grid0 *s0,
                     ST1 *s1,
                     ST2 *s2,
                     Particle *sp[NS+1],
                     Ghosts *ghosts,
                     HeckleBC *hbc,
                     int ipc);


//void (*collision)(Collision sc,
//                  STI si,
//                  STX sx,
//                  struct std *sd,
//                  Grid0 *s0,
//                  ST1 *s1,
//                  ST2 *s2,
//                  Particle *sp[NS+1],
//                  Ghosts *ghosts,
//                  HeckleBC *hbc,
//                  int ipc);


/*---------------------------------------------------------------------------
  setCollision()
  ---------------------------------------------------------------------------
  AIM :set collision structure sc
 ---------------------------------------------------------------------------*/
Collision setCollision(struct sti si,
                       struct stx sx,
                       struct std *sd,
                       char *dir);



/*---------------------------------------------------------------------------
  readCollision()
  ---------------------------------------------------------------------------
  AIM : read the collision.txt file
 ---------------------------------------------------------------------------*/
void readCollision(Collision *sc,
                   char *dir);


/*---------------------------------------------------------------------------
  collision1()
  ---------------------------------------------------------------------------
  AIM : wrapper for ion-ion collision module
 ---------------------------------------------------------------------------*/
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
                int ipc);

#endif // __COLLISION_H__


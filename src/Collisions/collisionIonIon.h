
#ifndef __COLLISION_ION_ION_H__
#define __COLLISION_ION_ION_H__


#include<structures.h>
#include<particle.h>
#include<particlebc.h>
#include <ghosts.h>
#include <hecklebc.h>


/*---------------------------------------------------------------------------
  coullog()
  ---------------------------------------------------------------------------
  AIM : calculate the coulomb logarithm
 ---------------------------------------------------------------------------*/
double coulombLog(STI si,
                  Collision sc,
                  ST2 sf,
                  double tempi[NS+1],
                  int s1,
                  int s2);



/*---------------------------------------------------------------------------
  collfrequency()
  ---------------------------------------------------------------------------
  AIM : calculate the collision frequency (Rutherford plasma)
 ---------------------------------------------------------------------------*/
double ionIonCollisionFrequency(STX sx,
                                Collision sc,
                                ST2 sf,
                                STI si,
                                struct std *sd,
                                double col,
                                Particle part1,
                                Particle part2,
                                double tempi[NS+1],
                                int s1,
                                int s2,
                                int ipc);



/*---------------------------------------------------------------------------
  takizuka()
  ---------------------------------------------------------------------------
  AIM : calculate the scattering angle (Takizuka and Abe)
 ---------------------------------------------------------------------------*/
void setDeviationAngles(STI si,
                        ST2 s2,
                        STX sx,
                        Collision sc,
                        struct std *sd,
                        double *sinthe,
                        double *costhe,
                        double *phi,
                        double col,
                        Particle part1,
                        Particle part2,
                        double dt,
                        double tempi[NS+1],
                        int specie1,
                        int specie2,
                        int ipc);



/*---------------------------------------------------------------------------
  vchange()
  ---------------------------------------------------------------------------
  AIM : change the velocities of a pair of particles
 ---------------------------------------------------------------------------*/
void velocityScattering(STI si,
                        STX sx,
                        Collision sc,
                        Particle *part1,
                        Particle *part2,
                        double sinthe,
                        double costhe,
                        double phi,
                        int specie1,
                        int specie2,
                        int ipc);



/*---------------------------------------------------------------------------
  collcell()
  ---------------------------------------------------------------------------
  AIM : compute the collisions over the domain
 ---------------------------------------------------------------------------*/
void ionIonCellCollisions(Collision sc,
                          Grid0 s0,
                          ST2 s2,
                          STI si,
                          STX sx,
                          struct std *sd,
                          Particle *sp[NS+1],
                          double tempi[NS+1],
                          int ipc);



/*---------------------------------------------------------------------------
  ionIonCollision()
  ---------------------------------------------------------------------------
  AIM : scatter ions because ion-ion collisions
 ---------------------------------------------------------------------------*/
void ionIonCollision(Collision sc,
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


/*---------------------------------------------------------------------------
  collision()
  ---------------------------------------------------------------------------
  AIM : wrapper for all kind of collisions
 ---------------------------------------------------------------------------*/
void (*collision)(Collision sc,
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

#endif // __COLLISION_ION_ION_H__


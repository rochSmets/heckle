
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "structures.h"
#include "defines.h"
#include "collision.h"
#include "ions.h"
#include "ghosts.h"
#include "hecklebc.h"

#define SWAP(type, a, b) do { type t=(a); (a)=(b); (b)=t; } while (0)


/* __ function for the coulomb logarithm __ */
double coulombLog(STI si,
               Collision sc,
               ST2 s2,
               double ionTemp[NS+1],
               int specie1,
               int specie2)
{
   double lambda;
   (void) s2;


   /* __ calculation of the log __ */
   lambda = sc.coulombLog - log(si.qs[specie1]*si.qs[specie2]
                              *(si.ms[specie1]+si.ms[specie2])
                              *sqrt((si.ns[specie1]*si.qs[specie1]*si.qs[specie1]/ionTemp[specie1])
                                   +(si.ns[specie2]*si.qs[specie2]*si.qs[specie2]/ionTemp[specie2]))
                              /(si.ms[specie1]*ionTemp[specie2]+si.ms[specie2]*ionTemp[specie1]));

   // truc a verifier !!!!!!!! la théorie n'est plus valide si lambda < 2
   // verifier aussi le : remplacer lambda par 2 si nan
   //if (lambda < 2 || lambda != lambda) {
   //   printf("warning : lambda < 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
   //   lambda = 2;
   //}

   return lambda;

}


/* __ return the collision frequency (Rutherford plasma) __ */
double ionIonCollisionFrequency(STX sx,
                                Collision sc,
                                ST2 s2,
                                STI si,
                                struct std *sd,
                                double nuIonIon,
                                Particle part1,
                                Particle part2,
                                double ionTemp[NS+1],
                                int specie1,
                                int specie2,
                                int ipc)
{
    double nu;      // collision frequency
    double u;       // relative velocity
    double mu;      // relative mass
    double lambda;  // coulomb logarithm
    double nMax;    // maximum density


    nMax = (s2.ns[specie1] < s2.ns[specie2]) ? s2.ns[specie1] : s2.ns[specie2];

    switch (ipc) {
        case 0 : /* __ predictor : relavite velocity depends on v __ */
           u = sqrt((part1.v[0]-part2.v[0])*(part1.v[0]-part2.v[0])
                   +(part1.v[1]-part2.v[1])*(part1.v[1]-part2.v[1])
                   +(part1.v[2]-part2.v[2])*(part1.v[2]-part2.v[2]));
        break;

        case 1 : /* __ corrector : relative velocity depends on w __ */
           u = sqrt((part1.w[0]-part2.w[0])*(part1.w[0]-part2.w[0])
                   +(part1.w[1]-part2.w[1])*(part1.w[1]-part2.w[1])
                   +(part1.w[2]-part2.w[2])*(part1.w[2]-part2.w[2]));
        break;

        /* __ no reason to get there __ */
        default :
           IAMDEAD(sx.r);
    }

    /* __ safegard for the value of u __ */
    u = (u > EPS6) ? u : EPS6;

    /* __ reduced mass __ */
    mu = (si.ms[specie1]*si.ms[specie2])/(si.ms[specie1]+si.ms[specie2]);

    lambda = coulombLog(si, sc, s2, ionTemp, specie1, specie2);

    nu = (nuIonIon*si.qs[specie1]*si.qs[specie1]*si.qs[specie2]*si.qs[specie2]*nMax*lambda)/(mu*mu*u*u*u);

    /* __ set lambda, min(nu) & max(nu) fir ions in sd struct __ */
    if (lambda < sd->lii) {
        sd->lii = lambda;
    }

    if (nu < sd->niimin) {
        sd->niimin = nu;
    }

    return nu;

}


/* __ function to calculate the scattering angle (Takizuka and Abe) __ */
void setDeviationAngles(STI si,
                        ST2 s2,
                        STX sx,
                        Collision sc,
                        struct std *sd,
                        double *sinTheta,
                        double *oneMinusCosTheta,
                        double *phi,
                        double nuIonIon,
                        Particle part1,
                        Particle part2,
                        double dt,
                        double ionTemp[NS+1],
                        int specie1,
                        int specie2,
                        int ipc)
{
    double nu;        // collision frequency
    double r1, r2;    // random numbers for the box and muller algo
    double delta;     // variance of the gaussian for the angle of diffusion
    double r3;


    /* __ calculation of the collision frequency __ */
    nu = ionIonCollisionFrequency(sx,
                                  sc,
                                  s2,
                                  si,
                                  sd,
                                  nuIonIon,
                                  part1,
                                  part2,
                                  ionTemp,
                                  specie1,
                                  specie2,
                                  ipc)*dt;

    /* __ computation of the delta parameter __ */
    r1 = RNM;
    r1 = (r1 > EPS8) ? r1 : EPS8;
    r1 = -2*log(r1);
    r2 = 2*PI*RNM;

    delta = sqrt(nu*r1)*cos(r2);

    /* __ scattering angles __ */
    if (nu < 1) {
        *sinTheta = 2*delta/(1+(delta*delta));
        *oneMinusCosTheta = 2*delta*delta/(1+(delta*delta));
    }

    else {
        r3 = 2*PI*RNM;

        *sinTheta = sin(r3);
        *oneMinusCosTheta = 1-cos(r3);
    }

    *phi = 2*PI*RNM;

}



/* __ function to change the velocities of a pair of particles __ */
void velocityScattering(STI si,
                        STX sx,
                        Collision sc,
                        Particle *part1,
                        Particle *part2,
                        double sinTheta,
                        double oneMinusCosTheta,
                        double phi,
                        int specie1,
                        int specie2,
                        int ipc)
{
    double u[3];               // 3 components of the relative velocity
    double uPerp;              // perp component of u (in xy plan
    double uModulus;           // perpendicular component of u, module of u
    double mu;                 // reduced mass
    double dx, dy, dz;         // variations of the relative velocity in each direction
    double r0;                 // random number
    int z1;                    // index for part 1 : 1 if collision, else 0
    int z2;                    // index for part 1 : 1 if collision, else 0


    /* __ reduced mass definition __ */
    mu = si.ms[specie1]*si.ms[specie2]/(si.ms[specie1]+si.ms[specie2]);

    /* __ see if the particles experience a diffusion depending on their weights __ */
    r0 = RNM;
    z1 = (r0 < si.ws[specie2]/si.ws[specie1]) ? 1 : 0;
    z2 = (r0 < si.ws[specie1]/si.ws[specie2]) ? 1 : 0;

    switch (ipc) {
        case 0 : /* __ predictor : change of velocities __ */
            /* __ parameters for the collision __ */
            u[0] = part1->v[0]-part2->v[0];
            u[1] = part1->v[1]-part2->v[1];
            u[2] = part1->v[2]-part2->v[2];

            uPerp = sqrt(u[0]*u[0]+u[1]*u[1]);
            uModulus = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

            if (uPerp == 0) {
                dx =  uModulus*sinTheta*cos(phi);
                dy =  uModulus*sinTheta*sin(phi);
                dz = -uModulus*oneMinusCosTheta;
            }

            else {
                dx = (u[0]/uPerp)*u[2]*sinTheta*cos(phi)-(u[1]/uPerp)*uModulus*sinTheta*sin(phi)-u[0]*(oneMinusCosTheta);
                dy = (u[1]/uPerp)*u[2]*sinTheta*cos(phi)+(u[0]/uPerp)*uModulus*sinTheta*sin(phi)-u[1]*(oneMinusCosTheta);
                dz = -uPerp*sinTheta*cos(phi)-u[2]*(oneMinusCosTheta);
            }

            /* __ change the velocities : x direction __ */
            part1->v[0] = part1->v[0]+z1*(mu/si.ms[specie1])*dx;
            part2->v[0] = part2->v[0]-z2*(mu/si.ms[specie2])*dx;

            /* __ change the velocities : y direction __ */
            part1->v[1] = part1->v[1]+z1*(mu/si.ms[specie1])*dy;
            part2->v[1] = part2->v[1]-z2*(mu/si.ms[specie2])*dy;

            /* __ change the velocities : z direction __ */
            part1->v[2] = part1->v[2]+z1*(mu/si.ms[specie1])*dz;
            part2->v[2] = part2->v[2]-z2*(mu/si.ms[specie2])*dz;
        break;

    case 1 : /* __ corrector : change of velocities __ */
        /* __ parameters for the collision __ */
        u[0] = part1->w[0]-part2->w[0];
        u[1] = part1->w[1]-part2->w[1];
        u[2] = part1->w[2]-part2->w[2];

        uPerp = sqrt(u[0]*u[0]+u[1]*u[1]);
        uModulus = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

        if (uPerp == 0) {
            dx=  uModulus*sinTheta*cos(phi);
            dy=  uModulus*sinTheta*sin(phi);
            dz= -uModulus*oneMinusCosTheta;
        }

        else {
            dx = (u[0]/uPerp)*u[2]*sinTheta*cos(phi)-(u[1]/uPerp)*uModulus*sinTheta*sin(phi)-u[0]*(oneMinusCosTheta);
            dy = (u[1]/uPerp)*u[2]*sinTheta*cos(phi)+(u[0]/uPerp)*uModulus*sinTheta*sin(phi)-u[1]*(oneMinusCosTheta);
            dz = -uPerp*sinTheta*cos(phi)-u[2]*(oneMinusCosTheta);
        }

        /* __ change the velocities : x direction __ */
        part1->w[0] = part1->w[0]+z1*(mu/si.ms[specie1])*dx;
        part2->w[0] = part2->w[0]-z2*(mu/si.ms[specie2])*dx;

        /* __ change the velocities : y direction __ */
        part1->w[1] = part1->w[1]+z1*(mu/si.ms[specie1])*dy;
        part2->w[1] = part2->w[1]-z2*(mu/si.ms[specie2])*dy;

        /* __ change the velocities : z direction __ */
        part1->w[2] = part1->w[2]+z1*(mu/si.ms[specie1])*dz;
        part2->w[2] = part2->w[2]-z2*(mu/si.ms[specie2])*dz;
        break;

    /* __ no reason to get there __ */
    default : IAMDEAD(sx.r);
    }

}


/* __ compute the ion-ion collisions in the given cell __ */
void ionIonCellCollisions(Collision sc,
                          Grid0 s0,
                          ST2 s2,
                          STI si,
                          STX sx,
                          struct std *sd,
                          Particle *sp[NS+1],
                          double ionTemp[NS+1],
                          int ipc)
{
    int specie1, specie2;      // indices for the species
    int pa, p, q, o;           // indices for the particules
    int r0;             // to randomize the particle indices array
    double sinTheta; // polar angle for scattering
    double oneMinusCosTheta;   // polar angle for scattering
    double phi;// azimuthal angle for scattering
    int i;                     // number of times a particle is selected for the inter-species collisions
    double r;
    double dtweit = 1.;        // correction of the frequency depending on the statistical weights


    /* __ randomize the particle array __ */
    for (specie1 = 1; specie1 < NS+1; specie1++) {

        int numSpecie1 = s0.npart[specie1];

        for (pa = 0; pa < numSpecie1; pa++) {
            do {
                r0 = RNM*numSpecie1;
            } while (r0 >= numSpecie1);

            SWAP(int, s0.ipart[specie1][r0], s0.ipart[specie1][pa]);
        }
    }

    /* __ collisions : loop overs the species for the first particle __ */
    for (specie1 = 1; specie1 < NS+1; specie1++) {

        /* __ INTRA-species collisions __ */
        int numSpecie1 = s0.npart[specie1];

        /* __ even number of particles __ */
        if (numSpecie1%2 == 0) {
            for (pa = 0; pa < numSpecie1/2; pa++) {
//if (2*pa >= numSpecie1 || (2*pa)+1 >= numSpecie1) printf("aaaaaaa\n");
                p = s0.ipart[specie1][ 2*pa   ];
                q = s0.ipart[specie1][(2*pa)+1];

                setDeviationAngles(si,
                                   s2,
                                   sx,
                                   sc,
                                   sd,
                                   &sinTheta,
                                   &oneMinusCosTheta,
                                   &phi,
                                   sc.nuIonIon,
                                   sp[specie1][p],
                                   sp[specie1][q],
                                   si.ts,
                                   ionTemp,
                                   specie1,
                                   specie1,
                                   ipc);

                velocityScattering(si,
                                   sx,
                                   sc,
                                   &sp[specie1][p],
                                   &sp[specie1][q],
                                   sinTheta,
                                   oneMinusCosTheta,
                                   phi,
                                   specie1,
                                   specie1,
                                   ipc);
            }
        }


        /* __ odd number of particles __ */
        else {
            for (pa = 0; pa < (numSpecie1/2)-1; pa++) {
//if (2*pa >= numSpecie1 || (2*pa)+1 >= numSpecie1) printf("bbbbbbb\n");
                p = s0.ipart[specie1][ 2*pa   ];
                q = s0.ipart[specie1][(2*pa)+1];

                setDeviationAngles(si,
                                   s2,
                                   sx,
                                   sc,
                                   sd,
                                   &sinTheta,
                                   &oneMinusCosTheta,
                                   &phi,
                                   sc.nuIonIon,
                                   sp[specie1][p],
                                   sp[specie1][q],
                                   si.ts,
                                   ionTemp,
                                   specie1,
                                   specie1,
                                   ipc);

                velocityScattering(si,
                                   sx,
                                   sc,
                                   &sp[specie1][p],
                                   &sp[specie1][q],
                                   sinTheta,
                                   oneMinusCosTheta,
                                   phi,
                                   specie1,
                                   specie1,
                                   ipc);
             }

            /* __ last three particles __ */

            /* __ pair 0 and 1 __ */
            p = numSpecie1-3;
            q = numSpecie1-2;
            setDeviationAngles(si,
                               s2,
                               sx,
                               sc,
                               sd,
                               &sinTheta,
                               &oneMinusCosTheta,
                               &phi,
                               0.5*sc.nuIonIon,
                               sp[specie1][p],
                               sp[specie1][q],
                               si.ts,
                               ionTemp,
                               specie1,
                               specie1,
                               ipc);

            velocityScattering(si,
                               sx,
                               sc,
                               &sp[specie1][p],
                               &sp[specie1][q],
                               sinTheta,
                               oneMinusCosTheta,
                               phi,
                               specie1,
                               specie1,
                               ipc);

            /* __ pair 0 and 2 __ */
            p = numSpecie1-3;
            q = numSpecie1-1;
            setDeviationAngles(si,
                               s2,
                               sx,
                               sc,
                               sd,
                               &sinTheta,
                               &oneMinusCosTheta,
                               &phi,
                               0.5*sc.nuIonIon,
                               sp[specie1][p],
                               sp[specie1][q],
                               si.ts,
                               ionTemp,
                               specie1,
                               specie1,
                               ipc);

            velocityScattering(si,
                               sx,
                               sc,
                               &sp[specie1][p],
                               &sp[specie1][q],
                               sinTheta,
                               oneMinusCosTheta,
                               phi,
                               specie1,
                               specie1,
                               ipc);

            /* __ pair 1 and 2 __ */
            p = numSpecie1-2;
            q = numSpecie1-1;
            setDeviationAngles(si,
                               s2,
                               sx,
                               sc,
                               sd,
                               &sinTheta,
                               &oneMinusCosTheta,
                               &phi,
                               0.5*sc.nuIonIon,
                               sp[specie1][p],
                               sp[specie1][q],
                               si.ts,
                               ionTemp,
                               specie1,
                               specie1,
                               ipc);

            velocityScattering(si,
                               sx,
                               sc,
                               &sp[specie1][p],
                               &sp[specie1][q],
                               sinTheta,
                               oneMinusCosTheta,
                               phi,
                               specie1,
                               specie1,
                               ipc);
        }

        /* __ INTER-species collisions __ */

        /* __ loop over the others species __ */
        for (specie2 = specie1+1; specie2 < NS+1; specie2++) {

            int numSpecie2 = s0.npart[specie2];

            /* __ case N_a > N_b __ */
            if (numSpecie1 >= numSpecie2) {

                /* __ max (W_1, W_2) / W_2 __ */
                dtweit = (si.ws[specie1] >= si.ws[specie2]) ? si.ws[specie1]/si.ws[specie2] : si.ws[specie2]/si.ws[specie2];

                // i is integer part & r the rest of ratio between number of species
                i = floor((double)numSpecie1/(double)numSpecie2);
                r =      ((double)numSpecie1/(double)numSpecie2)-i;

                /* __ first group of part of specie 2 __ */
                int firstGroupSpecie2 = (int)round(r*numSpecie2);

                /* __ first group __ */
                for (pa = 0; pa < firstGroupSpecie2; pa++) {
//if (pa >= numSpecie2) printf("ccccccc\n");
                    p = s0.ipart[specie2][pa];

                    /* __ each particle (specie2) selected i+1 times __ */
                    for (o = 0; o < i+1; o++) {
//if ((i+1)*pa+o >= numSpecie1) printf("ddddddd\n");
                    q = s0.ipart[specie1][(i+1)*pa+o];

                    setDeviationAngles(si,
                                       s2,
                                       sx,
                                       sc,
                                       sd,
                                       &sinTheta,
                                       &oneMinusCosTheta,
                                       &phi,
                                       dtweit*sc.nuIonIon,
                                       sp[specie2][p],
                                       sp[specie1][q],
                                       si.ts,
                                       ionTemp,
                                       specie2,
                                       specie1,
                                       ipc);

                    velocityScattering(si,
                                       sx,
                                       sc,
                                       &sp[specie2][p],
                                       &sp[specie1][q],
                                       sinTheta,
                                       oneMinusCosTheta,
                                       phi,
                                       specie2,
                                       specie1,
                                       ipc);
                    }
                }

                /* __ second group __ */
                for (pa = firstGroupSpecie2; pa < numSpecie2; pa++) {
//if (pa >= numSpecie2) printf("eeeeeee\n");
                    p = s0.ipart[specie2][pa];

                    /* __ each particle (specie2) selected i times __ */
                    for (o = 0; o < i; o++) {
//if (firstGroupSpecie2*(i+1)+(pa-firstGroupSpecie2)*i+o >= numSpecie1) printf("fffffff\n");
                        q = s0.ipart[specie1][firstGroupSpecie2*(i+1)+(pa-firstGroupSpecie2)*i+o];

                        setDeviationAngles(si,
                                           s2,
                                           sx,
                                           sc,
                                           sd,
                                           &sinTheta,
                                           &oneMinusCosTheta,
                                           &phi,
                                           dtweit*sc.nuIonIon,
                                           sp[specie2][p],
                                           sp[specie1][q],
                                           si.ts,
                                           ionTemp,
                                           specie2,
                                           specie1,
                                           ipc);

                        velocityScattering(si,
                                           sx,
                                           sc,
                                           &sp[specie2][p],
                                           &sp[specie1][q],
                                           sinTheta,
                                           oneMinusCosTheta,
                                           phi,
                                           specie2,
                                           specie1,
                                           ipc);
                    }
                }
            }

            else {
                /* __ max (W_1, W_2) / W_1 __ */
                dtweit = (si.ws[specie1] >= si.ws[specie2]) ? si.ws[specie1]/si.ws[specie1] : si.ws[specie2]/si.ws[specie1];

                i = floor((double)numSpecie2/(double)numSpecie1);
                r =      ((double)numSpecie2/(double)numSpecie1)-i;

                /* __ first group of part of specie 1 __ */
                int firstGroupSpecie1 = (int)round(r*numSpecie1);

                /* __ first group __ */
                for (pa = 0; pa < firstGroupSpecie1; pa++) {
//if (pa >= numSpecie1) printf("ggggggg\n");
                    p = s0.ipart[specie1][pa];

                    /* __ each particle (specie1) selected i+1 times __ */
                    for (o = 0; o < i+1; o++) {
//if ((i+1)*pa+o >= numSpecie2) printf("hhhhhhh\n");
                        q = s0.ipart[specie2][(i+1)*pa+o];

                        setDeviationAngles(si,
                                           s2,
                                           sx,
                                           sc,
                                           sd,
                                           &sinTheta,
                                           &oneMinusCosTheta,
                                           &phi,
                                           dtweit*sc.nuIonIon,
                                           sp[specie1][p],
                                           sp[specie2][q],
                                           si.ts,
                                           ionTemp,
                                           specie1,
                                           specie2,
                                           ipc);

                        velocityScattering(si,
                                           sx,
                                           sc,
                                           &sp[specie1][p],
                                           &sp[specie2][q],
                                           sinTheta,
                                           oneMinusCosTheta,
                                           phi,
                                           specie1,
                                           specie2,
                                           ipc);
                    }
                }
                /* __ second group __ */
                for (pa = firstGroupSpecie1; pa < numSpecie1; pa++) {
//if (pa >= numSpecie1) printf("iiiiiii\n");
                    p = s0.ipart[specie1][pa];

                    /* __ each particle selected i times __ */
                    for (o = 0; o < i; o++) {
//if (firstGroupSpecie1*(i+1)+(pa-firstGroupSpecie1)*i+o >= numSpecie2) printf("jjjjjjj\n");
                        q = s0.ipart[specie2][firstGroupSpecie1*(i+1)+(pa-firstGroupSpecie1)*i+o];

                        setDeviationAngles(si,
                                           s2,
                                           sx,
                                           sc,
                                           sd,
                                           &sinTheta,
                                           &oneMinusCosTheta,
                                           &phi,
                                           dtweit*sc.nuIonIon,
                                           sp[specie1][p],
                                           sp[specie2][q],
                                           si.ts,
                                           ionTemp,
                                           specie1,
                                           specie2,
                                           ipc);

                        velocityScattering(si,
                                           sx,
                                           sc,
                                           &sp[specie1][p],
                                           &sp[specie2][q],
                                           sinTheta,
                                           oneMinusCosTheta,
                                           phi,
                                           specie1,
                                           specie2,
                                           ipc);
                    }
                }
            }
        }
    }

}


/* __ function to sweep the domain and call the collisions __ */
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
                     int ipc)
{
    double xw, yw, zw;// particle position in the subdomain
    double ionTemp[NS+1]; // ion temperatures
    int i, j, k;      // indices for the grid points
    int s, m;         // incides for the specie and the particles
    int ijk0;         // index on Grid0
    int ijk;          // index for the grid point in 1D
    int ijk1;         // index for the grid point in the ST2 grid
    int *oldNumPerCell[NS+1];// former number of part per cell
    int n0;


    /* __ # of grid point on grid0 __ */
    n0 = sx.n[0]*sx.n[1]*sx.n[2];

    /* __ set lambda & nu from sd to init values __ */
    //sd->lii = sc.coulombLog;
    //sd->niimin = sc.nuIonIon;
    sd->lii = 100.0;
    sd->niimin = 1.0/si.ts;

    /* __ keep memory of the former num of part in each grid0 cell __ */
    for (s = 1; s < NS+1; s++) {
        oldNumPerCell[s] = malloc(n0*sizeof(int));

        for (ijk0 = 0; ijk0 < n0; ijk0++) {
            oldNumPerCell[s][ijk0] = s0[ijk0].npart[s];
            s0[ijk0].npart[s] = 0;

            /* __ keep the same amount of memory as at previous time step __ */
            s0[ijk0].ipart[s] = malloc(oldNumPerCell[s][ijk0]*sizeof(int));
        }
    }

    /* __ fill the list of particles for each cell (s0 structure) : loop on the part __ */
    for (s = 1; s < NS+1; s++) {

        /* __ loop on the part of specie s __ */
        for (m = 0; m < sx.ns[s]; m++) {

            /* __ corrector __ */
            if (ipc == 1) {
                /* __ keep memory of older part velocity __ */
                sp[s][m].w[0] = sp[s][m].v[0];
                sp[s][m].w[1] = sp[s][m].v[1];
                sp[s][m].w[2] = sp[s][m].v[2];
            }

            /* __ part "position" in the subdomain __ */
            xw = (sp[s][m].r[0])/si.dl[0]-sx.i0[0];
            yw = (sp[s][m].r[1])/si.dl[1]-sx.i0[1];
            zw = (sp[s][m].r[2])/si.dl[2]-sx.i0[2];

            /* __ index for the part "position" __ */
            i = (int)xw;
            j = (int)yw;
            k = (int)zw;

            /* __ index for the 1D array __ */
            ijk0 = IDX(i, j, k, sx.n[0], sx.n[1], sx.n[2]);

            /* __ beware to have enough space in s0[ijk0].ipart[s] pointer __ */
            if (s0[ijk0].npart[s] >= oldNumPerCell[s][ijk0]) {
                s0[ijk0].ipart[s] = realloc(s0[ijk0].ipart[s], (s0[ijk0].npart[s]+1)*sizeof(int));

                /* __ be sure the realloc is ok __ */
                if (s0[ijk0].ipart[s] == NULL) {
                    printf("the memory reallocation of s0.ipart went wrong !\n");
                }
            }

            /* __ update the list of particles in the cell __ */
            s0[ijk0].ipart[s][s0[ijk0].npart[s]] = m;
            s0[ijk0].npart[s]++;
        }
    }

    /* __ calculate the ion pressure tensor __ */
    ions(&si, &sx, s2, sp, ghosts, hbc, 0);

    /* __ loop over the cells __ */
    for (i = 0; i < sx.n[0]; i++) {
        for (j = 0; j < sx.n[1]; j++) {
            for (k = 0; k < sx.n[2]; k++) {

            ijk  = IDX(i  , j  , k  , sx.n[0], sx.n[1], sx.n[2]);
            ijk1 = IDX(i+1, j+1, k+1, sx.n[0]+2, sx.n[1]+2, sx.n[2]+2);

            /* __ calculate the ion temperatures __ */
            for (s = 1; s < NS+1; s++) {
                ionTemp[s] = (s2[ijk1].ps[s][0]
                             +s2[ijk1].ps[s][3]
                             +s2[ijk1].ps[s][5])
                            /(3*s2[ijk1].ns[s]);
            }

            /* __ compute ion-ion collision in the "ijk" cell __ */
            ionIonCellCollisions(sc, s0[ijk], s2[ijk1], si, sx, sd, sp, ionTemp, ipc);
            }
        }
    }

    /* __ then free the ipart pointer __ */
    for (ijk0 = 0; ijk0 < n0; ijk0++) {
        for (s = 1; s < NS+1; s++) {
            free(s0[ijk0].ipart[s]);
        }
    }

}

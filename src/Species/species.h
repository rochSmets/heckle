#ifndef SPECIES_H
#define SPECIES_H


#include <particlecomm.h>
#include <particle.h>



typedef struct s_species
{
    int id;
    int kind;
    double m;
    double q;

    Particle *parr;
    int parrsize;
    int npartloc;
    int npartglo;
    double ws;

} Species;



#endif // SPECIES_H


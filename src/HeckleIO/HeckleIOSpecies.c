
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <hdf5.h>
#include <stdint.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "iamdead.h"



// ___________________________________________________________________________
//
//                             PRIVATE FUNCTIONS
// ___________________________________________________________________________



typedef struct s_HeckleIOSpecies {

    hid_t file;            // identifier for the species file
    hid_t writePlist;      // property list for parallel writing
    hid_t memspace[NS+1];  // dataspace for subdomain
    hid_t filespace[NS+1]; // dataspace for domain (to write in the file)
} HeckleIOSpecies;



// ___________________________________________________________________________
//
// createSpeciesDataSpaces()
//
// AIM : creates the memory and file dataspaces for a given specie
// should be called at each write (# of part can change) & then closed after.
// ___________________________________________________________________________
//
void createSpeciesDataSpaces(HeckleIOSpecies *ioSpecie,
                             const STI* const si,
                             STX* const sx)
{

    hsize_t offset[NS+1];
    hsize_t nLoc[NS+1];
    hsize_t nGlob[NS+1];
    int s;
    int64_t *np[NS+1];
    int64_t ts[NS+1];
    int rank;


    // _____________________ MEMORY DATASPACE ________________________________

    // memspace size is always the # of part in subdomain
    for (s = 1; s < NS+1; s++) {
        nLoc[s] = (hsize_t) sx->ns[s];
        if (sx->ns[s] == 0) {
            printf("hoops : the memory dataspace for specie is going to be empty !\n");
            fflush(stdout);
        }
    }

    // create the dataspace in memory (memspace)
    for (s = 1; s < NS+1; s++) {
        ioSpecie->memspace[s] = H5Screate_simple(1, nLoc+s, NULL);
    }


    // _____________________ FILE DATASPACE __________________________________

    // filespace has the size of the total # of part from a given specie
    // the offset of node r is given by the total # of part in nodes < r
    for (s = 1; s < NS+1; s++) {
        nGlob[s]  = (hsize_t) si->ns[s];

        // this buffer records the # of parts of specie s in all nodes
        np[s] = malloc(sx->s * sizeof(*np[s]));

        // and collect these values from all nodes (for each nodes)
        MPI_Allgather(&sx->ns[s], 1, MPI_LONG_LONG_INT, np[s], 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

        // and this one is the total number of part in all nodes < r
        for (ts[s] = 0, rank = 0; rank < sx->r; rank++) {
            ts[s] += np[s][rank];
        }

        free(np[s]);

        // this value being the offset, so the part can be piled-up
        offset[s] = (hsize_t) ts[s];
    }


    // dataspace for the file
    for (s = 1; s < NS+1; s++) {
        ioSpecie->filespace[s] = H5Screate_simple(1, nGlob+s, NULL);

        H5Sselect_hyperslab(ioSpecie->filespace[s], H5S_SELECT_SET, offset+s, NULL, nLoc+s, NULL);
    }

}



// ___________________________________________________________________________
//
// closeSpeciesDataSpaces()
//
// AIM : close the memory and file dataspaces at the end of the run.
// ___________________________________________________________________________
//
void closeSpeciesDataSpaces(HeckleIOSpecies *ioSpecie) {

    int s;


    // close both dataspaces
    for (s = 1; s < NS+1; s++) {
        H5Sclose(ioSpecie->memspace[s]);
        H5Sclose(ioSpecie->filespace[s]);
    }

}



// ___________________________________________________________________________
//
// writeSpecieComponent()
//
// AIM : writes a single component of a given specie coordinate
// in the group groupe_id of the file species.h5
// ___________________________________________________________________________
//
void writeSpecieComponent(HeckleIOSpecies *ioSpecie, int s,
                          const float* const field,
                          const char* const fieldName,
                          hid_t group)
{
    hid_t dset;

    // create the dataset
    dset  = H5Dcreate(group, fieldName, H5T_NATIVE_FLOAT, ioSpecie->filespace[s],
                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // now the dataset is ready for us to write in it
    H5Dwrite(dset, H5T_NATIVE_FLOAT, ioSpecie->memspace[s],
             ioSpecie->filespace[s], ioSpecie->writePlist, field);

    // can close the dataset now
    H5Dclose(dset);
}



// ___________________________________________________________________________
//
// HeckleIOOpenSpeciesFile()
//
// AIM : Open the species file
// ___________________________________________________________________________
//
hid_t HeckleIOOpenSpeciesFile(HeckleIOSpecies *self)
{
    hid_t access;
    MPI_Info info = MPI_INFO_NULL;


    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    // then open the species.h5 file
    self->file = H5Fopen("species.h5", H5F_ACC_RDWR, access);

    // now close the plist since not needed
    H5Pclose(access);

    return self->file;

}






// ___________________________________________________________________________
//
//                             PUBLIC FUNCTIONS
// ___________________________________________________________________________



// ___________________________________________________________________________
//
// HeckleIOInitSpecies()
//
// AIM : Initialize the HeckleIOSpecies module
// ___________________________________________________________________________
//
HeckleIOSpecies *HeckleIOInitSpecies(const STI* const si,
                                     const STX* const sx)
{
    MPI_Info info = MPI_INFO_NULL;
    hid_t space;
    hid_t attrib;
    hid_t access;
    hsize_t nDim;
    HeckleIOSpecies *ioSpecie;


    // create the HeckleIOSpecies structure
    ioSpecie = malloc(sizeof *ioSpecie);
    if (ioSpecie == NULL) {
        printf("IOError - process %d could not allocate memory for Species\n", sx->r);
        return NULL;
    }

    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    // now create the HDF5 file
    ioSpecie->file = H5Fcreate("species.h5", H5F_ACC_TRUNC, H5P_DEFAULT, access);

    // now close the plist since not needed
    H5Pclose(access);

    // create the writing property list, to write in parallel
    ioSpecie->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(ioSpecie->writePlist, H5FD_MPIO_COLLECTIVE);

    // then we add some attribute : 'mass'
    nDim = NS+1;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioSpecie->file, "mass", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, si->ms);
    H5Aclose(attrib);
    H5Sclose(space);

    // attribute : 'charge'
    nDim = NS+1;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioSpecie->file, "charge", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, si->qs);
    H5Aclose(attrib);
    H5Sclose(space);

    // Attribute : 'weight'
    nDim = NS+1;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioSpecie->file, "weight", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, si->ws);
    H5Aclose(attrib);
    H5Sclose(space);

    // Attribute : 'numofpart'
    nDim = NS+1;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioSpecie->file, "nbrOfParticles", H5T_NATIVE_INT64, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_INT64, si->ns);
    H5Aclose(attrib);
    H5Sclose(space);

    // Attribute : dumpspecies
    space = H5Screate(H5S_SCALAR);
    attrib = H5Acreate(ioSpecie->file, "dumpSpecies", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT),
    H5Awrite(attrib, H5T_NATIVE_INT, &(si->tp));
    H5Aclose(attrib);
    H5Sclose(space);

    H5Fclose(ioSpecie->file);

    return ioSpecie;
}



// ___________________________________________________________________________
//
// HeckleIODeleteSpecies()
//
// AIM : close the write property list
// ___________________________________________________________________________
//
void HeckleIODeleteSpecies(HeckleIOSpecies *ioSpecie) {


    H5Pclose(ioSpecie->writePlist);       // close the write property list

}



// ___________________________________________________________________________
//
// writeSpecies()
//
// AIM : this routine writes the species (position, velocity & id)
// in a HDF5 file.
// ___________________________________________________________________________
//
void writeSpecies(HeckleIOSpecies *ioSpecie,
                  const STI* const si,
                  STX* const sx,
                  struct stp *sp[NS+1], double time) {

    int c, m, s;
    char timeName[800];
    char fieldName[800];
    char specieName[800];
    float *particlePosition[NS+1][3],
          *particleVelocity[NS+1][3],
          *particleGlobalIndex[NS+1];


    createSpeciesDataSpaces(ioSpecie, si, sx);

    // open the species.h5 file
    if (HeckleIOOpenSpeciesFile(ioSpecie) > 0) {
    }
    else {
        printf("your file species.h5 has a problem my man...\n");
    }

    /* __ print informations __ */
    if (sx->r == 0) {
        printf("________________ write species @ t = %11.5f _____\n\n", time);
    }

    for (s = 1; s < NS+1; s++)
    {
        for (c = 0; c < 3; c++)
        {
            particlePosition[s][c]  = malloc(sx->ns[s] * sizeof *particlePosition[s][c]);
            particleVelocity[s][c]  = malloc(sx->ns[s] * sizeof *particleVelocity[s][c]);
        }

        particleGlobalIndex[s] = malloc(sx->ns[s] * sizeof *particleGlobalIndex[s]);

        /* __ loop on the part. of specie "s" __ */
        for (m = 0; m < sx->ns[s]; m++)
        {
            particleGlobalIndex[s][m] = (float)sp[s][m].i;

            for (c = 0; c < 3; c++)
            {
                particlePosition[s][c][m] = (float)sp[s][m].r[c];
                particleVelocity[s][c][m] = (float)sp[s][m].v[c];
            }
        }

    }

    // now we have all the data, let's create the group
    sprintf(timeName, "time : %f", time);

    hid_t timegroup_id = H5Gcreate(ioSpecie->file, timeName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (s = 1; s < NS+1; s++) {
        // & another group under time group... one for each specie
        sprintf(specieName, "specie %d", s);
        hid_t group_id = H5Gcreate(timegroup_id, specieName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (c = 0; c < 3; c++) {
            sprintf(fieldName, "r[%d]", c);
            if (si->ns[s] != 0) writeSpecieComponent(ioSpecie, s, particlePosition[s][c], fieldName, group_id);
        }

        for (c = 0; c < 3; c++) {
            sprintf(fieldName, "v[%d]", c);
            if (si->ns[s] != 0) writeSpecieComponent(ioSpecie, s, particleVelocity[s][c], fieldName, group_id);
        }

        strcpy(fieldName, "index");
        if (si->ns[s] != 0) writeSpecieComponent(ioSpecie, s, particleGlobalIndex[s], fieldName, group_id);

        H5Gclose(group_id);
    }


    // now close the group
    H5Gclose(timegroup_id);
    H5Fflush(ioSpecie->file, H5F_SCOPE_GLOBAL);

    // now free the memory
    for (s = 1; s < NS+1; s++) {
        for (c = 0; c < 3; c++) {
            free(particlePosition[s][c]);
            free(particleVelocity[s][c]);
        }

        free(particleGlobalIndex[s]);
    }


    closeSpeciesDataSpaces(ioSpecie);
    H5Fclose(ioSpecie->file);
}



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <hdf5.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "iamdead.h"



// ___________________________________________________________________________
//
//                             PRIVATE FUNCTIONS
// ___________________________________________________________________________



typedef struct s_HeckleIOFields {
    hid_t file;       // identifier for the fields file
    hid_t writePlist; // property list for parallel writing
    hid_t memspace;   // dataspace for subdomain
    hid_t filespace;  // dataspace for domain (to write in the file)
    int n[3];         // # of points to write in each direction
} HeckleIOFields;



// ___________________________________________________________________________
//
// createFieldsDataSpaces()
//
// AIM : creates the memory and file dataspaces once and for all, based on
// the local and global size of the g1 grid and the offset of the local grid
// ___________________________________________________________________________
//
void createFieldsDataSpaces(HeckleIOFields *ioField,
                            const STI* const si,
                            const STX* const sx)
{

    hsize_t offset[3];
    hsize_t nLoc[3];
    hsize_t nGlob[3];
    int i;

    // _____________________ MEMORY DATASPACE ________________________________

    // for 1 process, memspace size == filespace size == whole grid
    if (sx->s == 1) {
        for (i = 0; i < 3; i++) {
            nLoc[i] = sx->n[i]+1;
        }
    }

    // now parallel case
    // the default behavior is to write only the first n points in each direction
    // i.e. only the points from 0 to sx->n[direction]-1 (included)
    // only the MPI domains reaching the border at xmax,ymax,zmax will write
    // n+1 points, i.e. from 0 to n[direction] (included) in the direction in
    // which they reach the border of the box.
    else {
        // default behavior : only writ the n-1 first points
        for (i = 0; i < 3; i++) {
            nLoc[i] = sx->n[i];
        }

        // for each direction look whether the current MPI domain
        // reaches the border
        for (i = 0; i < 3; i++) {
            // if we reach the border
            if (sx->i1[i] == si->n[i]) {
                nLoc[i]   = sx->n[i]+1;
            }
        }
    }


    // keep memory of # of points in subdomain for writeFields
    for (i = 0;  i < 3; i++) {
        ioField->n[i] = nLoc[i];
    }

    // create the dataspace in memory (memspace)
    ioField->memspace = H5Screate_simple(3, nLoc, NULL);



    // _____________________ FILE DATASPACE __________________________________

    // filespace has the size of the global grid g1
    // the offset of process P is given by the global index sx->i0[0..2]
    // for each direction
    for (i =0; i < 3; i++) {
        nGlob[i] = si->n[i]+1;
        offset[i] = sx->i0[i];
    }

    // dataspace for the file
    ioField->filespace = H5Screate_simple(3, nGlob, NULL);

    // we need to know where the local grid will be
    // in the file (offset) and how big it is (nLoc)
    H5Sselect_hyperslab(ioField->filespace, H5S_SELECT_SET, offset, NULL, nLoc, NULL);

}


// ___________________________________________________________________________
//
// writeSingleField()
//
// AIM : writes a single field array into the group 'group' of the file
// fields.h5
// ___________________________________________________________________________
//
void writeSingleField(HeckleIOFields *ioField,
                      const float* const field,
                      const char* const fieldName,
                      hid_t group)
{

    hid_t dset;


    // create the dataset
    dset  = H5Dcreate(group, fieldName, H5T_NATIVE_FLOAT, ioField->filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // now the dataset is ready for us to write in it
    H5Dwrite(dset, H5T_NATIVE_FLOAT, ioField->memspace, ioField->filespace, ioField->writePlist, field);

    // can close the dataset now
    H5Dclose(dset);

}



// ___________________________________________________________________________
//
// HeckleIOOpenFieldsFile()
//
// AIM : Open the field file
// ___________________________________________________________________________
//
hid_t HeckleIOOpenFieldsFile(HeckleIOFields *self)
{
    hid_t access;
    MPI_Info info = MPI_INFO_NULL;


    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    // then open the fields.h5 file
    self->file = H5Fopen("fields.h5", H5F_ACC_RDWR, access);

    // now close the property list since not needed
    H5Pclose(access);

    return self->file;
}



// ___________________________________________________________________________
//
//                             PUBLIC FUNCTIONS
// ___________________________________________________________________________



// ___________________________________________________________________________
//
// HeckleIOInitFields()
//
// AIM : Initialize the HeckleIOFields module
// ___________________________________________________________________________
//
HeckleIOFields *HeckleIOInitFields(const STI* const si,
                                   const STX* const sx)
{

    MPI_Info info = MPI_INFO_NULL;
    hid_t space;
    hid_t attrib;
    hid_t access;
    hsize_t nDim;
    HeckleIOFields *ioField;


    // create the HeckleIOFields structure
    ioField = malloc(sizeof(*ioField));
    if (ioField == NULL) {
        printf("IOError - process %d could not allocate memory for Fields\n", sx->r);
        return NULL;
    }

    // make the file access property list
    // and give it MPI/communicator info
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    // now create the HDF5 file
    ioField->file = H5Fcreate("fields.h5", H5F_ACC_TRUNC, H5P_DEFAULT, access);

    // now close the plist since not needed
    H5Pclose(access);

    // create the writing property list, to write in parallel
    ioField->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(ioField->writePlist, H5FD_MPIO_COLLECTIVE);


    // the file is now created...
    // we now create attributes to store useful data in the begining of the file

    // now some attribute : numofcells
    nDim = 3;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioField->file, "nbrOfCells", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_INT, &(si->n));
    H5Aclose(attrib);
    H5Sclose(space);

    // attribute : domainsize
    nDim = 3;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioField->file, "domainSize", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, &(si->l));
    H5Aclose(attrib);
    H5Sclose(space);

    // attribute : meshsize
    nDim = 3;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(ioField->file, "meshSize", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, &(si->dl));
    H5Aclose(attrib);
    H5Sclose(space);

    // attribute : dumpfields
    space = H5Screate(H5S_SCALAR);
    attrib = H5Acreate(ioField->file, "dumpFields", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT),
    H5Awrite(attrib, H5T_NATIVE_INT, &(si->tf));
    H5Aclose(attrib);
    H5Sclose(space);

    // attribute : hyperresistivity
    space = H5Screate(H5S_SCALAR);
    attrib = H5Acreate(ioField->file, "hyperResistivity", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT),
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, &(si->hyvi));
    H5Aclose(attrib);
    H5Sclose(space);

    // attribute : resistivity
    space = H5Screate(H5S_SCALAR);
    attrib = H5Acreate(ioField->file, "resistivity", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT),
    H5Awrite(attrib, H5T_NATIVE_DOUBLE, &(si->rsty));
    H5Aclose(attrib);
    H5Sclose(space);

    createFieldsDataSpaces(ioField, si, sx);

    H5Fclose(ioField->file);

    return ioField;

}


// ___________________________________________________________________________
//
// HeckleIODeleteFields()
//
// AIM : close the write property list & the dataspace
// ___________________________________________________________________________
//
void HeckleIODeleteFields(HeckleIOFields *ioField)
{


    // close write property list & dataspaces
    H5Pclose(ioField->writePlist);
    H5Sclose(ioField->memspace);
    H5Sclose(ioField->filespace);

}



// ___________________________________________________________________________
//
// writeFields()
//
// AIM : this routine writes the fields (electromagnetic and fluid moments)
// in a HDF5 file.
// ___________________________________________________________________________
//
void writeFields(HeckleIOFields *ioField,
                 const STX* const sx,
                 const ST1* const s1,
                 const ST2* const s2,
                 double time)
{

    int ijk, i, j, k;
    int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int ijkn;
    int numCellG1[3],numCellG2[3];
    int c, s, l;
    char timeName[800];
    char fieldName[800];
    int nbrOfNodes;
    hid_t group;
    float *B[3], *E[3], *J[3], *J_smoothed[3];
    float *Vi[3], *ns[NS+1],*Vs[NS+1][3], *Ps[NS+1][6];


    // open the fields.h5 file
    if (HeckleIOOpenFieldsFile(ioField) > 0) {
    }

    else {
        printf("your file fields.h5 has a problem my man...\n");
    }

    // print informations
    if (sx->r == 0) {
        printf("________________ write fields @ t = %11.5f ______\n\n", time);
    }

    // shortcut for # of grid points on g1 and g2 in each direction
    for (c = 0; c < 3; c++) {
        numCellG1[c] = sx->n[c]+1;
        numCellG2[c] = sx->n[c]+2;
    }

    // shortcuts for total # of points on g1, g2 & # of points to write (hio->n)
    nbrOfNodes   = ioField->n[0]*ioField->n[1]*ioField->n[2];

    for (c = 0; c < 3; c++) {
        B[c]  = malloc(nbrOfNodes * sizeof *B[c]);
        E[c]  = malloc(nbrOfNodes * sizeof *E[c]);
        J[c]  = malloc(nbrOfNodes * sizeof *J[c]);
        J_smoothed[c]  = malloc(nbrOfNodes * sizeof *J_smoothed[c]);
        Vi[c] = malloc(nbrOfNodes * sizeof *Vi[c]);
    }

    for (s = 0; s < NS+1; s++) {
        ns[s] = malloc(nbrOfNodes * sizeof *ns[s]);

        for (l=0; l < 3; l++) {
            Vs[s][l] = malloc(nbrOfNodes * sizeof *Vs[s][l]);
        }

        for (l=0; l < 6; l++) {
            Ps[s][l] = malloc(nbrOfNodes * sizeof *Ps[s][l]);
        }
    }


    // fill the arrays: nested loops on the subdomain
    // only loop over the points we actually want to write
    // hio->n excludes redundancy for interior domain boundary points
    for (i = 0; i < ioField->n[0]; i++) {
        for (j = 0; j < ioField->n[1]; j++) {
            for (k = 0; k < ioField->n[2]; k++) {
                /* __ set index on g1 __ */
                ijk = IDX(i, j, k, numCellG1[0], numCellG1[1], numCellG1[2]);

                // index for g1 sub-grid
                ijkn = IDX(i, j, k, ioField->n[0], ioField->n[1], ioField->n[2]);

                /* __ set indexes on g2 __ */
                ijk1 = IDX(i  , j  , k  , numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk2 = IDX(i+1, j  , k  , numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk3 = IDX(i  , j+1, k  , numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk4 = IDX(i+1, j+1, k  , numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk5 = IDX(i  , j  , k+1, numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk6 = IDX(i+1, j  , k+1, numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk7 = IDX(i  , j+1, k+1, numCellG2[0], numCellG2[1], numCellG2[2]);
                ijk8 = IDX(i+1, j+1, k+1, numCellG2[0], numCellG2[1], numCellG2[2]);

                /* __ for g1 grid __ */
                for (c = 0; c < 3; c++) {
                    B[c][ijkn] = (float)s1[ijk].b[c];
                }

                /* __ for g2 grid species independant quantities __ */
                for (c = 0; c < 3; c++) {
                    E[c][ijkn] = (float)0.125*(s2[ijk1].e[c]
                                              +s2[ijk2].e[c]
                                              +s2[ijk3].e[c]
                                              +s2[ijk4].e[c]
                                              +s2[ijk5].e[c]
                                              +s2[ijk6].e[c]
                                              +s2[ijk7].e[c]
                                              +s2[ijk8].e[c]);

                    J[c][ijkn] = (float)0.125*(s2[ijk1].j[c]
                                              +s2[ijk2].j[c]
                                              +s2[ijk3].j[c]
                                              +s2[ijk4].j[c]
                                              +s2[ijk5].j[c]
                                              +s2[ijk6].j[c]
                                              +s2[ijk7].j[c]
                                              +s2[ijk8].j[c]);

                    J_smoothed[c][ijkn] = (float)0.125*(s2[ijk1].j_smoothed[c]
                                              +s2[ijk2].j_smoothed[c]
                                              +s2[ijk3].j_smoothed[c]
                                              +s2[ijk4].j_smoothed[c]
                                              +s2[ijk5].j_smoothed[c]
                                              +s2[ijk6].j_smoothed[c]
                                              +s2[ijk7].j_smoothed[c]
                                              +s2[ijk8].j_smoothed[c]);

                    Vi[c][ijkn] = (float)0.125*(s2[ijk1].vi[c]
                                               +s2[ijk2].vi[c]
                                               +s2[ijk3].vi[c]
                                               +s2[ijk4].vi[c]
                                               +s2[ijk5].vi[c]
                                               +s2[ijk6].vi[c]
                                               +s2[ijk7].vi[c]
                                               +s2[ijk8].vi[c]);
                }


                /* __ for g2 grid species dependant quantities __ */
                for (s = 0; s < NS+1; s++) {
                    ns[s][ijkn] = (float)0.125*(s2[ijk1].ns[s]
                                               +s2[ijk2].ns[s]
                                               +s2[ijk3].ns[s]
                                               +s2[ijk4].ns[s]
                                               +s2[ijk5].ns[s]
                                               +s2[ijk6].ns[s]
                                               +s2[ijk7].ns[s]
                                               +s2[ijk8].ns[s]);
                    //if (s == 2) printf("%f\n", ns[s][ijkn]);

                    for (c = 0; c < 3; c++) {
                        Vs[s][c][ijkn] = (float)0.125*(s2[ijk1].vs[s][c]
                                                      +s2[ijk2].vs[s][c]
                                                      +s2[ijk3].vs[s][c]
                                                      +s2[ijk4].vs[s][c]
                                                      +s2[ijk5].vs[s][c]
                                                      +s2[ijk6].vs[s][c]
                                                      +s2[ijk7].vs[s][c]
                                                      +s2[ijk8].vs[s][c]);
                    }

                    for (c = 0; c < 6; c++) {
                        Ps[s][c][ijkn] = (float)0.125*(s2[ijk1].ps[s][c]
                                                      +s2[ijk2].ps[s][c]
                                                      +s2[ijk3].ps[s][c]
                                                      +s2[ijk4].ps[s][c]
                                                      +s2[ijk5].ps[s][c]
                                                      +s2[ijk6].ps[s][c]
                                                      +s2[ijk7].ps[s][c]
                                                      +s2[ijk8].ps[s][c]);
                    }
                }
            }
        }
    }


    // now we have all the data, let's create the group
    sprintf(timeName, "time : %f", time);
    group = H5Gcreate(ioField->file, timeName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // write the datasets
    writeSingleField(ioField, B[0], "Bx", group);
    writeSingleField(ioField, B[1], "By", group);
    writeSingleField(ioField, B[2], "Bz", group);

    writeSingleField(ioField, E[0], "Ex", group);
    writeSingleField(ioField, E[1], "Ey", group);
    writeSingleField(ioField, E[2], "Ez", group);

    writeSingleField(ioField, J[0], "Jx", group);
    writeSingleField(ioField, J[1], "Jy", group);
    writeSingleField(ioField, J[2], "Jz", group);

    writeSingleField(ioField, J_smoothed[0], "J_sm_x", group);
    writeSingleField(ioField, J_smoothed[1], "J_sm_y", group);
    writeSingleField(ioField, J_smoothed[2], "J_sm_z", group);

    writeSingleField(ioField, Vi[0], "Vix", group);
    writeSingleField(ioField, Vi[1], "Viy", group);
    writeSingleField(ioField, Vi[2], "Viz", group);


    // density of each species
    for (s = 0; s < NS+1; s++) {
        sprintf(fieldName, "n[%d]", s);
        writeSingleField(ioField, ns[s], fieldName, group);
    }

    // bulk velocity of each species
    for (s = 0; s < NS+1; s++) {
        sprintf(fieldName, "Vx[%d]", s);
        writeSingleField(ioField, Vs[s][0], fieldName, group);
        sprintf(fieldName, "Vy[%d]", s);
        writeSingleField(ioField, Vs[s][1], fieldName, group);
        sprintf(fieldName, "Vz[%d]", s);
        writeSingleField(ioField, Vs[s][2], fieldName, group);
    }


    // pressure tensor of each species
    for (s = 0; s < NS+1; s++) {
        sprintf(fieldName, "Pxx[%d]", s);
        writeSingleField(ioField, Ps[s][0], fieldName, group);

        sprintf(fieldName, "Pxy[%d]", s);
        writeSingleField(ioField, Ps[s][1], fieldName, group);

        sprintf(fieldName, "Pxz[%d]", s);
        writeSingleField(ioField, Ps[s][2], fieldName, group);

        sprintf(fieldName, "Pyy[%d]", s);
        writeSingleField(ioField, Ps[s][3], fieldName, group);

        sprintf(fieldName, "Pyz[%d]", s);
        writeSingleField(ioField, Ps[s][4], fieldName, group);

        sprintf(fieldName, "Pzz[%d]", s);
        writeSingleField(ioField, Ps[s][5], fieldName, group);
    }

    // close the group & flush the file
    H5Gclose(group);
    H5Fflush(ioField->file, H5F_SCOPE_GLOBAL);

   // free the memory
    for (c = 0; c < 3; c++) {
       free(B[c]);
       free(E[c]);
       free(J[c]);
       free(J_smoothed[c]);
       free(Vi[c]);
    }

    for (s = 0; s < NS+1; s++) {
        free(ns[s]);

        for(c = 0; c < 3; c++) {
            free(Vs[s][c]);
        }

        for(l = 0; l < 6; l++) {
            free(Ps[s][l]);
        }
    }

    // close the file
    H5Fclose(ioField->file);

}


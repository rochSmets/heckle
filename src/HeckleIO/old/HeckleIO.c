
#ifndef PYP
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <hdf5.h>

#include "structures.h"
#include "defines.h"
#include "misc.h"
#include "iamdead.h"

// ___________________________________________________________________________
//
//                             PRIVATE FUNCTIONS
// ___________________________________________________________________________



// structure for the associated h5 file
typedef struct s_HeckleIOFields
{
    hid_t FieldsFileID;     // identifier for the fields file
    hid_t writePlist;       // property list for parallel writing
    hid_t memspace;         // dataspace for subdomain
    hid_t filespace;        // dataspace for domain (to write in the file)
    int n[3];               // # of points to write in each direction
} HeckleIOFields;


typedef struct s_HeckleIOSpecies
{
    hid_t SpeciesFileID;    // identifier for the species file
    hid_t writePlist;       // property list for parallel writing
    hid_t memspace[NS+1];   // dataspace for subdomain
    hid_t filespace[NS+1];  // dataspace for domain (to write in the file)
} HeckleIOSpecies;


typedef struct s_HeckleIORestart
{
    hid_t RestartsFileID;   // identifier for the restarts file
    hid_t writePlist;       // property list for parallel writing
    hid_t *s1space;         // dataspace for st1 in subdomain
    hid_t *s2space;         // dataspace for st2 in for subdomain
    hid_t *spspace[NS+1];   // dataspace for istp[NS+1] in subdomain
    hid_t s1typeId;         // datatype for s1
    hid_t s2typeId;         // datatype for s2
    hid_t sptypeId;         // datatype for sp
    int numofcpu;
} HeckleIORestart;



// ___________________________________________________________________________
//
// createFieldsDataSpaces()
//
// AIM : creates the memory and file dataspaces once and for all, based on
// the local and global size of the g1 grid and the offset of the local grid.
// ___________________________________________________________________________
//
void createFieldsDataSpaces(HeckleIOFields *hiof, struct sti *si, struct stx *sx)
{
    hsize_t offset[3];
    hsize_t nloc[3];
    hsize_t nglob[3];
    int idim;

    // _____________________ MEMORY DATASPACE ________________________________

    // for 1 process, memspace size == filespace size == whole grid
    if (sx->s == 1) {
        for (idim = 0; idim < 3; idim++) {
            nloc[idim] = sx->n[idim]+1;
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
        for (idim = 0; idim < 3; idim++) {
            nloc[idim] = sx->n[idim];
        }

        // for each direction look whether the current MPI domain
        // reaches the border
        for (idim = 0; idim < 3; idim++) {
            // if we reach the border
            if (sx->i1[idim] == si->n[idim]) {
                nloc[idim]   = sx->n[idim]+1;
            }
        }
    }


    // keep memory of # of points in subdomain for writeFields
    for (idim = 0;  idim < 3; idim++) {
        hiof->n[idim] = nloc[idim];
    }

    // create the dataspace in memory (memspace)
    hiof->memspace = H5Screate_simple(3, nloc, NULL);



    // _____________________ FILE DATASPACE __________________________________

    // filespace has the size of the global grid g1
    // the offset of process P is given by the global index sx->i0[0..2]
    // for each direction
    for (idim =0; idim < 3; idim++) {
        nglob[idim] = si->n[idim]+1;
        offset[idim] = sx->i0[idim];
    }

    // dataspace for the file
    hiof->filespace = H5Screate_simple(3, nglob, NULL);

    // we need to know where the local grid will be
    // in the file (offset) and how big it is (nloc)
    H5Sselect_hyperslab(hiof->filespace,    // file dataspace
                        H5S_SELECT_SET,     // select
                        offset,             // offset
                        NULL,               // stride = 1
                        nloc,               // number of g1 points
                        NULL);              // block = 1 point

}


// ___________________________________________________________________________
//
// createSpeciesDataSpaces()
//
// AIM : creates the memory and file dataspaces for a given specie
// should be called at each write (# of part can change) & then closed after.
// ___________________________________________________________________________
//
void createSpeciesDataSpaces(HeckleIOSpecies *hios, struct sti *si, struct stx *sx)
{
    hsize_t offset[1];
    hsize_t nloc[NS+1];
    hsize_t nglob[NS+1];
    int ispe;
    int *np[NS+1];
    int ts[NS+1];
    int rank;


    // _____________________ MEMORY DATASPACE ________________________________

    // memspace size is always the # of part in subdomain
    for (ispe = 1; ispe <= NS; ispe++) {
        nloc[ispe]   = sx->ns[ispe];
    }

    // create the dataspace in memory (memspace)
    for (ispe = 1; ispe < NS+1; ispe++) {
        hios->memspace[ispe] = H5Screate_simple(1, nloc+ispe, NULL);
    }


    // _____________________ FILE DATASPACE __________________________________

    // filespace has the size of the total # of part from a given specie
    // the offset of node r is given by the total # of part in nodes < r
    for (ispe = 1; ispe < NS+1; ispe++) {
        nglob[ispe]  = si->ns[ispe];


        // this buffer records the # of parts of specie ispe in all nodes
        np[ispe] = malloc(sx->s * sizeof(*np[ispe]));

        // and collect these values from all nodes (for each nodes)
        MPI_Allgather(&sx->ns[ispe], 1, MPI_INT, np[ispe], 1, MPI_INT, MPI_COMM_WORLD);

        // and this one is the total number of part in all nodes < r
        for (ts[ispe] = 0, rank = 0; rank < sx->r; rank++) {
            ts[ispe] += np[ispe][rank];
        }

        free(np[ispe]);

        // this value being the offset, so the part can be piled-up
        offset[ispe] = ts[ispe];
    }


    // dataspace for the file
    for (ispe = 1; ispe < NS+1; ispe++) {
        hios->filespace[ispe] = H5Screate_simple(1, nglob+ispe, NULL);

        H5Sselect_hyperslab(hios->filespace[ispe], // file dataspace
                            H5S_SELECT_SET,        // select
                            offset+ispe,           // offset
                            NULL,                  // stride = 1
                            nloc+ispe,             // local # of particles
                            NULL);                 // block = 1 point
    }

}




// ___________________________________________________________________________
//
// createRestartDataSpaces()
//
// AIM : creates the memories and files dataspaces for s1, s2 & sp struct
// should be called at each write (# of part can change) & then closed after.
// ___________________________________________________________________________
//
void createRestartDataSpaces(HeckleIORestart *hior, struct sti *si, struct stx *sx)
{
    hsize_t nofcell[1];
    int i;
//  int ispe;
    int nx, ny, nz;
    int* sxn;


    // we need as many dataspace as nodes
    hior->s1space = malloc(hior->numofcpu * sizeof(hid_t));

    // and a buffer to collect the grid size & # of part for each nodes
    sxn = malloc(sx->s * 3 * sizeof(int));

    // gather all values of sx.n[3] for each node... on each node
    MPI_Allgather(sx->n, 3, MPI_INT, sxn, 3, MPI_INT, MPI_COMM_WORLD);

    // rebuild the size of the dataspace for each nodes
    for (i = 0; i < hior->numofcpu; i++) {
        nx = sxn[3*i];
        ny = sxn[3*i+1];
        nz = sxn[3*i+2];

        // set s1 dataspace for each node
        nofcell[0] = (nx+1)*(ny+1)*(nz+1);
        hior->s1space[i] = H5Screate_simple(1, nofcell, NULL);
    }


// all above has to be set for st2 & stp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  nofcell[0] = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);
//  hior->s2space = H5Screate_simple(1, nofcell, NULL);

    // memspace size is always the # of part in subdomain
//  for (ispe = 1; ispe < NS+1; ispe++) {
//      nofcell[0] = sx->ns[ispe];
//      hior->spspace[ispe] = H5Screate_simple(1, nofcell, NULL);
//  }

    free(sxn);

}



// ___________________________________________________________________________
//
// closeFieldDataSpaces()
//
// AIM : close the memory and file dataspaces at the end of the run.
// ___________________________________________________________________________
//
void closeFieldDataSpaces(HeckleIOFields *hiof)
{
    H5Sclose(hiof->memspace);
    H5Sclose(hiof->filespace);
}



// ___________________________________________________________________________
//
// closeSpeciesDataSpaces()
//
// AIM : close the memory and file dataspaces at the end of the run.
// ___________________________________________________________________________
//
void closeSpeciesDataSpaces(HeckleIOSpecies *hios)
{
    int ispe;


    for (ispe = 1; ispe < NS+1; ispe++) {
        H5Sclose(hios->memspace[ispe]);
        H5Sclose(hios->filespace[ispe]);
    }
}



// ___________________________________________________________________________
//
// closeRestartDataSpaces()
//
// AIM : close the dataspaces at the end of the run.
// ___________________________________________________________________________
//
void closeRestartDataSpaces(HeckleIORestart *hior)
{
    int i;
//  int ispe;


    for (i = 0; i < hior->numofcpu; i++) {
        H5Sclose(hior->s1space[i]);
//      H5Sclose(hior->s2space[i]);

//      for (ispe = 1; ispe < NS+1; ispe++)
//          H5Sclose(hior->spspace[ispe][i]);
    }

}



// ___________________________________________________________________________
//
//createRestartDataType()
//
// AIM : create the compound datatype to be used for restart (read & write)
// in restarts.h5
// ___________________________________________________________________________
//
void createRestartDataType(HeckleIORestart *hior, struct sti *si)
{
hid_t s1typeId;
//hid_t s2typeId;
//hid_t sptypeId;
hid_t typeFor1Double;
hid_t typeFor3Double;
hid_t typeFor6Double;
hid_t typeForNS1Double;
hid_t typeForNS1by3Double;
hid_t typeForNS1by6Double;
hid_t typeFor1Int;
hsize_t dim1[1] = {1};
hsize_t dim3[1] = {3};
hsize_t dim6[1] = {6};
hsize_t dimNS1[1] = {NS+1};
hsize_t dimNS1by3[2] = {NS+1, 3};
hsize_t dimNS1by6[2] = {NS+1, 6};


// create the needed types
typeFor1Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim1);
typeFor3Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim3);
typeFor6Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim6);
typeFor1Int = H5Tarray_create(H5T_NATIVE_INT, 1, dim1);
typeForNS1Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dimNS1);
typeForNS1by3Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, dimNS1by3);
typeForNS1by6Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, dimNS1by6);

// create the compound datatype for st1, st2 & stp
s1typeId = H5Tcreate(H5T_COMPOUND, sizeof(struct st1));
//s2typeId = H5Tcreate(H5T_COMPOUND, sizeof(struct st2));
//sptypeId = H5Tcreate(H5T_COMPOUND, sizeof(struct stp));

// create the compound datatype for st1
H5Tinsert(s1typeId, "b[3]", HOFFSET(struct st1, b), typeFor3Double);
H5Tinsert(s1typeId, "c[3]", HOFFSET(struct st1, c), typeFor3Double);
H5Tinsert(s1typeId, "ppc", HOFFSET(struct st1, ppc), typeFor1Int);

//// create the compound datatype for st2
//H5Tinsert(s2typeId, "e[3]", HOFFSET(struct st2, e), typeFor3Double);
//H5Tinsert(s2typeId, "f[3]", HOFFSET(struct st2, f), typeFor3Double);
//H5Tinsert(s2typeId, "r", HOFFSET(struct st2, r), typeFor1Double);
//H5Tinsert(s2typeId, "s", HOFFSET(struct st2, s), typeFor1Double);
//H5Tinsert(s2typeId, "te[6]", HOFFSET(struct st2, te), typeFor6Double);
//H5Tinsert(s2typeId, "vi[3]", HOFFSET(struct st2, vi), typeFor3Double);
//H5Tinsert(s2typeId, "j[3]", HOFFSET(struct st2, j), typeFor3Double);
//H5Tinsert(s2typeId, "ms[NS+1]", HOFFSET(struct st2, ms), typeForNS1Double);
//H5Tinsert(s2typeId, "ns[NS+1]", HOFFSET(struct st2, ns), typeForNS1Double);
//H5Tinsert(s2typeId, "os[NS+1]", HOFFSET(struct st2, os), typeForNS1Double);
//H5Tinsert(s2typeId, "vs[NS+1][3]", HOFFSET(struct st2, vs), typeForNS1by3Double);
//H5Tinsert(s2typeId, "ps[NS+1][6]", HOFFSET(struct st2, ps), typeForNS1by6Double);
//
//// create the compound datatype for stp
//H5Tinsert(sptypeId, "i", HOFFSET(struct stp, i), typeFor1Int);
//H5Tinsert(sptypeId, "r[3]", HOFFSET(struct stp, r), typeFor3Double);
//H5Tinsert(sptypeId, "s[3]", HOFFSET(struct stp, s), typeFor3Double);
//H5Tinsert(sptypeId, "v[3]", HOFFSET(struct stp, v), typeFor3Double);
//H5Tinsert(sptypeId, "w[3]", HOFFSET(struct stp, w), typeFor3Double);
//H5Tinsert(sptypeId, "b[3]", HOFFSET(struct stp, b), typeFor3Double);
//H5Tinsert(sptypeId, "ijk", HOFFSET(struct stp, ijk), typeFor1Int);

// and then keep them in hior struct
hior->s1typeId = H5Tcopy(s1typeId);
//hior->s2typeId = H5Tcopy(s2typeId);
//hior->sptypeId = H5Tcopy(sptypeId);

// close the intermediate datatype
H5Tclose(typeFor1Double);
H5Tclose(typeFor3Double);
H5Tclose(typeFor6Double);
H5Tclose(typeFor1Int);
H5Tclose(typeForNS1Double);
H5Tclose(typeForNS1by3Double);
H5Tclose(typeForNS1by6Double);

// close the st1, st2 & stp compound datatype
H5Tclose(s1typeId);
//H5Tclose(s2typeId);
//H5Tclose(sptypeId);

}



// ___________________________________________________________________________
//
// writeSingleField()
//
// AIM : writes a single field array into the group 'group_id' of the file
// fields.h5
// ___________________________________________________________________________
//
void writeSingleField(HeckleIOFields *hiof, float *field, char *fieldname, hid_t group_id)
{
    hid_t plist_id;
    hid_t dset_id;


    // might not be useful : try to use H5P-DEFAULT in the next H5Dcreate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // create the creation property list
    plist_id = H5Pcreate(H5P_DATASET_CREATE);

    // and then the dataset
    dset_id  = H5Dcreate(group_id,              // group handle
                         fieldname,             // name of the field
                         H5T_NATIVE_FLOAT,      // datatype
                         hiof->filespace,       // file dataspace
                         H5P_DEFAULT,           // link creation p. list
                         plist_id,              // creation p.list
                         H5P_DEFAULT);          // access p.list

    // don't need the plist anymore
    H5Pclose(plist_id);

    // now the dataset is ready for us to write in it
    H5Dwrite(dset_id,                           // dataset id
             H5T_NATIVE_FLOAT,                  // the datatype
             hiof->memspace,                    // memory space (on local grid)
             hiof->filespace,                   // file space (on global grid)
             hiof->writePlist,                  // the write property list (MPI)
             field);                            // and of course... the field to write

    // can close the dataset now
    H5Dclose(dset_id);
}



// ___________________________________________________________________________
//
// writeSpecieComponent()
//
// AIM : writes a single component of a given specie coordinate
// in the group groupe_id of the file species.h5
// ___________________________________________________________________________
//
void writeSpecieComponent(HeckleIOSpecies *hios, int ispe, float *field, char *fieldname, hid_t group_id)
{
    hid_t plist_id;
    hid_t dset_id;

    // might not be useful : try to use H5P-DEFAULT in the next H5Dcreate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // create the creation property list
    plist_id = H5Pcreate(H5P_DATASET_CREATE);

    // and then the dataset
    dset_id  = H5Dcreate(group_id,              // group handle
                         fieldname,             // name of the field
                         H5T_NATIVE_FLOAT,      // datatype
                         hios->filespace[ispe], // file dataspace
                         H5P_DEFAULT,           // link creation p. list
                         plist_id,              // creation p.list
                         H5P_DEFAULT);          // access p.list

    // don't need the plist anymore
    H5Pclose(plist_id);

    // now the dataset is ready for us to write in it
    H5Dwrite(dset_id,                           // dataset handle
             H5T_NATIVE_FLOAT,                  // datatype
             hios->memspace[ispe],              // memory dataspace (local # of part)
             hios->filespace[ispe],             // file dataspace (total # of part)
             hios->writePlist,                  // write Plist (MPI)
             field);                            // and the pointer on the part

    // can close the dataset now
    H5Dclose(dset_id);
}



// ___________________________________________________________________________
//
// HeckleIOCloseFieldsFile()
//
// AIM : close the field file
// ___________________________________________________________________________
//
void HeckleIOCloseFieldsFile(HeckleIOFields *self)
{


    H5Fclose(self->FieldsFileID);

}


// ___________________________________________________________________________
//
// HeckleIOCloseSpeciesFile()
//
// AIM : close the species file
// ___________________________________________________________________________
//
void HeckleIOCloseSpeciesFile(HeckleIOSpecies *self)
{


    H5Fclose(self->SpeciesFileID);

}



// ___________________________________________________________________________
//
// HeckleIOCloseRestartsFile()
//
// AIM : close the species file
// ___________________________________________________________________________
//
void HeckleIOCloseRestartsFile(HeckleIORestart *self)
{


    H5Fclose(self->RestartsFileID);

}



// ___________________________________________________________________________
//
// HeckleIOOpenFieldsFile()
//
// AIM : Open the field file
// ___________________________________________________________________________
//
void HeckleIOOpenFieldsFile(HeckleIOFields *self)
{
    hid_t fapl_id;
    MPI_Info info = MPI_INFO_NULL;


    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // then open the fields.h5 file
    self->FieldsFileID = H5Fopen("fields.h5", H5F_ACC_RDWR, fapl_id);

    // now close the plist since not needed
    H5Pclose(fapl_id);
}


// ___________________________________________________________________________
//
// HeckleIOOpenSpeciesFile()
//
// AIM : Open the species file
// ___________________________________________________________________________
//
void HeckleIOOpenSpeciesFile(HeckleIOSpecies *self)
{
    hid_t fapl_id;
    MPI_Info info = MPI_INFO_NULL;


    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // then open the species.h5 file
    self->SpeciesFileID = H5Fopen("species.h5", H5F_ACC_RDWR, fapl_id);

    // now close the plist since not needed
    H5Pclose(fapl_id);
}



// ___________________________________________________________________________
//
// HeckleIOOpenRestartsFile()
//
// AIM : Open the restarts file
// ___________________________________________________________________________
//
void HeckleIOOpenRestartsFile(HeckleIORestart *self)
{
    hid_t fapl_id;
    MPI_Info info = MPI_INFO_NULL;

    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // then open the restarts.h5 file
    self->RestartsFileID = H5Fopen("restarts.h5", H5F_ACC_RDWR, fapl_id);

    // now close the plist since not needed
    H5Pclose(fapl_id);
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
HeckleIOFields *HeckleIOInitFields(struct sti *si, struct stx *sx)
{
    MPI_Info info = MPI_INFO_NULL;
    hid_t sid;
    hid_t aid;
    hid_t fapl_id;
    hsize_t ndims;
    HeckleIOFields *hiof;


    // create the HeckleIOFields structure
    hiof = malloc(sizeof(*hiof));
    if (hiof == NULL) {
        printf("IOError - process %d could not allocate memory for Fields\n", sx->r);
        return NULL;
    }

    // make the file access property list
    // and give it MPI/communicator info
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // now create the HDF5 file
    hiof->FieldsFileID = H5Fcreate("fields.h5",  // filename
                                 H5F_ACC_TRUNC,  // overwrite if exists
                                 H5P_DEFAULT,    // default creat. p. list
                                 fapl_id);       // parallel access p. list

    // now close the plist since not needed
    H5Pclose(fapl_id);

    // create the writing property list, to write in parallel
    hiof->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(hiof->writePlist, H5FD_MPIO_COLLECTIVE);


    // the file is now created...
    // we now create attributes to store useful data in the begining of the file

    // Attribute 'ncells'
    ndims = 3;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hiof->FieldsFileID,          // file ID
                    "ncells",                    // attribute name
                    H5T_NATIVE_INT,              // datatype
                    sid,                         // dataspace
                    H5P_DEFAULT,                 // creat. p. list
                    hiof->writePlist);           // access p. list

    H5Awrite(aid, H5T_NATIVE_INT, &(si->n));
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'domsize'
    ndims = 3;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hiof->FieldsFileID,          // file ID
                    "domsize",                   // attribute name
                    H5T_NATIVE_DOUBLE,           // datatype
                    sid,                         // dataspace
                    H5P_DEFAULT,                 // creat. p. list
                    hiof->writePlist);           // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->l));
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'meshsize'
    ndims = 3;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hiof->FieldsFileID,          // file ID
                    "meshsize",                  // attribute name
                    H5T_NATIVE_DOUBLE,           // datatype
                    sid,                         // dataspace
                    H5P_DEFAULT,                 // creat. p. list
                    hiof->writePlist);           // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->dl));
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'fieldDtDump'
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(hiof->FieldsFileID,          // file ID
                    "fieldDtDump",               // attribute name
                    H5T_NATIVE_INT,              // datatype
                    sid,                         // dataspace
                    H5P_DEFAULT,                 // creat. p. list
                    hiof->writePlist),           // access p. list

    H5Awrite(aid, H5T_NATIVE_INT, &(si->tf));
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'hyperresistivity'
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(hiof->FieldsFileID,          // file ID
                    "hyperresistivity",          // attribute name
                    H5T_NATIVE_DOUBLE,           // datatype
                    sid,                         // dataspace
                    H5P_DEFAULT,                 // creat. p. list
                    hiof->writePlist),           // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->hyvi));
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'resistivity'
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(hiof->FieldsFileID,          // file ID
                    "resistivity",               // attribute name
                    H5T_NATIVE_DOUBLE,           // datatype
                    sid,                         // dataspace
                    H5P_DEFAULT,                 // creat. p. list
                    hiof->writePlist),           // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->rsty));
    H5Aclose(aid);
    H5Sclose(sid);


    createFieldsDataSpaces(hiof, si, sx);

    HeckleIOCloseFieldsFile(hiof);

    return hiof;

}



// ___________________________________________________________________________
//
// HeckleIOInitSpecies()
//
// AIM : Initialize the HeckleIOSpecies module
// ___________________________________________________________________________
//
HeckleIOSpecies *HeckleIOInitSpecies(struct sti *si, struct stx *sx)
{
    MPI_Info info = MPI_INFO_NULL;
    hid_t sid;
    hid_t aid;
    hid_t fapl_id;
    hsize_t ndims;
    HeckleIOSpecies *hios;


    // create the HeckleIOSpecies structure
    hios = malloc(sizeof *hios);
    if (hios == NULL) {
        printf("IOError - process %d could not allocate memory for Species\n", sx->r);
        return NULL;
    }

    // make the file access property list
    // and give it MPI/communicator info
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // now create the HDF5 file
    hios->SpeciesFileID = H5Fcreate("species.h5",      // filename
                                    H5F_ACC_TRUNC,     // overwrite if exists
                                    H5P_DEFAULT,       // default creat. p. list
                                    fapl_id);          // parallel access p. list

    // now close the plist since not needed
    H5Pclose(fapl_id);

    // create the writing property list, to write in parallel
    hios->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(hios->writePlist, H5FD_MPIO_COLLECTIVE);


    // the file is now created...
    // we now create attributes to store useful data in the begining of the file

    // attribute '# of part' need to be written !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // Attribute 'mass'
    ndims = NS+1;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hios->SpeciesFileID,               // file ID
                    "mass",                            // attribute name
                    H5T_NATIVE_DOUBLE,                 // datatype
                    sid,                               // dataspace
                    H5P_DEFAULT,                       // creat. p. list
                    hios->writePlist);                 // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, si->ms);
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'charge'
    ndims = NS+1;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hios->SpeciesFileID,               // file ID
                    "charge",                          // attribute name
                    H5T_NATIVE_DOUBLE,                 // datatype
                    sid,                               // dataspace
                    H5P_DEFAULT,                       // creat. p. list
                    hios->writePlist);                 // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, si->qs);
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'weight'
    ndims = NS+1;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hios->SpeciesFileID,               // file ID
                    "weight",                          // attribute name
                    H5T_NATIVE_DOUBLE,                 // datatype
                    sid,                               // dataspace
                    H5P_DEFAULT,                       // creat. p. list
                    hios->writePlist);                 // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, si->ws);
    H5Aclose(aid);
    H5Sclose(sid);

    // Attribute 'numofmacro'
    ndims = NS+1;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hios->SpeciesFileID,               // file ID
                    "numofmacro",                      // attribute name
                    H5T_NATIVE_INT,                    // datatype
                    sid,                               // dataspace
                    H5P_DEFAULT,                       // creat. p. list
                    hios->writePlist);                 // access p. list

    H5Awrite(aid, H5T_NATIVE_INT, si->ns);
    H5Aclose(aid);
    H5Sclose(sid);


//  createSpeciesDataSpaces(hios, si, sx);

    HeckleIOCloseSpeciesFile(hios);

    return hios;
}



// ___________________________________________________________________________
//
// HeckleIOInitRestarts()
//
// AIM : Initialize the HeckleIORestarts module
// ___________________________________________________________________________
//
HeckleIORestart *HeckleIOInitRestarts(struct sti *si, struct stx *sx)
{
    MPI_Info info = MPI_INFO_NULL;
    hid_t fapl_id;
    HeckleIORestart *hior;


    // create the HeckleIORestarts structure
    hior = malloc(sizeof *hior);
    if (hior == NULL) {
        printf("IOError - process %d could not allocate memory for Fields\n", sx->r);
        return NULL;
    }

    // make the file access property list
    // and give it MPI/communicator info
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // now create the HDF5 file
    hior->RestartsFileID = H5Fcreate("restarts.h5",    // filename
                                 H5F_ACC_TRUNC,        // overwrite if exists
                                 H5P_DEFAULT,          // default creat. p. list
                                 fapl_id);             // parallel access p. list

    // now close the plist since not needed
    H5Pclose(fapl_id);

    // create the writing property list, to write in parallel
    hior->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(hior->writePlist, H5FD_MPIO_COLLECTIVE);

    //createRestartDataSpaces(hior, si, sx); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! should be called at each time... so not in init

    hior->numofcpu = sx->s;

    createRestartDataType(hior, si);

    HeckleIOCloseRestartsFile(hior);

    return hior;
}



// ___________________________________________________________________________
//
// HeckleIODeleteFields()
//
// AIM : close the write property list & the dataspace
// ___________________________________________________________________________
//
void HeckleIODeleteFields(HeckleIOFields *hiof)
{


    H5Pclose(hiof->writePlist);       // close the write property list
    closeFieldDataSpaces(hiof);       // close memspace & filespace

}



// ___________________________________________________________________________
//
// HeckleIODeleteSpecies()
//
// AIM : close the write property list
// ___________________________________________________________________________
//
void HeckleIODeleteSpecies(HeckleIOSpecies *hios)
{


    H5Pclose(hios->writePlist);       // close the write property list
//  closeSpeciesDataSpaces(hios);     // close memspace & filespace !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

}



// ___________________________________________________________________________
//
// HeckleIODeleteRestarts()
//
// AIM : close the write property list
// ___________________________________________________________________________
//
void HeckleIODeleteRestarts(HeckleIORestart *hior)
{


    H5Pclose(hior->writePlist);       // close the write property list
//  closeRestartDataSpaces(hior);       // close dataspace really needed here ?????????????????????????!!!!!!!!!!!!!!!!!!!!!!!!!!!!

}



// ___________________________________________________________________________
//
// writeFields()
//
// AIM : this routine writes the fields (electromagnetic and fluid moments)
// in a HDF5 file.
// ___________________________________________________________________________
//
void writeFields(HeckleIOFields *hiof,
                 struct stx *sx,
                 struct st1 *s1,
                 struct st2 *s2,
                 double time)
{
    int ijk, i, j, k;
    int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int ijkn;
    int npg1[3],npg2[3];
    //int ng1, ng2;
    int c, ispe, l;
    char timestr[2000];
    int n;

    float *B[3], *E[3], *J[3];
    float *Vi[3], *ns[NS+1],*Vs[NS+1][3], *Ps[NS+1][6];


    // open the fields.h5 file
    HeckleIOOpenFieldsFile(hiof);

    // print informations
    if (sx->r == 0) {
        printf("________________ write fields @ t = %11.5f ______\n\n", time);
    }

    // shortcut for # of grid points on g1 and g2
    // in each direction
    for (c = 0; c < 3; c++) {
        npg1[c] = sx->n[c]+1;
        npg2[c] = sx->n[c]+2;
    }

    // shortcuts for total # of points on g1, g2
    // and # number of points to write (hio->n)
    //ng1 = npg1[0]*npg1[1]*npg1[2];
    //ng2 = npg2[0]*npg2[1]*npg2[2];
    n   = hiof->n[0]*hiof->n[1]*hiof->n[2];

    for (c = 0; c < 3; c++) {
        B[c]  = malloc(n * sizeof *B[c]);
        E[c]  = malloc(n * sizeof *E[c]);
        J[c]  = malloc(n * sizeof *J[c]);
        Vi[c] = malloc(n * sizeof *Vi[c]);
    }

    for (ispe = 0; ispe < NS+1; ispe++) {
        ns[ispe] = malloc(n * sizeof *ns[ispe]);

        for (l=0; l < 3; l++) {
            Vs[ispe][l] = malloc(n * sizeof *Vs[ispe][l]);
        }
        for (l=0; l < 6; l++) {
            Ps[ispe][l] = malloc(n * sizeof *Ps[ispe][l]);
        }
    }


    // fill the arrays: nested loops on the subdomain
    // only loop over the points we actually want to write
    // hio->n excludes redundancy for interior domain boundary points
    // (see createpattern)
    for (i = 0; i < hiof->n[0]; i++) {
        for (j = 0; j < hiof->n[1]; j++) {
            for (k = 0; k < hiof->n[2]; k++) {
                /* __ set index on g1 __ */
                ijk = IDX(i, j, k, npg1[0], npg1[1], npg1[2]);

                // index for g1 sub-grid
                ijkn = IDX(i, j, k, hiof->n[0], hiof->n[1], hiof->n[2]);

                /* __ set indexes on g2 __ */
                ijk1 = IDX(i  , j  , k  , npg2[0], npg2[1], npg2[2]);
                ijk2 = IDX(i+1, j  , k  , npg2[0], npg2[1], npg2[2]);
                ijk3 = IDX(i  , j+1, k  , npg2[0], npg2[1], npg2[2]);
                ijk4 = IDX(i+1, j+1, k  , npg2[0], npg2[1], npg2[2]);
                ijk5 = IDX(i  , j  , k+1, npg2[0], npg2[1], npg2[2]);
                ijk6 = IDX(i+1, j  , k+1, npg2[0], npg2[1], npg2[2]);
                ijk7 = IDX(i  , j+1, k+1, npg2[0], npg2[1], npg2[2]);
                ijk8 = IDX(i+1, j+1, k+1, npg2[0], npg2[1], npg2[2]);

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
                for (ispe = 0; ispe < NS+1; ispe++) {
                    ns[ispe][ijkn] = (float)0.125*(s2[ijk1].ns[ispe]
                                                  +s2[ijk2].ns[ispe]
                                                  +s2[ijk3].ns[ispe]
                                                  +s2[ijk4].ns[ispe]
                                                  +s2[ijk5].ns[ispe]
                                                  +s2[ijk6].ns[ispe]
                                                  +s2[ijk7].ns[ispe]
                                                  +s2[ijk8].ns[ispe]);

                    for (c = 0; c < 3; c++) {
                        Vs[ispe][c][ijkn] = (float)0.125*(s2[ijk1].vs[ispe][c]
                                                         +s2[ijk2].vs[ispe][c]
                                                         +s2[ijk3].vs[ispe][c]
                                                         +s2[ijk4].vs[ispe][c]
                                                         +s2[ijk5].vs[ispe][c]
                                                         +s2[ijk6].vs[ispe][c]
                                                         +s2[ijk7].vs[ispe][c]
                                                         +s2[ijk8].vs[ispe][c]);
                    }

                    for (c = 0; c < 6; c++) {
                        Ps[ispe][c][ijkn] = (float)0.125*(s2[ijk1].ps[ispe][c]
                                                         +s2[ijk2].ps[ispe][c]
                                                         +s2[ijk3].ps[ispe][c]
                                                         +s2[ijk4].ps[ispe][c]
                                                         +s2[ijk5].ps[ispe][c]
                                                         +s2[ijk6].ps[ispe][c]
                                                         +s2[ijk7].ps[ispe][c]
                                                         +s2[ijk8].ps[ispe][c]);
                    }
                }
            }
        }
    }


   // now we have all the data, let's create the group
    sprintf(timestr, "%010.5f", time); // name of the group

    hid_t group_id = H5Gcreate(hiof->FieldsFileID,   // group under root
                               timestr,              // time is the group name
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);

    writeSingleField(hiof, B[0], "Bx", group_id);
    writeSingleField(hiof, B[1], "By", group_id);
    writeSingleField(hiof, B[2], "Bz", group_id);

    writeSingleField(hiof, E[0], "Ex", group_id);
    writeSingleField(hiof, E[1], "Ey", group_id);
    writeSingleField(hiof, E[2], "Ez", group_id);

    writeSingleField(hiof, J[0], "Jx", group_id);
    writeSingleField(hiof, J[1], "Jy", group_id);
    writeSingleField(hiof, J[2], "Jz", group_id);

    writeSingleField(hiof, Vi[0], "Vix", group_id);
    writeSingleField(hiof, Vi[1], "Viy", group_id);
    writeSingleField(hiof, Vi[2], "Viz", group_id);


    // density of each species
    char fieldname[20];
    for (ispe = 0; ispe < NS+1; ispe++) {
        sprintf(fieldname, "n%02d", ispe);
        writeSingleField(hiof, ns[ispe], fieldname, group_id);
    }

    // change the group name for V0__, ... to be agree with convention (no x, y, z) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // bulk velocity of each species
    for (ispe = 0; ispe < NS+1; ispe++) {
        sprintf(fieldname, "Vx%02d", ispe);
        writeSingleField(hiof, Vs[ispe][0], fieldname, group_id);
        sprintf(fieldname, "Vy%02d", ispe);
        writeSingleField(hiof, Vs[ispe][1], fieldname, group_id);
        sprintf(fieldname, "Vz%02d", ispe);
        writeSingleField(hiof, Vs[ispe][2], fieldname, group_id);
    }


    // pressure tensor of each species
    for (ispe = 0; ispe < NS+1; ispe++) {
        sprintf(fieldname, "Pxx%02d", ispe);
        writeSingleField(hiof, Ps[ispe][0], fieldname, group_id);

        sprintf(fieldname, "Pxy%02d", ispe);
        writeSingleField(hiof, Ps[ispe][1], fieldname, group_id);

        sprintf(fieldname, "Pxz%02d", ispe);
        writeSingleField(hiof, Ps[ispe][2], fieldname, group_id);

        sprintf(fieldname, "Pyy%02d", ispe);
        writeSingleField(hiof, Ps[ispe][3], fieldname, group_id);

        sprintf(fieldname, "Pyz%02d", ispe);
        writeSingleField(hiof, Ps[ispe][4], fieldname, group_id);

        sprintf(fieldname, "Pzz%02d", ispe);
        writeSingleField(hiof, Ps[ispe][5], fieldname, group_id);
    }

    // now close the group
    H5Gclose(group_id);
    H5Fflush(hiof->FieldsFileID, H5F_SCOPE_GLOBAL);

   // now free the memory
    for (c = 0; c < 3; c++) {
       free(B[c]);
       free(E[c]);
       free(J[c]);
       free(Vi[c]);
    }

    for (ispe = 0; ispe < NS+1; ispe++) {
        free(ns[ispe]);

        for(c = 0; c < 3; c++) {
            free(Vs[ispe][c]);
        }

        for(l = 0; l < 6; l++) {
            free(Ps[ispe][l]);
        }
    }


    HeckleIOCloseFieldsFile(hiof);
}



// ___________________________________________________________________________
//
// writeSpecies()
//
// AIM : this routine writes the species (position, velocity & id)
// in a HDF5 file.
// ___________________________________________________________________________
//
void writeSpecies(HeckleIOSpecies *hios,
                 struct sti *si,
                 struct stx *sx,
                 struct stp *sp[NS+1],
                 double time)
{
    int c, m, s;
    char timestr[2000];
    char fieldname[2000];
    char partgrp[2000];
    float *R[NS+1][3], *V[NS+1][3], *I[NS+1];


    createSpeciesDataSpaces(hios, si, sx);

    HeckleIOOpenSpeciesFile(hios);

    /* __ print informations __ */
    if (sx->r == 0) {
        printf("________________ write species @ t = %11.5f _____\n\n", time);
    }

    for (s = 0; s < NS+1; s++) {
        for (c = 0; c < 3; c++) {
            R[s][c]  = malloc(sx->ns[s] * sizeof *R[s][c]);
            V[s][c]  = malloc(sx->ns[s] * sizeof *V[s][c]);
        }

        I[s] = malloc(sx->ns[s] * sizeof *I[s]);

        /* __ loop on the part. of specie "s" __ */
        for (m = 0; m < sx->ns[s]; m++) {
            I[s][m] = (float)sp[s][m].i;

            for (c = 0; c < 3; c++) {
                R[s][c][m] = (float)sp[s][m].r[c];
                V[s][c][m] = (float)sp[s][m].v[c];
            }
        }

    }

    // now we have all the data, let's create the group
    sprintf(timestr, "%010.5f", time); // name of the group

    hid_t timegroup_id = H5Gcreate(hios->SpeciesFileID,   // group under root
                                   timestr,               // time is the group name
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (s = 1; s < NS+1; s++) {
        // & another group under time group... one for each specie
        sprintf(partgrp, "specie%2d", s); // name of the field
        hid_t group_id = H5Gcreate(timegroup_id,          // group under timegroup
                                   partgrp,               // time is the group name
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (c = 0; c < 3; c++) {
            sprintf(fieldname, "r[%1d]", c); // name of the field
            if (si->ns[s] != 0) writeSpecieComponent(hios, s, R[s][c], fieldname, group_id);
        }

        for (c = 0; c < 3; c++) {
            sprintf(fieldname, "v[%1d]", c); // name of the field
            if (si->ns[s] != 0) writeSpecieComponent(hios, s, V[s][c], fieldname, group_id);
        }

        strcpy(fieldname, "index"); // name of the field
        if (si->ns[s] != 0) writeSpecieComponent(hios, s, I[s], fieldname, group_id);

        H5Gclose(group_id);
    }


    // now close the group
    H5Gclose(timegroup_id);
    H5Fflush(hios->SpeciesFileID, H5F_SCOPE_GLOBAL);

    // now free the memory
    for (s = 0; s < NS+1; s++) {
        for (c = 0; c < 3; c++) {
            free(R[s][c]);
            free(V[s][c]);
        }

        free(I[s]);
    }


    closeSpeciesDataSpaces(hios);

    HeckleIOCloseSpeciesFile(hios);
}



// ___________________________________________________________________________
//
// writeRestarts()
//
// AIM : this routine writes the restarts in a HDF5 file.
// ___________________________________________________________________________
//
void writeRestarts(HeckleIORestart *hior,
                 struct sti *si,
                 struct stx *sx,
                 struct st1 *s1,
                 struct st2 *s2,
                 struct stp *sp[NS+1],
                 struct std sd,
                 double time)
{
    int i, ispe;
    char grpstr[800];
    char timestr[800];
    char stpdsstr[800];
    hid_t timegrpId;
    hid_t *grpId;
    hid_t *st1dsId;
//  hid_t *st2dsId;
//  hid_t *stpdsId[NS+1];


    createRestartDataSpaces(hior, si, sx);

    HeckleIOOpenRestartsFile(hior);

    /* __ print informations __ */
    if (sx->r == 0) {
        printf("________________ write restarts @ t = %11.5f ____\n\n", time);
    }


    // let's create the groups, one for each rank
    grpId = malloc(sx->s * sizeof (hid_t));
    st1dsId = malloc(sx->s * sizeof (hid_t));
//  st2dsId = malloc(sx->s * sizeof (hid_t));
//  for (ispe = 1; ispe < NS+1; ispe++)
//      stpdsId[ispe] = malloc(sx->s * sizeof (hid_t));

    // name of the group
    sprintf(timestr, "%010.5f", time);

    // then create the "time" group
    timegrpId = H5Gcreate(hior->RestartsFileID,   // group under root
                          timestr,                // time is the group name
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    for (i = 0; i < sx->s; i++) {
        // name of the group
        sprintf(grpstr, "node%8d", i);

        grpId[i] = H5Gcreate(timegrpId,           // group under root
                             grpstr,              // index of node as group name
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        st1dsId[i] = H5Dcreate(grpId[i], "st1", hior->s1typeId, hior->s1space[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//      st2dsId[i] = H5Dcreate(grpId[i], "st2", hior->s2typeId, hior->s2space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (ispe = 1; ispe < NS+1; ispe++) {
            sprintf(stpdsstr, "stp[%2d]", ispe); // name of the group
//          stpdsId[ispe][i] = H5Dcreate(grpId[i], stpdsstr, hior->sptypeId, hior->spspace[ispe], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }

    }

    H5Dwrite(st1dsId[sx->r], hior->s1typeId, H5S_ALL, H5S_ALL, hior->writePlist, s1);
//  H5Dwrite(st2dsId[sx->r], hior->s2typeId, H5S_ALL, H5S_ALL, hior->writePlist, s2);

//  for (ispe = 1; ispe < NS+1; ispe++) {
//      H5Dwrite(stpdsId[ispe][sx->r], hior->sptypeId, H5S_ALL, H5S_ALL, hior->writePlist, sp);
//  }


    // flush the file
//  H5Fflush(hior->RestartsFileID, H5F_SCOPE_GLOBAL);

    // close the root group associated to the time
    H5Gclose(timegrpId);

    // now free dataset and group, 1 for each node
    for (i = 0; i < sx->s; i++) {
        H5Gclose(grpId[i]);
        H5Dclose(st1dsId[i]);
//      H5Dclose(st2dsId[i]);
//      for (ispe = 1; ispe < NS+1; ispe++)
//          H5Gclose(stpdsId[ispe][i]);
    }

    free(grpId);
    free(st1dsId);
//  free(st2dsId);
//  for (ispe = 1; ispe < NS+1; ispe++)
//      free(stpdsId[ispe]);

    closeRestartDataSpaces(hior);

    HeckleIOCloseRestartsFile(hior);

}

































/* __ write the xplosed part. _______________________________________________ */
void writepart(struct sti si, struct stx sx, struct stp *sp[NS+1], int it)
{
    float *xw[NS+1], *yw[NS+1], *zw[NS+1];
    float *uw[NS+1], *vw[NS+1], *ww[NS+1];
    int *iw[NS+1];
    int m, s;
    int wi;
    char nfile[16];
    FILE *fp;
    fp = NULL;

    /* __ set current ts __ */
    wi = it/si.tp;

    /* __ print informations __ */
    if (sx.r == 0)
    {
        printf("________________ write p-dump # %4i for node %3i ____\n", wi, sx.r);
        printf("\n");
    }

    /* __ make the file name __ */
    sprintf(nfile, "xp%03i-%04i.dat", sx.r, wi);

    /* __ open the file __ */
    fp = fopen(nfile, "wb");
    if (fp == NULL) printf("problem in opening file %s\n", nfile);

    else
    {
        /* __ loop on the part. __ */
        for (s = 1; s < NS+1; s++)
        {
            /* __ memory allocation __ */
            if (sx.ns[s] != 0)
            {
                iw[s] = (int *)malloc(sx.ns[s]*sizeof(int));

                xw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
                yw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
                zw[s] = (float *)malloc(sx.ns[s]*sizeof(float));

                uw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
                vw[s] = (float *)malloc(sx.ns[s]*sizeof(float));
                ww[s] = (float *)malloc(sx.ns[s]*sizeof(float));
            }
        }

        /* __ fill the arrays : loop on the species __ */
        for (s = 1; s < NS+1; s++)
        {
            /* __ loop on the part. of specie "s" __ */
            for (m = 0; m < sx.ns[s]; m++)
            {
                /* __ fill the buffers __ */
                iw[s][m] = sp[s][m].i;

                xw[s][m] = (float)sp[s][m].r[0];
                yw[s][m] = (float)sp[s][m].r[1];
                zw[s][m] = (float)sp[s][m].r[2];

                uw[s][m] = (float)sp[s][m].v[0];
                vw[s][m] = (float)sp[s][m].v[1];
                ww[s][m] = (float)sp[s][m].v[2];
            }
        }

        /* __ write the file __ */
        fwrite(sx.ns, sizeof(int), NS+1, fp);
        fwrite(si.ns, sizeof(int), NS+1, fp);

        /* __ loop on the part. __ */
        for (s = 1; s < NS+1; s++)
        {
            /* __ write the buffers in the file __ */
            if (sx.ns[s] != 0)
            {
                fwrite(iw[s], sizeof(int), sx.ns[s], fp);

                fwrite(xw[s], sizeof(float), sx.ns[s], fp);
                fwrite(yw[s], sizeof(float), sx.ns[s], fp);
                fwrite(zw[s], sizeof(float), sx.ns[s], fp);

                fwrite(uw[s], sizeof(float), sx.ns[s], fp);
                fwrite(vw[s], sizeof(float), sx.ns[s], fp);
                fwrite(ww[s], sizeof(float), sx.ns[s], fp);
            }
        }

        /* __ loop on the species __ */
        for (s = 1; s < NS+1; s++)
        {
            /* __ clean-up the pointers __ */
            if (sx.ns[s] != 0)
            {
                free(iw[s]);

                free(xw[s]);
                free(yw[s]);
                free(zw[s]);

                free(uw[s]);
                free(vw[s]);
                free(ww[s]);
            }
        }

        /* __ close the file __ */

        fclose(fp);
    }

}
/*===========================================================================*/









/* __ write the time dump ___________________________________________________ */
void writedump(struct sti si, struct stx sx, struct std sd, int it)
{
float tt;
float db, dw, ma, pb,pa, pe, ab, aa, ae, ea, ee, fx, fy, fz;
FILE *fp;


/* __ only on node 0 __ */
if (sx.r == 0)
   {
   tt = (float)it*si.ts;
   db = (float)sd.db;
   dw = (float)sd.dw;
   ma = (float)sd.ma;
   pb = (float)sd.pb;
   pa = (float)sd.pa;
   pe = (float)sd.pe;
   ab = (float)sd.ab;
   aa = (float)sd.aa;
   ae = (float)sd.ae;
   ea = (float)sd.ea;
   ee = (float)sd.ee;
   fx = (float)sd.fx;
   fy = (float)sd.fy;
   fz = (float)sd.fz;

   /* __ print informations __ */
   printf("________________ time step : %8d ________________\n", it);
   printf("\n");

   /* __ remove old dump __ */
   if (it == 0) remove("hdump.dat");

   /* __ open the file __ */
   fp = fopen("hdump.dat", "a+b");
   if (fp == NULL) printf("problem in opening file hdump.dat\n");

   /* __ write the file __ */
   fwrite(&tt, sizeof(float), 1, fp);
   fwrite(&db, sizeof(float), 1, fp);
   fwrite(&dw, sizeof(float), 1, fp);
   fwrite(&ma, sizeof(float), 1, fp);
   fwrite(&pb, sizeof(float), 1, fp);
   fwrite(&pa, sizeof(float), 1, fp);
   fwrite(&pe, sizeof(float), 1, fp);
   fwrite(&ab, sizeof(float), 1, fp);
   fwrite(&aa, sizeof(float), 1, fp);
   fwrite(&ae, sizeof(float), 1, fp);
   fwrite(&ea, sizeof(float), 1, fp);
   fwrite(&ee, sizeof(float), 1, fp);
   fwrite(&fx, sizeof(float), 1, fp);
   fwrite(&fy, sizeof(float), 1, fp);
   fwrite(&fz, sizeof(float), 1, fp);

   /* __ close the file __ */
   fclose(fp);
   }

}






/* __ write the xplosed orbits ______________________________________________ */
void writeorbit(struct sti si, struct stx *sx, struct st1 *s1, struct st2 *s2,
                struct stp *sp[NS+1], struct sto *so, int it)
{
float xw, yw, zw;
float ro[3], wo[3], bo[3], eo[3], jo[3], io[3];
float no[NS+1], vo[NS+1][3], po[NS+1][6];
float lx, ly, lz;
float w1, w2, w3, w4, w5, w6, w7, w8;
int sw, iw;
int wi;
int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
int m0, m1, m2, n0, n1, n2;
int i, j, k, l, m, n, s;
char nfile[16];
FILE *fp;


/* __ # of grid points __ */
m0 = sx->n[0]+1;
m1 = sx->n[1]+1;
m2 = sx->n[2]+1;
n0 = sx->n[0]+2;
n1 = sx->n[1]+2;
n2 = sx->n[2]+2;

/* __ set the so.no value __ */
so->no = 0;

/* __ identify the orbits to follow (with s & m) : loop on the species __ */
for (s = 1; s < NS+1; s++)
    {
    /* __ loop on all the part. __ */
    for (m = 0; m < sx->ns[s]; m++)
        {
        /* __ loop on all the orbits to follow __ */
        for (n = 0; n < so->wo; n++)
            {
            /* __ if the part. is an orbit __ */
            if (s == so->so[n] && sp[s][m].i == so->io[n])
               {
               /* __ keep memory of s & m __ */
               so->s[so->no] = s;
               so->m[so->no] = m;

               /* __ increase the # of orbits on the node __ */
               so->no++;
               }
            }
        }
    }

/* __ set current ts __ */
wi = it/si.tt;

/* __ print informations __ */
if (sx->r == 0)
   {
   printf("________________ write o-dump # %4i for node %3i ____\n", wi, sx->r);
   printf("\n");
   }

/* __ make the file name __ */
sprintf(nfile, "xo%03i-%04i.dat", sx->r, wi);

/* __ open the file __ */
fp = fopen(nfile, "a+b");
if (fp == NULL) printf("problem in opening file %s\n", nfile);

/* __ write the # of orbits & the current ts __ */
fwrite(&so->no, sizeof(int), 1, fp);
fwrite(&wi, sizeof(int), 1, fp);

/* __ loop on the particles __ */
for (n = 0; n < so->no; n++)
    {
    /* __ part orbit __ */
    sw = so->s[n];
    iw = sp[so->s[n]][so->m[n]].i;

    /* __ interpolate b field __ */
    xw = (sp[so->s[n]][so->m[n]].r[0]/si.dl[0]-sx->i0[0] < sx->n[0]) ?
          sp[so->s[n]][so->m[n]].r[0]/si.dl[0]-sx->i0[0] : sx->n[0]-EPS4;
    yw = (sp[so->s[n]][so->m[n]].r[1]/si.dl[1]-sx->i0[1] < sx->n[1]) ?
          sp[so->s[n]][so->m[n]].r[1]/si.dl[1]-sx->i0[1] : sx->n[1]-EPS4;
    zw = (sp[so->s[n]][so->m[n]].r[2]/si.dl[2]-sx->i0[2] < sx->n[2]) ?
          sp[so->s[n]][so->m[n]].r[2]/si.dl[2]-sx->i0[2] : sx->n[2]-EPS4;

    /* __ index for the part. "position" __ */
    i = (int)floor(xw);
    j = (int)floor(yw);
    k = (int)floor(zw);

    #ifdef BUG
    if (i < 0 || i >= sx->n[0])
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (j < 0 || j >= sx->n[1])
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (k < 0 || k >= sx->n[2])
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    #endif

//  #ifdef BUG
//  if (i < 0 || i >= sx->n[0]) shit(sx->r);
//  if (j < 0 || j >= sx->n[1]) shit(sx->r);
//  if (k < 0 || k >= sx->n[2]) shit(sx->r);
//  #endif

    /* __ part. location in the cell __ */
    lx = xw-i;
    ly = yw-j;
    lz = zw-k;

    /* __ indexes of the rounding grid points on g1 __ */
    ijk1 = IDX(i  , j  , k  , m0, m1, m2);
    ijk2 = IDX(i+1, j  , k  , m0, m1, m2);
    ijk3 = IDX(i  , j+1, k  , m0, m1, m2);
    ijk4 = IDX(i+1, j+1, k  , m0, m1, m2);
    ijk5 = IDX(i  , j  , k+1, m0, m1, m2);
    ijk6 = IDX(i+1, j  , k+1, m0, m1, m2);
    ijk7 = IDX(i  , j+1, k+1, m0, m1, m2);
    ijk8 = IDX(i+1, j+1, k+1, m0, m1, m2);

    /* __ weight for each vertices of the rounding grid points __ */
    w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
    w2 = (    lx)*(1.0-ly)*(1.0-lz);
    w3 = (1.0-lx)*(    ly)*(1.0-lz);
    w4 = (    lx)*(    ly)*(1.0-lz);
    w5 = (1.0-lx)*(1.0-ly)*(    lz);
    w6 = (    lx)*(1.0-ly)*(    lz);
    w7 = (1.0-lx)*(    ly)*(    lz);
    w8 = (    lx)*(    ly)*(    lz);

    /* __ loop in the 3 directions __ */
    for (l = 0; l < 3; l++)
        {
        /* __ set b field seen by the part. __ */
        bo[l] = w1*s1[ijk1].b[l]
               +w2*s1[ijk2].b[l]
               +w3*s1[ijk3].b[l]
               +w4*s1[ijk4].b[l]
               +w5*s1[ijk5].b[l]
               +w6*s1[ijk6].b[l]
               +w7*s1[ijk7].b[l]
               +w8*s1[ijk8].b[l];
        }

    /* __ interpolate e, n, j, v & pe fields __ */
    xw = sp[so->s[n]][so->m[n]].r[0]/si.dl[0]-sx->i0[0]+0.5;
    yw = sp[so->s[n]][so->m[n]].r[1]/si.dl[1]-sx->i0[1]+0.5;
    zw = sp[so->s[n]][so->m[n]].r[2]/si.dl[2]-sx->i0[2]+0.5;

    /* __ index for the part. "position" __ */
    i = (int)floor(xw);
    j = (int)floor(yw);
    k = (int)floor(zw);

    #ifdef BUG
    if (i < 0 || i >= sx->n[0]+1)
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (j < 0 || j >= sx->n[1]+1)
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    if (k < 0 || k >= sx->n[2]+1)
       {
       deadpart(si, sx, s1, s2, sp, s, m, it, __FILE__, __LINE__);
       }
    #endif

//  #ifdef BUG
//  if (i < 0 || i >= sx->n[0]+1) shit(sx->r);
//  if (j < 0 || j >= sx->n[1]+1) shit(sx->r);
//  if (k < 0 || k >= sx->n[2]+1) shit(sx->r);
//  #endif

    /* __ part. location in the cell __ */
    lx = xw-i;
    ly = yw-j;
    lz = zw-k;

    /* __ indexes of the rounding grid points on g2 __ */
    ijk1 = IDX(i  , j  , k  , n0, n1, n2);
    ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
    ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
    ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
    ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
    ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
    ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
    ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

    /* __ weight for each vertices of the rounding grid points __ */
    w1 = (1.0-lx)*(1.0-ly)*(1.0-lz);
    w2 = (    lx)*(1.0-ly)*(1.0-lz);
    w3 = (1.0-lx)*(    ly)*(1.0-lz);
    w4 = (    lx)*(    ly)*(1.0-lz);
    w5 = (1.0-lx)*(1.0-ly)*(    lz);
    w6 = (    lx)*(1.0-ly)*(    lz);
    w7 = (1.0-lx)*(    ly)*(    lz);
    w8 = (    lx)*(    ly)*(    lz);

    /* __ loop in the 3 directions __ */
    for (l = 0; l < 3; l++)
        {
        /* __ set e field seen by the part. __ */
        eo[l] = w1*s2[ijk1].e[l]
               +w2*s2[ijk2].e[l]
               +w3*s2[ijk3].e[l]
               +w4*s2[ijk4].e[l]
               +w5*s2[ijk5].e[l]
               +w6*s2[ijk6].e[l]
               +w7*s2[ijk7].e[l]
               +w8*s2[ijk8].e[l];

        /* __ set e field seen by the part. __ */
        jo[l] = w1*s2[ijk1].j[l]
               +w2*s2[ijk2].j[l]
               +w3*s2[ijk3].j[l]
               +w4*s2[ijk4].j[l]
               +w5*s2[ijk5].j[l]
               +w6*s2[ijk6].j[l]
               +w7*s2[ijk7].j[l]
               +w8*s2[ijk8].j[l];

        /* __ set e field seen by the part. __ */
        io[l] = w1*s2[ijk1].vi[l]
               +w2*s2[ijk2].vi[l]
               +w3*s2[ijk3].vi[l]
               +w4*s2[ijk4].vi[l]
               +w5*s2[ijk5].vi[l]
               +w6*s2[ijk6].vi[l]
               +w7*s2[ijk7].vi[l]
               +w8*s2[ijk8].vi[l];

        /* __ set position of the part. __ */
        ro[l] = (float) sp[so->s[n]][so->m[n]].r[l]+sp[so->s[n]][so->m[n]].b[l]*si.l[l];

        /* __ set velocity of the part. __ */
        wo[l] = (float) sp[so->s[n]][so->m[n]].v[l];
        }

    /* __ loop on the species __ */
    for (s = 0; s < NS+1; s++)
        {
        no[s] = w1*s2[ijk1].ns[s]
               +w2*s2[ijk2].ns[s]
               +w3*s2[ijk3].ns[s]
               +w4*s2[ijk4].ns[s]
               +w5*s2[ijk5].ns[s]
               +w6*s2[ijk6].ns[s]
               +w7*s2[ijk7].ns[s]
               +w8*s2[ijk8].ns[s];

        for (l = 0; l < 3; l++)
            {
            vo[s][l] = w1*s2[ijk1].vs[s][l]
                      +w2*s2[ijk2].vs[s][l]
                      +w3*s2[ijk3].vs[s][l]
                      +w4*s2[ijk4].vs[s][l]
                      +w5*s2[ijk5].vs[s][l]
                      +w6*s2[ijk6].vs[s][l]
                      +w7*s2[ijk7].vs[s][l]
                      +w8*s2[ijk8].vs[s][l];
            }

        for (l = 0; l < 6; l++)
            {
            po[s][l] = w1*s2[ijk1].ps[s][l]
                      +w2*s2[ijk2].ps[s][l]
                      +w3*s2[ijk3].ps[s][l]
                      +w4*s2[ijk4].ps[s][l]
                      +w5*s2[ijk5].ps[s][l]
                      +w6*s2[ijk6].ps[s][l]
                      +w7*s2[ijk7].ps[s][l]
                      +w8*s2[ijk8].ps[s][l];
            }
        }


    /* __ write the file __ */
    fwrite(&sw, sizeof(int), 1, fp);
    fwrite(&iw, sizeof(int), 1, fp);
    fwrite(ro, sizeof(float), 3, fp);
    fwrite(wo, sizeof(float), 3, fp);
    fwrite(bo, sizeof(float), 3, fp);
    fwrite(eo, sizeof(float), 3, fp);
    fwrite(jo, sizeof(float), 3, fp);
    fwrite(io, sizeof(float), 3, fp);
    for (s = 0; s < NS+1; s++) fwrite(&no[s], sizeof(float), 1, fp);
    for (s = 0; s < NS+1; s++) fwrite(vo[s], sizeof(float), 3, fp);
    for (s = 0; s < NS+1; s++) fwrite(po[s], sizeof(float), 6, fp);
    }

/* __ close the file __ */
fclose(fp);

/* __ last time step __ */
if (it == si.nt)
   {
   /* __ clean-up the pointers __ */
   free(so->so);
   free(so->io);
   free(so->s);
   free(so->m);
   }

}


/* _____ write the xplosed restart __________________________________________ */
void writerestart(struct sti si, struct stx sx,
                  struct st1 *s1, struct st2 *s2,
                  struct stp *sp[NS+1], struct std sd,
                  int it)
{
    int nn1, nn2;
    int s;
    char nfile[10];
    FILE *fp;


    for (int s=0; s < NS; s++)
    {
        /* __ total # of part. __ */
        MPI_Allreduce(&(sx.ns[s]), &(si.ns[s]), 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        /* __ max # of part. in domain __ */
        MPI_Allreduce(&(sx.ns[s]), &(sx.nm[s]), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }

    /* __ # of grid points on g1 & g2 __ */
    nn1 = (sx.n[0]+1)*(sx.n[1]+1)*(sx.n[2]+1);
    nn2 = (sx.n[0]+2)*(sx.n[1]+2)*(sx.n[2]+2);

    /* __ print informations __ */
    if (sx.r == 0)
    {
        printf("________________ write restart file for node  %3i ____\n", sx.r);
        printf("\n");
    }

    /* __ make the file name __ */
    if ((it/si.tr) % 2 == 0) sprintf(nfile, "hr%03i.dat", sx.r);
    if ((it/si.tr) % 2 == 1) sprintf(nfile, "hR%03i.dat", sx.r);

    /* __ open the file __ */
    fp = fopen(nfile, "wb");
    if (fp == NULL) printf("problem in opening file %s\n", nfile);

    /* __ write the time step __ */
    fwrite(&it, sizeof(int), 1, fp);

    ///* __ write the sti structure __ */
    //fwrite(&si, sizeof(struct sti), 1, fp);

    /* __ write ns & ws __ */
    for (s = 1; s < NS+1; s++) {
        fwrite(&(si.ns[s]), sizeof(int), 1, fp);
        fwrite(&(si.ws[s]), sizeof(double), 1, fp);
    }

    /* __ write the stx structure __ */
    fwrite(&sx, sizeof(struct stx), 1, fp);

    /* __ write the st1 structure __ */
    fwrite(s1, sizeof(struct st1), nn1, fp);

    /* __ write the st2 structure __ */
    fwrite(s2, sizeof(struct st2), nn2, fp);

    /* __ write the stp structure __ */
    for (s = 1; s < NS+1; s++) {
        fwrite(sp[s], sizeof(struct stp), sx.ns[s], fp);
    }

    /* __ write the initial total energy __ */
    fwrite(&sd.e0, sizeof(double), 1, fp);

    /* __ close the file __ */
    fclose(fp);

}


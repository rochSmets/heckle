
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


#define RUN_FROM_SCRATCH  0
#define RUN_RESTARTING    1


// ___________________________________________________________________________
//
//                             PRIVATE FUNCTIONS
// ___________________________________________________________________________



typedef struct s_HeckleIORestart {

    hid_t fileR;         // identifier for the restarts file to read (reload.h5)
    hid_t fileW;         // identifier for the restarts file to write (restarts.h5)
    hid_t writePlist;    // property list for parallel writing
    hid_t *space0;       // dataspace for s0 in subdomain
    hid_t *space1;       // dataspace for s1 in subdomain
    hid_t *space2;       // dataspace for s2 in subdomain
    hid_t *spacep[NS+1]; // dataspace for sp[NS+1] in subdomain
    hid_t *spacex;       // dataspace for sx in subdomain
    hid_t type0;         // datatype for s0
    hid_t type1;         // datatype for s1
    hid_t type2;         // datatype for s2
    hid_t typep;         // datatype for sp
    hid_t typex;         // datatype for sx
} HeckleIORestart;



// ___________________________________________________________________________
//
// createRestartDataSpaces()
//
// AIM : creates dataspaces for s1, s2, sp a sx struct
// called at each write (# of part can change) & then closed after
// ___________________________________________________________________________
//
void createRestartDataSpaces(HeckleIORestart *ioRestart, const STI* const si,
                             STX* const sx)
{

    hsize_t numCellF0[1];
    hsize_t numCellF1[1];
    hsize_t numCellF2[1];
    hsize_t numCellSp[1];
    hsize_t sizeOne[1] = {1};
    int i;
    int s;
    int nx, ny, nz;
    const int ns1 = NS+1;
    int* sxn;
    int64_t* sxs;


    // we need as many dataspace as ranks
    ioRestart->space0 = malloc(sx->s*sizeof(hid_t));
    ioRestart->space1 = malloc(sx->s*sizeof(hid_t));
    ioRestart->space2 = malloc(sx->s*sizeof(hid_t));

    for (s = 1; s < NS+1; s++) {
        ioRestart->spacep[s] = malloc(sx->s*sizeof(hid_t));
    }

    ioRestart->spacex = malloc(sx->s*sizeof(hid_t));

    // and a buffer to collect the grid size & # of part for each ranks
    sxn = malloc(sx->s*3*sizeof(int64_t));
    sxs = malloc(sx->s*ns1*sizeof(int64_t));

    // gather all values of sx.n & sx.ns for each rank... on each rank
    MPI_Allgather(sx->n, 3, MPI_INT, sxn, 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(sx->ns, ns1, MPI_LONG_LONG_INT, sxs, ns1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);

    // rebuild the size of the dataspace for each ranks
    for (i = 0; i < sx->s; i++) {
        nx = sxn[3*i  ];
        ny = sxn[3*i+1];
        nz = sxn[3*i+2];

        // set s0, s1 & s2 dataspace for each rank
        numCellF0[0] = (nx  )*(ny  )*(nz  );
        numCellF1[0] = (nx+1)*(ny+1)*(nz+1);
        numCellF2[0] = (nx+2)*(ny+2)*(nz+2);

        ioRestart->space0[i] = H5Screate_simple(1, numCellF0, NULL);
        ioRestart->space1[i] = H5Screate_simple(1, numCellF1, NULL);
        ioRestart->space2[i] = H5Screate_simple(1, numCellF2, NULL);

        // set sp dataspace for each rank
        for (s = 1; s < NS+1; s++) {
            numCellSp[0] = sxs[ns1*i+s];

            // not very nice... but write an empty dataset is not acceped !
            // just need to be aware of the fact that sp must be allocated
            if ((int)numCellSp[0] == 0) numCellSp[0] = 1;

            ioRestart->spacep[s][i] = H5Screate_simple(1, numCellSp, NULL);
        }

        // set sx dataspace for each rank
        ioRestart->spacex[i] = H5Screate_simple(1, sizeOne, NULL);
    }

    free(sxn);
    free(sxs);

}



// ___________________________________________________________________________
//
// closeRestartDataSpaces()
//
// AIM : close the dataspaces at the end of each writes
// ___________________________________________________________________________
//
void closeRestartDataSpaces(HeckleIORestart *ioRestart,
                            const STX* const sx)
{
    int i, s;

    // close all dataspaces
    for (i = 0; i < sx->s; i++) {
        H5Sclose(ioRestart->space0[i]);
        H5Sclose(ioRestart->space1[i]);
        H5Sclose(ioRestart->space2[i]);

        for (s = 1; s < NS+1; s++) {
            H5Sclose(ioRestart->spacep[s][i]);
        }

        H5Sclose(ioRestart->spacex[i]);
    }

    // then free the associated pointers
    free(ioRestart->space0);
    free(ioRestart->space1);
    free(ioRestart->space2);

    for (s = 1; s < NS+1; s++) {
        free(ioRestart->spacep[s]);
    }

    free(ioRestart->spacex);

}



// ___________________________________________________________________________
//
// createRestartDataType()
//
// AIM : create once the compound datatype to be used for read & write restarts
// ___________________________________________________________________________
//
void createRestartDataType(HeckleIORestart *ioRestart,
                           const STI* const si)
{

    hid_t type1;
    hid_t type2;
    hid_t typep;
    hid_t typex;
    hid_t typeFor1Double;
    hid_t typeFor3Double;
    hid_t typeFor6Double;
    hid_t typeForNS1Double;
    hid_t typeForNS1by3Double;
    hid_t typeForNS1by6Double;
    hid_t typeFor1Int;
    hid_t typeFor3Int;
    hid_t typeForNS1Int;
    hid_t typeFor27Int;
    hsize_t dim1[1] = {1};
    hsize_t dim3[1] = {3};
    hsize_t dim6[1] = {6};
    hsize_t dimNS1[1] = {NS+1};
    hsize_t dim27[1] = {27};
    hsize_t dimNS1by3[2] = {NS+1, 3};
    hsize_t dimNS1by6[2] = {NS+1, 6};


    // create the needed types
    typeFor1Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim1);
    typeFor3Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim3);
    typeFor6Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim6);
    typeForNS1Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dimNS1);
    typeFor1Int = H5Tarray_create(H5T_NATIVE_INT, 1, dim1);
    typeFor3Int = H5Tarray_create(H5T_NATIVE_INT, 1, dim3);
    typeForNS1Int = H5Tarray_create(H5T_NATIVE_INT64, 1, dimNS1);
    typeFor27Int = H5Tarray_create(H5T_NATIVE_INT, 1, dim27);
    typeForNS1by3Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, dimNS1by3);
    typeForNS1by6Double = H5Tarray_create(H5T_NATIVE_DOUBLE, 2, dimNS1by6);

    // create the compound datatype for st1, st2, stp & stx
    type1 = H5Tcreate(H5T_COMPOUND, sizeof(struct st1));
    type2 = H5Tcreate(H5T_COMPOUND, sizeof(struct st2));
    typep = H5Tcreate(H5T_COMPOUND, sizeof(struct stp));
    typex = H5Tcreate(H5T_COMPOUND, sizeof(struct stx));

    // create the compound datatype for st1
    H5Tinsert(type1, "b[3]", HOFFSET(struct st1, b), typeFor3Double);
    H5Tinsert(type1, "c[3]", HOFFSET(struct st1, c), typeFor3Double);
    #ifdef FULLP
    H5Tinsert(type1, "bKept[3]", HOFFSET(struct st1, bKept), typeFor3Double);
    #endif

    // create the compound datatype for st2
    H5Tinsert(type2, "e[3]", HOFFSET(struct st2, e), typeFor3Double);
    H5Tinsert(type2, "f[3]", HOFFSET(struct st2, f), typeFor3Double);
    H5Tinsert(type2, "r", HOFFSET(struct st2, r), typeFor1Double);
    H5Tinsert(type2, "s", HOFFSET(struct st2, s), typeFor1Double);
    H5Tinsert(type2, "te[6]", HOFFSET(struct st2, te), typeFor6Double);
    H5Tinsert(type2, "vi[3]", HOFFSET(struct st2, vi), typeFor3Double);
    H5Tinsert(type2, "j[3]", HOFFSET(struct st2, j), typeFor3Double);
    H5Tinsert(type2, "ms[NS+1]", HOFFSET(struct st2, ms), typeForNS1Double);
    H5Tinsert(type2, "ns[NS+1]", HOFFSET(struct st2, ns), typeForNS1Double);
    H5Tinsert(type2, "os[NS+1]", HOFFSET(struct st2, os), typeForNS1Double);
    H5Tinsert(type2, "vs[NS+1][3]", HOFFSET(struct st2, vs), typeForNS1by3Double);
    H5Tinsert(type2, "ps[NS+1][6]", HOFFSET(struct st2, ps), typeForNS1by6Double);
    #ifdef FULLP
    H5Tinsert(type2, "pKept[6]", HOFFSET(struct st2, pKept), typeFor6Double);
    H5Tinsert(type2, "dHalf[6]", HOFFSET(struct st2, dHalf), typeFor6Double);
    H5Tinsert(type2, "dFull[6]", HOFFSET(struct st2, dFull), typeFor6Double);
    #endif

    // create the compound datatype for stp
    H5Tinsert(typep, "i", HOFFSET(struct stp, i), typeFor1Int);
    H5Tinsert(typep, "r[3]", HOFFSET(struct stp, r), typeFor3Double);
    H5Tinsert(typep, "s[3]", HOFFSET(struct stp, s), typeFor3Double);
    H5Tinsert(typep, "v[3]", HOFFSET(struct stp, v), typeFor3Double);
    H5Tinsert(typep, "w[3]", HOFFSET(struct stp, w), typeFor3Double);
    H5Tinsert(typep, "b[3]", HOFFSET(struct stp, b), typeFor3Int);
    H5Tinsert(typep, "ijk", HOFFSET(struct stp, ijk), typeFor1Int);

    // create the compound datatype for stx
    H5Tinsert(typex, "s", HOFFSET(struct stx, s), typeFor1Int);
    H5Tinsert(typex, "r", HOFFSET(struct stx, r), typeFor1Int);
    H5Tinsert(typex, "d[3]", HOFFSET(struct stx, d), typeFor3Int);
    H5Tinsert(typex, "ns[NS+1]", HOFFSET(struct stx, ns), typeForNS1Int);
    H5Tinsert(typex, "i0[3]", HOFFSET(struct stx, i0), typeFor3Int);
    H5Tinsert(typex, "i1[3]", HOFFSET(struct stx, i1), typeFor3Int);
    H5Tinsert(typex, "n[3]", HOFFSET(struct stx, n), typeFor3Int);
    H5Tinsert(typex, "l[3]", HOFFSET(struct stx, l), typeFor3Double);
    H5Tinsert(typex, "nt[27]", HOFFSET(struct stx, nt), typeFor27Int);
    H5Tinsert(typex, "nf[27]", HOFFSET(struct stx, nf), typeFor27Int);
    H5Tinsert(typex, "nm[NS+1]", HOFFSET(struct stx, nm), typeForNS1Int);
    H5Tinsert(typex, "irun", HOFFSET(struct stx, irun), typeFor1Int);

    // and then keep them in ioRestart struct
    ioRestart->type0 = H5Tcopy(H5T_NATIVE_INT);
    ioRestart->type1 = H5Tcopy(type1);
    ioRestart->type2 = H5Tcopy(type2);
    ioRestart->typep = H5Tcopy(typep);
    ioRestart->typex = H5Tcopy(typex);

    // close the intermediate datatype
    H5Tclose(typeFor1Double);
    H5Tclose(typeFor3Double);
    H5Tclose(typeFor6Double);
    H5Tclose(typeForNS1Double);
    H5Tclose(typeForNS1by3Double);
    H5Tclose(typeForNS1by6Double);
    H5Tclose(typeFor1Int);
    H5Tclose(typeFor3Int);
    H5Tclose(typeForNS1Int);
    H5Tclose(typeFor27Int);

    // close the st0, st1, st2 & stp compound datatype
    H5Tclose(type1);
    H5Tclose(type2);
    H5Tclose(typep);
    H5Tclose(typex);

}



// ___________________________________________________________________________
//
// HeckleIOOpenRestartsFile()
//
// AIM : Open the restarts file for parallel writing
// ___________________________________________________________________________
//
hid_t HeckleIOOpenRestartsFile(HeckleIORestart *self, char* sfr) {

    hid_t access;
    MPI_Info info = MPI_INFO_NULL;


    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    // then open the restarts.h5 file
    self->fileW = H5Fopen(sfr, H5F_ACC_RDWR, access);//fonction appellee seulement poour s'assurer qu'elle retourne une valeur positive... mais pas gardee en memeoire

    // now close the property list since not needed
    H5Pclose(access);

    return self->fileW;

}






// ___________________________________________________________________________
//
// HeckleIOOpenReloadFile()
//
// AIM : Open the reload file for parallel reading
// ___________________________________________________________________________
//
hid_t HeckleIOOpenReloadFile(HeckleIORestart *self, char* sfr) {

    hid_t access;
    MPI_Info info = MPI_INFO_NULL;


    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    // then open the reload.h5 file
    self->fileR = H5Fopen(sfr, H5F_ACC_RDONLY, access);//fonction appellee seulement poour s'assurer qu'elle retourne une valeur positive... mais pas gardee en memeoire

    // now close the property list since not needed
    H5Pclose(access);

    return self->fileR;

}






// ___________________________________________________________________________
//
//                             PUBLIC FUNCTIONS
// ___________________________________________________________________________



// ___________________________________________________________________________
//
// HeckleIOInitRestarts()
//
// AIM : Initialize the HeckleIORestarts module
// ___________________________________________________________________________
//
HeckleIORestart *HeckleIOInitRestarts(const STI* const si,
                                      const STX* const sx, char *dir)
{
    char sfr[80];
    MPI_Info info = MPI_INFO_NULL;
    hid_t access;
    HeckleIORestart *ioRestart;


    // create the HeckleIORestarts structure
    ioRestart = malloc(sizeof *ioRestart);
    if (ioRestart == NULL) {
        printf("IOError - process %d could not allocate memory for Fields\n", sx->r);
        return NULL;
    }

    // needed because mpi... and default communicator (MPI_COMM_WORLD)
    access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);


    /* __ build the file name __ */
    strcpy(sfr, dir);
    strcat(sfr, "restarts.h5");

    // now create the HDF5 file whatever restart or not
    // this file will be at the same location as all other h5 files
    ioRestart->fileW = H5Fcreate("restarts.h5", H5F_ACC_TRUNC, H5P_DEFAULT, access);

    // for restarted run, need to read reload.h5
    if (si->rst != RUN_FROM_SCRATCH)
        // the restarts read from restarts has to be in the dir provided as entry when starting the code
        ioRestart->fileR = H5Fopen(sfr, H5F_ACC_RDONLY, H5P_DEFAULT);

    // now close the property list since not needed
    H5Pclose(access);

    // create the writing property list, to write in parallel
    ioRestart->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(ioRestart->writePlist, H5FD_MPIO_COLLECTIVE);

    // create the compound datatype
    createRestartDataType(ioRestart, si);

    // close the created file
    H5Fclose(ioRestart->fileW);

    return ioRestart;

}



// ___________________________________________________________________________
//
// HeckleIODeleteRestarts()
//
// AIM : close the write property list
// ___________________________________________________________________________
//
void HeckleIODeleteRestarts(HeckleIORestart *ioRestart) {


    // close the write property mist
    H5Pclose(ioRestart->writePlist);

    // close the restart datatype
    H5Tclose(ioRestart->type0);
    H5Tclose(ioRestart->type1);
    H5Tclose(ioRestart->type2);

    H5Tclose(ioRestart->typep);

    H5Tclose(ioRestart->typex);
}



// ___________________________________________________________________________
//
// writeRestarts()
//
// AIM : this routine writes the restarts in a HDF5 file.
// ___________________________________________________________________________
//
void writeRestarts(HeckleIORestart *ioRestart,
                   STI* si,
                   STX* sx,
                   Grid0 * s0,
                   const ST1 * const s1,
                   const ST2 * const s2,
                   struct stp *sp[NS+1], struct std sd, int it, char *dir)
{

    int i, s;
    char rankName[800];
    char timeString[800];
    char spDsetName[800];
    char sfr[80];
    hid_t timeGroup;
    hid_t *rankGroup;
    hid_t *s0Dset;
    hid_t *s1Dset;
    hid_t *s2Dset;
    hid_t *spDset[NS+1];
    hid_t *sxDset;
    hsize_t nDim;
    hid_t space;
    hid_t attrib;
    int *myS0, n0, ijk;


    // only keep informations for s0.ppc
    n0 = (sx->n[0])*(sx->n[1])*(sx->n[2]);
    myS0 = malloc(n0*sizeof(int));

    for (ijk = 0; ijk < n0; ijk++) {
        myS0[ijk] = s0[ijk].ppc;
    }

    // create the dataspace
    createRestartDataSpaces(ioRestart, si, sx);

    /* __ build the file name __ */
    strcpy(sfr, "restarts.h5");

    // open the "restarts.h5" file
    if (HeckleIOOpenRestartsFile(ioRestart, sfr) > 0) {
    }

    else {
         printf("open restarts.h5 impossible, file not found... : %s\n", sfr);
    }

    /* __ print informations __ */
    if (sx->r == 0) {
        printf("________________ write restarts @ t = %11.5f ____\n\n", it*si->ts);
    }

    // reduce si.ns & sx.nm
    for (s = 0; s < NS+1; s++) {
        if (s == 0) {
            si->ns[0] = 0;
            sx->nm[0] = 0;
        }

        else {
            /* __ total # of part. __ */
            MPI_Allreduce(&(sx->ns[s]), &(si->ns[s]), 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

            /* __ max # of part. in domain __ */
            MPI_Allreduce(&(sx->ns[s]), &(sx->nm[s]), 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
        }
    }

    // let's create the groups, one for each rank
    rankGroup = malloc(sx->s * sizeof (hid_t));

    s0Dset = malloc(sx->s * sizeof (hid_t));
    s1Dset = malloc(sx->s * sizeof (hid_t));
    s2Dset = malloc(sx->s * sizeof (hid_t));

    for (s = 1; s < NS+1; s++)
        spDset[s] = malloc(sx->s * sizeof (hid_t));

    sxDset = malloc(sx->s * sizeof (hid_t));

    // name of the "time" group
    sprintf(timeString, "time_%f", it*si->ts);

    // then create the "time" group
    timeGroup = H5Gcreate(ioRestart->fileW, timeString, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // attribute : "sd.e0" (has to be written just once)
    if (H5Aexists(ioRestart->fileW, "sde0") == 0) {
        nDim = 1;
        space = H5Screate_simple(1, &nDim, NULL);
        attrib = H5Acreate(ioRestart->fileW, "sde0", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attrib, H5T_NATIVE_DOUBLE, &(sd.e0));
        H5Aclose(attrib);
        H5Sclose(space);
    }

    // attribute : "it"
    nDim = 1;
    space = H5Screate_simple(1, &nDim, NULL);
    attrib = H5Acreate(timeGroup, "it", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attrib, H5T_NATIVE_INT, &(it));
    H5Aclose(attrib);
    H5Sclose(space);

    // loop needed to make a group for each cpu
    for (i = 0; i < sx->s; i++) {

        // name of the group
        sprintf(rankName, "rank_%d", i);
        rankGroup[i] = H5Gcreate(timeGroup, rankName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // attribute : "si.ns"
        nDim = NS+1;
        space = H5Screate_simple(1, &nDim, NULL);
        attrib = H5Acreate(rankGroup[i], "si_ns", H5T_NATIVE_INT64, space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attrib, H5T_NATIVE_INT64, si->ns);
        H5Aclose(attrib);
        H5Sclose(space);

        // attribute : "si.ws"
        nDim = NS+1;
        space = H5Screate_simple(1, &nDim, NULL);
        attrib = H5Acreate(rankGroup[i], "si_ws", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attrib, H5T_NATIVE_DOUBLE, si->ws);
        H5Aclose(attrib);
        H5Sclose(space);

        // create dataset for s0, s1, s2, sp[NS+1] & sx
        s0Dset[i] = H5Dcreate(rankGroup[i], "s0", ioRestart->type0, ioRestart->space0[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        s1Dset[i] = H5Dcreate(rankGroup[i], "s1", ioRestart->type1, ioRestart->space1[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        s2Dset[i] = H5Dcreate(rankGroup[i], "s2", ioRestart->type2, ioRestart->space2[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (s = 1; s < NS+1; s++) {
            sprintf(spDsetName, "sp_%d", s);
            spDset[s][i] = H5Dcreate(rankGroup[i], spDsetName, ioRestart->typep, ioRestart->spacep[s][i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        }

        sxDset[i] = H5Dcreate(rankGroup[i], "sx", ioRestart->typex, ioRestart->spacex[i], H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    }

    // then write the datasets
    H5Dwrite(s0Dset[sx->r], ioRestart->type0, H5S_ALL, H5S_ALL, ioRestart->writePlist, myS0);
    H5Dwrite(s1Dset[sx->r], ioRestart->type1, H5S_ALL, H5S_ALL, ioRestart->writePlist, s1);
    H5Dwrite(s2Dset[sx->r], ioRestart->type2, H5S_ALL, H5S_ALL, ioRestart->writePlist, s2);

    for (s = 1; s < NS+1; s++) {
        H5Dwrite(spDset[s][sx->r], ioRestart->typep, H5S_ALL, H5S_ALL, ioRestart->writePlist, sp[s]);
    }

    H5Dwrite(sxDset[sx->r], ioRestart->typex, H5S_ALL, H5S_ALL, ioRestart->writePlist, sx);

    // flush the file
    H5Fflush(ioRestart->fileW, H5F_SCOPE_GLOBAL);

    // close the root group associated to the time
    H5Gclose(timeGroup);

    // now free dataset and group, 1 for each rank
    for (i = 0; i < sx->s; i++) {
        H5Gclose(rankGroup[i]);

        H5Dclose(s0Dset[i]);
        H5Dclose(s1Dset[i]);
        H5Dclose(s2Dset[i]);

        for (s = 1; s < NS+1; s++)
            H5Dclose(spDset[s][i]);

        H5Dclose(sxDset[i]);
    }

    // free the pointer to dataset (for each rank)
    free(rankGroup);

    free(s0Dset);
    free(s1Dset);
    free(s2Dset);

    for (s = 1; s < NS+1; s++)
        free(spDset[s]);

    free(sxDset);

    // close the dataspace
    closeRestartDataSpaces(ioRestart, sx);

    // close the "restarts.h5" file
    H5Fclose(ioRestart->fileW);

    free(myS0);

}



// ___________________________________________________________________________
//
// readRestarts()
//
// AIM : this routine read the restarts in a HDF5 file.
// ___________________________________________________________________________
//
void readRestarts(HeckleIORestart *ioRestart, struct sti *si, struct stx *sx, Grid0 **s0, struct st1 **s1, struct st2 **s2, struct stp *(*sp)[NS+1], struct std *sd, int *it, char *dir) {

    char sfr[80];
    char timeName[800];
    char rankName[800];
    char spName[800];
    int i, s;
    hid_t timeGroup;
    hid_t *rankGroup;
    hid_t attrib;
    hid_t* s0Dset;
    hid_t* s1Dset;
    hid_t* s2Dset;
    hid_t* spDset[NS+1];
    hid_t* sxDset;
    int n1, n2;
    int *myS0, n0, ijk;


    /* __ build the file name __ */
    strcpy(sfr, dir);
    strcat(sfr, "reload.h5");

    // print informations
    printf("________________ read restart file for rank   %3i ____\n\n", sx->r);

    // open the "reload.h5" file
    if (HeckleIOOpenReloadFile(ioRestart, sfr) > 0) {
    }

    else {
         printf("open reload.h5 impossible, file not found...\n");
    }


    // read attribute "sd.e0"
    attrib = H5Aopen(ioRestart->fileR, "sde0", H5P_DEFAULT);
    H5Aread(attrib, H5T_NATIVE_DOUBLE, &((*sd).e0));
    H5Aclose(attrib);

    // name of the group
    sprintf(timeName, "time_%f", si->time4rst);

    //create time group
    timeGroup = H5Gopen(ioRestart->fileR, timeName, H5P_DEFAULT);

    // read attribute "it"
    attrib = H5Aopen(timeGroup, "it", H5P_DEFAULT);
    H5Aread(attrib, H5T_NATIVE_INT, it);
    H5Aclose(attrib);

    // create groups for each rank
    rankGroup = malloc(sx->s*sizeof(hid_t));
    for (i = 0; i < sx->s; i++) {
        sprintf(rankName, "rank_%d", i);
        rankGroup[i] = H5Gopen(timeGroup, rankName, H5P_DEFAULT);
    }

    // read attribute "si.ns[NS+1]"
    attrib = H5Aopen(rankGroup[sx->r], "si_ns", H5P_DEFAULT);
    H5Aread(attrib, H5T_NATIVE_INT64, &(si->ns[0]));
    H5Aclose(attrib);

    // read attribute "si.ws[NS+1]"
    attrib = H5Aopen(rankGroup[sx->r], "si_ws", H5P_DEFAULT);
    H5Aread(attrib, H5T_NATIVE_DOUBLE, &(si->ws[0]));
    H5Aclose(attrib);

    // memory allocation for each dataset (1 for each rank)
    s0Dset =  malloc(sx->s*sizeof(hid_t));
    s1Dset =  malloc(sx->s*sizeof(hid_t));
    s2Dset =  malloc(sx->s*sizeof(hid_t));

    for (s = 1; s < NS+1; s++) {
        spDset[s] = (hid_t *) malloc(sx->s*sizeof(hid_t));
    }

    sxDset = malloc(sx->s*sizeof(hid_t));

    // open each dataset : s1, s2, sp[NS+1] & sx
    for (i = 0; i < sx->s; i++) {
        s0Dset[i] = H5Dopen(rankGroup[i], "s0", H5P_DEFAULT);
        s1Dset[i] = H5Dopen(rankGroup[i], "s1", H5P_DEFAULT);
        s2Dset[i] = H5Dopen(rankGroup[i], "s2", H5P_DEFAULT);

        for (s = 1; s < NS+1; s++) {
            sprintf(spName, "sp_%d", s);
            spDset[s][i] = H5Dopen(rankGroup[i], spName, H5P_DEFAULT);
        }

        sxDset[i] = H5Dopen(rankGroup[i], "sx", H5P_DEFAULT);
    }

    // read sx dataset
    H5Dread(sxDset[sx->r], ioRestart->typex, H5S_ALL, H5S_ALL, H5P_DEFAULT, sx);

    // set the size of memory allocation of s0, s1, s2 & sp[NS+1]
    n0 = (sx->n[0]  )*(sx->n[1]  )*(sx->n[2]  );
    n1 = (sx->n[0]+1)*(sx->n[1]+1)*(sx->n[2]+1);
    n2 = (sx->n[0]+2)*(sx->n[1]+2)*(sx->n[2]+2);

    // only keep informations for s0.ppc
    myS0 = malloc(n0*sizeof(int));

    *s0 =  malloc(n0 * sizeof(struct st0));
    *s1 =  malloc(n1 * sizeof(struct st1));
    *s2 =  malloc(n2 * sizeof(struct st2));

    for (s = 1; s < NS+1; s++) {
        (*sp)[s] =  malloc(si->nm * sizeof(struct stp));
    }

    // read s1, s2 & sp[NS+1] dataset
    H5Dread(s0Dset[sx->r], ioRestart->type0, H5S_ALL, H5S_ALL, H5P_DEFAULT, myS0);
    H5Dread(s1Dset[sx->r], ioRestart->type1, H5S_ALL, H5S_ALL, H5P_DEFAULT, *s1);
    H5Dread(s2Dset[sx->r], ioRestart->type2, H5S_ALL, H5S_ALL, H5P_DEFAULT, *s2);

    for (s = 1; s < NS+1; s++) {
        H5Dread(spDset[s][sx->r], ioRestart->typep, H5S_ALL, H5S_ALL, H5P_DEFAULT, (*sp)[s]);
    }

    for (ijk = 0; ijk < n0; ijk++) {
        (*s0)[ijk].ppc = myS0[ijk];
    }

    // close & release resources
    for (i = 0; i < sx->s; i++) {
        H5Dclose(s0Dset[i]);
        H5Dclose(s1Dset[i]);
        H5Dclose(s2Dset[i]);

        for (s = 1; s < NS+1; s++) {
            H5Dclose(spDset[s][i]);
        }

        H5Dclose(sxDset[i]);
    }

    for (i = 0; i < sx->s; i++) {
        H5Gclose(rankGroup[i]);
    }

    free(rankGroup);

    H5Gclose(timeGroup);
    H5Fclose(ioRestart->fileR);

    // needed to be sure "reload.h5" is not removed by a rank while another one still need it
    MPI_Barrier(MPI_COMM_WORLD);

     
    MPI_Info info= MPI_INFO_NULL;
    hid_t access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(access, MPI_COMM_WORLD, info);

    free(myS0);
}


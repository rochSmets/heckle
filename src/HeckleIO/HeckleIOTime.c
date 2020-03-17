
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

#define NBR_OF_DIAGS (10 + 3*NS)


// ___________________________________________________________________________
//
//                             PRIVATE FUNCTIONS
// ___________________________________________________________________________



typedef struct s_HeckleIOTime {

    hid_t file;                      // identifier for the fields file
    hid_t creationPlist;             // dataset creation property list (to enable chunking)
    hid_t space[NBR_OF_DIAGS];       // space for 1 float
    hid_t dataset[NBR_OF_DIAGS];     // dataset for 10+3*NS float
    hsize_t siz[NBR_OF_DIAGS];       // number of already written fields of time dumps
    char dsetName[NBR_OF_DIAGS][48]; // name of the dataset
} HeckleIOTime;




// ___________________________________________________________________________
//
// createTimeDataSpace()
//
// AIM : creates the space. it has to be extendable, as it will be
// appended at each si.tt. create also the associated dataset
// ___________________________________________________________________________
//
void createTimeDataSpace(HeckleIOTime *ioTime,
                         const STI* const si,
                         const STX* const sx)
{

    hsize_t sizeOne[1] = {1};
    hsize_t sizeMax[1] = {H5S_UNLIMITED};
    hsize_t sizeChunk[1] = {1};
    //const float fillValue = 0.0;
    int i;


    // create a resizeable dataspace
    for (i = 0; i < NBR_OF_DIAGS; i++) {
        ioTime->space[i] = H5Screate_simple(1, sizeOne, sizeMax);
    }

    // modify the dataset creation property list (enable chunking)
    ioTime->creationPlist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(ioTime->creationPlist, 1, sizeChunk);

}



// ___________________________________________________________________________
//
// closeTimeDataSpaces()
//
// AIM : close the dataspace at the end of the run.
// ___________________________________________________________________________
//
void closeTimeDataSpaces(HeckleIOTime *ioTime) {

    int i;


    for (i = 0; i < NBR_OF_DIAGS; i++) {
        H5Sclose(ioTime->space[i]);
    }

}



// ___________________________________________________________________________
//
// HeckleIOOpenTimeFile()
//
// AIM : Open the time.h5 file
// ___________________________________________________________________________
//
hid_t HeckleIOOpenTimeFile(HeckleIOTime *self)
{
    // open the time.h5 file
    self->file = H5Fopen("time.h5", H5F_ACC_RDWR, H5P_DEFAULT);

    return self->file;

}



// ___________________________________________________________________________
//
//                             PUBLIC FUNCTIONS
// ___________________________________________________________________________



// ___________________________________________________________________________
//
// HeckleIOInitTime()
//
// AIM : Initialize the HeckleIOTime module
// ___________________________________________________________________________
//
HeckleIOTime *HeckleIOInitTime(const STI* const si,
                               const STX* const sx)
{

    int i, s;
    HeckleIOTime *ioTime = NULL;

    // only process 0 works on this module
    if (sx->r != 0)
    {
        return NULL;
    }

    else
    {
        ioTime = malloc(sizeof(*ioTime));
        if (ioTime == NULL) {
            printf("IOError - process %d could not allocate memory for Time\n", sx->r);
            return NULL;
        }


        // now create the HDF5 file
        ioTime->file = H5Fcreate("time.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        // create the dataspace
        createTimeDataSpace(ioTime, si, sx);

        // create the dataset names
        strcpy(ioTime->dsetName[0], "Time");
        strcpy(ioTime->dsetName[1], "divB");
        strcpy(ioTime->dsetName[2], "pseudo divB");
        strcpy(ioTime->dsetName[3], "Mag. energy");
        strcpy(ioTime->dsetName[4], "X Elec. Fluctuations");
        strcpy(ioTime->dsetName[5], "Y Elec. Fluctuations");
        strcpy(ioTime->dsetName[6], "Z Elec. Fluctuations");



        for (s = 0; s < NS+1; s++)
        {
            sprintf(ioTime->dsetName[7+3*s], "Bulk Energy sp = %d", s);
            sprintf(ioTime->dsetName[8+3*s], "Parallel Energy sp = %d", s);
            sprintf(ioTime->dsetName[9+3*s], "Perp Energy sp = %d", s);
        }


        // create the dataset
        for (i = 0; i < NBR_OF_DIAGS; i++)
        {
            // create the dataset with the customized dataspace
           ioTime->dataset[i] = H5Dcreate(ioTime->file,
                                           ioTime->dsetName[i],
                                           H5T_NATIVE_FLOAT,
                                           ioTime->space[i],
                                           H5P_DEFAULT,
                                           ioTime->creationPlist,
                                           H5P_DEFAULT);

            // and close the dataset
            H5Dclose(ioTime->dataset[i]);

            // init the # of time dmp written in time.h5
            ioTime->siz[i] = 0;
        }


        H5Fclose(ioTime->file);
        return ioTime;
        }

}



// ___________________________________________________________________________
//
// HeckleIODeleteTime()
//
// AIM : close the write property list & the dataspace
// ___________________________________________________________________________
//
void HeckleIODeleteTime(HeckleIOTime *ioTime)
{
    if (ioTime)
    {
        closeTimeDataSpaces(ioTime);
    }
}



// ___________________________________________________________________________
//
// writeFields()
//
// AIM : this routine writes the fields (electromagnetic and fluid moments)
// in a HDF5 file.
// ___________________________________________________________________________
//
void writeTime(HeckleIOTime *ioTime,
               const STI* const si,
               const STX* const sx,
               const struct std* const sd, int it)
{
    if (ioTime)
    {
        float floatfield[NBR_OF_DIAGS];
        hid_t space[NBR_OF_DIAGS];
        hsize_t offset[NBR_OF_DIAGS];
        hsize_t sizeOne[1] = {1};
        hsize_t sizeAll[1];
        int i, s;


        // open the time.h5 file
        if (HeckleIOOpenTimeFile(ioTime) > 0) {
        }

        else {
            printf("your file time.h5 has a problem my man...\n");
        }

        // print informations
        if (sx->r == 0) {
              printf("________________ time step : %8d ________________\n", it);
              printf("\n");
        }

        // fill the float fields
        floatfield[0] = (float)it*si->ts;
        floatfield[1] = (float)sd->db;
        floatfield[2] = (float)sd->dw;
        floatfield[3] = (float)sd->ma;
        floatfield[4] = (float)sd->fx;
        floatfield[5] = (float)sd->fy;
        floatfield[6] = (float)sd->fz;

        for (s = 0; s < NS+1; s++) {
            floatfield[7+3*s] = (float)sd->pb[s];
            floatfield[8+3*s] = (float)sd->ta[s];
            floatfield[9+3*s] = (float)sd->te[s];
        }


        // the open, extend & write the datasets
        for (i = 0; i < NBR_OF_DIAGS; i++) {

            // open the existing dataset (created with customized datapace)
            ioTime->dataset[i] = H5Dopen(ioTime->file, ioTime->dsetName[i], H5P_DEFAULT);

            // size of the dataset to write : add 1 float to the existing siz[0] ones
            sizeAll[0] = ioTime->siz[i]+1;

            // be sure the dataset has the size to write a single float
            H5Dset_extent(ioTime->dataset[i], sizeAll);

            // get the file dataspace from the dataset
            space[i] = H5Dget_space(ioTime->dataset[i]);

            // set the offset
            offset[i] = ioTime->siz[i];

            // get an hyperslab in dataspace
            H5Sselect_hyperslab(space[i], H5S_SELECT_SET, &offset[i], NULL, sizeOne, NULL);

            // then write the hyperslab
            H5Dwrite(ioTime->dataset[i],
                     H5T_NATIVE_FLOAT,
                     ioTime->space[i], space[i],
                     H5P_DEFAULT,
                     &floatfield[i]);

            // close the dataspace
            H5Sclose(space[i]);

            // update the # of floats already written
            ioTime->siz[i]++;

            // and close the dataset
            H5Dclose(ioTime->dataset[i]);

        }

        H5Fflush(ioTime->file, H5F_SCOPE_GLOBAL);
        H5Fclose(ioTime->file);
    }
}


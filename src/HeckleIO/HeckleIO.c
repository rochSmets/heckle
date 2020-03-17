
#ifndef PYP
#include <mpi.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "structures.h"
//#include "plasma.h"
#include "defines.h"
#include "misc.h"
#include "iamdead.h"
#include <hdf5.h>

/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                          PRIVATE FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */



typedef struct s_HeckleIO
{
    hid_t FieldFileID;		/* hdf5 identifier for the field file      */
    hid_t writePlist;		/* hdf5 property list for parallel writing */
    hid_t memspace;
    hid_t filespace;
    int n[3];               /* # of points to write in each direction  */
} HeckleIO;







/*---------------------------------------------------------------------------
  createpattern()
  ---------------------------------------------------------------------------
  AIM : creates the memory and file dataspaces once and for all, based on
  the local and global size of the g1 grid and the offset of the local grid.
 ---------------------------------------------------------------------------*/
void createpattern(HeckleIO *hio, struct sti *si, struct stx *sx)
{
    hsize_t offset[3];
    hsize_t nloc[3];
    hsize_t nglob[3];
    int idim;

    /* --------------------------------------------------------------------- */
    /*                       MEMORY DATASPACE                                */
    /* --------------------------------------------------------------------- */

     // for 1 process, memspace size == filespace size == whole grid
    if (sx->s == 1)
    {
        for (idim= 0;  idim < 3; idim++)
        {
            nloc[idim]   = sx->n[idim]+1;
            hio->n[idim] = sx->n[idim]+1;
        }
    }

    // now parallel case
    // the default behavior is to write only the first n points in each direction
    // i.e. only the points from 0 to sx->n[direction]-1 (included)
    // only the MPI domains reaching the border at xmax,ymax,zmax will write
    // n+1 points, i.e. from 0 to n[direction] (included) in the direction in
    // which they reach the border of the box.
    else
    {
       // default behavior : only writ the n-1 first points
        hio->n[0] = sx->n[0];
        hio->n[1] = sx->n[1];
        hio->n[2] = sx->n[2];

        nloc[0]   = sx->n[0];
        nloc[1]   = sx->n[1];
        nloc[2]   = sx->n[2];


        // for each direction look whether the current MPI domain
        // reaches the border
        for (idim=0; idim < 3; idim++)
        {
            // if we reach the border
            if (sx->i1[idim] == si->n[idim])
            {
                hio->n[idim] = sx->n[idim]+1;
                nloc[idim]   = sx->n[idim]+1;
            }
        }
    } //end parallel case




    // create the dataspace in memory (memspace)
    hio->memspace = H5Screate_simple(3, nloc, NULL);



    /* --------------------------------------------------------------------- */
    /*                         FILE DATASPACE                                */
    /* --------------------------------------------------------------------- */
    // filespace has the size of the global grid g1
    // the offset of process P is given by the global index sx->i0[0..2]
    // for each direction
    for (idim=0; idim < 3; idim++)
    {
        nglob[idim]  = si->n[idim]+1;
        offset[idim] = sx->i0[idim];
    }


    // dataspace for the file
    hio->filespace = H5Screate_simple(3, nglob, NULL);


    // we need to know where the local grid will be
    // in the file (offset) and how big it is (nloc)
    H5Sselect_hyperslab(hio->filespace,     // file dataspace
                        H5S_SELECT_SET,     // select
                        offset,             // offset
                        NULL,               // stride = 1
                        nloc,               // number of g1 points
                        NULL);              // block = 1 point
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  writeSingleField()
  ---------------------------------------------------------------------------
  AIM : writes a single field array into the group 'group_id' of the file
  fields.h5
 ---------------------------------------------------------------------------*/
void writeSingleField(HeckleIO *hio, float *field, char *fieldname, hid_t group_id)
{
    hid_t plist_id;
    hid_t dset_id;

    // create the dataset (and before its creation property list

    plist_id = H5Pcreate(H5P_DATASET_CREATE);

    dset_id  = H5Dcreate(group_id,              // group handle
                         fieldname,             // name of the field
                         H5T_NATIVE_FLOAT,      // datatype
                         hio->filespace,        // file dataspace
                         H5P_DEFAULT,           // link creation p. list
                         plist_id,              // creation p.list
                         H5P_DEFAULT);          // access p.list

    // don't need the plist anymore
    H5Pclose(plist_id);

    // now the dataset is ready for us to write in it

    H5Dwrite(dset_id,
             H5T_NATIVE_FLOAT,
             hio->memspace,
             hio->filespace,
             hio->writePlist,
             field);

    // can close the dataset now
    H5Dclose(dset_id);
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  HeckleIOCloseFieldFile()
  ---------------------------------------------------------------------------
  AIM : close the field file
 ---------------------------------------------------------------------------*/
void HeckleIOCloseFieldFile(HeckleIO *self)
{
    H5Fclose(self->FieldFileID);
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  HeckleIOOpenFieldFile()
  ---------------------------------------------------------------------------
  AIM : Open the field file
 ---------------------------------------------------------------------------*/
void HeckleIOOpenFieldFile(HeckleIO *self)
{
    hid_t fapl_id;
    MPI_Info info = MPI_INFO_NULL;

    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    self->FieldFileID = H5Fopen("fields.h5", H5F_ACC_RDWR, fapl_id);

    // now close the plist since not needed
    H5Pclose(fapl_id);
}
/*===========================================================================*/







/* ///////////////////////////////////////////////////////////////////////////
   //                                                                       //
   //                           PUBLIC FUNCTIONS                            //
   //                                                                       //
   /////////////////////////////////////////////////////////////////////////// */




/*---------------------------------------------------------------------------
  ioInit()
  ---------------------------------------------------------------------------
  AIM : Initialize the heckleIO module
 ---------------------------------------------------------------------------*/
HeckleIO *HeckleIOInit(struct sti *si, struct stx *sx)
{
    MPI_Info info = MPI_INFO_NULL;
    hid_t sid;
    hid_t aid;
    hid_t fapl_id;
    hsize_t ndims;
    HeckleIO *hio;

    hio = malloc(sizeof *hio);
    if (hio == NULL)
    {
        printf("IOError - process %d could not allocate memory\n", sx->r);
        return NULL;
    }

    // make the file access property list
    // and give it MPI/communicator info
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);

    // now create the HDF5 file

    hio->FieldFileID = H5Fcreate("fields.h5",        // filename
                                 H5F_ACC_TRUNC,     // overwrite if exists
                                 H5P_DEFAULT,       // default creat. p. list
                                 fapl_id);          // parallel access p. list

    // now close the plist since not needed
    H5Pclose(fapl_id);

    // create the writing property list, to write in parallel
    hio->writePlist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(hio->writePlist, H5FD_MPIO_COLLECTIVE);


    /* ok the file is created
       we now create attributes to store useful data in the begining of the file */

    // Attribute 'ncells'
    ndims = 3;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hio->FieldFileID,            // file ID
                    "ncells",                   // attribute name
                    H5T_NATIVE_INT,          // datatype
                    sid,                        // dataspace
                    H5P_DEFAULT,                // creat. p. list
                    H5P_DEFAULT);
                    //hio->writePlist);            // access p. list


    H5Awrite(aid, H5T_NATIVE_INT, &(si->n));
    H5Aclose(aid); // close attribute
    H5Sclose(sid); // close dataspace



    // Attribute 'domsize'
    ndims = 3;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hio->FieldFileID,            // file ID
                    "domsize",                  // attribute name
                    H5T_NATIVE_DOUBLE,          // datatype
                    sid,                        // dataspace
                    H5P_DEFAULT,                // creat. p. list
                                        H5P_DEFAULT);
                    //hio->writePlist);            // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->l));
    H5Aclose(aid); // close attribute
    H5Sclose(sid); // close dataspace



    // Attribute 'meshsize'
    ndims = 3;
    sid = H5Screate_simple(1, &ndims, NULL);
    aid = H5Acreate(hio->FieldFileID,            // file ID
                    "meshsize",                 // attribute name
                    H5T_NATIVE_DOUBLE,          // datatype
                    sid,                        // dataspace
                    H5P_DEFAULT,                // creat. p. list
                                        H5P_DEFAULT);
                    //hio->writePlist);            // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->dl));
    H5Aclose(aid); // close attribute
    H5Sclose(sid); // close dataspace


    // Attribute 'fieldDtDump'
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(hio->FieldFileID,           // file ID
                    "fieldDtDump",              // attribute name
                    H5T_NATIVE_INT,             // datatype
                    sid,                        // dataspace
                    H5P_DEFAULT,                // creat. p. list
                                        H5P_DEFAULT);
                    //hio->writePlist),           // access p. list

    H5Awrite(aid, H5T_NATIVE_INT, &(si->tf));
    H5Aclose(aid); // close attribute
    H5Sclose(sid); // close dataspace



    // Attribute 'hyperresistivity'
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(hio->FieldFileID,           // file ID
                    "hyperresistivity",         // attribute name
                    H5T_NATIVE_DOUBLE,          // datatype
                    sid,                        // dataspace
                    H5P_DEFAULT,                // creat. p. list
                                        H5P_DEFAULT);
                    //hio->writePlist),           // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->hyvi));
    H5Aclose(aid); // close attribute
    H5Sclose(sid); // close dataspace




    // Attribute 'resistivity'
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(hio->FieldFileID,           // file ID
                    "resistivity",              // attribute name
                    H5T_NATIVE_DOUBLE,          // datatype
                    sid,                        // dataspace
                    H5P_DEFAULT,                // creat. p. list
                                        H5P_DEFAULT);
                    //hio->writePlist),           // access p. list

    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(si->rsty));
    H5Aclose(aid); // close attribute
    H5Sclose(sid); // close dataspace



    createpattern(hio, si, sx);

    HeckleIOCloseFieldFile(hio);

    return hio;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  ioDelete()
  ---------------------------------------------------------------------------
  AIM :
 ---------------------------------------------------------------------------*/
void HeckleIODelete(HeckleIO *io)
{
    H5Pclose(io->writePlist);       // close the write property list
    //H5Fclose(io->FieldFileID);      // close the file
}
/*===========================================================================*/










/*---------------------------------------------------------------------------
  writeFields()
  ---------------------------------------------------------------------------
  AIM : this routine writes the fields (electromagnetic and fluid moments)
  in a HDF5 file.
 ---------------------------------------------------------------------------*/
void writeFields(HeckleIO *hio,
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



    HeckleIOOpenFieldFile(hio);


    /* __ print informations __ */
    if (sx->r == 0)
    {
        printf("________________ write f-dump t = %010.5f for node %3i ____\n\n",
               time, sx->r);
    }

    // shortcut for # of grid points on g1 and g2
    // in each direction
    for (c=0; c < 3; c++)
    {
        npg1[c] = sx->n[c]+1;
        npg2[c] = sx->n[c]+2;
    }

    // shortcuts for total # of points on g1, g2
    // and # number of points to write (hio->n)
    //ng1 = npg1[0]*npg1[1]*npg1[2];
    //ng2 = npg2[0]*npg2[1]*npg2[2];
    n   = hio->n[0]*hio->n[1]*hio->n[2];

    for (c=0; c < 3; c++)
    {
        B[c]  = malloc(n * sizeof *B[c]);
        E[c]  = malloc(n * sizeof *E[c]);
        J[c]  = malloc(n * sizeof *J[c]);
        Vi[c] = malloc(n * sizeof *Vi[c]);
    }

    for (ispe=0; ispe < NS+1; ispe++)
    {
        ns[ispe] = malloc(n * sizeof *ns[ispe]);

        for (l=0; l < 3; l++)
        {
            Vs[ispe][l] = malloc(n * sizeof *Vs[ispe][l]);
        }
        for (l=0; l < 6; l++)
        {
            Ps[ispe][l] = malloc(n * sizeof *Ps[ispe][l]);
        }
    }


    /* fill the arrays: nested loops on the subdomain
       only loop over the points we actually want to write
       hio->n excludes redundancy for interior domain boundary points
       (see createpattern) */
    for (i = 0; i < hio->n[0]; i++)
    {
        for (j = 0; j < hio->n[1]; j++)
        {
            for (k = 0; k < hio->n[2]; k++)
            {
                /* __ set index on g1 __ */
                ijk = IDX(i, j, k, npg1[0], npg1[1], npg1[2]);

                // index for g1 sub-grid
                ijkn = IDX(i, j, k, hio->n[0], hio->n[1], hio->n[2]);

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
                for (c = 0; c < 3; c++)
                {
                    B[c][ijkn] = (float)s1[ijk].b[c];
                }

                /* __ for g2 grid species independant quantities __ */
                for (c = 0; c < 3; c++)
                {
                    E[c][ijkn] = (float)0.125*(  s2[ijk1].e[c]
                                               +s2[ijk2].e[c]
                                               +s2[ijk3].e[c]
                                               +s2[ijk4].e[c]
                                               +s2[ijk5].e[c]
                                               +s2[ijk6].e[c]
                                               +s2[ijk7].e[c]
                                               +s2[ijk8].e[c]);

                     J[c][ijkn] = (float)0.125*( s2[ijk1].j[c]
                                               +s2[ijk2].j[c]
                                               +s2[ijk3].j[c]
                                               +s2[ijk4].j[c]
                                               +s2[ijk5].j[c]
                                               +s2[ijk6].j[c]
                                               +s2[ijk7].j[c]
                                               +s2[ijk8].j[c]);

                    Vi[c][ijkn] = (float)0.125*( s2[ijk1].vi[c]
                                               +s2[ijk2].vi[c]
                                               +s2[ijk3].vi[c]
                                               +s2[ijk4].vi[c]
                                               +s2[ijk5].vi[c]
                                               +s2[ijk6].vi[c]
                                               +s2[ijk7].vi[c]
                                               +s2[ijk8].vi[c]);
                }


                /* __ for g2 grid species dependant quantities __ */
                for (ispe = 0; ispe < NS+1; ispe++)
                {
                    ns[ispe][ijkn] = (float)0.125*(  s2[ijk1].ns[ispe]
                                                   +s2[ijk2].ns[ispe]
                                                   +s2[ijk3].ns[ispe]
                                                   +s2[ijk4].ns[ispe]
                                                   +s2[ijk5].ns[ispe]
                                                   +s2[ijk6].ns[ispe]
                                                   +s2[ijk7].ns[ispe]
                                                   +s2[ijk8].ns[ispe]);

                    for (c = 0; c < 3; c++)
                    {
                        Vs[ispe][c][ijkn] = (float)0.125*( s2[ijk1].vs[ispe][c]
                                                         +s2[ijk2].vs[ispe][c]
                                                         +s2[ijk3].vs[ispe][c]
                                                         +s2[ijk4].vs[ispe][c]
                                                         +s2[ijk5].vs[ispe][c]
                                                         +s2[ijk6].vs[ispe][c]
                                                         +s2[ijk7].vs[ispe][c]
                                                         +s2[ijk8].vs[ispe][c]);
                    }

                    for (c = 0; c < 6; c++)
                    {
                        Ps[ispe][c][ijkn] = (float)0.125*( s2[ijk1].ps[ispe][c]
                                                         +s2[ijk2].ps[ispe][c]
                                                         +s2[ijk3].ps[ispe][c]
                                                         +s2[ijk4].ps[ispe][c]
                                                         +s2[ijk5].ps[ispe][c]
                                                         +s2[ijk6].ps[ispe][c]
                                                         +s2[ijk7].ps[ispe][c]
                                                         +s2[ijk8].ps[ispe][c]);
                    }
                } // end loop on species
            } // end loop on k
        } // end loop on j
    } //end loop on i




   // now we have all the data, let's create the group
    sprintf(timestr, "%010.5f", time); // name of the group

    hid_t group_id = H5Gcreate(hio->FieldFileID,     // group under root
                               timestr,             // time is the group name
                               H5P_DEFAULT,
                               H5P_DEFAULT,
                               H5P_DEFAULT);



    writeSingleField(hio, B[0], "Bx", group_id);
    writeSingleField(hio, B[1], "By", group_id);
    writeSingleField(hio, B[2], "Bz", group_id);

    writeSingleField(hio, E[0], "Ex", group_id);
    writeSingleField(hio, E[1], "Ey", group_id);
    writeSingleField(hio, E[2], "Ez", group_id);

    writeSingleField(hio, J[0], "Jx", group_id);
    writeSingleField(hio, J[1], "Jy", group_id);
    writeSingleField(hio, J[2], "Jz", group_id);

    writeSingleField(hio, Vi[0], "Vix", group_id);
    writeSingleField(hio, Vi[1], "Viy", group_id);
    writeSingleField(hio, Vi[2], "Viz", group_id);


    // density of each species
    char fieldname[20];
    for (ispe=0; ispe < NS+1; ispe++)
    {
        sprintf(fieldname, "n%02d", ispe);
        writeSingleField(hio, ns[ispe], fieldname, group_id);
    }

    // bulk velocity of each species
    for (ispe=0; ispe < NS+1; ispe++)
    {
        sprintf(fieldname, "Vx%02d", ispe);
        writeSingleField(hio, Vs[ispe][0], fieldname, group_id);
        sprintf(fieldname, "Vy%02d", ispe);
        writeSingleField(hio, Vs[ispe][1], fieldname, group_id);
        sprintf(fieldname, "Vz%02d", ispe);
        writeSingleField(hio, Vs[ispe][2], fieldname, group_id);
    }


    // pressure tensor of each species
    for (ispe=0; ispe < NS+1; ispe++)
    {
        sprintf(fieldname, "Pxx%02d", ispe);
        writeSingleField(hio, Ps[ispe][0], fieldname, group_id);

        sprintf(fieldname, "Pxy%02d", ispe);
        writeSingleField(hio, Ps[ispe][1], fieldname, group_id);

        sprintf(fieldname, "Pxz%02d", ispe);
        writeSingleField(hio, Ps[ispe][2], fieldname, group_id);

        sprintf(fieldname, "Pyy%02d", ispe);
        writeSingleField(hio, Ps[ispe][3], fieldname, group_id);

        sprintf(fieldname, "Pyz%02d", ispe);
        writeSingleField(hio, Ps[ispe][4], fieldname, group_id);

        sprintf(fieldname, "Pzz%02d", ispe);
        writeSingleField(hio, Ps[ispe][5], fieldname, group_id);
    }

    // now close the group
    H5Gclose(group_id);
    H5Fflush(hio->FieldFileID, H5F_SCOPE_GLOBAL);

   // now free the memory
    for (c=0; c < 3; c++)
    {
       free(B[c]);
       free(E[c]);
       free(J[c]);
       free(Vi[c]);
    }

    for (ispe=0; ispe < NS+1; ispe++)
    {
        free(ns[ispe]);

        for(c=0; c < 3; c++)
        {
            free(Vs[ispe][c]);
        }

        for(l=0; l < 6; l++)
        {
            free(Ps[ispe][l]);
        }
    }


    HeckleIOCloseFieldFile(hio);
}
/*===========================================================================*/













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







/* --------------------------------------------------------------------- */
/*                           OBSOLETE STUFF                              */
/* --------------------------------------------------------------------- */




/*---------------------------------------------------------------------------
  writefield()
  ---------------------------------------------------------------------------
  AIM : obsolete function... version of writefield before HDF5
 ---------------------------------------------------------------------------*/
void obsolete_writefield(struct sti si,
                struct stx sx,
                struct st1 *s1,
                struct st2 *s2,
                int it)
{
    float *bw[3], *ew[3], *jw[3], *iw[3];
    float *nw[NS+1], *vw[NS+1][3], *pw[NS+1][6];
    int i, j, k, l, s;
    int ijk;
    int ijk1, ijk2, ijk3, ijk4, ijk5, ijk6, ijk7, ijk8;
    int m0, m1, m2, n0, n1, n2;
    int wi;
    int nn1;
    char nfile[16];
    FILE *fp;


    /* __ #of grid points __ */
    m0 = sx.n[0]+1;
    m1 = sx.n[1]+1;
    m2 = sx.n[2]+1;
    n0 = sx.n[0]+2;
    n1 = sx.n[1]+2;
    n2 = sx.n[2]+2;

    /* __ set current ts __ */
    wi = it/si.tf;

    /* __ # of grid points on g1 __ */
    nn1 = m0*m1*m2;

    /* __ print informations __ */
    if (sx.r == 0)
    {
        printf("________________ write f-dump # %4i for node %3i ____\n\n",
               wi, sx.r);
    }

    /* __ make the file name __ */
    sprintf(nfile, "xf%03i-%04i.dat", sx.r, wi);

    /* __ open the file __ */
    fp = fopen(nfile, "wb");
    if (fp == NULL) printf("problem in opening file %s\n", nfile);



    /* __ memory allocation __ */
    for (l = 0; l < 3; l++)
    {
        bw[l] = malloc(nn1 * sizeof(float));
        ew[l] = malloc(nn1 * sizeof(float));
        jw[l] = malloc(nn1 * sizeof(float));
        iw[l] = malloc(nn1 * sizeof(float));
    }

    for (s = 0; s < NS+1; s++)
    {
        nw[s] = malloc(nn1*sizeof(float));

        for (l = 0; l < 3; l++)
        {
            vw[s][l] = malloc(nn1*sizeof(float));
        }

        for (l = 0; l < 6; l++)
        {
            pw[s][l] = malloc(nn1*sizeof(float));
        }
    }



    /* __ fill the arrays : nested loops on the subdomain __ */
    for (i = 0; i < m0; i++)
    {
        for (j = 0; j < m1; j++)
        {
            for (k = 0; k < m2; k++)
            {
                /* __ set index on g1 __ */
                ijk = IDX(i, j, k, m0, m1, m2);

                /* __ set indexes on g2 __ */
                ijk1 = IDX(i  , j  , k  , n0, n1, n2);
                ijk2 = IDX(i+1, j  , k  , n0, n1, n2);
                ijk3 = IDX(i  , j+1, k  , n0, n1, n2);
                ijk4 = IDX(i+1, j+1, k  , n0, n1, n2);
                ijk5 = IDX(i  , j  , k+1, n0, n1, n2);
                ijk6 = IDX(i+1, j  , k+1, n0, n1, n2);
                ijk7 = IDX(i  , j+1, k+1, n0, n1, n2);
                ijk8 = IDX(i+1, j+1, k+1, n0, n1, n2);

                /* __ for g1 grid __ */
                for (l = 0; l < 3; l++)
                {
                    bw[l][ijk] = (float)s1[ijk].b[l];
                }

                /* __ for g2 grid species independant quantities __ */
                for (l = 0; l < 3; l++)
                {
                    ew[l][ijk] = (float)0.125*( s2[ijk1].e[l]
                                               +s2[ijk2].e[l]
                                               +s2[ijk3].e[l]
                                               +s2[ijk4].e[l]
                                               +s2[ijk5].e[l]
                                               +s2[ijk6].e[l]
                                               +s2[ijk7].e[l]
                                               +s2[ijk8].e[l]);

                    jw[l][ijk] = (float)0.125*( s2[ijk1].j[l]
                                               +s2[ijk2].j[l]
                                               +s2[ijk3].j[l]
                                               +s2[ijk4].j[l]
                                               +s2[ijk5].j[l]
                                               +s2[ijk6].j[l]
                                               +s2[ijk7].j[l]
                                               +s2[ijk8].j[l]);

                    iw[l][ijk] = (float)0.125*( s2[ijk1].vi[l]
                                               +s2[ijk2].vi[l]
                                               +s2[ijk3].vi[l]
                                               +s2[ijk4].vi[l]
                                               +s2[ijk5].vi[l]
                                               +s2[ijk6].vi[l]
                                               +s2[ijk7].vi[l]
                                               +s2[ijk8].vi[l]);
                }


                /* __ for g2 grid species dependant quantities __ */
                for (s = 0; s < NS+1; s++)
                {
                    nw[s][ijk] = (float)0.125*( s2[ijk1].ns[s]
                                               +s2[ijk2].ns[s]
                                               +s2[ijk3].ns[s]
                                               +s2[ijk4].ns[s]
                                               +s2[ijk5].ns[s]
                                               +s2[ijk6].ns[s]
                                               +s2[ijk7].ns[s]
                                               +s2[ijk8].ns[s]);

                    for (l = 0; l < 3; l++)
                    {
                        vw[s][l][ijk] = (float)0.125*( s2[ijk1].vs[s][l]
                                                      +s2[ijk2].vs[s][l]
                                                      +s2[ijk3].vs[s][l]
                                                      +s2[ijk4].vs[s][l]
                                                      +s2[ijk5].vs[s][l]
                                                      +s2[ijk6].vs[s][l]
                                                      +s2[ijk7].vs[s][l]
                                                      +s2[ijk8].vs[s][l]);
                    }

                    for (l = 0; l < 6; l++)
                    {
                        pw[s][l][ijk] = (float)0.125*( s2[ijk1].ps[s][l]
                                                      +s2[ijk2].ps[s][l]
                                                      +s2[ijk3].ps[s][l]
                                                      +s2[ijk4].ps[s][l]
                                                      +s2[ijk5].ps[s][l]
                                                      +s2[ijk6].ps[s][l]
                                                      +s2[ijk7].ps[s][l]
                                                      +s2[ijk8].ps[s][l]);
                    }
                } // end loop on species
            } // end loop on k
        } // end loop on j
    } //end loop on i



    /* __ write the file __ */
    fwrite(sx.i0, sizeof(int), 3, fp);
    fwrite(sx.n, sizeof(int), 3, fp);

    for (l = 0; l < 3; l++)
    {
        fwrite(bw[l], sizeof(float), nn1, fp);
    }

    for (l = 0; l < 3; l++)
    {
        fwrite(ew[l], sizeof(float), nn1, fp);
    }

    for (l = 0; l < 3; l++)
    {
        fwrite(jw[l], sizeof(float), nn1, fp);
    }

    for (l = 0; l < 3; l++)
    {
        fwrite(iw[l], sizeof(float), nn1, fp);
    }

    for (s = 0; s < NS+1; s++)
    {
        fwrite(nw[s], sizeof(float), nn1, fp);
    }

    for (s = 0; s < NS+1; s++)
    {
        for (l = 0; l < 3; l++)
        {
            fwrite(vw[s][l], sizeof(float), nn1, fp);
        }
    }

    for (s = 0; s < NS+1; s++)
    {
        for (l = 0; l < 6; l++)
        {
            fwrite(pw[s][l], sizeof(float), nn1, fp);
        }
    }


    /* __ clean-up the pointers __ */
    for (l = 0; l < 3; l++)
    {
        free(bw[l]);
        free(ew[l]);
        free(jw[l]);
        free(iw[l]);
    }

    for (s = 0; s < NS+1; s++)
    {
        free(nw[s]);

        for (l = 0; l < 3; l++)
        {
            free(vw[s][l]);
        }

        for (l = 0; l < 6; l++)
        {
            free(pw[s][l]);
        }
    }



    /* __ close the file __ */
    fclose(fp);

}

/*===========================================================================*/






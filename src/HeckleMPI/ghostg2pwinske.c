#ifdef FULLP

#include <stdlib.h>
#include <stdio.h>

#include <structures.h>
#include <defines.h>

#include <mpi.h>

#define GHOSTFIELD_P 0
#define GHOSTFIELD_DRIVER 1

static MPI_Datatype mpighost;

/* ///////////////////////////////////////////////////////////////////////////
 //                                                                       //
 //                          PRIVATE FUNCTIONS                            //
 //                                                                       //
 /////////////////////////////////////////////////////////////////////////// */

/* this structure holds the indices of ghost points
 i,j,k represent the dimensions while 0 and 1 mean
 the first and last points on the border to be transfered*/
typedef struct ghostindices_s {
    int i0, i1;
    int j0, j1;
    int k0, k1;
} GhostIndices;

/* each ghost node has the electron pressure tensor or driver term
 at full timestep.
 */
typedef struct ghostpoint_s {
    double pxx;
    double pyy;
    double pzz;
    double pxy;
    double pxz;
    double pyz;
} GhostPoint;

typedef struct GhostG2winske_s {

    GhostPoint *sendbuff[27];
    GhostPoint *recvbuff[27];

    // number of ghost points to be exchanged
    int counter[27];

    // indices of ghost points
    GhostIndices isend[27];
    GhostIndices irecv[27];

    void (*getVariable)(const ST2 * const s2, const STX * const sx, int i,
        int j, int k, double *pxx, double *pxy, double *pxz, double *pyy,
        double *pyz, double *pzz);

    void (*setVariable)(double pxx, double pxy, double pxz, double pyy,
        double pyz, double pzz, int i, int j, int k, const STX* sx, ST2 *s2);

} GhostG2Winske;

/*---------------------------------------------------------------------------
 findGhostIndices()
 ---------------------------------------------------------------------------
 AIM :
 ---------------------------------------------------------------------------*/
static void findGhostIndices(GhostG2Winske *self, const STX* const sx) {

	int a, b, c;
	int t;

	/* First define the indices of ghost nodes for each neighbor */

	for (a = -1; a <= 1; a++) {
		for (b = -1; b <= 1; b++) {
			for (c = -1; c <= 1; c++) {
				t = (1 + c) + 3 * ((1 + b) + 3 * (1 + a));

				switch (a) {
				case -1:
					self->isend[t].i0 = 1;
					self->isend[t].i1 = 1;
					self->irecv[t].i0 = sx->n[0] + 1;
					self->irecv[t].i1 = sx->n[0] + 1;
					break;

				case 0:
					self->isend[t].i0 = 1;
					self->isend[t].i1 = sx->n[0];
					self->irecv[t].i0 = 1;
					self->irecv[t].i1 = sx->n[0];
					break;

				case 1:
					self->isend[t].i0 = sx->n[0];
					self->isend[t].i1 = sx->n[0];
					self->irecv[t].i0 = 0;
					self->irecv[t].i1 = 0;
					break;
				} // end switch a

				switch (b) {
				case -1:
					self->isend[t].j0 = 1;
					self->isend[t].j1 = 1;
					self->irecv[t].j0 = sx->n[1] + 1;
					self->irecv[t].j1 = sx->n[1] + 1;
					break;

				case 0:
					self->isend[t].j0 = 1;
					self->isend[t].j1 = sx->n[1];
					self->irecv[t].j0 = 1;
					self->irecv[t].j1 = sx->n[1];
					break;

				case 1:
					self->isend[t].j0 = sx->n[1];
					self->isend[t].j1 = sx->n[1];
					self->irecv[t].j0 = 0;
					self->irecv[t].j1 = 0;
					break;
				} // end switch b

				switch (c) {
				case -1:
					self->isend[t].k0 = 1;
					self->isend[t].k1 = 1;
					self->irecv[t].k0 = sx->n[2] + 1;
					self->irecv[t].k1 = sx->n[2] + 1;
					break;

				case 0:
					self->isend[t].k0 = 1;
					self->isend[t].k1 = sx->n[2];
					self->irecv[t].k0 = 1;
					self->irecv[t].k1 = sx->n[2];
					break;

				case 1:
					self->isend[t].k0 = sx->n[2];
					self->isend[t].k1 = sx->n[2];
					self->irecv[t].k0 = 0;
					self->irecv[t].k1 = 0;
					break;
				} // end switch c

			} // end loop on c
		} // end loop on b
	} // end loop on a

	// now count how many elements we need to exhange with each neighbor

	// first initialize counter to 0
	for (t = 0; t < 27; t++)
		self->counter[t] = 0;

	// then count...
	for (t = 0; t < 27; t++) {
		self->counter[t] = (self->isend[t].i1 - self->isend[t].i0 + 1)
				* (self->isend[t].j1 - self->isend[t].j0 + 1)
				* (self->isend[t].k1 - self->isend[t].k0 + 1);
	}

}
/*===========================================================================*/

void getPressure(const ST2 * const s2, const STX * const sx, int i, int j,
		int k, double *pxx, double *pxy, double *pxz, double *pyy, double *pyz,
		double *pzz) {

	int ijk, nxg2, nyg2, nzg2;
	nxg2 = sx->n[0] + 2;
	nyg2 = sx->n[1] + 2;
	nzg2 = sx->n[2] + 2;
	ijk = IDX(i, j, k, nxg2, nyg2, nzg2);
	*pxx = s2[ijk].ps[0][0];
	*pxy = s2[ijk].ps[0][1];
	*pxz = s2[ijk].ps[0][2];
	*pyy = s2[ijk].ps[0][3];
	*pyz = s2[ijk].ps[0][4];
	*pzz = s2[ijk].ps[0][5];
}

void getDriver(const ST2 * const s2, const STX * const sx, int i, int j, int k,
		double *pxx, double *pxy, double *pxz, double *pyy, double *pyz,
		double *pzz) {
	int ijk, nxg2, nyg2, nzg2;
	nxg2 = sx->n[0] + 2;
	nyg2 = sx->n[1] + 2;
	nzg2 = sx->n[2] + 2;
	ijk = IDX(i, j, k, nxg2, nyg2, nzg2);
	*pxx = s2[ijk].dFull[0];
	*pxy = s2[ijk].dFull[1];
	*pxz = s2[ijk].dFull[2];
	*pyy = s2[ijk].dFull[3];
	*pyz = s2[ijk].dFull[4];
	*pzz = s2[ijk].dFull[5];
}

void setPressure(double pxx, double pxy, double pxz, double pyy, double pyz,
		double pzz, int i, int j, int k, const STX* sx, ST2 *s2) {
	int ijk, nxg2, nyg2, nzg2;
	nxg2 = sx->n[0] + 2;
	nyg2 = sx->n[1] + 2;
	nzg2 = sx->n[2] + 2;
	ijk = IDX(i, j, k, nxg2, nyg2, nzg2);
	s2[ijk].ps[0][0] = pxx;
	s2[ijk].ps[0][1] = pxy;
	s2[ijk].ps[0][2] = pxz;
	s2[ijk].ps[0][3] = pyy;
	s2[ijk].ps[0][4] = pyz;
	s2[ijk].ps[0][5] = pzz;
}

void setDriver(double pxx, double pxy, double pxz, double pyy, double pyz,
		double pzz, int i, int j, int k, const STX* sx, ST2 *s2) {
	int ijk, nxg2, nyg2, nzg2;
	nxg2 = sx->n[0] + 2;
	nyg2 = sx->n[1] + 2;
	nzg2 = sx->n[2] + 2;
	ijk = IDX(i, j, k, nxg2, nyg2, nzg2);
	s2[ijk].dFull[0] = pxx;
	s2[ijk].dFull[1] = pxy;
	s2[ijk].dFull[2] = pxz;
	s2[ijk].dFull[3] = pyy;
	s2[ijk].dFull[4] = pyz;
	s2[ijk].dFull[5] = pzz;
}

/* ///////////////////////////////////////////////////////////////////////////
 //                                                                       //
 //                           PUBLIC FUNCTIONS                            //
 //                                                                       //
 /////////////////////////////////////////////////////////////////////////// */

/*---------------------------------------------------------------------------
 GhostG2WinskeInit()
 ---------------------------------------------------------------------------
 AIM :  Initialize the module for ghost communications (P)
 ---------------------------------------------------------------------------*/
GhostG2Winske* GhostG2WinskeInit(const STX * const sx, int fieldID) {
	int t;
	GhostG2Winske* self = malloc(sizeof *self);

	// determine the indices of ghost nodes
	findGhostIndices(self, sx);

	// at this point, we know ghost node indices and how many of them we have
	// so we can allocate memory for buffers for incoming/outgoing data

	for (t = 0; t < 27; t++) {
		self->sendbuff[t] = malloc(
				self->counter[t] * sizeof *self->sendbuff[t]);
		self->recvbuff[t] = malloc(
				self->counter[t] * sizeof *self->recvbuff[t]);
	}

	switch (fieldID) {

	case GHOSTFIELD_P:
		self->getVariable = getPressure;
		self->setVariable = setPressure;
		break;

	case GHOSTFIELD_DRIVER:
		self->getVariable = getDriver;
		self->setVariable = setDriver;
		break;

	default:
		printf("error - ghostfield ID %d invalid\n", fieldID);
		MPI_Abort(MPI_COMM_WORLD, -1);
		exit(-1);
	}

	// now define the MPI Datatype
	int count;
	int blocklength[6];
	MPI_Aint displacements[6];
	MPI_Datatype types[6];
	MPI_Aint startaddress, tmpaddress;
	GhostPoint gptmp;

	count = 6;

	blocklength[0] = 1; //  pxx
	blocklength[1] = 1; //  pyy
	blocklength[2] = 1; //  pzz
	blocklength[3] = 1; //  pxy
	blocklength[4] = 1; //  pxz
	blocklength[5] = 1; //  pyz

	MPI_Get_address(&gptmp, &startaddress);
	MPI_Get_address(&gptmp.pxx, &tmpaddress);
	displacements[0] = tmpaddress - startaddress;

	MPI_Get_address(&gptmp.pyy, &tmpaddress);
	displacements[1] = tmpaddress - startaddress;

	MPI_Get_address(&gptmp.pzz, &tmpaddress);
	displacements[2] = tmpaddress - startaddress;

	MPI_Get_address(&gptmp.pxy, &tmpaddress);
	displacements[3] = tmpaddress - startaddress;

	MPI_Get_address(&gptmp.pxz, &tmpaddress);
	displacements[4] = tmpaddress - startaddress;

	MPI_Get_address(&gptmp.pyz, &tmpaddress);
	displacements[5] = tmpaddress - startaddress;

	types[0] = MPI_DOUBLE;
	types[1] = MPI_DOUBLE;
	types[2] = MPI_DOUBLE;
	types[3] = MPI_DOUBLE;
	types[4] = MPI_DOUBLE;
	types[5] = MPI_DOUBLE;

	MPI_Type_create_struct(count, blocklength, displacements, types, &mpighost);
	MPI_Type_commit(&mpighost);

	return self;
}
/*===========================================================================*/

/*---------------------------------------------------------------------------
 GhostG2WinskeDelete()
 ---------------------------------------------------------------------------
 AIM : Delete the ghost communication module
 ---------------------------------------------------------------------------*/
void GhostG2WinskeDelete(GhostG2Winske* self) {
	int t;

	if (self) {
		for (t = 0; t < 27; t++) {
			free(self->sendbuff[t]);
			free(self->recvbuff[t]);
		}
		free(self);
	}
}
/*===========================================================================*/

/*---------------------------------------------------------------------------
 GhostG2WinskeSendRecv()
 ---------------------------------------------------------------------------
 AIM : Send and receive ghost points
 ---------------------------------------------------------------------------*/
void GhostG2WinskeSendRecv(GhostG2Winske* self, const STX * const sx, ST2 * s2) {
	int t;
	int i, j, k;
	int is0, is1, js0, js1, ks0, ks1;
	int ir0, ir1, jr0, jr1, kr0, kr1;
	int cpt = 0;
	MPI_Status st;

	for (t = 0; t < 27; t++) {
		cpt = 0;

		is0 = self->isend[t].i0;
		is1 = self->isend[t].i1;
		js0 = self->isend[t].j0;
		js1 = self->isend[t].j1;
		ks0 = self->isend[t].k0;
		ks1 = self->isend[t].k1;
		// fill the buffer
		for (i = is0; i <= is1; i++) {
			for (j = js0; j <= js1; j++) {
				for (k = ks0; k <= ks1; k++) {
					self->getVariable(s2, sx, i, j, k,
							&self->sendbuff[t][cpt].pxx,
							&self->sendbuff[t][cpt].pxy,
							&self->sendbuff[t][cpt].pxz,
							&self->sendbuff[t][cpt].pyy,
							&self->sendbuff[t][cpt].pyz,
							&self->sendbuff[t][cpt].pzz);
					cpt++; // one more ghost point stored
				} // end k loop
			} // end j loop
		} // end i loop
	} // end neighbor loop

	// ok now all ghost nodes have been stored in the outgoing buffers
	// proceed to sending to neighbors

	for (t = 0; t < 27; t++) {
		if (t != 13) // I'm not sending to myself
				{
			MPI_Sendrecv(self->sendbuff[t],         // send buffer
					self->counter[t],          // # of ghost nodes to send
					mpighost,                  // MPI datatype
					sx->nt[t],                 // destination process,
					t,                         // message send tag,
					self->recvbuff[t],         // receiving buffer
					self->counter[t],          // # of ghost nodes to recv
					mpighost,                  // recv datatype
					sx->nf[t],                 // reception process,
					t,                         // message recv tag
					MPI_COMM_WORLD, &st);

		}
	} //end neighbor loop

	// ok now all the received data is stored in recvbuff
	// we need to unpack it in

	// I've not received anything from myself, so avoid unpacking t=13 recvbuffs
	// and also do not unpack data if receiving neighbor is MPI_PROC_NULL
	// because that is the border of the simulation domain and those points
	// will be fixed by boundary conditions
	for (t = 0; t < 27; t++) {
		if (t != 13 && sx->nf[t] != MPI_PROC_NULL) {
			cpt = 0;

			ir0 = self->irecv[t].i0;
			ir1 = self->irecv[t].i1;
			jr0 = self->irecv[t].j0;
			jr1 = self->irecv[t].j1;
			kr0 = self->irecv[t].k0;
			kr1 = self->irecv[t].k1;

			for (i = ir0; i <= ir1; i++) {
				for (j = jr0; j <= jr1; j++) {
					for (k = kr0; k <= kr1; k++) {
						self->setVariable(self->recvbuff[t][cpt].pxx,
								self->recvbuff[t][cpt].pxy,
								self->recvbuff[t][cpt].pxz,
								self->recvbuff[t][cpt].pyy,
								self->recvbuff[t][cpt].pyz,
								self->recvbuff[t][cpt].pzz, i, j, k, sx, s2);
						cpt++; // oe more ghost node unpacked
					} // end k loop
				} // end j loop
			} // end i loop
		} // end of not 13 and not proc_null
	} // end neighbor loop
}
/*===========================================================================*/
#endif


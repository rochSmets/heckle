/* --------------------------------------------------------------------- */
/*  MODULE             : MEMORY                                          */
/*  Short Description  : prints an estimation of the memory consumption  */
/* --------------------------------------------------------------------- */


#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>

#include <mpi.h>
#include <structures.h>
#include <particle.h>
#include <particlecomm.h>




typedef struct s_Memory
{
    size_t bytes;
    size_t KB;
    size_t MB;
    size_t GB;
}Memory;



static Memory *PerProcMem;
static Memory *GlobMem;






/*---------------------------------------------------------------------------
  memory_init()
  ---------------------------------------------------------------------------
  AIM : initialize the memory module
 ---------------------------------------------------------------------------*/
Memory *memory_init(void)
{
    Memory *mem;
    mem = malloc(sizeof*mem);
    mem->bytes = 0;
    mem->KB    = 0;
    mem->MB    = 0;
    mem->GB    = 0;

    return mem;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  memory_delete()
  ---------------------------------------------------------------------------
  AIM : deletes the memory object
 ---------------------------------------------------------------------------*/
void memory_delete(Memory* mem)
{
    free(mem);
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  memory_start()
  ---------------------------------------------------------------------------
  AIM : starts the memory tracking
 ---------------------------------------------------------------------------*/
void memory_start(void)
{
    PerProcMem = memory_init();
    GlobMem    = memory_init();
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  memory_stop()
  ---------------------------------------------------------------------------
  AIM : stops the memory tracking
 ---------------------------------------------------------------------------*/
void memory_stop(void)
{
    memory_delete(PerProcMem);
    memory_delete(GlobMem);
}
/*===========================================================================*/









/*---------------------------------------------------------------------------
  memory_sum()
  ---------------------------------------------------------------------------
  AIM : add operator for memory objects
 ---------------------------------------------------------------------------*/
void memory_sum(const Memory* const mem1, const Memory* const mem2, Memory *res)
{
    res->bytes = mem1->bytes + mem2->bytes;
    res->KB    = mem1->KB + mem2->KB;
    res->MB    = mem1->MB + mem2->MB;
    res->GB    = mem1->GB + mem2->GB;
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  memory_add()
  ---------------------------------------------------------------------------
  AIM : add X bytes to the memory
 ---------------------------------------------------------------------------*/
void memory_add(Memory *mem, size_t bytes)
{
    mem->bytes += bytes;
    mem->KB     = mem->bytes/1024.;
    mem->MB     = mem->KB/1024.;
    mem->GB     = mem->MB/1024.;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  memory_remove()
  ---------------------------------------------------------------------------
  AIM : remove X bytes from the memory
 ---------------------------------------------------------------------------*/
void memory_remove(Memory *mem, size_t bytes)
{
    mem->bytes -= bytes;
    mem->KB     = mem->bytes/1024.;
    mem->MB     = mem->KB/1024.;
    mem->GB     = mem->MB/1024.;
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  memory_allocate()
  ---------------------------------------------------------------------------
  AIM : allocate memory
 ---------------------------------------------------------------------------*/
void *memory_allocate(size_t bytes)
{
    void *ptr = malloc(bytes);

    if (ptr != NULL)
    {
        memory_add(PerProcMem, bytes);
    }

    return ptr;
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  memory_free()
  ---------------------------------------------------------------------------
  AIM : free memory
 ---------------------------------------------------------------------------*/
void memory_free(void *ptr, size_t bytes)
{
    if (ptr != NULL)
    {
        free(ptr);
        memory_remove(PerProcMem, bytes);
    }
}
/*===========================================================================*/






/*---------------------------------------------------------------------------
  memory_print()
  ---------------------------------------------------------------------------
  AIM : prints the memory consumption
 ---------------------------------------------------------------------------*/
//void memory_print(const Memory *const mem)
void memory_print(Memory * mem)
{
    int myrank;
//  int irank;
//  int nproc;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int64_t bytes;



    /* __ find the size of the bigest rank __ */
    MPI_Allreduce(&(mem->bytes), &(bytes), 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);

    if (myrank == 0)
        {
        printf("\n");
        printf("________________ memory status _______________________\n");
//      printf("divergence b     :\n");
        printf("RAM needed (Mb)  :%12zu\n", bytes/(1024*1024));
//      printf("\n");
//      printf("______________________________________________________\n");
        fflush(stdout);
        }

//    if (myrank == 0)
//    {
//        printf("________________ memory status _______________________\n");
//        fflush(stdout);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    for (irank = 0; irank < nproc; irank++)
//    {
//        MPI_Barrier(MPI_COMM_WORLD);
//        if (irank == myrank)
//        {
//            printf("cpu # %6d     : %16ld (GB)\n", myrank, mem->GB);
//            printf("cpu # %6d     : %16ld (MB)\n", myrank, mem->MB);
//            printf("cpu # %6d     : %16ld (KB)\n", myrank, mem->KB);
//            printf("cpu # %6d     : %16ld  (B)\n", myrank, mem->bytes);
//
//            fflush(stdout);
//            MPI_Barrier(MPI_COMM_WORLD);
//        }
//    }
//
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (myrank == 0)
//    {
//        printf("______________________________________________________\n");
//        fflush(stdout);
//    }
}
/*===========================================================================*/







/*---------------------------------------------------------------------------
  memory_status()
  ---------------------------------------------------------------------------
  AIM : prints the actual memory consumption
 ---------------------------------------------------------------------------*/
void memory_status(void)
{
    memory_print(PerProcMem);
}
/*===========================================================================*/








/*---------------------------------------------------------------------------
  memory_estimate()
  ---------------------------------------------------------------------------
  AIM : estimates the max memory consumption per process
 ---------------------------------------------------------------------------*/
void memory_estimate(const STI* const si, const STX* const sx)
{
    Memory *mem;
    int ngp1, ngp2, ngp4;
    int npart_bc_storage = 50000;
    int npart_ppack = 50000;

    ngp1 = (sx->n[0]+1) * (sx->n[1]+1) * (sx->n[2]+1);
    ngp2 = (sx->n[0]+2) * (sx->n[1]+2) * (sx->n[2]+2);
    ngp4 = (sx->n[0]+4) * (sx->n[1]+4) * (sx->n[2]+4);

    mem = memory_init();

    memory_add(mem, si->nm * NS * sizeof(Particle));
    memory_add(mem, ngp1 * sizeof(ST1));
    memory_add(mem, ngp2 * sizeof(ST2));
    memory_add(mem, ngp4 * sizeof(ST4));
    memory_add(mem, npart_bc_storage * sizeof (int) * 2);
    memory_add(mem, npart_ppack * 2 * 27 *sizeof (Ppack));
    memory_print(mem);


    memory_delete(mem);
}
/*===========================================================================*/













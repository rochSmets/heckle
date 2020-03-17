//Routines pour les versions "securisees" de malloc et free

// Copyright Joel Stienlet

//void *smalloc(char *fct_appelante,unsigned int input_size)
//void sfree(char *fct_appelante,void *Ptr)

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

//#include <sys/types.h>
//#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <errno.h>
#include <time.h>

#include <signal.h>

#include <string.h>

int was_in_malloc_func;

int nber_allocated_blocs; // compteur de blocs en cours d'allocation

#define qword long long int




/*

   fonction pour allouer une page

*/
void *mmap_alloc(size_t size)
{
void *res;
int pagesize;
pagesize = getpagesize();
int nber_pages;
nber_pages = ((unsigned long long int)size - 1) / pagesize + 1;



int fd;
fd = open("/dev/zero",O_RDWR);//O_RDWR O_RDONLY
if(fd == -1)
        {
        printf("ERROR: could not open /dev/zero\n");
        exit(-1);
        }
res = mmap(NULL,nber_pages * pagesize,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);

if(res == MAP_FAILED)
        {
        printf("ERROR: ");

        if(errno == EBADF)   printf("ebadf: The file descriptor is not a valid open file");
        if(errno == EACCES)  printf("eacces");
        if(errno == EINVAL)  printf("Einval");
        if(errno == ETXTBSY) printf("etxtbsy");
        if(errno == EAGAIN)  printf("Eagain: the mapping could not be locked into memory");
        if(errno == ENOMEM)  printf("enomem: no memory");
        if(errno == EMFILE)  printf("emfile: too many mapped regions");
        if(errno == ENODEV)  printf("enodev: file type not supported by mmap");
        if(errno == ENOTSUP) printf("enotsup");
        if(errno == ENXIO)   printf("ENXIO");
        if(errno == EOVERFLOW) printf("overflow");

        printf("\n");
        exit(-1);
        }

close(fd);

//printf("%d pages allouees de %d octets\n",nber_pages,pagesize);

return res;
}

void mmap_free(void *area, size_t size)
{
int pagesize;
pagesize = getpagesize();
int nber_pages;
nber_pages = ((unsigned long long int)size - 1) / pagesize + 1;

munmap(area,nber_pages * pagesize);
}

void *smalloc(char *fct_appelante, size_t size)
{
//printf("Secure malloc appele\n");
was_in_malloc_func = 1;
int pagesize;
pagesize = getpagesize();

int nber_pages;
nber_pages = ((unsigned long long int)size - 1) / pagesize + 1;

int size_not_used;
size_not_used = (nber_pages * pagesize) - size;

// on essaie d'allouer 3 plages consecutives
void *area1;
void *area2;
void *area3;

void *last_junk_alloc = NULL;
void **tmp_junk;

int is_echec;
int is_echec1;
int is_echec2;
int reverse;

do
        {
        is_echec = 1;
        reverse = 0;
        area1 = mmap_alloc(pagesize); // page de debut
        area2 = mmap_alloc(nber_pages * pagesize); // page utile
        area3 = mmap_alloc(pagesize); // page de fin

        is_echec1 = 1;
        if(((qword)area2 - (qword)area1) == pagesize) is_echec1 = 0;
        if(((qword)area2 - (qword)area1) == -nber_pages *pagesize)
                {
                 is_echec1 = 0;
                reverse = 1;
                }

        is_echec2 = 1;

        if((((qword)area3 - (qword)area2) == nber_pages * pagesize) && (reverse == 0)) is_echec2 = 0;
        if((((qword)area3 - (qword)area2) == -pagesize) && (reverse == 1)) is_echec2 = 0;

        if(!is_echec1 && !is_echec2) is_echec = 0;


        if(is_echec)
                {
                //printf("ECHEC: is_echec1 = %d, is_echec2 = %d\n",is_echec1,is_echec2);

                // free areas:
                mmap_free(area1,pagesize);
                mmap_free(area2,nber_pages * pagesize);
                mmap_free(area3,pagesize);

                tmp_junk = mmap_alloc(pagesize);
                *tmp_junk = last_junk_alloc;
                last_junk_alloc = tmp_junk;
                }
        else
                {
                // free junk allocs:
                if(last_junk_alloc != NULL)
                        {
                        do
                                {
                                tmp_junk = last_junk_alloc;
                                last_junk_alloc = *tmp_junk;
                                mmap_free(tmp_junk, pagesize);
                                }while(last_junk_alloc != NULL);

                        }
                }
}while(is_echec);

void *low_area;

if(reverse)
        {
        low_area = area3;
        }
else
        {
        low_area = area1;
        }

// on ecrit les donnees concernant l'allocation dans low_area:
size_t *taille;
taille = low_area;
*taille = nber_pages * pagesize;

// reste a modifier les droits des pages area1 et area2:
mprotect(area1, pagesize, PROT_NONE);
mprotect(area3, pagesize, PROT_NONE);
// en passant la protection "no access" nous assure que les donnees dans low_area ne seront pas modifiees!

if(drand48() > 0.5)
        {
        //printf("--------------------------- DECALAGE !!!!!!!!!!!!!!!!! ------------------\n");
        area2 = area2 + size_not_used;
        }
else
        {
        //printf("--------------XXXXXXXXXXXXXXXXXxx  DECALAGE !!!!!!!!!!!!!!!!! ------------------\n");
        }

was_in_malloc_func = 0;

nber_allocated_blocs ++;

return area2;
}

// -------------------------------------------------------------------
//
//    pour liberer une zone memoire allouee avec smalloc
//
// -------------------------------------------------------------------
void sfree(char *fct_appelante, void **area_ptr)
{
//printf("Secure free appele\n");
was_in_malloc_func = 1;
if(area_ptr == NULL)
        {
        printf("ERROR: pointeur NULL passe en argument de sfree\n");
        exit(-1);
        }
was_in_malloc_func = 2;

void *area;
area = *area_ptr;

if(area == NULL)
        {
        printf("ERROR: pointeur vers NULL passe en argument de sfree\n");
        exit(-1);
        }

was_in_malloc_func = 3;

// -------- a partir d'ici on sait que le pointeur transmis est bon en principe

int pagesize;
pagesize = getpagesize();

// on cherche la plus petite adresse alignee sur une page:
void *real_area;
real_area = (void *)(((qword)area / pagesize) * pagesize);

void *data_area;
data_area = real_area - pagesize;

mprotect(data_area, pagesize, PROT_READ);

was_in_malloc_func = 4;

size_t *taille;
taille = data_area;

void *upper_area;
upper_area = real_area + *taille;

was_in_malloc_func = 5;

mmap_free(upper_area, pagesize);
mmap_free(real_area, *taille);
mmap_free(data_area, pagesize);

was_in_malloc_func = 6;

*area_ptr = NULL;
was_in_malloc_func = 0;

nber_allocated_blocs --;
}


//
//
//      pour reallouer une zone memoire
//
//
void *srealloc(char *fct_appelante, void *area_ptr, size_t new_size)
{
printf("Secure realloc appele\n");
was_in_malloc_func = 1;
if(area_ptr == NULL)
        {
        //printf("ERROR: pointeur NULL passe en argument de srealloc\n");
        //exit(-1);
        return smalloc(fct_appelante,new_size);
        }

void *area;
area = area_ptr;
/*
if(area == NULL)
        {
        printf("ERROR: pointeur vers NULL passe en argument de sfree\n");
        exit(-1);
        }
*/
int pagesize;
pagesize = getpagesize();

// on cherche la plus petite adresse alignee sur une page:
void *real_area;
real_area = (void *)(((qword)area / pagesize) * pagesize);

void *data_area;
data_area = real_area - pagesize;

mprotect(data_area, pagesize, PROT_READ);

size_t *taille;
taille = data_area;

void *new_ptr;
new_ptr = smalloc("srealloc",new_size);

memcpy(new_ptr,area,*taille);

// desalloue les pages precedentes:

void *upper_area;
upper_area = real_area + *taille;

mmap_free(upper_area, pagesize);
mmap_free(real_area, *taille);
mmap_free(data_area, pagesize);

//*area_ptr = NULL;
was_in_malloc_func = 0;

nber_allocated_blocs --;

return new_ptr;
}


//--------------------------------------------------------------------------------

void sigsegv_handler(int n)
{
printf("ERREUR: Erreur de segmentation!\n");
if(was_in_malloc_func > 0) printf("cause: erreur interne de secure malloc %d\n",was_in_malloc_func);
exit(-1);
}


int secure_malloc_main(int nber_args, char *args[], int (*smain)(int, char **))
{
int ret_main_norm;
srand48(time(NULL));

signal(SIGSEGV,sigsegv_handler);

nber_allocated_blocs = 0;

ret_main_norm = smain(nber_args,args);

if(nber_allocated_blocs != 0)
        {
        printf("----------------------------------------------------\n");
        printf("   ERREUR: %d blocs n'ont pas ete liberes!     \n",nber_allocated_blocs);
        printf("----------------------------------------------------\n");
        }

return ret_main_norm;
}





/*------------------------------------------
 * Function to free a n dimensional
 *  table.
 *
 * INPUT:
        - the table to free
        - the number of dimensions
        - the sizes of each dimension
 *------------------------------------------
 */
void sFree_Tab_nD(char *caller,void *Tab,int n, int *sizes)
{
int i;

/* check if all arguments are OK: */

if(Tab == NULL)
        {
        printf("ERROR: in Free_Tab_nD: tried to free the NULL pointer !\n");
        exit(-1);
        }

if(n < 0)
        {
        printf("ERROR: in Free_Tab_nD: the number of dimensions should be > 0 !\n");
        exit(-1);
        }

for(i = 0; i < n; i++)
        {
        if(sizes[i] < 0)
                {
                printf("ERROR: in Free_Tab_nD: all dimensions should have a >0 size !\n");
                exit(-1);
                }
        }

/* free allocations recursively: */
if(n > 1) sFree_Tab_nD(caller,*((void **)Tab),n-1,sizes + 1);
sfree(caller, &Tab);
}




/*------------------------------------------
 * Function to allocate a n dimensional
 *  table.
 *
 * INPUT:
        - the number of dimensions
        - the size of each dimension
        - the size of the variable
 *------------------------------------------
 */
void *sAlloc_Tab_nD(char *caller,int n, int *sizes, size_t var_size)
{
int i,j;
void *tab = NULL;
void **tmp_tab;
long long int total_size = 1;
unsigned int max_size;
max_size = -1;

/* printf("%lg\n",(double)max_size); */

int pagesize;
pagesize = getpagesize();

/* check if all arguments are OK: */

if(n < 0)
        {
        printf("ERROR: in Alloc_Tab_nD: the number of dimensions should be > 0 !\n");
        exit(-1);
        }

for(i = 0; i < n; i++)
        {
        if(sizes[i] < 0)
                {
                printf("ERROR: in Alloc_Tab_nD: all dimensions should have a >0 size !\n");
                exit(-1);
                }
        }

for(i=0;i<n;i++)
        {
        total_size *= sizes[i];
        if(total_size > max_size)
                {
                printf("ERROR: in Alloc_Tab_nD: you allocate too much memory !\n");
                exit(-1);
                }
        }

/* allocation: */

for(i = n-1; i >= 0; i--)
        {
        if(i == n-1)
                {
                /* last table */
                tab = smalloc(caller,total_size * var_size);
                if(tab == NULL)
                        {
                        printf("ERROR: in Alloc_Tab_nD: malloc 1 failed \n");exit(-1); /* terminate */
                        }
                //printf("last table\n");
                }
        else
                {
                /* not the last table */
                tmp_tab = smalloc(caller,total_size * sizeof(void *));
                if(tmp_tab == NULL)
                        {
                        printf("ERROR: in Alloc_Tab_nD: malloc 2 failed \n");exit(-1); /* terminate */
                        }

                if(i== n-2)
                        {/* pointing into the last table */
                        for(j=0;j<total_size;j++) tmp_tab[j] = tab + j * var_size * sizes[i+1];
                        //printf("intermediate table: last table\n");
                        }
                else
                        { /* pointing into a pointer table */
                        for(j=0;j<total_size;j++) tmp_tab[j] = tab + j * sizeof(void *) * sizes[i+1];
                        //printf("intermediate table: pointer table\n");
                        }
                void *real_area;
                real_area = (void *)(((qword)tmp_tab / pagesize) * pagesize);
                mprotect(real_area, total_size * sizeof(void *), PROT_READ);
                tab = tmp_tab;
                }
        total_size /= sizes[i];
        }

return tab;
}




// /*ajout nico*/
// void *nicomalloc(size_t size)
// {
//     return smalloc(__func__, size);
// }
// 
// void nicofree(void *area_ptr)
// {
//     sfree(__func__, area_ptr);
// }








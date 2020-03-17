/*

Fichier contenant les definitions de smalloc

*/ 


#include <stdlib.h>
#include <time.h>

#define malloc(a) nicomalloc(a)
#define free(a) nicofree(a)

void *nicomalloc(size_t size);
void *nicofree(void *area_ptr);
void *smalloc(char *fct_appelante, size_t size);
void *srealloc(char *fct_appelante, void *area_ptr, size_t new_size);
void *scalloc(char *fct_appelante, int n, size_t size);
//void *calloc(size_t nelem, size_t elsize);

void sfree(char *fct_appelante, void *area_ptr);

int secure_malloc_main(int nber_args, char *args[], int (*smain)(int, char **));

void sfree_compat(void *ptr);
void *smalloc_compat(size_t size);

// pour la compatibilite avec les malloc et free de glibc:
void frjo(void *ptr);
void *joaloc(size_t size);


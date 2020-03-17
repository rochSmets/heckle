
#ifndef DEFINES
#define DEFINES


#define PI   3.141592653589
#define EPS2 1.0E-2
#define EPS4 1.0E-4
#define EPS6 1.0E-6
#define EPS8 1.0E-8
#define GAMMA 1.666666666666
#define RNM  ((double)rand()/(double)RAND_MAX)
#define PSMOOTH 8


//#define  idx(i, j, k, n0, n1, n2) ((k)+(n2)*((j)+(n1)*(i)))
#define  IDX(i, j, k, n0, n1, n2) ((k)+(n2)*((j)+(n1)*((i)+(n0)*(0))))

//#define ijk2k(ijk, n0, n1, n2) (ijk%n2)
//#define ijk2j(ijk, n0, n1, n2) (((ijk-ijk2k(ijk, n0, n1, n2))/n2)%n1)
//#define ijk2i(ijk, n0, n1, n2) ((((ijk-ijk2k(ijk, n0, n1, n2))/n2)-ijk2j(ijk, n0, n1, n2))/n1)

//#define shit(x) (printf("%s @ %d on %d\n", __FILE__, __LINE__, x), exit(0))
//#define CRASH(x) (printf("%s @ %d on %d\n", __FILE__, __LINE__, x), exit(0))
//#define SCRUNCH(x) (printf("%s @ %d on %d\n", __FILE__, __LINE__, x), exit(0))
#define IAMDEAD(x) (printf("%s @ %d on %d\n", __FILE__, __LINE__, x), exit(0))

#define NEIGHBOR_LEFT   4
#define NEIGHBOR_RIGHT  22

#define NEIGHBOR_TOP    16
#define NEIGHBOR_BOTTOM 10

#define NEIGHBOR_BACK   12
#define NEIGHBOR_FRONT  14




#endif


#include "hyb.h"


extern int nx, nz;
extern int tspt;

extern float dx, dz;
extern float me;
extern float ki;

extern float **bx, **by, **bz;
extern float **dn;
extern float **pxx, **pxy, **pxz;
extern float **pyx, **pyy, **pyz;
extern float **pzx, **pzy, **pzz;
extern float **qxx, **qxy, **qxz;
extern float **qyx, **qyy, **qyz;
extern float **qzx, **qzy, **qzz;
extern float **vix, **viy, **viz;
extern float **wix, **wiy, **wiz;

extern float temp_e;
extern float dt;

extern int it;
extern int ipc;
extern int gyrotrop;
extern int istart;



/* ____ 2 D version of the full pression tensor _____________________________ */
void ram()
{

float dbydx, dbzdx, dbxdz, dbydz;
float dlx, dlz;
float bxav, byav, bzav, btav;
float omega;
float FI, KP;
float bunx, buny, bunz;
float uv[3][3];
float cxp, cxc, cxm;
float czp, czc, czm;
float dvxdx, dvydx, dvzdx;
float dvxdz, dvydz, dvzdz;
float buff1, buff2, buff3, buff4;
float b_p[3][3];
float b_q[3][3];
float b_r[3][3];
float b_s[3][3];

float **vx, **vy, **vz;

int i, j, k, l, m, n;


/* _____ memory allocation _____ */
vx = (float **) malloc((nx+2) * sizeof(float *));
vy = (float **) malloc((nx+2) * sizeof(float *));
vz = (float **) malloc((nx+2) * sizeof(float *));

for (i = 0; i < nx+2; i++)
    {
    vx[i] = (float *) malloc((nz+2) * sizeof(float));
    vy[i] = (float *) malloc((nz+2) * sizeof(float));
    vz[i] = (float *) malloc((nz+2) * sizeof(float));
    }

/* _____ initial conditions _____ */
if (it == 0 && ipc == 0)
   {
   for (i = 0; i < nx+2; i++)
       {
       for (j = 0; j < nz+2; j++)
           {
           /* _____ at even # of half time step _____ */
           pxx[i][j] = dn[i][j]*temp_e;
           pyy[i][j] = dn[i][j]*temp_e;
           pzz[i][j] = dn[i][j]*temp_e;
           pxy[i][j] = 0.0;
           pxz[i][j] = 0.0;
           pyx[i][j] = 0.0;
           pyz[i][j] = 0.0;
           pzx[i][j] = 0.0;
           pzy[i][j] = 0.0;

           /* _____ at full time step _____ */
           qxx[i][j] = dn[i][j]*temp_e;
           qyy[i][j] = dn[i][j]*temp_e;
           qzz[i][j] = dn[i][j]*temp_e;
           qxy[i][j] = 0.0;
           qxz[i][j] = 0.0;
           qyx[i][j] = 0.0;
           qyz[i][j] = 0.0;
           qzx[i][j] = 0.0;
           qzy[i][j] = 0.0;
           }
       }
   }

/* _____ gyrotropic case _____ */
if (gyrotrop == 1)
   {
   for (i = 0; i < nx+2; i++)
       {
       for (j = 0; j < nz+2; j++)
           {
           pxx[i][j] = dn[i][j]*temp_e;
           pyy[i][j] = dn[i][j]*temp_e;
           pzz[i][j] = dn[i][j]*temp_e;
           pxy[i][j] = 0.0;
           pxz[i][j] = 0.0;
           pyx[i][j] = 0.0;
           pyz[i][j] = 0.0;
           pzx[i][j] = 0.0;
           pzy[i][j] = 0.0;
           }
       }
   }
 
/* _____ non-gyrotropic case _____ */
if (gyrotrop == 0)
   {
   /* _____ electron velocities : Ve = Vi - J / n e _____ */
   for (i = 1; i < nx+1; i++)
       {
       for (j = 1; j < nz+1; j++)
           {
           dlx       = 0.5/(dn[i][j]*dx);
           dlz       = 0.5/(dn[i][j]*dz);

           dbydx    = by[i][j]-by[i-1][j]+by[i][j-1]-by[i-1][j-1];
           dbzdx    = bz[i][j]-bz[i-1][j]+bz[i][j-1]-bz[i-1][j-1];
           dbxdz    = bx[i][j]+bx[i-1][j]-bx[i][j-1]-bx[i-1][j-1];
           dbydz    = by[i][j]+by[i-1][j]-by[i][j-1]-by[i-1][j-1];

           vx[i][j] = 0.5*(vix[i][j]+wix[i][j])-dlz*(-dbydz);
           vy[i][j] = 0.5*(viy[i][j]+wiy[i][j])-dlz*(+dbxdz)-dlx*(-dbzdx);
           vz[i][j] = 0.5*(viz[i][j]+wiz[i][j])-dlx*(+dbydx);
           }
       }

   /* _____ boundary conditions _____ */
   fill(vx, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
   fill(vy, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
   fill(vz, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);

   /* _____ dt P = -V.G P-P G.V-P.G V-(P.G V)T-q/M [P X B+(P X B)T] _____ */
   for (i = 1; i < nx+1; i++)
       {
       for (j = 1; j < nz+1; j++)
           {
           bxav  = 0.25*(bx[i][j]+bx[i][j-1]+bx[i-1][j]+bx[i-1][j-1]);
           byav  = 0.25*(by[i][j]+by[i][j-1]+by[i-1][j]+by[i-1][j-1]);
           bzav  = 0.25*(bz[i][j]+bz[i][j-1]+bz[i-1][j]+bz[i-1][j-1]);

           omega = sqrt(bxav*bxav+byav*byav+bzav*bzav)/me;
           FI    = (1.0-ki)*omega;
           KP    = ki*dt*omega;

           /* _____ create unit vectors _____ */
           btav = sqrt(bxav*bxav+byav*byav+bzav*bzav);

           bunx = bxav/btav;
           buny = byav/btav;
           bunz = bzav/btav;

           /* _____ uv[0][*] is along B _____ */
           uv[0][0] = bunx;
           uv[0][1] = buny;
           uv[0][2] = bunz;

           /* _____ uv[1][*] is in the X - Z plan _____ */
           uv[1][0] = +bunz/sqrt(bunx*bunx+bunz*bunz);
           uv[1][1] = 0.0;
           uv[1][2] = -bunx/sqrt(bunx*bunx+bunz*bunz);

           /* _____ uv[2][*] is as uv[2][*] = uv[0][*] X uv[1][*] _____ */
           uv[2][0] = -bunx*buny/sqrt(bunx*bunx+bunz*bunz);
           uv[2][1] = (bunx*bunx+bunz*bunz)/sqrt(bunx*bunx+bunz*bunz);
           uv[2][2] = -bunz*buny/sqrt(bunx*bunx+bunz*bunz);
 
           /* _____ coefficients for donor cells _____ */
           cxp = +0.5*(1.0-sign(vx[i][j]))/dlx;
           cxm = -0.5*(1.0+sign(vx[i][j]))/dlx;
           czp = +0.5*(1.0-sign(vz[i][j]))/dlz;
           czm = -0.5*(1.0+sign(vz[i][j]))/dlz;

           cxc = sign(vx[i][j])/dlx;
           czc = sign(vz[i][j])/dlz;

           dvxdx = (vx[i+1][j]-vx[i-1][j])/(2.0*dlx);
           dvydx = (vy[i+1][j]-vy[i-1][j])/(2.0*dlx);
           dvzdx = (vz[i+1][j]-vz[i-1][j])/(2.0*dlx);
           dvxdz = (vx[i][j+1]-vx[i][j-1])/(2.0*dlz);
           dvydz = (vy[i][j+1]-vy[i][j-1])/(2.0*dlz);
           dvzdz = (vz[i][j+1]-vz[i][j-1])/(2.0*dlz);

           /* _____ upwind scheme for the gyro-kinetic term _____ */

           /* _____ xx component _____ */
           buff1 = 2.0*(qxx[i][j]*dvxdx+qxz[i][j]*dvxdz);
           buff2 = qxx[i][j]*(dvxdx+dvzdz);
           buff3 = vx[i][j]*(cxp*qxx[i+1][j]+cxc*qxx[i][j]+cxm*qxx[i-1][j])+
                   vz[i][j]*(czp*qxx[i][j+1]+czc*qxx[i][j]+czm*qxx[i][j-1]);
           buff4 = FI*2.0*(pxy[i][j]*bunz-pxz[i][j]*buny);

           b_p[0][0] = pxx[i][j]-dt*(buff1+buff2+buff3+buff4);

           /* _____ xy component _____ */
           buff1 = qxx[i][j]*dvydx+qxz[i][j]*dvydz+
                   qyx[i][j]*dvxdx+qyz[i][j]*dvxdz;
           buff2 = qxy[i][j]*(dvxdx+dvzdz);
           buff3 = vx[i][j]*(cxp*qxy[i+1][j]+cxc*qxy[i][j]+cxm*qxy[i-1][j])+
                   vz[i][j]*(czp*qxy[i][j+1]+czc*qxy[i][j]+czm*qxy[i][j-1]);
           buff4 = FI*(pxz[i][j]*bunx-pxx[i][j]*bunz+
                       pyy[i][j]*bunz-pyz[i][j]*buny);

           b_p[0][1] = pxy[i][j]-dt*(buff1+buff2+buff3+buff4);

           /* _____ xz component _____ */
           buff1 = qxx[i][j]*dvzdx+qxz[i][j]*dvzdz+
                   qzx[i][j]*dvxdx+qzz[i][j]*dvxdz;
           buff2 = qxz[i][j]*(dvxdx+dvzdz);
           buff3 = vx[i][j]*(cxp*qxz[i+1][j]+cxc*qxz[i][j]+cxm*qxz[i-1][j])+
                   vz[i][j]*(czp*qxz[i][j+1]+czc*qxz[i][j]+czm*qxz[i][j-1]);
           buff4 = FI*(pxx[i][j]*buny-pxy[i][j]*bunx+
                       pzy[i][j]*bunz-pzz[i][j]*buny);

           b_p[0][2] = pxz[i][j]-dt*(buff1+buff2+buff3+buff4);

           /* _____ yy component _____ */
           buff1 = 2.0*(qyx[i][j]*dvydx+qyz[i][j]*dvydz);
           buff2 = qyy[i][j]*(dvxdx+dvzdz);
           buff3 = vx[i][j]*(cxp*qyy[i+1][j]+cxc*qyy[i][j]+cxm*qyy[i-1][j])+
                   vz[i][j]*(czp*qyy[i][j+1]+czc*qyy[i][j]+czm*qyy[i][j-1]);
           buff4 = FI*2.0*(pyz[i][j]*bunx-pyx[i][j]*bunz);

           b_p[1][1] = pyy[i][j]-dt*(buff1+buff2+buff3+buff4);

           /* _____ yz component _____ */
           buff1 = qyx[i][j]*dvzdx+qyz[i][j]*dvzdz+
                   qzx[i][j]*dvydx+qzz[i][j]*dvydz;
           buff2 = qyz[i][j]*(dvxdx+dvzdz);
           buff3 = vx[i][j]*(cxp*qyz[i+1][j]+cxc*qyz[i][j]+cxm*qyz[i-1][j])+
                   vz[i][j]*(czp*qyz[i][j+1]+czc*qyz[i][j]+czm*qyz[i][j-1]);
           buff4 = FI*(pyx[i][j]*buny-pyy[i][j]*bunx+
                       pzz[i][j]*bunx-pzx[i][j]*bunz);

           b_p[1][2] = pyz[i][j]-dt*(buff1+buff2+buff3+buff4);

           /* _____ zz component _____ */
           buff1 = 2.0*(qzx[i][j]*dvzdx+qzz[i][j]*dvzdz);
           buff2 = qzz[i][j]*(dvxdx+dvzdz);
           buff3 = vx[i][j]*(cxp*qzz[i+1][j]+cxc*qzz[i][j]+cxm*qzz[i-1][j])+
                   vz[i][j]*(czp*qzz[i][j+1]+czc*qzz[i][j]+czm*qzz[i][j-1]);
           buff4 = FI*2.0*(pzx[i][j]*buny-pzy[i][j]*bunx);

           b_p[2][2] = pzz[i][j]-dt*(buff1+buff2+buff3+buff4);

           /* _____ use the symetry of the pressure tensor _____ */
           b_p[1][0] = b_p[0][1];
           b_p[2][0] = b_p[0][2];
           b_p[2][1] = b_p[1][2];

           /* _____ transform the driver in uv coordinates _____ */
           for (k = 0; k < 3; k++)
               {
               for (l = 0; l < 3; l++)
                   {
                   b_q[k][l] = 0.0;

                   for (m = 0; m < 3; m++)
                       {
                       for (n = 0; n < 3; n++)
                           {
                           b_q[k][l] += uv[k][m]*b_p[m][n]*uv[l][n];
                           }
                       }
                   }
               }

           /* _____ solve for the transformed pressure _____ */
           b_r[0][0] = b_q[0][0];
           b_r[0][1] = (b_q[0][1]-KP*b_q[0][2])/(1.0+KP*KP);
           b_r[0][2] = (b_q[0][2]+KP*b_q[0][1])/(1.0+KP*KP);
           b_r[1][0] = b_r[0][1];
           b_r[1][1] = (b_q[1][1]*(1.0+2.0*KP*KP)-2.0*KP*b_q[1][2]+
                        2.0*KP*KP*b_q[2][2])/(1.0+4.0*KP*KP);
           b_r[1][2] = (KP*b_q[1][1]+b_q[1][2]-KP*b_q[2][2])/(1.0+4.0*KP*KP);
           b_r[2][0] = b_r[0][2];
           b_r[2][1] = b_r[1][2];
           b_r[2][2] = (2.0*KP*KP*b_q[1][1]+2.0*KP*b_q[1][2]+
                        (1.0+2.0*KP*KP)*b_q[2][2])/(1.0+4.0*KP*KP);

           /* _____ transform back in lab. coordinates _____ */
           for (k = 0; k < 3; k++)
               {
               for (l = 0; l < 3; l++)
                   {
                   b_s[k][l] = 0.0;

                   for (m = 0; m < 3; m++)
                       {
                       for (n = 0; n < 3; n++)
                           {
                           b_s[k][l] += uv[m][k]*b_r[m][n]*uv[n][l];
                           }
                       }
                   }
               }

           pxx[i][j] = b_s[0][0];
           pxy[i][j] = b_s[0][1];
           pxz[i][j] = b_s[0][2];
           pyx[i][j] = b_s[1][0];
           pyy[i][j] = b_s[1][1];
           pyz[i][j] = b_s[1][2];
           pzx[i][j] = b_s[2][0];
           pzy[i][j] = b_s[2][1];
           pzz[i][j] = b_s[2][2];

           }
       }

   /* _____ set the boundaries conditions _____ */
   fill(pxx, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
   fill(pyy, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
   fill(pzz, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);

   if (it%tspt == 0)
      {
      fill(pxy, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      fill(pxz, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      fill(pyx, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      fill(pyz, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      fill(pzx, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      fill(pzy, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      }

   if (it%tspt != 0)
      {
      fill(pxy, nx+2, nz+2, 0.0, 0.0, 0.0, 0.0);
      fill(pxz, nx+2, nz+2, 0.0, 0.0, 0.0, 0.0);
      fill(pyx, nx+2, nz+2, 0.0, 0.0, 0.0, 0.0);
      fill(pyz, nx+2, nz+2, 0.0, 0.0, 0.0, 0.0);
      fill(pzx, nx+2, nz+2, 0.0, 0.0, 0.0, 0.0);
      fill(pzy, nx+2, nz+2, 0.0, 0.0, 0.0, 0.0);
      }

   /* _____ smooth somewhat the pressure tensor _____ */
   if (it%tspt == 0)
      {
      smooth(pxx, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      smooth(pyy, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);
      smooth(pzz, nx+2, nz+2, +1.0, +1.0, +1.0, +1.0);

      smooth(pxy, nx+2, nz+2, +0.0, +0.0, +0.0, +0.0);
      smooth(pxz, nx+2, nz+2, +0.0, +0.0, +0.0, +0.0);
      smooth(pyx, nx+2, nz+2, +0.0, +0.0, +0.0, +0.0);
      smooth(pyz, nx+2, nz+2, +0.0, +0.0, +0.0, +0.0);
      smooth(pzx, nx+2, nz+2, +0.0, +0.0, +0.0, +0.0);
      smooth(pzy, nx+2, nz+2, +0.0, +0.0, +0.0, +0.0);
      }
   }

/* _____ free the pointers _____ */
free(vx);
free(vy);
free(vz);

}


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <omp.h>
#define NX 1000
#define NY 100000
////#include "mpi.h"

int main (int argc, char** argv)
{

float * A;
float ain;
int i,j;
int nsize,nx,ny;
double ti;
 //,tf;
double tt;
//int nt,ierr;
int index;

int setA(float*, int, int, int*, int);


nx = NX;
ny = NY;
nsize = NX*NY;

A = malloc(sizeof(float)*nsize);

#pragma omp target data map(A[0:nsize])
{

#pragma omp target teams distribute parallel for
for (i=1;i<nsize;i++)A[i]=1.0;

//printf("Total Devices: %d\n", omp_get_num_devices());

ti = omp_get_wtime();

////#pragma omp parallel for collapse(2)
#pragma omp target teams distribute parallel for collapse(2) private(index,i,j,ain)
 for (i=0;i<nx;i++)
   ////#pragma omp parallel for 
   for (j=0;j<ny;j++)
   {
      setA(&ain,i,j,&index,ny);
      A[index] = ain;
   }

// end target data map
}

tt = omp_get_wtime() - ti;

//for (i=0;i<nsize;i++)
//printf("%f\n",A[i]);
printf("Loop time: %f\n",tt);


}

DECLARE_TARGET
int setA(float *ain, int i, int j, int *iii, int ny)
{

      int index;
      double xx,yy;
      float zz;
      index = ny*i+j;
      xx = (double)(i-j);
      yy = xx*xx/(double)(i+j+1);
      zz = (float)sin(yy);
      *ain = zz;
      *iii = index;
      return 0;
}

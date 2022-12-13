#include <omp.h>
#include <stdio.h>
#include "../offload_acc_omp.h"
DECLARE_TARGET
int go_do_something(void);
DECLARE_TARGET_END

int main ( void ){

  int n;

/*   printf("Doing something!\n"); */

/* #pragma omp target device(0) */
/*   n=go_do_something(); */


/*   printf("Did something (%d)!\n",n); */
  return 0;
}

int go_do_something(void){

  /*  const int m = 1000;
  int i,j,k[m];

#pragma omp parallel for private(j)
  for (i=0;i<m;i++) {
    j=i*2;
    k[i] = j;
  }
  */
  return 0;//k[m-1];
}

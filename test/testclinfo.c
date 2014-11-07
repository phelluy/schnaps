#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  
  // unit tests
    
  int resu=TestCLInfo();
	 
  if (resu) printf("CLInfo test OK !\n");
  else printf("CLInfo test failed !\n");

  return !resu;
} 




int TestCLInfo(void){

  int test = (1==1);

  CLInfo cli;

  InitCLInfo(&cli,_CL_PLATFORM,_CL_DEVICE);
  PrintCLInfo(&cli);

  char prog[]="#pragma OPENCL EXTENSION cl_khr_fp64: enable \n"
    "__kernel  \n "
    "void testadd(__global double* a,__global double* b){ \n"
    " int i = get_global_id(0); \n"
    "  a[i] += b[i]; \n"
    "}\n";

  BuildKernels(&cli,prog);

  char* s;
  ReadFile("src/field.cl",&s);

  //printf("%s\n\n\n",s);

  BuildKernels(&cli,s);
  



  return test;

}

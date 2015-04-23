#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int TestCLInfo(void) {
  int test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  CLInfo cli;

  InitCLInfo(&cli, nplatform_cl, ndevice_cl);
  PrintCLInfo(&cli);

  char prog[]=
    "__kernel  \n "
    "void testadd(__global float* a,__global float* b){ \n"
    " int i = get_global_id(0); \n"
    "  a[i] += b[i]; \n"
    "}\n";

  BuildKernels(&cli, prog, NULL);


  GetOpenCLCode();

  char *s;
  ReadFile("schnaps.cl", &s);

  // Pass the value of m to the OpenCL code via the preprocessor
  char buildoptions[1000];
  sprintf(buildoptions, "-D _M=%d", 1);
  sprintf(numflux_cl_name, "%s", "NumFlux");
  strcat(buildoptions," -D NUMFLUX=");
  strcat(buildoptions, numflux_cl_name);
  
  BuildKernels(&cli, s, buildoptions);

  return test;
}

int main(void) {
  int resu = TestCLInfo();
  if (resu) 
    printf("CLInfo test OK !\n");
  else 
    printf("CLInfo test failed !\n");
  return !resu;
} 

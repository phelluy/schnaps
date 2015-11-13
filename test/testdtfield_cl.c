#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestDtfield_CL()
{
  int retval = 0;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  field f;
  init_empty_field(&f);
  
  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  char *mshname = "../test/disque2d.msh";
  //char *mshname = "test/unit-cube.msh";
  
  ReadMacroMesh(&(f.macromesh), mshname);
  Detect2DMacroMesh(&f.macromesh);
  BuildConnectivity(&f.macromesh);

#if 1
  // 2D version
  assert(f.macromesh.is2d);

  f.model.cfl = 0.05;
  f.model.m = 1;
  m = f.model.m;

  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

#else
  // 3D version
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 2; // z direction degree
  f.interp.interp_param[4] = 3; // x direction refinement
  f.interp.interp_param[5] = 3; // y direction refinement
  f.interp.interp_param[6] = 3; // z direction refinement
#endif

  set_global_m(f.model.m);
  set_source_CL(&f, "OneSource");
  Initfield(&f);
  CopyfieldtoGPU(&f);
  
  real tnow = 0.0;
  dtfield_CL(&f, tnow, f.wn_cl, 0, NULL, NULL);

  CopyfieldtoCPU(&f);

  // Displayfield(&f);
  show_cl_timing(&f);

  real *saveptr = f.dtwn;
  f.dtwn = calloc(f.wsize, sizeof(real));

  f.model.Source = OneSource;
  dtfield(&f, tnow, f.wn, f.dtwn);
 
  real maxerr = 0;
  for(int i = 0; i < f.wsize; i++) {
    real error = f.dtwn[i] - saveptr[i];
    //printf("error= \t%f\t%f\t%f\n", error, f.dtwn[i], saveptr[i]);
    maxerr = fmax(fabs(error), maxerr);
  }
  printf("max error: %f\n", maxerr);

  real tolerance;
  if(sizeof(real) == sizeof(double))
    tolerance = 1e-8;
  else
    tolerance = 1e-4;
  if(maxerr > tolerance) {
    retval += 1;
  }

  return retval;
}

int main()
{
  int retval = TestDtfield_CL();

  if(retval == 0) 
    printf("Dtfield_CL test OK !\n");
  else 
    printf("Dtfield_CL test failed !\n");

  return retval;
} 

#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestDtfield_CL(void);

int main(void) {
  int resu = TestDtfield_CL();

  if (resu) 
    printf("Dtfield_CL test OK !\n");
  else 
    printf("Dtfield_CL test failed !\n");

  return !resu;
} 

int TestDtfield_CL(void){
  bool test = true;

  field f;
  
  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  char *mshname =  "test/disque2d.msh";
  
  ReadMacroMesh(&(f.macromesh), mshname);
  Detect2DMacroMesh(&(f.macromesh));
  BuildConnectivity(&(f.macromesh));

#if 1
  // 2D version
  assert(f.macromesh.is2d);

  f.model.cfl = 0.05;
  f.model.m = 1;

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

  Initfield(&f);

  dtfield_CL(&f);
  CopyfieldtoCPU(&f);

  // Displayfield(&f);

  double *saveptr = f.dtwn;
  f.dtwn = calloc(f.wsize, sizeof(double));

  dtfield(&f, f.wn, f.dtwn);
 
  double maxerr = 0;
  for(int i = 0; i < f.wsize; i++) {
    double error = f.dtwn[i] - saveptr[i];
    //printf("error= \t%f\t%f\t%f\n", error, f.dtwn[i], saveptr[i]);
    maxerr = fmax(fabs(error), maxerr);
  }
  printf("max error: %f\n", maxerr);

  test = (maxerr < 1e-8);

  return test;
}

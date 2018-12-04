#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int TestDtfield3D_CL(void){  
  bool test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  
  MacroMesh mesh;
  //char *mshname = "../test/disque2d.msh";
  char *mshname = "../test/unit-cube.msh";
  
  //ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  ReadMacroMesh(&mesh, mshname);
  //Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

  // 3D version
  assert(!mesh.is2d);
  
  model.cfl = 0.05;
  model.m = 1;
  m = model.m;

  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TransBoundaryFlux;
  model.InitData = TransInitData;
  model.ImposedData = TransImposedData;

  int deg[]={2, 2, 2};
  int raf[]={3, 3, 3};
  
  Simulation simu;
  EmptySimulation(&simu);

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TransBoundaryFlux");
  strcat(cl_buildoptions, buf);


  set_global_m(model.m);

  
  model.Source = OneSource;
  set_source_CL(&simu, "OneSource");
  //model.Source = NULL;


  assert(simu.use_source_cl);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  cl_int status;
  cl_event clv_dtfield = clCreateUserEvent(simu.cli.context, &status);
  
  dtfield_CL(&simu, &simu.w_cl, 0, NULL, &clv_dtfield);
  clWaitForEvents(1, &clv_dtfield);
  CopyfieldtoCPU(&simu);

  // Displayfield(&f);
  show_cl_timing(&simu);

  schnaps_real *saveptr = simu.dtw;
  simu.dtw = calloc(simu.wsize, sizeof(schnaps_real));

  DtFields_old(&simu, simu.w, simu.dtw);
 
  schnaps_real maxerr = 0;
  for(int i = 0; i < simu.wsize; i++) {
    schnaps_real error = simu.dtw[i] - saveptr[i];
    //printf("error= \t%f\t%f\t%f\n", error, f.dtwn[i], saveptr[i]);
    maxerr = fmax(fabs(error), maxerr);
  }
  printf("max error: %f\n", maxerr);

  test = (maxerr < _SMALL * 10);

  return test;

}

int main(void) {
  int resu = TestDtfield3D_CL();

  if (resu) 
    printf("3D Dtfield_CL test OK !\n");
  else 
    printf("3D Dtfield_CL test failed !\n");

  return !resu; 
} 

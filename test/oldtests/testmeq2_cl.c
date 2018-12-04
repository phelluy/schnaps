#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "model.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#define _XOPEN_SOURCE 700

schnaps_real maxerr(schnaps_real *a, schnaps_real *b, int n) 
{
  schnaps_real err = 0.0;
  for(int i = 0; i < n; ++i) {
    err = fmax(fabs(a[i] - b[i]), err);
  }
  return err;
}

int TestmEq2(void)  { 
  bool test = true;

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  char *mshname =  "../test/testcube.msh";
  
  MacroMesh mesh;
  //ReadMacroMesh(&mesh,"../test/testcube.msh");
  ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

  model.m = 2;

  model.NumFlux = VecTransNumFlux2d;
  //f.model.NumFlux = Maxwell2DNumFlux_centered;
  model.BoundaryFlux = VecTransBoundaryFlux2d;
  model.InitData = VecTransInitData2d;
  model.ImposedData = VecTransImposedData2d;
  model.Source = NULL;


  int deg[]={3, 3, 0};
  int raf[]={2, 2, 1};

  assert(mesh.is2d);

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif
  
  CheckMacroMesh(&mesh, deg, raf);



  Simulation simu;
  EmptySimulation(&simu);

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  set_source_CL(&simu, "ZeroSource");
  sprintf(numflux_cl_name, "%s", "VecTransNumFlux2d");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "VecTransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);



  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = .5;
  simu.cfl=0.05;
  simu.vmax=1;

#if 0
  // C version
  RK2(&simu, tmax);
#else
  // OpenCL version
  schnaps_real dt = 0;
  RK4_CL(&simu, tmax, dt, 0, 0, 0);

  CopyfieldtoCPU(&simu); 
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&simu);
  printf("\n");
#endif


  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%f\n", dd);

  schnaps_real tolerance = 0.009;

  test = dd < tolerance;
  
 
  return test;

} ;

int main(void) {
  int resu = TestmEq2();
  if (resu) 
    printf("OpenCL m greater than 1 test OK!\n");
  else 
    printf("OpenCL m greater than 1 test failed !\n");
  return !resu;
}

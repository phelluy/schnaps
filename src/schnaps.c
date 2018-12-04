#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "getopt.h"
#ifdef _WITH_OPENCL
#include "clutils.h"
#include "clinfo.h"
#endif
int main(int argc, char *argv[]) 
{
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

  model.m = 7;

  model.NumFlux = Maxwell2DNumFlux_upwind;
  //f.model.NumFlux = Maxwell2DNumFlux_centered;
  model.BoundaryFlux = Maxwell2DBoundaryFlux_upwind;
  model.InitData = Maxwell2DInitData;
  model.ImposedData = Maxwell2DImposedData;
  model.Source = Maxwell2DSource;
  //model.Source = NULL;


  int deg[]={3, 3, 0};
  int raf[]={4, 4, 1};

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

#ifdef _WITH_OPENCL

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  set_source_CL(&simu, "Maxwell2DSource");
  sprintf(numflux_cl_name, "%s", "Maxwell2DNumFlux_upwind");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell2DBoundaryFlux_upwind");
  strcat(cl_buildoptions, buf);
#endif


  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = .5;
  simu.cfl=0.2;
  simu.vmax=1;

#if 1
  // C version
  RK4(&simu, tmax);
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

  schnaps_real tolerance = 0.0025;

  test = dd < tolerance;

  FreeMacroMesh(&mesh);

  int exit_status = !test; // because 0 = success...
  
  return exit_status;
}

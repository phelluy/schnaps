#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include <string.h>

int TestMaxwell3D()
{
  int retval = 0;

  MacroMesh mesh;
  //ReadMacroMesh(&mesh, "../test/testdisque.msh");
  ReadMacroMesh(&mesh, "../test/testcube.msh");
  BuildConnectivity(&mesh);

  schnaps_real dt = 0.0;
  schnaps_real cfl = 0.1;
  schnaps_real vmax = 1;
  schnaps_real tmax = 1;

  int thedeg = 3;
  int theraf = 2;
  int deg[] = {thedeg, thedeg, thedeg};
  int raf[] = {theraf, theraf, theraf};
  CheckMacroMesh(&mesh, deg, raf);
  
  Simulation simu;
  EmptySimulation(&simu);

  Model model;

  model.m = 8;

  model.InitData = Maxwell3DInitData;
  model.ImposedData = Maxwell3DImposedData;
  
  model.NumFlux = Maxwell3DNumFluxClean_upwind;
  model.BoundaryFlux = Maxwell3DBoundaryFluxClean_upwind;
  model.Source = NULL;

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  //set_source_CL(&simu, "Maxwell3DSource");
  sprintf(numflux_cl_name, "%s", "Maxwell3DNumFluxClean_upwind");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);


  //sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell3DBoundaryFluxClean_upwind");
  
  sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell3DBoundaryFluxClean_upwind");
  strcat(cl_buildoptions, buf);

  
  InitSimulation(&simu, &mesh, deg, raf, &model);

  simu.cfl = cfl;
  simu.vmax = vmax;

#if 0
  // C version
  RK4(&simu, tmax);
#else
  // OpenCL version
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return 0;
  }
  

  RK4_CL(&simu, tmax, dt, 0, 0, 0);

  CopyfieldtoCPU(&simu);
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&simu);
  printf("\n");
#endif

  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = L2error(&simu);
  printf("erreur L2=%f\n", dd);

  schnaps_real tolerance = 0.08;
  if(dd >= tolerance)
    retval += 1;

  return retval;
}

int main()
{
  int retval = TestMaxwell3D();
  if(retval == 0)
    printf("Maxwell3D test OK!\n");
  else
    printf("Maxwell3D failed !\n");
  return retval;
}

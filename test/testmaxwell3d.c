#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include <string.h>

int TestMaxwell3D() 
{
  int retval = 0;
  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;  
  f.model.m = 8; // num of conservative variables

  f.model.NumFlux = Maxwell3DCleanNumFlux_upwind;
  f.model.BoundaryFlux = Maxwell3DCleanBoundaryFlux_upwind;
  f.model.InitData = Maxwell3DCleanInitData;
  f.model.ImposedData = Maxwell3DCleanImposedData;
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 2; // z direction refinement

  ReadMacroMesh(&f.macromesh, "../test/testcube.msh");

  BuildConnectivity(&f.macromesh);

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  // Source is for rho and J, which are zero here.
  //set_source_CL(&f, "Maxwell3DSource");
  sprintf(numflux_cl_name, "%s", "Maxwell3DCleanNumFlux_upwind");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell3DCleanBoundaryFlux_upwind");
  strcat(cl_buildoptions, buf);

  Initfield(&f);
  f.tnow = 0.0;
  
  CheckMacroMesh(&f.macromesh, f.interp.interp_param + 1);

  real tmax = 0.1;
  f.vmax = 1.0;
  real dt = set_dt(&f);
  //  dt = 1e-4;

#if 0
  // C version
  RK4(&f, tmax, dt);
#else
  // OpenCL version
  CopyfieldtoGPU(&f);
  RK4_CL(&f, tmax, dt, 0, 0, 0);
  CopyfieldtoCPU(&f);
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&f);
  printf("\n");
#endif

  // Save the results and the error
  /* Plotfield(0, false, &f, NULL, "dgvisu.msh"); */
  /* Plotfield(0, true, &f, "error", "dgerror.msh"); */

  real dd = L2error(&f);
  real tolerance = 0.08;
  if(dd > tolerance)
    retval += 1;
  
  printf("L2 error: %f\n", dd);

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

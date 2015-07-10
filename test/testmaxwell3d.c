#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include <string.h>

int TestMaxwell3D() 
{
  bool test = true;
  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;  
  f.model.m = 8; // num of conservative variables

  f.model.NumFlux = Maxwell3DNumFluxClean_uncentered;
  f.model.BoundaryFlux = Maxwell3DBoundaryFlux_uncentered;
  f.model.InitData = Maxwell3DInitData;
  f.model.ImposedData = Maxwell3DImposedData;
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 2; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");

  // FIXME: temp
  /* Detect2DMacroMesh(&(f.macromesh)); */
  /* assert(f.macromesh.is2d); */

  BuildConnectivity(&(f.macromesh));

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  // Source is for rho and J, which are zero here.
  //set_source_CL(&f, "Maxwell3DSource");
  sprintf(numflux_cl_name, "%s", "Maxwell3DNumFluxClean_uncentered");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "Maxwell3DBoundaryFlux_uncentered");
  strcat(cl_buildoptions, buf);

  Initfield(&f);
  
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  real tmax = 0.1;
  f.vmax = 1.0;
  real dt = set_dt(&f);
  //  dt = 1e-4;

#if 0
  // C version
  RK2(&f, tmax, dt);
#else
  // OpenCL version
  RK4_CL(&f, tmax, dt, 0, 0, 0);
  CopyfieldtoCPU(&f);
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&f);
  printf("\n");
#endif

  // Save the results and the error
  Plotfield(0, false, &f, NULL, "dgvisu.msh");
  Plotfield(0, true, &f, "error", "dgerror.msh");

  real dd = L2error(&f);
  real tolerance = 0.08;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
}

int main() 
{
  int resu = TestMaxwell3D();
  if (resu) 
    printf("Maxwell3D test OK!\n");
  else 
    printf("Maxwell3D failed !\n");
  return !resu;
}

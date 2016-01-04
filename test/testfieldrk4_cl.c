#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldRK4_CL()
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
  // test/testmacromesh.msh
  // test/unit-cube.msh

  char *mshname =  "../test/disque2d.msh";
  //char *mshname =  "../test/testdisque2d.msh";
  //char *mshname =  "../test/unit-cube.msh";
  
  ReadMacroMesh(&f.macromesh, mshname);
  Detect2DMacroMesh(&f.macromesh);
  BuildConnectivity(&f.macromesh);

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

  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);
 
  real tmax = 0.01;
  f.vmax=1;
  real dt = 0.0;

#if 1
  CopyfieldtoGPU(&f);
  clFinish(f.cli.commandqueue);
  RK4_CL(&f, tmax, dt, 0, NULL, NULL);
  clFinish(f.cli.commandqueue);
  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);
  show_cl_timing(&f);
#else
  RK4(&f, tmax, dt);
#endif
  
  /* Plotfield(0, false, &f, NULL, "dgvisu.msh"); */
  /* Plotfield(0, true , &f, "error", "dgerror.msh"); */

  real dd = L2error(&f);

  printf("L2 error: %f\n", dd);

  real tolerance = 0.002;
  if(dd > tolerance)
    retval += 1;
  
  return retval;
}

int main() {
  int retval = TestfieldRK4_CL();

  if(retval == 0) 
    printf("field RK4_CL test OK !\n");
  else 
    printf("field RK4_CL test failed !\n");

  return retval;
} 

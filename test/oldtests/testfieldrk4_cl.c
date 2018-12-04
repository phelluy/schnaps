#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int TestfieldRK4_CL(void){
  int test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  Simulation simu;
  EmptySimulation(&simu);


  MacroMesh mesh;
  char *mshname =  "../test/disque2d.msh";
  
  ReadMacroMesh(&mesh, mshname);
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);
  int deg[]={3, 3, 0};
  int raf[]={2, 2, 2};
  CheckMacroMesh(&mesh, deg, raf);

  /* f.model.cfl = 0.05; */
  /* f.model.m = 1; */
  /* f.model.NumFlux = TransNumFlux; */
  /* f.model.BoundaryFlux = TestTransBoundaryFlux; */
  /* f.model.InitData = TestTransInitData; */
  /* f.model.ImposedData = TestTransImposedData; */
  /* f.varindex = GenericVarindex; */

  /* f.interp.interp_param[0] = f.model.m; */
  /* f.interp.interp_param[1] = 3; // x direction degree */
  /* f.interp.interp_param[2] = 3; // y direction degree */
  /* f.interp.interp_param[3] = 3; // z direction degree */
  /* f.interp.interp_param[4] = 1; // x direction refinement */
  /* f.interp.interp_param[5] = 1; // y direction refinement */
  /* f.interp.interp_param[6] = 1; // z direction refinement */

  Model model;

#if 1
  // 2D version
  assert(mesh.is2d);

  simu.cfl = 0.04;
  model.m = 1;

  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  model.Source = NULL;

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux2d");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);

#else
  // 3D version
  simu.cfl = 0.05;
  model.m = 1;
  
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;
#endif

  //AffineMapMacroMesh(&(simu.macromesh));
  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = .1;
  simu.vmax = 1;
  schnaps_real dt = 0;
  RK4_CL(&simu, tmax, dt,  0, NULL, NULL);
  
  CopyfieldtoCPU(&simu);
 
  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  //Plotfield(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = L2error(&simu);

  printf("L2 error: %f\n", dd);

  show_cl_timing(&simu);

  schnaps_real tolerance = 0.007;

  test = dd < tolerance;
  
  return test;
};

int main(void) {
  int resu = TestfieldRK4_CL();

  if(resu) 
    printf("field RK4_CL test OK !\n");
  else 
    printf("field RK4_CL test failed !\n");

  return !resu;
} 

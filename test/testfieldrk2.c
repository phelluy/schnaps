#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldRK2(void){
  bool test = true;

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  char *mshname =  "../test/disque2d.msh";

  MacroMesh mesh;
  ReadMacroMesh(&mesh,mshname);
  //ReadMacroMesh(&mesh,"../test/testdisque.msh");
  //ReadMacroMesh(&mesh,"../test/testcube2.msh");
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

#if 1
  // 2D version
  model.cfl = 0.05;
  model.m = 1;

  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  model.Source = NULL;

  int deg[]={3, 3, 0};
  int raf[]={3, 3, 1};

  assert(mesh.is2d);
#else
  // 3D version
  model.m = 1;
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;
  model.Source = NULL;

  int deg[]={2, 1, 1};
  //int raf[]={2, 1, 1};
  int raf[]={4, 4, 4};

#endif

  // 2015-01-19: the below parameters fail with testmacrocellinterface
  // but pass the test here (perhaps because the error is hidden by
  // the RK error?)
  /*
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement
  */

  //AffineMapMacroMesh(&(f.macromesh));

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  schnaps_real tmax = 0.25;
  //tmax = 0.006848;
  simu.cfl=0.2;
  simu.vmax=1;
  RK2(&simu,tmax);

  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%f\n", dd);

  schnaps_real tolerance = 0.003;

  test = dd < tolerance;

  FreeMacroMesh(&mesh);

  return test;
}

int main(void) {
  int resu = TestfieldRK2();
  if(resu)
    printf("field RK2 test OK !\n");
  else
    printf("field RK2 test failed !\n");
  return !resu;
}

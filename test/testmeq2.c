#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestmEq2(void){
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
  ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  //ReadMacroMesh(&mesh,"../test/testcube2.msh");
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

  model.cfl = 0.05;
  model.m = 2;

  model.NumFlux = VecTransNumFlux2d;
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

  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = 0.5;
  simu.cfl=0.05;
  simu.vmax=1;
  RK2(&simu,tmax);
 
  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%f\n", dd);

  schnaps_real tolerance = 0.015;

  test = dd < tolerance;

  FreeMacroMesh(&mesh);
  
  return test;
}

int main(void) {
  int resu = TestmEq2();
  if(resu) 
    printf("m greater than 1 test OK !\n");
  else 
    printf("m greater than 1 test failed !\n");
  return !resu;
} 

#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldRK4(void){
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

  int deg[]={4, 4, 0};
  int raf[]={4, 4, 1};

  assert(mesh.is2d);
#else
  // 3D version
  model.m = 1;
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;

  int deg[]={4, 4, 2};
  int raf[]={8, 8, 3};

#endif


#ifdef _WITH_OPENCL
  // if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
  //  printf("OpenCL device not acceptable.\n");
  //  return true;
  // }
#endif
  
  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = 0.02;
  simu.cfl=0.1;
  simu.vmax=1;
  RK4(&simu,tmax);
 
  //PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  //PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%.8e\n", dd);

  schnaps_real tolerance = 0.007;

  test = dd < tolerance;

  FreeMacroMesh(&mesh);
  
  return test;
}

int main(void) {
  int resu = TestfieldRK4();
  if(resu) 
    printf("field RK4 test OK !\n");
  else 
    printf("field RK4 test failed !\n");
  return !resu;
} 

#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldRK2_2D(void) {
  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testdisque2d.msh");
  Detect2DMacroMesh(&mesh);
  // require a 2d computation
  assert(mesh.is2d);
  BuildConnectivity(&mesh);

  Model model;

  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  model.Source = NULL;

  int deg[]={2, 2, 0};
  int raf[]={1, 1, 1};



  CheckMacroMesh(&mesh,deg,raf);

  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  printf("cfl param =%f\n",simu.hmin);
 
  schnaps_real tmax = 0.2;
  simu.cfl=0.2;
  simu.vmax=1;
  RK2(&simu,tmax);
 
  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%f\n", dd);

  schnaps_real tolerance = 0.07;

  test = dd < tolerance;

  return test;
}

int main(void) {
  int resu = TestfieldRK2_2D();
  if (resu) 
    printf("field RK2 2D test OK !\n");
  else 
    printf("field RK2 2D test failed !\n");
  return !resu;
} 

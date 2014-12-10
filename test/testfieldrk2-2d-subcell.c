#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  // unit tests
  int resu=TestFieldRK2_2D_SubCell();
  if (resu) printf("Field RK2 2D Subcell test OK !\n");
  else printf("Field RK2 2D Subcell test failed !\n");
  return !resu;
} 

int TestFieldRK2_2D_SubCell(void) {
  bool test = true;

  Field f;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testmacromesh.msh");
  bool is2d = Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  InitField(&f);
  f.is2d = true;

  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param: %f\n", f.hmin);

  assert(f.is2d);

  RK2(&f, 0.2);
 
  PlotField(0, false, &f, NULL, "dgvisu.msh");
  PlotField(0, true, &f, "error", "dgerror.msh");

  double dd = L2error(&f);

  printf("erreur L2=%f\n", dd);

  test = test && (dd < 0.006);
  return test;
};

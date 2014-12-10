#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"

int main(void) {
  // Unit tests
  int resu = TestmEq2();
 
  if (resu) printf("multiple velocity transport test OK !\n");
  else printf("multiple velocity transport test failed !\n");
  return !resu;
}

int TestmEq2(void) {
  bool test = true;
  Field f;
  
  f.model.m = 9; // num of conservative variables
  f.model.mx = 3;
  f.model.my = 3;
  f.model.mz = 1;
  assert(f.model.m == f.model.mx * f.model.my * f.model.mz);
  f.model.vmax = 1.0;
  f.model.NumFlux = vlaTransNumFlux2d;
  f.model.BoundaryFlux = vlaTransBoundaryFlux2d;
  f.model.InitData = vlaTransInitData2d;
  f.model.ImposedData = vlaTransImposedData2d;
  f.varindex = GenericVarindex;

  // Set the global parameters for the Vlasov equation
  set_vlasov_params(f.model.mx, f.model.my, f.model.mz, f.model.vmax); 
    
  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "geo/square.msh");
  // Try to detect a 2d mesh
  bool is2d = Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);

  // Mesh preparation
  BuildConnectivity(&(f.macromesh));
 
  // Prepare the initial fields
  InitField(&f);
  f.is2d = true;

  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param: %f\n", f.hmin);

  dtField(&f);

  double tmax = 0.1;
  RK2(&f, tmax);
 
  // Save the results and the error

  int mplot = f.model.m / 2; // m / 2 gives (vx, vy) = (0, 0) 
  printf("mplot: %d\n", mplot);
  PlotField(mplot, false, &f, "dgvisu.msh");
  /* PlotField(mplot, true, &f, "dgerror.msh"); */

  /* double dd = L2error(&f); */

  /* printf("erreur L2=%f\n", dd); */
  /* test = test && (dd < 1e-7); */
  return test;
};

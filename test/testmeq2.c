#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"

int TestmEq2(void) {
  bool test = true;
  field f;
  init_empty_field(&f);
 
  int vec = 2;

  f.model.cfl = 0.05;  
  if(vec == 2) {
    f.model.m = 2; // num of conservative variables
  } else {
    f.model.m = 1; // num of conservative variables
  }
  f.model.NumFlux = VecTransNumFlux2d;
  f.model.BoundaryFlux = VecTransBoundaryFlux2d;
  f.model.InitData = VecTransInitData2d;
  f.model.ImposedData = VecTransImposedData2d;
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");
  // Try to detect a 2d mesh
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  // Mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // Prepare the initial fields
  
  Initfield(&f);
  //f.dt = 1e-3;
  
  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  // time derivative
  //dtfield(&f);
  //Displayfield(&f);
 
  real tmax = 0.1;
  f.vmax=1;
  real dt = set_dt(&f);
  RK2(&f, tmax, dt);
 
  // Save the results and the error
  Plotfield(0, false, &f, NULL, "dgvisu.msh");
  Plotfield(0, true, &f, "error", "dgerror.msh");

  real dd = L2error(&f);
  real tolerance = 1e-4;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
};

int main(void) {
  int resu = TestmEq2();
  if (resu) 
    printf("m greater than 1 test OK!\n");
  else 
    printf("m greater than 1 test failed !\n");
  return !resu;
}

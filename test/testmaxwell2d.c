#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"

int TestMaxwell2D(void) {
  bool test = true;
  field f;

  f.model.cfl = 0.05;  
  f.model.m = 4; // num of conservative variables

  f.model.NumFlux = Maxwell2DNumFlux;
  f.model.BoundaryFlux = Maxwell2DBoundaryFlux;
  f.model.InitData = Maxwell2DInitData;
  f.model.ImposedData = Maxwell2DImposedData;
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
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
  f.is2d = true;
  //f.dt = 1e-3;
  
  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  // time derivative
  //dtfield(&f);
  //Displayfield(&f);
 
  double tmax = 1.0;
  RK2(&f, tmax);
 
  // Save the results and the error
  Plotfield(0, false, &f, NULL, "dgvisu.msh");
  Plotfield(0, true, &f, "error", "dgerror.msh");

  double dd = L2error(&f);
  double tolerance = 8e-3;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
};

int main(void) {
  int resu = TestMaxwell2D();
  if (resu) 
    printf("Maxwell2D test OK!\n");
  else 
    printf("Maxwell2D failed !\n");
  return !resu;
}

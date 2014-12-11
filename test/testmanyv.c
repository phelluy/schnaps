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
  
  f.model.mx = 5;
  f.model.my = 5;
  f.model.mz = 1;
  f.model.m = f.model.mx * f.model.my * f.model.mz;
  f.model.vmax = 0.5;
  f.model.NumFlux = vlaTransNumFlux2d;
  f.model.BoundaryFlux = vlaTransBoundaryFlux2d;
  f.model.InitData = vlaTransInitData2d;
  f.model.ImposedData = vlaTransImposedData2d;
  f.varindex = GenericVarindex;

  // Set the global parameters for the Vlasov equation
  set_vlasov_params(&(f.model));
  
  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
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

  { // Check error wrt initial conditions
    double dd = L2error(&f);
    printf("L2 error: %f\n", dd);
  }

  double tmax = 1.0;
  RK2(&f, tmax);
 
  // Save the results and the error
  for(int ix = 0; ix < f.model.mx; ++ix) {
    for(int iy = 0; iy < f.model.my; ++iy) {
      int mplot = ix * f.model.my + iy; 
      printf("mplot: %d\n", mplot);
      double vx = vlasov_vel(ix, f.model.mx, f.model.vmax);
      double vy = vlasov_vel(iy, f.model.my, f.model.vmax);
      char fieldname[100];
      sprintf(fieldname, "output field has v = (%f,%f)", vx, vy);
      printf("%s\n", fieldname);
      
      char filename[100];
      sprintf(filename, "dgvisuix%diy%d.msh", ix, iy);
      PlotField(mplot, false, &f, fieldname, filename);
    }
  }
  /* PlotField(mplot, true, &f, "dgerror.msh"); */

  double dd = L2error(&f);
  printf("L2 error: %f\n", dd);
  /* test = test && (dd < 1e-7); */

  return test;
};

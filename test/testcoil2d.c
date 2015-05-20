#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "maxwell.h"





void Coil2DImposedData(real x[3], real t,real w[]) {
  w[0]=0;
  w[1]=0;
  w[2]=1;
  w[3]=0;
  w[4]=0;
  w[5]=0;
  w[6]=0;

  real r = x[0] * x[0] + x[1] * x[1];

  if (r > 1) w[2]=0;

}

void Coil2DBoundaryFlux(real x[3], real t, real wL[], real *vnorm,
		       real *flux) {
  real wR[7];
  Coil2DImposedData(x, t, wR);
  Maxwell2DNumFlux(wL, wR, vnorm, flux);
}


void Coil2DInitData(real x[3], real w[]) {
  real t = 0;
  Coil2DImposedData(x, t, w);
}


int TestCoil2D(void) {
  bool test = true;
  field f;

  f.model.cfl = 0.05;  
  f.model.m = 7; // num of conservative variables

  f.model.NumFlux = Maxwell2DNumFlux;
  f.model.BoundaryFlux = Coil2DBoundaryFlux;
  f.model.InitData = Coil2DInitData;
  f.model.ImposedData = Coil2DImposedData;
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 8; // x direction refinement
  f.interp.interp_param[5] = 8; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testmacromesh.msh");
  // Try to detect a 2d mesh
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  // Mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // Prepare the initial fields
  Initfield(&f);
  f.model.Source = Maxwell2DSource;
  f.is2d = true;
  //f.dt = 1e-3;
  
  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  // time derivative
  //dtfield(&f);
  //Displayfield(&f);
 
  // init the particles on a circle
  PIC pic;
  InitPIC(&pic,1000);
  CreateCoil2DParticles(&pic,&(f.macromesh));
  PlotParticles(&pic,&(f.macromesh));

  // time evolution
  real tmax = 0.3;
  f.vmax=1;
  RK2(&f, tmax);
 
  // Save the results and the error
  Plotfield(2, false, &f, NULL, "dgvisu.msh");
  Plotfield(2, true, &f, "error", "dgerror.msh");

  real dd = L2error(&f);
  real tolerance = 8e-3;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
};

int main(void) {
  int resu = TestCoil2D();
  if (resu) 
    printf("Coil2D test OK!\n");
  else 
    printf("Coil2D failed !\n");
  return !resu;
}

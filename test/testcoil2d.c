#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "maxwell.h"

void Coil2DImposedData(const real x[3],const real t,real w[])
{
  real r = x[0] * x[0] + x[1] * x[1];
  w[0] = 0;
  w[1] = 0;
  w[2] = r > 1 ? 0 : 1;
  w[3] = 0;
  w[4] = 0;
  w[5] = 0;
  w[6] = 0;
}

void coil_pre_dtfield(void *f, real *w);

void coil_pre_dtfield(void *fv, real *w){
  AccumulateParticles(fv,w);
}


void Coil2DSource(const real *x, const real t, const real *w, real *source)
{
  // w: (Ex, Ey, Hz, Hz, \lambda, rho, Jx, Jy)
  
  // FIXME add documentation
  
  const real khi = 1.0;
  source[0] = -w[4];
  source[1] = -w[5];
  source[2] = 0;
  source[3] = 0;//instead of khi * w[6]: we want div E =0
  source[4] = 0;
  source[5] = 0;
  source[6] = 0;
}




void Coil2DBoundaryFlux(real x[3], real t, real wL[], real *vnorm,
			real *flux)
{
  real wR[7];
  Coil2DImposedData(x, t, wR);
  Maxwell2DNumFlux_uncentered(wL, wR, vnorm, flux);
}

void Coil2DInitData(real x[3], real w[])
{
  real t = 0;
  Coil2DImposedData(x, t, w);
}

int TestCoil2D(void)
{
  bool test = true;
  field f;
  init_empty_field(&f);

  init_empty_field(&f);

  f.model.cfl = 0.2;  
  f.model.m = 7; // num of conservative variables

  f.model.NumFlux = Maxwell2DNumFlux_uncentered;
  f.model.BoundaryFlux = Coil2DBoundaryFlux;
  f.model.InitData = Coil2DInitData;
  f.model.ImposedData = Coil2DImposedData;
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
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
  f.model.Source = Coil2DSource;
  f.pre_dtfield = coil_pre_dtfield;
  //f.dt = 1e-3;
  
  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  // time derivative
  //dtfield(&f);
  //Displayfield(&f);
 
  // init the particles on a circle
  PIC pic;
  InitPIC(&pic, 100);
  CreateCoil2DParticles(&pic, &f.macromesh);
  PlotParticles(&pic, &f.macromesh);

  f.pic = &pic;

  // time evolution
  real tmax = 0.1;
  f.vmax = 1;
  real dt = set_dt(&f);
  RK2(&f, tmax, dt);
 
  // Save the results and the error
  Plotfield(2, false, &f, NULL, "dgvisu.msh");
  Plotfield(2, true, &f, "error", "dgerror.msh");

  real dd = L2error(&f);
  real tolerance = 0.3;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
}

int main(void) {
  int resu = TestCoil2D();
  if (resu) 
    printf("Coil2D test OK!\n");
  else 
    printf("Coil2D failed !\n");
  return !resu;
}

#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"

void Test_TransportVP_ImposedData(const real *x, const real t, real *w);
void Test_TransportVP_InitData(real *x, real *w);
real TransportVP_ImposedKinetic_Data(const real *x, const real t, real v);
void Test_TransportVP_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux);

void UpdateVlasovPoisson(void *field, real *w);
void PlotVlasovPoisson(void *vf, real *w);

int main(void) {
  
  // unit tests
    
  int resu = Test_TransportVP();
	 
  if (resu) printf("poisson test OK !\n");
  else printf("poisson test failed !\n");

  return !resu;
} 

int Test_TransportVP(void) {

  bool test = true;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=_INDEX_MAX; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  f.model.NumFlux=VlasovP_Lagrangian_NumFlux;
 
  //f.model.Source = NULL;
 
  f.model.InitData = Test_TransportVP_InitData;
  f.model.ImposedData = Test_TransportVP_ImposedData;
  f.model.BoundaryFlux = Test_TransportVP_BoundaryFlux;

  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;  // _M
  f.interp.interp_param[1] = 3;  // x direction degree
  f.interp.interp_param[2] = 0;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 32;  // x direction refinement
  f.interp.interp_param[5] = 1;  // y direction refinement
  f.interp.interp_param[6] = 1;  // z direction refinement
 // read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");
  // try to detect a 2d mesh
  Detect1DMacroMesh(&(f.macromesh));
  bool is1d = f.macromesh.is1d;
  assert(is1d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  f.model.cfl = 0.05;
  Initfield(&f);
  f.vmax = _VMAX; // maximal wave speed
  f.nb_diags = 3;
  f.pre_dtfield = UpdateVlasovPoisson;
  f.update_after_rk = PlotVlasovPoisson;
  f.model.Source = VlasovP_Lagrangian_Source;
  // prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  real tmax = 0.03;
  real dt = set_dt(&f);
  RK2(&f, tmax, dt);
  //RK2(&f,0.03,0.05);

   // save the results and the error
  int iel = 2 * _NB_ELEM_V / 3;
  int iloc = _DEG_V;
  printf("Trace vi=%f\n", -_VMAX + iel * _DV + _DV * glop(_DEG_V, iloc));
  Plotfield(iloc + iel * _DEG_V, false, &f, "sol","dgvisu.msh");
  Plotfield(iloc + iel * _DEG_V, true, &f, "error","dgerror.msh");
  Plot_Energies(&f, dt);

  real dd_Kinetic = L2_Kinetic_error(&f);
  
  printf("erreur kinetic L2=%lf\n", dd_Kinetic);
  test= test && (dd_Kinetic < 1e-2);

  return test;
}

void Test_TransportVP_ImposedData(const real *x, const real t, real *w) {

  for(int i = 0; i <_INDEX_MAX_KIN + 1; i++) {
    int j = i % _DEG_V; // local connectivity put in function
    int nel = i / _DEG_V; // element num (TODO : function)

    real vi = (-_VMAX + nel * _DV + _DV * glop(_DEG_V, j));
 
    w[i] = TransportVP_ImposedKinetic_Data(x, t, vi);
  }
  // exact value of the potential and electric field
  w[_INDEX_PHI] = x[0];
  w[_INDEX_EX] = 1;
  w[_INDEX_RHO] = 0.; //rho init
  w[_INDEX_VELOCITY] = 0; // u init
  w[_INDEX_PRESSURE] = 0; // p init
  w[_INDEX_TEMP] = 0; // e ou T init
}

void Test_TransportVP_InitData(real *x, real *w) {
  real t = 0;
  Test_TransportVP_ImposedData(x, t, w);
}

real TransportVP_ImposedKinetic_Data(const real *x, const real t, real v) {
  real f;
  real pi = 4.0 * atan(1.0);
  real xnew = 0, vnew = 0;
  f = exp(-(v - t) * (v - t)) *
    exp(-36 * ((x[0] - v * t + 0.5 * t * t) - 0.5)
	* ((x[0] - v * t + 0.5 * t * t) - 0.5));
  return f;
}

void Test_TransportVP_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[_MV + 6];
  Test_TransportVP_ImposedData(x , t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
}

void UpdateVlasovPoisson(void *vf, real *w) {
  field *f = vf;
  
  int type_bc = 1;
  real bc_l = 0;
  real bc_r = 1;
    
  // Computation_charge_density(f,w);
  
  // Solving poisson
  SolvePoisson1D(f,w,type_bc,bc_l,bc_r);    
  
}

void PlotVlasovPoisson(void *vf, real *w) {
  real k_energy = 0, e_energy = 0, t_energy = 0;
  
  field *f = vf;
  
  Energies(f, w, k_energy, e_energy, t_energy);
  vf = f;
}

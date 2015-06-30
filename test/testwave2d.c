#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"


void TestPeriodic_Wave_ImposedData(const real *x, const real t, real *w);
void TestPeriodic_Wave_InitData(real *x, real *w);
void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);
void Wave_Upwind_NumFlux(real wL[],real wR[],real* vnorm,real* flux);

#define _SPEED_WAVE (1)

int main(void) {
  
  // unit tests
    
  int resu = Test_Wave_Periodic();
	 
  if (resu) printf("poisson test OK !\n");
  else printf("poisson test failed !\n");

  return !resu;
} 

int Test_Wave_Periodic(void) {

  bool test = true;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=3; 
  f.model.NumFlux=Wave_Upwind_NumFlux;
 
  //f.model.Source = NULL;
 
  f.model.InitData = TestPeriodic_Wave_InitData;
  f.model.ImposedData = TestPeriodic_Wave_ImposedData;
  f.model.BoundaryFlux = Wave_Upwind_BoundaryFlux;

  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;  // _M
  f.interp.interp_param[1] = 2;  // x direction degree
  f.interp.interp_param[2] = 2;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 8;  // x direction refinement
  f.interp.interp_param[5] = 8;  // y direction refinement
  f.interp.interp_param[6] = 1;  // z direction refinement
 // read the gmsh file

  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");

  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  real A[3][3] = {{2.0, 0, 0}, {0, 2.0, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&(f.macromesh),A,x0);

  f.macromesh.period[0]=2.0;
  f.macromesh.period[1]=2.0;
  BuildConnectivity(&(f.macromesh));

  // prepare the initial fields
  f.model.cfl = 0.5;
  Initfield(&f);
  f.vmax = _SPEED_WAVE; // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = NULL;
  // prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  real tmax = 0.1;
  real dt = set_dt(&f);
  RK2(&f, tmax, dt);
  //RK2(&f,0.03,0.05);

   // Save the results and the error
  Plotfield(0, false, &f, NULL, "dgvisu.msh");
  Plotfield(0, true, &f, "error", "dgerror.msh");

  real dd = L2error(&f);
  real tolerance = 9e-3;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  return test;
}

void Wave_Upwind_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real flux_temp=0;
  
  flux[0]=0.5*((wL[1]+wR[1])*vnorm[0] + (wL[2]+wR[2])*vnorm[1])+0.5*(wL[0]-wR[0]);
  
  flux_temp=0.5*(wL[0]+wR[0]) + 0.5*((wL[1]-wR[1])*vnorm[0] + (wL[2]-wR[2])*vnorm[1]);
  flux[1]=flux_temp*vnorm[0];
  flux[2]=flux_temp*vnorm[1];
 

  flux[0]=_SPEED_WAVE*flux[0];
  flux[1]=_SPEED_WAVE*flux[0];
  flux[2]=_SPEED_WAVE*flux[0];
  
};


void TestPeriodic_Wave_ImposedData(const real *x, const real t, real *w) {
  real pi=4.0*atan(1.0);

  w[0] = -pi*sqrt(2.0)*sin(pi*sqrt(2.0)*t)*cos(pi*x[0])*cos(pi*x[1]);
  w[1] = pi*cos(pi*sqrt(2.0)*t)*sin(pi*x[0])*cos(pi*x[1]);
  w[2] = pi*cos(pi*sqrt(2.0)*t)*cos(pi*x[0])*sin(pi*x[1]);
  

}

void TestPeriodic_Wave_InitData(real *x, real *w) {
  real t = 0;
  TestPeriodic_Wave_ImposedData(x, t, w);
}



void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[_MV + 6];
  TestPeriodic_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}



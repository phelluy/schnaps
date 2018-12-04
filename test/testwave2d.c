#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"

void TestSteady_Wave_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);
void TestSteady_Wave_InitData(schnaps_real *x, schnaps_real *w);
void TestSteady_Wave_Source(const schnaps_real *xy, const schnaps_real t, const schnaps_real *w, schnaps_real *S);
void Wave_Steady_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
			      schnaps_real *flux);
void Wave_Periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				schnaps_real *flux);

int Test_Wave_Periodic(void);
int Test_Wave_Steady(void);

int main(void) {
  
  // unit tests
    
  int resu1 = 0;
  int resu2 = 0;
  resu1=Test_Wave_Periodic();
  resu2=Test_Wave_Steady();
	 
  if (resu1 + resu2 >1) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu1;
} 

int Test_Wave_Periodic(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m = 3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = TestPeriodic_Wave_InitData;
  model.ImposedData = TestPeriodic_Wave_ImposedData;
  model.BoundaryFlux = Wave_Periodic_BoundaryFlux;
  model.Source = NULL;

  int deg[]={4, 4, 0};
  int raf[]={4, 4, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  schnaps_real tmax = 0.025;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%.12e\n", dd);

  schnaps_real tolerance = 0.002;

  test = test && (dd < tolerance);

  PlotFields(0,false, &simu, "p", "dgvisu_exp.msh");
  PlotFields(1,false, &simu, "u", "dgvisu_exu.msh");
  PlotFields(2,false, &simu, "v", "dgvisu_exv.msh");

  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  ThetaTimeScheme(&simu2, tmax, simu.dt);

  
  dd = L2error(&simu2);

  printf("erreur implicit L2=%.12e\n", dd);

  test = test && (dd < tolerance);

  PlotFields(0,false, &simu2, "p", "dgvisu_imp.msh");
  PlotFields(1,false, &simu2, "u", "dgvisu_imu.msh");
  PlotFields(2,false, &simu2, "v", "dgvisu_imv.msh");


  FreeMacroMesh(&mesh);

  return test;
}

int Test_Wave_Steady(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m=3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = TestSteady_Wave_InitData; 
  model.ImposedData = TestSteady_Wave_ImposedData; 
  model.BoundaryFlux = Wave_Steady_BoundaryFlux; 
  model.Source = TestSteady_Wave_Source; 

  int deg[]={3, 3, 0};
  int raf[]={2, 2, 1};
  
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);

  schnaps_real tmax = 0.01;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);
 
  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur explicit L2=%.12e\n", dd);

  PlotFields(0,false, &simu, "p", "dgvisu_exp.msh");
  PlotFields(1,false, &simu, "u", "dgvisu_exu.msh");
  PlotFields(2,false, &simu, "v", "dgvisu_exv.msh");

  schnaps_real tolerance = 0.0000001;

  test = test && (dd < tolerance);

  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  ThetaTimeScheme(&simu2, tmax, simu.dt);
  
  dd = L2error(&simu2);

  printf("erreur implicit L2=%.12e\n", dd);

  PlotFields(0,false, &simu2, "p", "dgvisu_imp.msh");
  PlotFields(1,false, &simu2, "u", "dgvisu_imu.msh");
  PlotFields(2,false, &simu2, "v", "dgvisu_imv.msh");

  test = test && (dd < tolerance);

  FreeMacroMesh(&mesh);
  
  return test;
}

void TestPeriodic_Wave_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  schnaps_real pi=4.0*atan(1.0);
  schnaps_real L=_LENGTH_DOMAIN;
  schnaps_real Coef=(2.0*pi)/L;
  schnaps_real a=_SPEED_WAVE;

  w[0] = -a*Coef*sqrt(2.0)*sin(a*Coef*sqrt(2.0)*t)*cos(Coef*x[0])*cos(Coef*x[1]);
  w[1] = a*Coef*cos(Coef*a*sqrt(2.0)*t)*sin(Coef*x[0])*cos(Coef*x[1]);
  w[2] = a*Coef*cos(Coef*a*sqrt(2.0)*t)*cos(Coef*x[0])*sin(Coef*x[1]);
  

}

void TestPeriodic_Wave_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestPeriodic_Wave_ImposedData(x, t, w);
}


void Wave_Periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[3];
  TestPeriodic_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}


void TestSteady_Wave_ImposedData(const schnaps_real *xy, const schnaps_real t, schnaps_real *w) {

  schnaps_real x=xy[0];
  schnaps_real y=xy[1];


  w[0] = x*(1-x)*y*(1-y)+1;
  w[1] = 2*x*(1-x)*y*(1-y)+2;
  w[2] = 3*x*(1-x)*y*(1-y)+3;


}


void TestSteady_Wave_Source(const schnaps_real *xy, const schnaps_real t, const schnaps_real *w, schnaps_real *S){
  
  schnaps_real x=xy[0];
  schnaps_real y=xy[1];

  S[0] = 2*(1-2*x)*(y*(1-y))+3*(1-2*y)*(x*(1-x));
  S[1] = (1-2*x)*(y*(1-y));
  S[2] = (1-2*y)*(x*(1-x));

  S[0] *= _SPEED_WAVE;
  S[1] *= _SPEED_WAVE;
  S[2] *= _SPEED_WAVE;

}

void TestSteady_Wave_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSteady_Wave_ImposedData(x, t, w);
}


void Wave_Steady_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[3];
  TestSteady_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}

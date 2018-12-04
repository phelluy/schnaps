#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"


//! \brief boundary rusanov flux based for a Steady State with velocity imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_U_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary HLL flux based for a Steady State with velocity imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_U_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary Roe flux based for a Steady State with velocity imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_U_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary rusanov flux based for a Steady State with pressure imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_P_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary HLL flux based for a Steady State with pressure imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_P_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary Roe flux based for a Steady State with pressure imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_P_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);


//! \brief boundary rusanov flux based for periodic Shallow Water case
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Periodic_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary HLL flux based for periodic Shallow Water case
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Periodic_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary Roe flux based for periodic Shallow Water case
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Periodic_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary rusanov flux based for equilibirum Shallow Water case
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Equilibrium_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary HLL flux based for equilibrium Shallow Water case
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Equilibrium_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

//! \brief boundary Roe flux based for equilibrium Shallow Water case
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Equilibrium_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);

int Test_SH_equilibrium(void);
int Test_SH_periodic(void);
int Test_SH_SteadyState_UImposed(void);
int Test_SH_SteadyState_PImposed(void);

int main(void) {
  
  // unit tests
    
  int resu1 =0;
  int resu2 =0;
  int resu3 =0;
  int resu4 =0;
  printf("/---------/ Test equilibrium /-----------/");
  resu1=Test_SH_equilibrium();
  printf("/----------------------------------------/");
  printf("/------------/ Test periodic /-----------/");
  resu2=Test_SH_periodic();
  printf("/----------------------------------------/");
  printf("/--------/ Test SteadyState u /----------/");
  resu3=Test_SH_SteadyState_UImposed();
  printf("/----------------------------------------/");
  printf("/--------/ Test SteadyState p /----------/");
  resu4=Test_SH_SteadyState_PImposed();	 
  if (resu1 + resu2 + resu3 + resu4 >3) printf("shallow water  test OK !\n");
  else printf("shallow water test failed !\n");

  return !resu1;
} 

int Test_SH_equilibrium(void) {
  int test1 = 0,test2 = 0,test3=0,test = 0;
  schnaps_real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  BuildConnectivity(&mesh);

  Model model;

  int deg[]={2, 2, 0};
  int raf[]={8, 8, 1};

  assert(mesh.is2d);
  CheckMacroMesh(&mesh, deg, raf);
  

  /******* Test for Rusanov ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_equilibrium_InitData;
  model.ImposedData = TestSH_equilibrium_ImposedData;
  model.BoundaryFlux = Equilibrium_Rusanov_BoundaryFlux;
  model.Source = ShallowWater_classical_SourceTerm;

  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu.vmax = _SPEED_WAVE;
  simu.cfl = 0.05;
  RK4(&simu, tmax);

  dd = L2error(&simu);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd);
  
  /******* Test for hLL ******/
  model.m=6; 
  model.NumFlux=ShallowWater_HLL_NumFlux;
  model.InitData = TestSH_equilibrium_InitData;
  model.ImposedData = TestSH_equilibrium_ImposedData;
  model.BoundaryFlux = Equilibrium_HLL_BoundaryFlux;
  model.Source = ShallowWater_classical_SourceTerm;

  Simulation simu2;
  EmptySimulation(&simu2);

  InitSimulation(&simu2, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu2.vmax = _SPEED_WAVE;
  simu2.cfl = 0.1;
  RK4(&simu2, tmax);

  dd = L2error(&simu2);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  /******* Test for Roe ******/
   model.m=6; 
  model.NumFlux=ShallowWater_Roe_NumFlux;
  model.InitData = TestSH_equilibrium_InitData;
  model.ImposedData = TestSH_equilibrium_ImposedData;
  model.BoundaryFlux = Equilibrium_Roe_BoundaryFlux;
  model.Source = ShallowWater_classical_SourceTerm;

  Simulation simu3;
  EmptySimulation(&simu3);

  InitSimulation(&simu3, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu3.vmax = _SPEED_WAVE;
  simu3.cfl = 0.1;
  RK4(&simu3, tmax);

  dd = L2error(&simu3);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test3=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  /*PlotFields(0,false, &simu3, "h", "dgvisu_h.msh");
  PlotFields(1,false, &simu3, "u1", "dgvisu_u1.msh");
  PlotFields(2,false, &simu3, "u2", "dgvisu_u2.msh");*/

  if(test1 +test2+test3 > 2){
    test=1;     
  }

    FreeMacroMesh(&mesh);
  
  return test;
}


int Test_SH_periodic(void) {
  int test1 = 0,test2 = 0,test3=0,test = 0;
  schnaps_real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0,dd2=0,dd3=0;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  BuildConnectivity(&mesh);

  Model model;

  int deg[]={2, 2, 0};
  int raf[]={8, 8, 1};

  assert(mesh.is2d);
  CheckMacroMesh(&mesh, deg, raf);
  

  /******* Test for Rusanov ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = Periodic_Rusanov_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu.vmax = 1;//_SPEED_WAVE;
  simu.cfl = 0.025;
  RK2(&simu, tmax);

  dd = L2error_onefield(&simu,0);
  dd2 = L2error_onefield(&simu,1);
  dd3 = L2error_onefield(&simu,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd+dd2+dd3);
  
  /******* Test for hLL ******/
  model.m=6; 
  model.NumFlux=ShallowWater_HLL_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = Periodic_HLL_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu2.vmax = 1;//_SPEED_WAVE;
  simu2.cfl = 0.025;
  RK2(&simu2, tmax);

  dd = L2error_onefield(&simu2,0);
  dd2 = L2error_onefield(&simu2,1);
  dd3 = L2error_onefield(&simu2,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd+dd2+dd3);

  /******* Test for Roe ******/
   model.m=6; 
  model.NumFlux=ShallowWater_Roe_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = Periodic_Roe_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  Simulation simu3;
  EmptySimulation(&simu3);
  InitSimulation(&simu3, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu3.vmax = 1;//_SPEED_WAVE;
  simu3.cfl = 0.025;
  RK2(&simu3, tmax);

  dd = L2error_onefield(&simu3,0);
  dd2 = L2error_onefield(&simu3,1);
  dd3 = L2error_onefield(&simu3,2);
  tolerance = 5e-3;

  schnaps_real dd4 = L2error_onefield(&simu2,3);
  schnaps_real dd5 = L2error_onefield(&simu2,4);
  schnaps_real dd6 = L2error_onefield(&simu2,5);

  printf("erreur h L2=%.12e\n", dd);
  printf("erreur u1 L2=%.12e\n", dd2);
  printf("erreur u2 L2=%.12e\n", dd3);
  printf("erreur b L2=%.12e\n", dd4);
  printf("erreur bx L2=%.12e\n", dd5);
  printf("erreur by L2=%.12e\n", dd6);
  
  if(dd+dd2+dd3 < tolerance){
    test3=1;     
  }
  printf("L2 Roe error %.8e\n", dd+dd2+dd3);

  if(test1 +test2+test3 > 2){
    test=1;     
  }

    FreeMacroMesh(&mesh);

  return test;
}

int Test_SH_SteadyState_UImposed(void) {
  int test1 = 0,test2 = 0,test3=0,test = 0;
  schnaps_real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0,dd2=0,dd3=0;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  BuildConnectivity(&mesh);

  Model model;

  int deg[]={4, 4, 0};
  int raf[]={4, 4, 1};

  assert(mesh.is2d);
  CheckMacroMesh(&mesh, deg, raf);
  

  /******* Test for Rusanov ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_SteadyState_U_InitData;
  model.ImposedData = TestSH_SteadyState_U_ImposedData;
  model.BoundaryFlux = SteadyState_U_Rusanov_BoundaryFlux;
  model.Source = ShallowWater_SteadyState_U_SourceTerm;

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  tmax = 0.002;
  simu.vmax = 1;//_SPEED_WAVE;
  simu.cfl = 0.025;
  RK2(&simu, tmax);

  dd = L2error_onefield(&simu,0);
  dd2 = L2error_onefield(&simu,1);
  dd3 = L2error_onefield(&simu,2);
  tolerance = 5e-8;
  
  if(dd+dd2+dd3 < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd+dd2+dd3);
  
  /******* Test for hLL ******/
  model.m=6; 
  model.NumFlux=ShallowWater_HLL_NumFlux;
  model.InitData = TestSH_SteadyState_U_InitData;
  model.ImposedData = TestSH_SteadyState_U_ImposedData;
  model.BoundaryFlux = SteadyState_U_HLL_BoundaryFlux;
  model.Source = ShallowWater_SteadyState_U_SourceTerm;

  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  
  tmax = 0.002;
  simu2.vmax = 1;//_SPEED_WAVE;
  simu2.cfl = 0.025;
  RK2(&simu2, tmax);

  dd = L2error_onefield(&simu2,0);
  dd2 = L2error_onefield(&simu2,1);
  dd3 = L2error_onefield(&simu2,2);
  tolerance = 5e-8;
  
  if(dd+dd2+dd3 < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd+dd2+dd3);

  /******* Test for Roe ******/
   model.m=6; 
  model.NumFlux=ShallowWater_Roe_NumFlux;
  model.InitData = TestSH_SteadyState_U_InitData;
  model.ImposedData = TestSH_SteadyState_U_ImposedData;
  model.BoundaryFlux = SteadyState_U_Roe_BoundaryFlux;
  model.Source = ShallowWater_SteadyState_U_SourceTerm;

  Simulation simu3;
  EmptySimulation(&simu3);
  InitSimulation(&simu3, &mesh, deg, raf, &model);
  
  tmax = 0.002;
  simu3.vmax = 1;//_SPEED_WAVE;
  simu3.cfl = 0.025;
  RK2(&simu3, tmax);

  dd = L2error_onefield(&simu3,0);
  dd2 = L2error_onefield(&simu3,1);
  dd3 = L2error_onefield(&simu3,2);
  tolerance = 5e-8;
  
  printf("erreur h L2=%.12e\n", dd);
  printf("erreur u1 L2=%.12e\n", dd2);
  printf("erreur u2 L2=%.12e\n", dd3);
  
  if(dd+dd2+dd3 < tolerance){
    test3=1;     
  }
  printf("L2 Roe error %.8e\n", dd+dd2+dd3);

  if(test1 +test2+test3 > 2){
    test=1;     
  }
    FreeMacroMesh(&mesh);
  return test;
}

int Test_SH_SteadyState_PImposed(void) {
  int test1 = 0,test2 = 0,test3=0,test = 0;
  schnaps_real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0,dd2=0,dd3=0;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  BuildConnectivity(&mesh);

  Model model;
  int deg[]={4, 4, 0};
  int raf[]={4, 4, 1};
  assert(mesh.is2d);
  CheckMacroMesh(&mesh, deg, raf);
  

  /******* Test for Rusanov ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_SteadyState_P_InitData;
  model.ImposedData = TestSH_SteadyState_P_ImposedData;
  model.BoundaryFlux = SteadyState_P_Rusanov_BoundaryFlux;
  model.Source = ShallowWater_SteadyState_P_SourceTerm;

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  tmax = 0.002;
  simu.vmax = 1;//_SPEED_WAVE;
  simu.cfl = 0.025;
  RK2(&simu, tmax);

  dd = L2error_onefield(&simu,0);
  dd2 = L2error_onefield(&simu,1);
  dd3 = L2error_onefield(&simu,2);
  tolerance = 5e-8;
  
  if(dd+dd2+dd3 < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd+dd2+dd3);
  
  /******* Test for hLL ******/
  model.m=6; 
  model.NumFlux=ShallowWater_HLL_NumFlux;
  model.InitData = TestSH_SteadyState_P_InitData;
  model.ImposedData = TestSH_SteadyState_P_ImposedData;
  model.BoundaryFlux = SteadyState_P_HLL_BoundaryFlux;
  model.Source = ShallowWater_SteadyState_P_SourceTerm;

  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  
  tmax = 0.002;
  simu2.vmax = 1;//_SPEED_WAVE;
  simu2.cfl = 0.025;
  RK2(&simu2, tmax);

  dd = L2error_onefield(&simu2,0);
  dd2 = L2error_onefield(&simu2,1);
  dd3 = L2error_onefield(&simu2,2);
  tolerance = 5e-8;
  
  if(dd+dd2+dd3 < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd+dd2+dd3);

  /******* Test for Roe ******/
   model.m=6; 
  model.NumFlux=ShallowWater_Roe_NumFlux;
  model.InitData = TestSH_SteadyState_P_InitData;
  model.ImposedData = TestSH_SteadyState_P_ImposedData;
  model.BoundaryFlux = SteadyState_P_Roe_BoundaryFlux;
  model.Source = ShallowWater_SteadyState_P_SourceTerm;

  Simulation simu3;
  EmptySimulation(&simu3);
  InitSimulation(&simu3, &mesh, deg, raf, &model);
  
  tmax = 0.002;
  simu3.vmax = 1;//_SPEED_WAVE;
  simu3.cfl = 0.025;
  RK2(&simu3, tmax);

  dd = L2error_onefield(&simu3,0);
  dd2 = L2error_onefield(&simu3,1);
  dd3 = L2error_onefield(&simu3,2);
  tolerance = 5e-8;

  printf("erreur h L2=%.12e\n", dd);
  printf("erreur u1 L2=%.12e\n", dd2);
  printf("erreur u2 L2=%.12e\n", dd3);
  
  if(dd+dd2+dd3 < tolerance){
    test3=1;     
  }
  printf("L2 Roe error %.8e\n", dd+dd2+dd3);

  if(test1 +test2+test3 > 2){
    test=1;     
  }
    FreeMacroMesh(&mesh);

  return test;
}

void TestSH_equilibrium_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {

  w[0] = 1.0-0.8*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[1] = 0.0;
  w[2] = 0.0;
  w[3] = 0.8*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[4] = -0.8*50*(2.0*x[0]-_LENGTH_DOMAIN)*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[5] = -0.8*50*(2.0*x[1]-_LENGTH_DOMAIN)*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
}

void TestSH_equilibrium_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSH_equilibrium_ImposedData(x, t, w);
}


void Equilibrium_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void Equilibrium_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void Equilibrium_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}

void TestSH_periodic_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  schnaps_real u0=2;
  schnaps_real v0=2;
  schnaps_real pi=4.0*atan(1.0);
  
  w[0] = 1.0+0.2*sin(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t));
  w[1] = u0*w[0];
  w[2] = v0*w[0];
  w[3] = 1.0-(1.0+0.2*sin(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t)));
  w[4] = -0.2*((2.*pi)/_LENGTH_DOMAIN)*cos(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t));
  w[5] = -0.2*((2.*pi)/_LENGTH_DOMAIN)*cos(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t));
}

void TestSH_periodic_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSH_periodic_ImposedData(x, t, w);
}

void ShallowWater_periodic_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source){
  schnaps_real g=_GRAVITY;
  schnaps_real hL=0, hR=0, uL=0, uR=0, vL=0, vR=0;
  schnaps_real S=0,b=0,bx=0,by=0;
  schnaps_real wexact[6]; 
  hL = w[0];
  uL=w[1]/w[0];
  vL=w[2]/w[0];

  TestSH_periodic_ImposedData(x , t, wexact);

  source[0]= 0.0;
  source[1]= -g*hL*wexact[4];
  source[2]= -g*hL*wexact[5];
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
};

void Periodic_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void Periodic_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void Periodic_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}


void TestSH_SteadyState_U_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  schnaps_real alpha=1.0; 
  w[0] = 1.0+alpha*((x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]));
  w[1] = ((x[0]-x[0]*x[0])*(1.0-2*x[1]))*w[0];
  w[2] = -((x[1]-x[1]*x[1])*(1.0-2*x[0]))*w[0];
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = 0.0;
}

void TestSH_SteadyState_U_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSH_SteadyState_U_ImposedData(x, t, w);
}

void ShallowWater_SteadyState_U_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source){
  schnaps_real g=_GRAVITY;
  schnaps_real alpha=1.0;
  schnaps_real S_11, S_12,S_13, S_21, S_22, S_23, S_factor;
  schnaps_real wexact[6];

  S_factor=alpha*(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])+1.0;

  S_11=(x[0]-x[0]*x[0])*(1.0-2.0*x[1])*(1.0-2.0*x[1]);
  S_12=2.0*(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
  S_13=g*alpha*(x[1]-x[1]*x[1]);

  S_21=(x[1]-x[1]*x[1])*(1.0-2.0*x[0])*(1.0-2.0*x[0]);
  S_22=2.0*(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
  S_23=g*alpha*(x[0]-x[0]*x[0]);

  source[0]= 0.0;
  source[1]= S_factor*(1.0-2.0*x[0])*(S_11+S_12+S_13);
  source[2]= S_factor*(1.0-2.0*x[1])*(S_21+S_22+S_23);
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
};

void SteadyState_U_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_SteadyState_U_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void SteadyState_U_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_SteadyState_U_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}

void SteadyState_U_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_SteadyState_U_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}



void TestSH_SteadyState_P_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  schnaps_real g=_GRAVITY;
  schnaps_real u01=2.0,u02=2.0; 
  
  w[0] = sqrt(2.0/g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])));
  w[1] = u01*w[0];
  w[2] = u02*w[0];
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = 0.0;
}

void TestSH_SteadyState_P_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSH_SteadyState_P_ImposedData(x, t, w);
}


void ShallowWater_SteadyState_P_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source){
  schnaps_real g=_GRAVITY;
  schnaps_real S_h1, S_h2;
  schnaps_real wexact[6];
  schnaps_real u01=2.0,u02=2.0; 

  S_h1=1.0/sqrt(2.0*g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])));
  S_h2=u01*(1-2.0*x[0])*(x[1]-x[1]*x[1])+u02*(1-2.0*x[1])*(x[0]-x[0]*x[0]);
  
  source[0]= S_h1*S_h2;
  source[1]= (1-2.0*x[0])*(x[1]-x[1]*x[1])+u01*source[0];
  source[2]= (1-2.0*x[1])*(x[0]-x[0]*x[0])+u02*source[0];
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
};

void SteadyState_P_Rusanov_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_SteadyState_P_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void SteadyState_P_HLL_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_SteadyState_P_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void SteadyState_P_Roe_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm, schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_SteadyState_P_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}

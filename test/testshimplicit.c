#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"
#include "implicit.h"

int Test_SH_equilibrium_Implicit(void);
void ShallowWater_Rusanov_equilibrium_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
						   schnaps_real *flux);
void ShallowWater_Rusanov_periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
						   schnaps_real *flux);

int main(void) {
  
  // unit tests

#ifdef _DOUBLE_PRECISION
    
  int resu = Test_SH_equilibrium_Implicit();
	 
  if (resu ==1) printf("shallow water equilibrium  test OK !\n");
  else printf("shallow water equilibrium test failed !\n");

  return !resu;

#else
  printf("shallow water equilibrium  test not available in single precision !\n");
 return 0;
#endif
} 

int Test_SH_equilibrium_Implicit(void) {

  int test1 = 0,test2 = 0,test = 0;
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
  int raf[]={6, 6, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  

  /******* Test equilibrium  ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_equilibrium_InitData;
  model.ImposedData = TestSH_equilibrium_ImposedData;
  model.BoundaryFlux = ShallowWater_Rusanov_equilibrium_BoundaryFlux;
  model.Source = ShallowWater_classical_SourceTerm;


  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);

  
  tmax = 0.01;
  simu.vmax = 1;
  simu.cfl = 0.2;
  simu.dt = Get_Dt_RK(&simu);
  
  ThetaTimeScheme_WithJF(&simu,tmax,simu.dt);
  
  dd = L2error(&simu);
  printf("erreur equilibirum L2=%.12e\n", dd);

  PlotFields(0,false, &simu, "h", "dgvisu_h.msh");
  PlotFields(1,false, &simu, "u1", "dgvisu_u1.msh");
  PlotFields(2,false, &simu, "u2", "dgvisu_u2.msh");

  PlotFields(0,true, &simu, "h", "dgerror_h.msh");
  PlotFields(1,true, &simu, "u1", "dgerror_u1.msh");
  PlotFields(2,true, &simu, "u2", "dgerror_u2.msh");
  
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test1=1;     
  }


  /******* Test periodic  ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = ShallowWater_Rusanov_periodic_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  int deg2[]={2, 2, 0};
  int raf2[]={8, 8, 1};

  Simulation simu2;
  EmptySimulation(&simu2);

  InitSimulation(&simu2, &mesh, deg2, raf2, &model);

  
  tmax=0.01;
  simu2.vmax = 1;
  simu2.cfl = 0.025;
  simu2.dt = Get_Dt_RK(&simu2);
  
  ThetaTimeScheme_WithJF(&simu2,tmax,simu2.dt);
  
  schnaps_real dd1 = L2error_onefield(&simu2,0);
  schnaps_real dd2 = L2error_onefield(&simu2,1);
  schnaps_real dd3 = L2error_onefield(&simu2,2);
  
  printf("erreur periodic L2=%.12e\n", dd1+dd2+dd3);

  PlotFields(0,false, &simu2, "h", "dgvisu_periodic_h.msh");
  PlotFields(1,false, &simu2, "u1", "dgvisu_periodic_u1.msh");
  PlotFields(2,false, &simu2, "u2", "dgvisu_periodic_u2.msh");

  tolerance = 5e-3;
  
  if(dd1+dd2+dd3 < tolerance){
    test2=1;     
  }
 

  if(test1 + test2 > 1){
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


void ShallowWater_Rusanov_equilibrium_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
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

void ShallowWater_Rusanov_periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}





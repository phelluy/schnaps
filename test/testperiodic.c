#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"

int TestPeriodic(void);


void TestPeriodic_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);
void TestPeriodic_InitData(schnaps_real *x, schnaps_real *w);
void TestPeriodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
			       schnaps_real *flux);

int main(void) {
  // unit tests
    
  int resu = TestPeriodic();
	 
  if (resu) printf("periodic test OK !\n");
  else printf("periodic test failed !\n");

  return !resu;
} 

int TestPeriodic(void) {
  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect1DMacroMesh(&mesh);
  assert(mesh.is1d);
  // periodic mesh
  mesh.period[0] = 1;
  BuildConnectivity(&mesh);
  int deg[]={2, 0, 0};
  int raf[]={16, 1, 1};

  CheckMacroMesh(&mesh,deg,raf);
  PrintMacroMesh(&mesh);


  Model model;
  schnaps_real degV=2;
  schnaps_real nbEV=24;
  KineticData * kd=&schnaps_kinetic_data;

  InitKineticData(kd,nbEV,degV);
 
  
  model.m=kd->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.Source = NULL;
  model.BoundaryFlux = TestPeriodic_BoundaryFlux;
  model.InitData = TestPeriodic_InitData;
  model.ImposedData = TestPeriodic_ImposedData;

  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  printf("cfl param =%f\n",simu.hmin);

  simu.vmax = kd->vmax; // maximal wave speed 
  simu.cfl = 0.05;
  schnaps_real tmax = 0.4;
 
  RK4(&simu,tmax);

  // save the results and the error
  PlotFields(0, false, &simu, "sol","dgvisu.msh");
  PlotFields(0, true, &simu, "error","dgerror.msh");

  schnaps_real dd = L2error(&simu);

  printf("erreur L2: %lf\n", dd);
  test = test && (dd<2e-1);

  FreeMacroMesh(&mesh);

  return test;
}

void TestPeriodic_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[])
{
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real pi = 4 * atan(1.0);
  for(int i = 0; i < kd->index_max_kin + 1 ; ++i) {
    int j = i % kd->deg_v; // local connectivity put in function
    int nel = i / kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax + nel * kd->dv + kd->dv * glop(kd->deg_v, j));

    w[i] = cos(2 * pi * ( x[0] - vi * t) );
  }
  // exact value of the potential and electric field
  w[kd->index_phi] = 0;
  w[kd->index_ex] = 0;
  w[kd->index_ey] = 0;
  w[kd->index_ez] = 0;
  w[kd->index_rho] = 2.0; //rho init
  w[kd->index_u] = 0; // u init
  w[kd->index_P] = 0; // p init
  w[kd->index_T] = 0; // e ou T init
}

void TestPeriodic_InitData(schnaps_real x[3], schnaps_real w[])
{
  schnaps_real t = 0;
  TestPeriodic_ImposedData(x, t, w);
}

void TestPeriodic_BoundaryFlux(schnaps_real x[3], schnaps_real t, schnaps_real wL[], schnaps_real *vnorm,
			       schnaps_real* flux)
{
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  TestPeriodic_ImposedData(x, t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
  assert(false);
}



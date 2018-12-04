#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "solverpoisson.h"
#include "linear_solver.h"

int TestPoisson2d(void) ;
void TestPoisson_ImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]);
void TestPoisson_InitData(schnaps_real x[3],schnaps_real w[]);
void TestPoisson_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
			      schnaps_real* flux);

int main(void) 
{
  
  // unit tests
    
  int resu = TestPoisson2d();
	 
  if (resu) printf("2d poisson test OK !\n");
  else printf("2d poisson test failed !\n");

  return !resu;
}

int TestPoisson2d(void) 
{
  bool test = true;

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testdisque2d.msh");
  Detect2DMacroMesh(&mesh);
  bool is2d=mesh.is2d;
  assert(is2d);
  BuildConnectivity(&mesh);

  Model model;


  schnaps_real degV=2;
  schnaps_real nbEV=12;
  KineticData * kd=&schnaps_kinetic_data;

  InitKineticData(kd,nbEV,degV);
 
  model.m=kd->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.Source = VlasovP_Lagrangian_Source;
  
  model.BoundaryFlux = TestPoisson_BoundaryFlux;
  model.InitData = TestPoisson_InitData;
  model.ImposedData = TestPoisson_ImposedData;
 
  int deg[]={2, 2, 0};
  int raf[]={1, 1, 1};
   
 
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=kd->index_phi;
  
  InitContinuousSolver(&ps,&simu,1,nb_var,listvar);

  ps.matrix_assembly=ContinuousOperator_Poisson2D;
  ps.rhs_assembly=RHSPoisson_Continuous;
  ps.bc_assembly= ExactDirichletContinuousMatrix;
  ps.postcomputation_assembly=Computation_ElectricField_Poisson;

  ps.lsol.solver_type = GMRES;
  ps.lsol.pc_type=NONE;


  SolveContinuous2D(&ps);

  schnaps_real errl2 = L2error(&simu);

  printf("Erreur L2=%f\n",errl2);

  test = test && (errl2 < 2e-2);

  printf("Plot...\n");


  PlotFields(kd->index_rho, false, &simu, NULL, "dgvisu.msh");
  PlotFields(kd->index_ex, false, &simu, NULL, "dgex.msh");

#ifdef PARALUTION 
  paralution_end();
#endif

  FreeMacroMesh(&mesh);


  return test;
}

void TestPoisson_ImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData * kd=&schnaps_kinetic_data;
  for(int i = 0; i < kd->index_max+1; i++){
    w[i] = 0;
  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi] = (x[0] * x[0] + x[1] * x[1])/4;
  w[kd->index_ex] =  -x[0]/2;
  w[kd->index_ey] =  -x[1]/2;
  w[kd->index_ez] =  0;
  w[kd->index_rho] = -1; //rho init
  /* w[_INDEX_PHI] = x[0] ; */
  /* w[_INDEX_EX] =  -1; */
  /* w[_INDEX_RHO] = 0; //rho init */
}

void TestPoisson_InitData(schnaps_real x[3], schnaps_real w[])
{
  schnaps_real t = 0;
  TestPoisson_ImposedData(x, t, w);
}

void TestPoisson_BoundaryFlux(schnaps_real x[3], schnaps_real t, schnaps_real wL[], schnaps_real *vnorm, 
			      schnaps_real *flux)
{
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  TestPoisson_ImposedData(x, t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
}



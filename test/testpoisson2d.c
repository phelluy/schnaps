#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "solverpoisson.h"


void TestPoisson_ImposedData(const real x[3],const real t,real w[]);
void TestPoisson_InitData(real x[3],real w[]);
void TestPoisson_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
			      real* flux);

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

  field f;
  init_empty_field(&f);

  int vec=1;

  // num of conservative variables f(vi) for each vi, phi, E, rho, u,
  // p, e (ou T)
  f.model.m=_INDEX_MAX; 
  f.vmax = _VMAX; // maximal wave speed
  f.model.NumFlux = VlasovP_Lagrangian_NumFlux;
  f.model.Source = VlasovP_Lagrangian_Source;
   //f.model.Source = NULL;
  
  f.model.BoundaryFlux = TestPoisson_BoundaryFlux;
  f.model.InitData = TestPoisson_InitData;
  f.model.ImposedData = TestPoisson_ImposedData;
  f.model.Source = NULL;
 
  f.varindex = GenericVarindex;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL; 
   
    
  f.interp.interp_param[0] = f.model.m;  // _M
  f.interp.interp_param[1] = 3;  // x direction degree
  f.interp.interp_param[2] = 3;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 2;  // x direction refinement
  f.interp.interp_param[5] = 2;  // y direction refinement
  f.interp.interp_param[6] = 1;  // z direction refinement
  // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testdisque2d.msh");
  //ReadMacroMesh(&(f.macromesh),"geo/square.msh");
  // try to detect a 2d mesh
  //bool is1d=Detect1DMacroMesh(&(f.macromesh));
  Detect2DMacroMesh(&(f.macromesh));
  bool is2d=f.macromesh.is2d;
  assert(is2d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));
  //assert(1==2);
  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  Initfield(&f);
  f.nb_diags=0;
  
  // prudence...
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);


  PoissonSolver ps;

  InitPoissonSolver(&ps,&f,_INDEX_PHI);

  SolvePoisson2D(&ps,_Dirichlet_Poisson_BC);

  real errl2 = L2error(&f);

  printf("Erreur L2=%f\n",errl2);

  test = test && (errl2 < 4e-4);

  printf("Plot...\n");


  Plotfield(_INDEX_PHI, false, &f, NULL, "dgvisu.msh");
  Plotfield(_INDEX_EX, false, &f, NULL, "dgex.msh");


  return test;
}

void TestPoisson_ImposedData(const real x[3], const real t, real w[])
{
  for(int i = 0; i < _INDEX_MAX; i++){
    w[i] = 0;
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI] = (x[0] * x[0] + x[1] * x[1])/4;
  w[_INDEX_EX] =  -x[0]/2;
  w[_INDEX_RHO] = -1; //rho init
  /* w[_INDEX_PHI] = x[0] ; */
  /* w[_INDEX_EX] =  -1; */
  /* w[_INDEX_RHO] = 0; //rho init */
}

void TestPoisson_InitData(real x[3], real w[])
{
  real t = 0;
  TestPoisson_ImposedData(x, t, w);
}

void TestPoisson_BoundaryFlux(real x[3], real t, real wL[], real *vnorm, 
			      real *flux)
{
  real wR[_INDEX_MAX];
  TestPoisson_ImposedData(x, t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
}



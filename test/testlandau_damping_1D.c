#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"


void Test_Landau_Damping_ImposedData(const real x[3], const real t,real w[]);
void Test_Landau_Damping_InitData(real x[3],real w[]);
void Test_Landau_Damping_BoundaryFlux(real x[3],real t,real wL[],real* vnorm, real* flux);

void UpdateVlasovPoisson(void* field, real *w);
void PlotVlasovPoisson(void* vf, real * w);

int main(void) {
  
  // unit tests
    
  int resu=TestLandau_Damping_1D();
	 
  if (resu) printf("landau test OK !\n");
  else printf("landau test failed !\n");

  return !resu;
} 


int TestLandau_Damping_1D(void) {

  bool test=true;

  field f;
  init_empty_field(&f);

  int vec=1;
  real k=0.5;
  real pi=4.0*atan(1.0);
  
  f.model.m=_INDEX_MAX; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  f.model.NumFlux=VlasovP_Lagrangian_NumFlux;
 
  //f.model.Source = NULL;
 
  f.model.InitData=Test_Landau_Damping_InitData;
  f.model.ImposedData=Test_Landau_Damping_ImposedData;
  f.model.BoundaryFlux=Test_Landau_Damping_BoundaryFlux;

  f.varindex=GenericVarindex;
    
  f.interp.interp_param[0]=f.model.m;  // _M
  f.interp.interp_param[1]=3;  // x direction degree
  f.interp.interp_param[2]=0;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=32;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement
 // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  
  real A[3][3] = {{2*pi/k, 0, 0}, {0, 1, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&(f.macromesh),A,x0);
  
  // try to detect a 2d mesh
  Detect1DMacroMesh(&(f.macromesh));
  bool is1d=f.macromesh.is1d;
  assert(is1d);

  // mesh preparation
  f.macromesh.period[0]=2.0*pi/k;
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  f.model.cfl=0.05;
  Initfield(&f);
  f.vmax = _VMAX; // maximal wave speed
  f.macromesh.is1d=true;
  //f.macromesh.is1d=true;
  f.nb_diags=3;
  f.pre_dtfield=UpdateVlasovPoisson;
  f.update_after_rk=PlotVlasovPoisson;
  f.model.Source = VlasovP_Lagrangian_Source;
  // prudence...
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);

  real dt = set_dt(&f);
  RK2(&f,0.1, dt);
  //RK2(&f,0.03,0.05);

   // save the results and the error
  int iel=_NB_ELEM_V/2;
  int iloc=0;//_DEG_V;
  printf("Trace vi=%f\n",-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc));
  Plotfield(iloc+iel*_DEG_V,(1==0),&f,"sol","dgvisu.msh");
  Plotfield(iloc+iel*_DEG_V,(1==1),&f,"error","dgerror.msh");
  Plot_Energies(&f, dt);

  test= 1;

  return test;

}

void Test_Landau_Damping_ImposedData(const real x[3], const real t, real w[])
{
  //parameters of the case
  
  real k=0.5;
  real eps = 0.005;
  real my_pi= 4.0*atan(1.0);
  
  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    real vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));
 
    w[i]=(1.0+eps*cos(k*x[0]))*(1.0/sqrt(2.0*my_pi))*exp(-(vi*vi)/2.0);

  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI]=-eps*cos(k*x[0]);
  w[_INDEX_EX]=(eps/k)*sin(k*x[0]);
  w[_INDEX_RHO]=1.; //rho init
  w[_INDEX_VELOCITY]=0; // u init
  w[_INDEX_PRESSURE]=0; // p init
  w[_INDEX_TEMP]=0; // e ou T init

};

void Test_Landau_Damping_InitData(real x[3],real w[]){

  real t=0;
  Test_Landau_Damping_ImposedData(x,t,w);

};



void Test_Landau_Damping_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
				       real* flux){
  real wR[_MV+6];
  Test_Landau_Damping_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
  assert(1==2);
};


void UpdateVlasovPoisson(void* vf, real * w){
  int type_bc;
  real bc_l, bc_r;
  int itermax;
  field* f=vf;
  type_bc=2;
  bc_l=0;
  bc_r=0;
    
  Computation_charge_density(f,w);
  
  SolvePoisson1D(f,w,type_bc,bc_l,bc_r);    
  
}


void PlotVlasovPoisson(void* vf, real * w){
  real k_energy=0,e_energy=0,t_energy=0;
  
  field* f=vf;
  
  Energies(f,w,k_energy,e_energy,t_energy);
  vf=f;
}

#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "solverpoisson.h"


void TestPoisson_ImposedData(double x[3],double t,double w[]);
void TestPoisson_InitData(double x[3],double w[]);
void TestPoisson_BoundaryFlux(double x[3],double t,double wL[],double* vnorm, double* flux);

int main(void) {
  
  // unit tests
    
  int resu=TestPoisson();
	 
  if (resu) printf("poisson test OK !\n");
  else printf("poisson test failed !\n");

  return !resu;
} 


int TestPoisson(void) {

  bool test=true;

  field f;

  int vec=1;
  
  f.model.m=_MV+6; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  f.vmax = _VMAX; // maximal wave speed
  f.model.NumFlux=VlasovP_Lagrangian_NumFlux;
  f.model.Source = VlasovP_Lagrangian_Source;
   //f.model.Source = NULL;
  
  f.model.BoundaryFlux=TestPoisson_BoundaryFlux;
  f.model.InitData=TestPoisson_InitData;
  f.model.ImposedData=TestPoisson_ImposedData;
 
  
  f.varindex=GenericVarindex;
  f.update_before_rk=NULL;
  f.update_after_rk=NULL; 
    
    
  f.interp.interp_param[0]=f.model.m;  // _M
  f.interp.interp_param[1]=3;  // x direction degree
  f.interp.interp_param[2]=0;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=32;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement
  // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  // try to detect a 2d mesh
  //bool is1d=Detect1DMacroMesh(&(f.macromesh));
  Detect1DMacroMesh(&(f.macromesh));
  bool is1d=f.macromesh.is1d;
  assert(is1d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));


  PrintMacroMesh(&(f.macromesh));
  //assert(1==2);
  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  Initfield(&f);
  f.macromesh.is1d=true;
  f.is1d=true;
  f.nb_diags=0;
  
  // prudence...
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);

  // time derivative
  //dtField(&f);
  //DisplayField(&f);
  //assert(1==2);
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
 
  /*Compute_electric_field(&f);

  // check the gradient on every glop
  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    printf("elem %d\n",ie);
    for(int ipg=0;ipg<NPG(f.interp_param+1);ipg++){
      double xref[3],wpg;
      ref_pg_vol(f.interp_param+1,ipg,xref,&wpg,NULL);
      printf("Gauss point %d %f %f %f \n",ipg,xref[0],xref[1],xref[2]);
      int imem=f.varindex(f.interp_param,ie,ipg,_MV+1);
      printf("gradphi exact=%f gradphinum=%f\n",1-2*xref[0],f.wn[imem]);
      test=test && (fabs(f.wn[imem]-(1-2*xref[0]))<1e-10);
    }
    }*/

  //Computation_charge_density(f);
  
  SolvePoisson(&f,1,0.0,0.0);

  // check the gradient given by the poisson solver
  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    printf("elem %d\n",ie);
    for(int ipg=0;ipg<NPG(f.interp_param+1);ipg++){
      double xref[3],wpg;
      ref_pg_vol(f.interp_param+1,ipg,xref,&wpg,NULL);
      printf("Gauss point %d %f %f %f \n",ipg,xref[0],xref[1],xref[2]);
      int imem=f.varindex(f.interp_param,ie,ipg,_MV+1);
      printf("gradphi exact=%f gradphinum=%f rap=%f\n",
	     1-2*xref[0],f.wn[imem],(1-2*xref[0])/f.wn[imem]);
      test=test && (fabs(f.wn[imem]-(1-2*xref[0]))<1e-10);
    }
  }

  return test;

};

void TestPoisson_ImposedData(double x[3],double t,double w[]){

  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));

    w[i]=1./_VMAX;
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI]=x[0]*(1-x[0]);
  w[_INDEX_EX]=1.-2.*x[0];
  w[_INDEX_RHO]=2.; //rho init
  w[_INDEX_VELOCITY]=0; // u init
  w[_INDEX_PRESSURE]=0; // p init
  w[_INDEX_TEMP]=0; // e ou T init

};

void TestPoisson_InitData(double x[3],double w[]){

  double t=0;
  TestPoisson_ImposedData(x,t,w);

};


void TestPoisson_BoundaryFlux(double x[3],double t,double wL[],double* vnorm,
				       double* flux){
  double wR[_MV+6];
  TestPoisson_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};



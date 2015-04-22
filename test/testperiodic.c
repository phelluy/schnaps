#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"


void TestPeriodic_ImposedData(double x[3],double t,double w[]);
void TestPeriodic_InitData(double x[3],double w[]);
void TestPeriodic_BoundaryFlux(double x[3],double t,double wL[],double* vnorm, double* flux);

int main(void) {
  
  // unit tests
    
  int resu=TestPeriodic();
	 
  if (resu) printf("periodic test OK !\n");
  else printf("periodic test failed !\n");

  return !resu;
} 


int TestPeriodic(void) {

  bool test=true;

 #ifndef _PERIOD
  printf("peridicity disabled\n");
  return test;
#endif

  field f;

  int vec=1;
  
  f.model.m=_MV+1; // num of conservative variables
  f.vmax = _VMAX; // maximal wave speed 
  f.model.NumFlux=VlasovP_Lagrangian_NumFlux;
   f.model.Source = NULL;
  
  f.model.BoundaryFlux=TestPeriodic_BoundaryFlux;
  f.model.InitData=TestPeriodic_InitData;
  f.model.ImposedData=TestPeriodic_ImposedData;
 
  f.varindex=GenericVarindex;
  f.update_before_rk=NULL;
  f.update_after_rk=NULL; 
    
    
  f.interp.interp_param[0]=_MV+1;  // _M
  f.interp.interp_param[1]=3;  // x direction degree
  f.interp.interp_param[2]=0;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=10;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement
  // read the gmsh file
  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  // try to detect a 2d mesh
  Detect1DMacroMesh(&(f.macromesh));
  bool is1d=f.macromesh.is1d;
  assert(is1d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

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
  //RK2(&f,0.5,0.1);
  RK2(&f,0.5);
 
  // save the results and the error
  Plotfield(0,(1==0),&f,"sol","dgvisu.msh");
  Plotfield(0,(1==1),&f,"error","dgerror.msh");

  double dd=L2error(&f);
  double dd_Kinetic=L2_Kinetic_error(&f);
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  test= test && (dd<3e-3);


  //SolvePoisson(&f);

  return test;

};


void TestPeriodic_ImposedData(double x[3],double t,double w[]){
  double pi=4*atan(1.);
  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    double vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));

    w[i]=cos(2*pi*(x[0]-vi*t));
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI]=0;
  w[_INDEX_EX]=0;
  w[_INDEX_RHO]=2.; //rho init
  w[_INDEX_VELOCITY]=0; // u init
  w[_INDEX_PRESSURE]=0; // p init
  w[_INDEX_TEMP]=0; // e ou T init

};

void TestPeriodic_InitData(double x[3],double w[]){

  double t=0;
  TestPeriodic_ImposedData(x,t,w);

};


void TestPeriodic_BoundaryFlux(double x[3],double t,double wL[],double* vnorm,
				       double* flux){
  double wR[_MV+6];
  TestPeriodic_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};



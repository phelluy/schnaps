#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"


int TestLandau_Damping_1D(void);

void Test_Landau_Damping_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[]);
void Test_Landau_Damping_InitData(schnaps_real x[3],schnaps_real w[]);
void Test_Landau_Damping_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm, schnaps_real* flux);

void UpdateVlasovPoisson(void* field);
void PlotVlasovPoisson(void* vf, schnaps_real * w);

int main(void) {
  
  // unit tests
    
  int resu=TestLandau_Damping_1D();
	 
  if (resu) printf("landau test OK !\n");
  else printf("landau test failed !\n");

  return !resu;
} 


int TestLandau_Damping_1D(void) {
  
  int test=0;
  int vec=1;
  schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  schnaps_real A[3][3] = {{2.0*pi/k, 0, 0}, {0, 1, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  // try to detect a 1d mesh
  Detect1DMacroMesh(&mesh);
  bool is1d=mesh.is1d;
  assert(is1d);

  // mesh preparation
  mesh.period[0]=2.0*pi/k;
  BuildConnectivity(&mesh);

  
  Model model;
  schnaps_real degV=2;
  schnaps_real nbEV=24;
  KineticData * kd=&schnaps_kinetic_data;

  InitKineticData(kd,nbEV,degV);
  
  model.m=kd->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.InitData =Test_Landau_Damping_InitData;
  model.ImposedData = Test_Landau_Damping_ImposedData;
  model.BoundaryFlux = Test_Landau_Damping_BoundaryFlux;
  model.Source = VlasovP_Lagrangian_Source;

  
  int deg[]={3, 0, 0};
  int raf[]={24, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = kd->vmax; // maximal wave speed
  simu.cfl=0.25;
  simu.nb_diags = 4;
  simu.pre_dtfields = UpdateVlasovPoisson;
  simu.post_dtfields=NULL;
  simu.update_after_rk = PlotVlasovPoisson;
 
  
  schnaps_real tmax = 30;
  RK2(&simu, tmax);

    // save the results and the error
  int iel = 2 * kd->nb_elem_v / 3;
  int iloc = kd->deg_v;
  printf("Trace vi=%f\n", -kd->vmax + iel * kd->dv + kd->dv * glop(kd->deg_v, iloc));
  PlotFields(iloc+iel*kd->deg_v,false,&simu,"sol f ","dgvisu.msh");
  PlotFields(kd->index_ex,false,&simu,"sol","dgvisuEx.msh");
  PlotFields(kd->index_phi,false,&simu,"sol","dgvisuPhi.msh");
  PlotFields(kd->index_rho,false,&simu,"sol","dgvisuRho.msh");
  Plot_Energies(&simu, simu.dt);
 

  test= 1;

  return test; 

}

void Test_Landau_Damping_ImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  //parameters of the case
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real k=0.5;
  schnaps_real eps = 0.001;
  schnaps_real my_pi= 4.0*atan(1.0);
  
  for(int i=0;i<kd->index_max_kin+1;i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax+nel*kd->dv +
		 kd->dv* glop(kd->deg_v,j));
 
    w[i]=(1.0+eps*cos(k*x[0]))*(1.0/sqrt(2.0*my_pi))*exp(-(vi*vi)/2.0);

  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi]=-(eps/(k*k))*cos(k*x[0]);
  w[kd->index_ex]=(eps/k)*sin(k*x[0]);
  w[kd->index_rho]=1; //rho init
  w[kd->index_u]=0; // u init
  w[kd->index_P]=0; // p init
  w[kd->index_T]=0; // e ou T init

};

void Test_Landau_Damping_InitData(schnaps_real x[3],schnaps_real w[]){

  schnaps_real t=0;
  Test_Landau_Damping_ImposedData(x,t,w);

};



void Test_Landau_Damping_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  Test_Landau_Damping_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
  assert(1==2);
};


void UpdateVlasovPoisson(void *si) {
  Simulation *simu = si;
  KineticData * kd=&schnaps_kinetic_data;
  int type_bc = 1;

  Computation_charge_density(simu);
  static ContinuousSolver ps;
  static bool is_init = false;

  if (!is_init){
    is_init = true;
    int nb_var=1;
    int * listvar= malloc(nb_var * sizeof(int));
    listvar[0]=kd->index_phi;
    InitContinuousSolver(&ps,simu,1,nb_var,listvar);
    
    ps.matrix_assembly=ContinuousOperator_Poisson1D;
    ps.rhs_assembly=RHSPoisson_Continuous;
    ps.bc_assembly=Periodic_BoundaryCondition_Poisson1D;
    ps.postcomputation_assembly=Computation_ElectricField_Poisson;
    
    ps.lsol.solver_type = LU;
    ps.lsol.pc_type=NONE;
  }
  
  SolveContinuous2D(&ps);
  //freeContinuousSolver(&ps);
}

void PlotVlasovPoisson(void *si, schnaps_real *w) {
  schnaps_real k_energy = 0, e_energy = 0, t_energy = 0, t_charge=0;
  
  Simulation *simu = si;
  
  Energies(simu, w, k_energy, e_energy, t_energy,1);
  Charge_total(simu,w,t_charge,4);
  si = simu; 
}








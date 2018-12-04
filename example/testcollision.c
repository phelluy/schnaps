#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"


int TestCollision(void);

void Equilibrium_VelocityPerturbation_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[]);
void Equilibrium_VelocityPerturbation_InitData(schnaps_real x[3],schnaps_real w[]);
void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm, schnaps_real* flux);
void Equilibrium_VelocityPerturbation_TotalSource(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, schnaps_real* source);

void Equilibrium_SpacePerturbation_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[]);
void Equilibrium_SpacePerturbation_InitData(schnaps_real x[3],schnaps_real w[]);
void Equilibrium_SpacePerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm, schnaps_real* flux);
void Equilibrium_SpacePerturbation_TotalSource(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, schnaps_real* source);

void Collision_VlasovPoissonI(void* field);
void Collision_VlasovPoissonII(void* field);
void PlotVlasovPoisson(void* vf, schnaps_real * w);

int main(void) {
  
  // unit tests
    
  int resu=TestCollision();
	 
  if (resu) printf("landau test OK !\n");
  else printf("landau test failed !\n");

  return !resu;
} 


int TestCollision(void) {
  
  int test=0;
  int vec=1;
  schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  // try to detect a 1d mesh
  Detect1DMacroMesh(&mesh);
  bool is1d=mesh.is1d;
  assert(is1d);

  // mesh preparation
  //mesh.period[0]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  schnaps_real degV=4;
  schnaps_real nbEV=32;
  KineticData * kd=&schnaps_kinetic_data;

  InitKineticData(kd,nbEV,degV);
  
  model.m=kd->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  
  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  //model.InitData =Equilibrium_SpacePerturbation_InitData;
  //model.ImposedData = Equilibrium_SpacePerturbation_ImposedData;
  //model.BoundaryFlux = Equilibrium_SpacePerturbation_BoundaryFlux;
  //model.Source = Equilibrium_SpacePerturbation_TotalSource;

  model.InitData =Equilibrium_VelocityPerturbation_InitData;
  model.ImposedData = Equilibrium_VelocityPerturbation_ImposedData;
  model.BoundaryFlux = Equilibrium_VelocityPerturbation_BoundaryFlux;
  model.Source = Equilibrium_VelocityPerturbation_TotalSource;

  
  int deg[]={4, 0, 0};
  int raf[]={32, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = kd->vmax; // maximal wave speed
  simu.cfl=0.5;
  simu.nb_diags = 4;
  simu.pre_dtfields = Collision_VlasovPoissonI;
  simu.update_after_rk = PlotVlasovPoisson;
 
  schnaps_real tmax = 0.1;
  if(kd->time_order == 1){
    RK1(&simu, tmax);
    simu.post_dtfields= NULL;
  }
  else {
    simu.post_dtfields= Collision_VlasovPoissonII;
    RK2(&simu, tmax);
  }
  

    // save the results and the error
  int iel = 2 * kd->nb_elem_v / 3;
  int iloc = kd->deg_v;
  printf("Trace vi=%f\n", -kd->vmax + iel * kd->dv + kd->dv * glop(kd->deg_v, iloc));
  PlotFields(iloc + iel * kd->deg_v, false, &simu, "sol","dgvisu_kin.msh");
  PlotFields(iloc + iel * kd->deg_v, true, &simu, "error","dgerror_kin.msh");
  
  Plot_Energies(&simu, simu.dt);

  schnaps_real dd_Kinetic = L2_Kinetic_error(&simu);
  schnaps_real ddrho = L2error_onefield(&simu,kd->index_rho);
  schnaps_real ddt = L2error_onefield(&simu,kd->index_T);
  schnaps_real dde = L2error_onefield(&simu,kd->index_ex);
  schnaps_real ddp = L2error_onefield(&simu,kd->index_phi);

  
  printf("erreur kinetic L2=%.5e, erreur Rho L2=%.5e, erreur T L2=%.5e,  erreur Ex L2=%.5e, erreur Phi L2=%.5e  \n", dd_Kinetic,ddrho,ddt,dde,ddp);
  test= test && (dd_Kinetic < 1e-2);


  test= 1;

  return test; 

}


void Equilibrium_SpacePerturbation_ImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real rho=4*my_pi*my_pi*(1.+0.5*(sin(2*my_pi*x[0])));
  
  for(int i=0;i<kd->index_max_kin+1;i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax+nel*kd->dv +
		 kd->dv* glop(kd->deg_v,j));
 
    w[i]=(rho/sqrt(2.0*my_pi))*exp(-(vi*vi)/2.0);
  }

  w[kd->index_phi]=0.5*sin(2.0*my_pi*x[0]);
  w[kd->index_ex]=-my_pi*cos(2.0*my_pi*x[0]);
  w[kd->index_ey]=0.0;
  w[kd->index_ez]=0.0;
  w[kd->index_rho]=rho;  //rho init
  w[kd->index_u]=0; // u init
  w[kd->index_P]=0.; // p init
  w[kd->index_T]=1.0; // e ou T init
};

void Equilibrium_VelocityPerturbation_ImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real my_pi= 4.0*atan(1.0), eps=0.01;
  
  for(int i=0;i<kd->index_max_kin+1;i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax+nel*kd->dv +
		 kd->dv* glop(kd->deg_v,j));
 
    w[i]=(1.0/sqrt(2.0*my_pi))*exp(-(vi*vi)/2.0) + eps*(kd->vmax*kd->vmax - vi*vi);
  }

  w[kd->index_phi]=0.0;
  w[kd->index_ex]=0.0;
  w[kd->index_ey]=0.0;
  w[kd->index_ez]=0.0;
  w[kd->index_rho]=1.0+(4./3.)*eps*pow(kd->vmax,3.0);  //rho init
  w[kd->index_u]=0; // u init
  w[kd->index_P]=0.5*(kd->gamma-1.)*(1.0+(4.0/15.0)*eps*pow(kd->vmax,5.0)); // p init
  w[kd->index_T]=(1.0+(4.0/15.0)*eps*pow(kd->vmax,5.0))/w[kd->index_rho]; // e ou T init
};


void Equilibrium_SpacePerturbation_TotalSource(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, 
			       schnaps_real* source) {
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real Transport_source[kd->index_max];
  schnaps_real M=0;
  schnaps_real my_pi= 4.0*atan(1.0);

  VlasovP_Lagrangian_Source(x,t,w,Transport_source);

   for(int iv=0;iv< kd->index_max_kin + 1;iv++){
      int j=iv%kd->deg_v; // local connectivity put in function
      int nel=iv/kd->deg_v; // element num (TODO : function)
   
      schnaps_real vn = (-kd->vmax+nel*kd->dv +
			 kd->dv* glop(kd->deg_v,j));
      
      M=(1.0/sqrt(2.0*my_pi))*exp(-(vn*vn)/2.0);	
      source[iv] = Transport_source[iv];
      source[iv] += 0.5*(vn+(1+0.5*sin(2.0*my_pi*x[0]))*vn)*8.0*pow(my_pi,3.0)*cos(2.0*my_pi*x[0])*M;
   }
  source[kd->index_phi]=0.0;
  source[kd->index_ex]=0.0;
  source[kd->index_ey]=0.0;
  source[kd->index_ez]=0.0;
  source[kd->index_rho]=0.0;  //rho init
  source[kd->index_u]=0; // u init
  source[kd->index_P]=0.; // p init
  source[kd->index_T]=0.0;
};

void Equilibrium_VelocityPerturbation_TotalSource(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, 
			       schnaps_real* source) {
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real Transport_source[kd->index_max];
  schnaps_real M=0,M_0=0,my_pi= 4.0*atan(1.0),eps=0.01;
  schnaps_real rho=1.0+(4./3.)*eps*pow(kd->vmax,3.0);  //rho init
  schnaps_real T=(1.0+(4.0/15.0)*eps*pow(kd->vmax,5.0))/rho; // e ou T init

  VlasovP_Lagrangian_Source(x,t,w,Transport_source);

   for(int iv=0;iv< kd->index_max_kin + 1;iv++){
      int j=iv%kd->deg_v; // local connectivity put in function
      int nel=iv/kd->deg_v; // element num (TODO : function)
   
      schnaps_real vn = (-kd->vmax+nel*kd->dv +
			 kd->dv* glop(kd->deg_v,j));
      
      M=(rho/sqrt(2.0*T*my_pi))*exp(-(vn*vn)/(2.0*T));
      M_0=(1.0/sqrt(2.0*my_pi))*exp(-(vn*vn)/2.0);
      source[iv] = Transport_source[iv]-kd->knud*(M - M_0 - eps*(kd->vmax*kd->vmax - vn*vn));
   }
  source[kd->index_phi]=0.0;
  source[kd->index_ex]=0.0;
  source[kd->index_ey]=0.0;
  source[kd->index_ez]=0.0;
  source[kd->index_rho]=0.0;  //rho init
  source[kd->index_u]=0; // u init
  source[kd->index_P]=0.; // p init
  source[kd->index_T]=0.0;
};


void Equilibrium_SpacePerturbation_InitData(schnaps_real x[3],schnaps_real w[]){
  schnaps_real t=0;
  Equilibrium_SpacePerturbation_ImposedData(x,t,w);
};

void Equilibrium_VelocityPerturbation_InitData(schnaps_real x[3],schnaps_real w[]){
  schnaps_real t=0;
  Equilibrium_VelocityPerturbation_ImposedData(x,t,w);
};


void Equilibrium_SpacePerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  Equilibrium_SpacePerturbation_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};

void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  Equilibrium_VelocityPerturbation_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};



void Collision_VlasovPoissonI(void *si) {
  Simulation *simu = si;
  KineticData * kd=&schnaps_kinetic_data;

  Collision_Source(simu);
  
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
    ps.bc_assembly=Periodic_BoundaryCondition_Poisson1D;//;
    ps.postcomputation_assembly=Computation_ElectricField_Poisson;
    
    ps.lsol.solver_type = LU;
    ps.lsol.pc_type=NONE;
  }
  
  SolveContinuous2D(&ps);
  //freeContinuousSolver(&ps);
}

void Collision_VlasovPoissonII(void *si) {
  Simulation *simu = si;
  KineticData * kd=&schnaps_kinetic_data;

  Collision_Source(simu);
  
}


void PlotVlasovPoisson(void *si, schnaps_real *w) {
  schnaps_real k_energy = 0, e_energy = 0, t_energy = 0, t_charge=0;
  
  Simulation *simu = si;
  
  Energies(simu, w, k_energy, e_energy, t_energy,1);
  Charge_total(simu,w,t_charge,4);
  si = simu; 
}








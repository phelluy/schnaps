#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "quantities_vp.h"
#include "gyro.h"
#include "solverpoisson.h"

int TestSliceMiniPoisson(void);
void SliceInitData(schnaps_real x[3],schnaps_real w[]);
void SliceImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[]);
void SliceBoundaryFlux(schnaps_real x[3],schnaps_real t,
		      schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux);

int main(void) {
  
  // unit tests
    
  int resu=TestSliceMiniPoisson();
	 
  if (resu) printf("slice poisson test OK !\n");
  else printf("slice poisson test failed !\n");

  return !resu;
} 


int TestSliceMiniPoisson(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh");

  mesh.period[2]=1;
  BuildConnectivity(&mesh);

  int vec=1;
  
    
  int deg[]={2, 2, 2};
  int raf[]={2, 2, 2};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  int nbelemv = 2;
  int deg_v = 2;
  
  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  kd->solve_quasineutrality = true;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=GyroUpwindNumFlux;
  //model.NumFlux=GyroZeroNumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=GyroBoundaryFlux;
  model.InitData=SliceInitData;
  model.ImposedData=SliceImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;


  InitSimulation(&simu, &mesh, deg, raf, &model);

  

  //simu.pre_dtfields = UpdateGyroPoisson;
   simu.vmax = kd->vmax; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0;
  schnaps_real tmax = 0;
  //RK4(&simu,tmax);
  //Computation_charge_density(&simu);
  // save the results and the error
  //PlotFields(1,(1==0),&simu,"sol","dgvisu.msh");

  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=kd->index_phi;

  int type_bc = 1;
  
  InitContinuousSolver(&ps,&simu,type_bc,nb_var,listvar);

  ps.matrix_assembly=ContinuousOperator_Poisson2D;
  ps.rhs_assembly=RHSPoisson_Continuous;
  ps.bc_assembly= ExactDirichletContinuousMatrix;
  ps.postcomputation_assembly=Computation_ElectricField_Poisson;
  //ps.postcomputation_assembly=NULL;

#undef PARALUTION
#ifdef PARALUTION
  ps.lsol.solver_type = PAR_LU;
  ps.lsol.pc_type=NONE;
#else
  ps.lsol.solver_type = LU;
  ps.lsol.pc_type=NONE;
#endif

  SolveContinuous2D(&ps);
  /* DisplayLinearSolver(&ps.lsol); */
  /* assert(1==2); */
  
  
  PlotFields(kd->index_phi,(1==0),&simu,"sol","dgvisu.msh");
  PlotFields(kd->index_phi,(1==1),&simu,"error","dgerror.msh");

  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&f);
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd < 0.03);


  return test; 

};


void SliceInitData(schnaps_real x[3],schnaps_real w[]){
  schnaps_real t=0;
  SliceImposedData(x,t,w);
}

void SliceImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData *kd = &schnaps_kinetic_data;
  for(int i = 0; i <kd->index_max_kin + 1; i++){
    w[i] = 0;
  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi]=(x[0] * x[0] + x[1] * x[1])/4;
  //w[kd->index_phi]= 0;
  w[kd->index_rho] = -1;
  w[kd->index_ex]=-x[0]/2;
  w[kd->index_ey]=-x[1]/2;
  //w[kd->index_ex]=0;
  //w[kd->index_ey]=0;
  w[kd->index_ez]=0;
  w[kd->index_u] = 0; // u init
  w[kd->index_P] = 0; // p init
  w[kd->index_T] = 0; // e ou T init
}

void SliceBoundaryFlux(schnaps_real x[3],schnaps_real t,
		      schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux)
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  SliceImposedData(x,t,wR);
  GyroUpwindNumFlux(wL,wR,vnorm,flux);
}


#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "global.h"
#include "lattice.h"
#include "implicit.h"
//#include "field.h"
//
typedef struct LbmSimuParams{
  int degx;
  int degy;
  int rafx;
  int rafy;
  schnaps_real cfl;
  schnaps_real tmax;
  schnaps_real cref;
  schnaps_real tau;
  schnaps_real diag_2d_period;
} LbmSimuParams;
typedef struct DoubleShearKHParams{
  schnaps_real kappa;
  schnaps_real delta;
  schnaps_real uref;
} DoubleShearKHParams;
typedef struct Linear2DWaveParams{
  int nkx;
  int nky;
  schnaps_real offset;
} Linear2DWaveParams;
//
int TestLattice_isothermal_DoubleShearKH(void);
void DoubleShearKH_InitData(schnaps_real x[3], schnaps_real w[]);
void DoubleshearKH_Plot_Fields(void *s,schnaps_real *w);
//
int TestLattice_isothermal_Linear2DWave(void);
int TestLattice_isothermal_Linear2DWaveImplicit(void);
void Linear2DWave_InitData(schnaps_real x[3],schnaps_real w[]);
void Linear2DWave_ImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]);
void Linear2DWave_ImposedData_OneNode(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]);
void Linear2DWave_Periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
               schnaps_real *flux);
void Linear2DWave_Periodic_BoundaryFlux_OneNode(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
               schnaps_real *flux);
void Linear2DWave_Plot_Fields(void *s,schnaps_real *w);
void Linear2DWave_CollectDiags(void *simu,schnaps_real *diag_vals);
//
void Relaxation(void* s);
void Moments(void * s);
// global parameters with default values
LbmSimuParams SimParams={
  .degx=4,.degy=4,
  .rafx=8,.rafy=8,
  .cfl=1.0,.tmax=1.0,
  .cref=1.0,.tau=0.00001,
  .diag_2d_period=1.0};
DoubleShearKHParams DKHParams={.kappa=80.0,.delta=0.05,.uref=0.05};
Linear2DWaveParams  LW2DParams={.nkx=1,.nky=0,.offset=0.0};
//
// global dat for taggig diags (this should be in a structure somewhere simu ?)
char simutag[4]="TAG";
//
int main(void) {
#define LBMDOUBLESHEARKH 0
#define LBMWAVE2D 1
//
#define LBMTESTCASE 1
#if (LBMTESTCASE==LBMWAVE2D)
  printf(" LBM - 2D - Isothermal - Linear2DWave\n");
  //
  SimParams.degx=4;
  SimParams.degy=4;
  SimParams.rafx=4;
  SimParams.rafy=4;
  
  SimParams.cfl=1.0;
  SimParams.tmax=0.001;
  SimParams.tau=0.0000001;
  SimParams.diag_2d_period=0.1;
  //
  LW2DParams.nkx=1;
  LW2DParams.nky=1;
  LW2DParams.offset=0.0;
  // 
  //int resu= TestLattice_isothermal_Linear2DWave();
  int resu= TestLattice_isothermal_Linear2DWaveImplicit();
#endif
#if (LBMTESTCASE==LBMDOUBLESHEARKH)
  printf(" LBM - 2D - Isothermal - Double sheared flow  KH\n"); 
  int resu=TestLattice_isothermal_DoubleShearKH();
#endif
  if (resu) printf("lattice test OK !\n");
  else printf("lattice test failed !\n");
  return !resu;
} 
//
// Basic routines common to all tests
void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  
  LatticeData * ld=&schnaps_lattice_data;

};


void Relaxation(void* s){
  Simulation * simu =s;
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real w_eq[simu->wsize];
  Compute_distribution_eq(simu,w_eq);
  Compute_relaxation(simu,w_eq);
}

void Moments(void* s){
  Simulation * simu =s;
  
  Compute_moments(simu);
}
/* ************************************************************************** */
/* KH Double shear test routines*/
/* ************************************************************************** */
int TestLattice_isothermal_DoubleShearKH(void) {
  
  int test=0;
  //int vec=1;
  //schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_D2Q9;
  schnaps_real cref= SimParams.cref;
  InitLatticeData(ld,2,9,0,cref);
  ld->diag_2d_period=SimParams.diag_2d_period;
  ld->tau=SimParams.tau;
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  model.InitData = DoubleShearKH_InitData;
  model.ImposedData = NULL;
  model.BoundaryFlux = NULL;
  model.Source = NULL;
  
  int deg[]={SimParams.degx, SimParams.degy, 0};
  int raf[]={SimParams.rafx,SimParams.rafy, 1};
  //
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  //
  InitSimulation(&simu, &mesh, deg, raf, &model);
  Moments(&simu);
  //simu.vmax = 2*ld->c; 
  simu.vmax = ld->c *sqrt(2.0); 
  simu.cfl=SimParams.cfl;
  simu.nb_diags = 1;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk = DoubleshearKH_Plot_Fields;
  schnaps_real tmax = SimParams.tmax;
  //
  DoubleshearKH_Plot_Fields(&simu, NULL);
  //
  RK2(&simu, tmax);
  //
  test= 1;
  return test; 
}
/*KH Double shear Init*/
void DoubleShearKH_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real delta = DKHParams.delta, kappa=DKHParams.kappa,uref=DKHParams.uref;
  //
  schnaps_real rho = 1.0;
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld->c)/3.0;
  schnaps_real p =1.0;
  //
  //printf(" c0 :%f \t temp:%f\n",ld->c, temp);
  if(x[1]<0.5){
    ux=uref * tanh(kappa*(x[1]-0.25));
  }
  else{
    ux=uref * tanh(kappa*(0.75-x[1]));
  }    
  uy = uref * delta * sin(2.0 * my_pi*(x[0]+0.25));
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
};
/* ********* KH Double shear data dumps ********* */
void DoubleshearKH_Plot_Fields(void *s,schnaps_real *w){
  Simulation *simu=s;
  LatticeData *ld=&schnaps_lattice_data;
  schnaps_real period=ld->diag_2d_period;
  schnaps_real dt=simu->dt;
  int diagperiod = (int) (period/dt);
  int istep=simu->iter_time_rk;
  schnaps_real t=simu->tnow;
  int tmax=simu->tmax;
  int create_file=0;
  if (istep==0){
    create_file=1;
  }
  else
  {
  create_file = 0;
  }
  if (diagperiod || create_file){
  if ((istep%diagperiod ==0) || (t== tmax)){
    printf("Dumping fields at it=%i\n",istep);
    PlotFieldsBinSparseMultitime(ld->index_rho,false,simu,"rho","lbm_KH2shear_rho.msh",create_file,t,istep);
    PlotFieldsBinSparseMultitime(ld->index_ux,false,simu,"ux","lbm_KH2shear_ux.msh",create_file,t,istep);
    PlotFieldsBinSparseMultitime(ld->index_uy,false,simu,"uy","lbm_KH2shear_uy.msh",create_file,t,istep);
    Compute_and_dump_Vorticity_2d(simu,"lbm_KH2shear_uvort.msh",create_file,t,istep);
  };
  };
};
/* ****************************************************************************/
/* **************************** 2D Wave equation linear LBM *******************/
int TestLattice_isothermal_Linear2DWave(void) {
  
  int test=0;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);
  //
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_linearwave_D2Q9;
  ld->collect_diags=&Linear2DWave_CollectDiags;
  //
  //schnaps_real csound= sqrt(3.0/2.0);
  InitLatticeData(ld,2,9,0,SimParams.cref);
  ld->diag_2d_period=SimParams.diag_2d_period;
  ld->tau=SimParams.tau;
  //  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  model.InitData = Linear2DWave_InitData;
  model.ImposedData = Linear2DWave_ImposedData;
  //model.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux;
  model.BoundaryFlux =NULL;
  model.Source = NULL;
  //
  int deg[]={SimParams.degx,SimParams.degy, 0};
  int raf[]={SimParams.rafx,SimParams.rafy,1};
  //
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  //
  InitSimulation(&simu, &mesh, deg, raf, &model);
  Moments(&simu);
  simu.vmax = ld->c *sqrt(2.0); 
  simu.cfl=SimParams.cfl;
  simu.nb_diags = 3;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk = Linear2DWave_Plot_Fields;
  schnaps_real tmax = SimParams.tmax;
  //
  //Linear2DWave_Plot_Fields(&simu,NULL);
  sprintf(simutag,"RK2");
  RK2(&simu, tmax);
  //
  schnaps_real end_diags[simu.nb_diags];
  Linear2DWave_CollectDiags(&simu,end_diags);
  printf("RK2 end error:t=%f \tabs : %f \tmean : %f \trel:%f\n",simu.tnow,end_diags[0],end_diags[1],end_diags[2]);
  //
  Dump_Lattice_Diagnostics(&simu,simutag);
  test= 1;
  return test; 
}
int TestLattice_isothermal_Linear2DWaveImplicit(void) {
  
  int test=0;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);
  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);
  //
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_linearwave_D2Q9;
  ld->collect_diags=&Linear2DWave_CollectDiags;
  //
  //schnaps_real csound= sqrt(3.0/2.0);
  InitLatticeData(ld,2,9,0,SimParams.cref);
  ld->diag_2d_period=SimParams.diag_2d_period;
  ld->tau=SimParams.tau;
  //  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  model.InitData = Linear2DWave_InitData;
  model.ImposedData = Linear2DWave_ImposedData;
  model.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux;
  model.Source = NULL;
  //
  int deg[]={SimParams.degx,SimParams.degy, 0};
  int raf[]={SimParams.rafx,SimParams.rafy,1};
  //
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
  //
  InitSimulation(&simu, &mesh, deg, raf, &model);
  Moments(&simu);
  simu.vmax = ld->c *sqrt(2.0); 
  simu.cfl=SimParams.cfl;
  simu.nb_diags = 3;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk=Linear2DWave_Plot_Fields;
  schnaps_real tmax = SimParams.tmax;
  //
  // basic explicit RK2
  sprintf(simutag,"RK2");
  RK2(&simu, tmax);
  schnaps_real end_diags[simu.nb_diags];
  Linear2DWave_CollectDiags(&simu,end_diags);
  printf("RK2 end error: abs : %f \tmean : %f \trel:%f\n",end_diags[0],end_diags[1],end_diags[2]);
  Dump_Lattice_Diagnostics(&simu,simutag);
  //
  Simulation simu2; // implicit scheme simulation
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  Moments(&simu2);
  simu2.vmax = ld->c *sqrt(2.0); 
  simu2.cfl=SimParams.cfl; // no impact we impose dt
  simu2.nb_diags = 3;
  simu2.pre_dtfields = Relaxation;
  simu2.post_dtfields = Moments;
  simu2.update_after_rk=Linear2DWave_Plot_Fields;
  //simu2.update_after_rk=NULL;
  tmax = SimParams.tmax;
  schnaps_real dt= simu.dt; // we use the timestep of previous RK2 run
  printf("Testing Implicit scheme with dt=%f\n",dt);
  //LatticeThetaTimeScheme_MultiSim(&simu2,simu_advection,tmax,dt);
  sprintf(simutag,"IMP");
  Model model_advec; // model for advection of a single velocity node
  model_advec.m=1;
  model_advec.InitData =Lattice_Dummy_InitData;
  model_advec.NumFlux=Lattice_OneNodeNumFlux;
  model_advec.ImposedData = Linear2DWave_ImposedData_OneNode;
  model_advec.BoundaryFlux = Linear2DWave_Periodic_BoundaryFlux_OneNode;
  model_advec.Source = NULL;
  LatticeThetaTimeScheme(&simu2,&model_advec,tmax,dt);
  /* ************************************************************* */
  // test implicit
  schnaps_real end_diags2[simu2.nb_diags];
  Linear2DWave_CollectDiags(&simu2,end_diags2);
  //
  printf("Summary\n");
  printf("RK2  :t=%f \t abs : %f \tmean : %f \trel:%f\n",simu.tnow,end_diags[0],end_diags[1],end_diags[2]);
  printf("Imp  :t=%f \t abs : %f \tmean : %f \trel:%f\n",simu2.tnow,end_diags2[0],end_diags2[1],end_diags2[2]);
  //
  Dump_Lattice_Diagnostics(&simu2,simutag);
  //
  test= 1;
  return test; 
}
void Linear2DWave_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  // wave mode numbers in half integer units
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset=LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky= 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix= kx * x[0];
  schnaps_real phiy= ky * x[1];
  //
  schnaps_real rho = offset+cos(phix) * cos(phiy);
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
}
void Linear2DWave_ImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]){
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset=LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky= 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix= kx * x[0];
  schnaps_real phiy= ky * x[1];
  //
  schnaps_real k = sqrt(kx * kx + ky * ky);
  schnaps_real c0 = ld->c * sqrt(1/3.0);
  schnaps_real omega= k* c0;
  schnaps_real phit= omega *t;
  //
  //schnaps_real rho = 1.0 + pert * cos(phix);
  schnaps_real rho = offset+cos(phix) * cos(phiy) * cos(phit);
  //schnaps_real ux =uref * pert * sin(phiy);
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
}
void Linear2DWave_ImposedData_OneNode(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]){
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset=LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky= 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix= kx * x[0];
  schnaps_real phiy= ky * x[1];
  //
  schnaps_real k = sqrt(kx * kx + ky * ky);
  schnaps_real c0 = ld->c * sqrt(1/3.0);
  schnaps_real omega= k* c0;
  schnaps_real phit= omega *t;
  //
  //schnaps_real rho = 1.0 + pert * cos(phix);
  schnaps_real rho = offset+cos(phix) * cos(phiy) * cos(phit);
  //schnaps_real ux =uref * pert * sin(phiy);
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  int i_node= ld->current_node_index;
  w[0]= ld->feq(i_node,ld,rho,ux,uy,uz,temp,p);
}
/* ************************************************************************** */
void Linear2DWave_Periodic_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
	LatticeData *ld=&schnaps_lattice_data;
  schnaps_real wR[ld->index_max];
  Linear2DWave_ImposedData(x,t,wR);
  Lattice_NumFlux(wL,wR,vnorm,flux);
}
void Linear2DWave_Periodic_BoundaryFlux_OneNode(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux){ 
	LatticeData *ld=&schnaps_lattice_data;
	int i_node=ld->current_node_index;
  schnaps_real wR[1];
  Linear2DWave_ImposedData_OneNode(x,t,wR);
  Lattice_OneNodeNumFlux(wL,wR,vnorm,flux);
}
/* ************************************************************************** */
void Linear2DWave_Plot_Fields(void *s,schnaps_real *w){
  Simulation *simu=s;
  LatticeData *ld=&schnaps_lattice_data;
  schnaps_real period=ld->diag_2d_period;
  schnaps_real dt=simu->dt;
  int diagperiod = (int) (period/dt);
  int istep=simu->iter_time_rk;
  schnaps_real t=simu->tnow;
  int tmax=simu->tmax;
  int create_file=0;
  Store_Lattice_diags(simu);
  if (istep==0){
    create_file=1;
  }
  else
  {
  create_file = 0;
  }
  if (diagperiod || create_file){
  if ((istep%diagperiod ==0) || (t== tmax)){
    istep=istep+1;
    printf("Dumping fields at it=%i (period %i)\n",istep,diagperiod);
    int raf=simu->fd[0].raf[0];
    schnaps_real cfl=simu->cfl;
    char filename_rho[sizeof("lbm2DWave_rho_TAG_raf000_cfl0.000.msh")];
    sprintf(filename_rho,"lbm_2DWave_rho_%s_raf%03d_cfl%1.3f.msh",simutag,raf,cfl);
    //char filename_rho_error[sizeof("lbm2DWave_rho_000.msh")];
    //sprintf(filename_rho_error,"lbm_2DWave_rho_error_%03d.msh",raf);
    PlotFieldsBinSparseMultitime(ld->index_rho,false,simu,"rho",filename_rho,create_file,t,istep);
    PlotFieldsBinSparseMultitime(ld->index_rho,true,simu,"rho_error",filename_rho,0,t,istep);
  };
  }
};
/* ************************************************************************** */
void Linear2DWave_CollectDiags(void *s,schnaps_real *diag_vals){
  Simulation *simu=s;
  schnaps_real error = 0;
  schnaps_real mean = 0;
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field *f = simu->fd + ie;
    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      schnaps_real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      w[iv] = f->wn[imem];
      }
      schnaps_real wex[f->model.m];
      schnaps_real wpg, det;
      // Compute wpg, det, and the exact solution
      schnaps_real xphy[3], xpgref[3];
      schnaps_real dtau[3][3], codtau[3][3];
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
      xpgref, // xref
      NULL, -1, // dpsiref, ifa
      xphy, dtau, // xphy, dtau
      codtau, NULL, NULL); // codtau, dpsi, vnds
      det = dot_product(dtau[0], codtau[0]);
      // Get the exact value
      f->model.ImposedData(xphy, simu->tnow, wex);
    	schnaps_real diff = w[0] - wex[0];
      error += diff * diff * wpg * det;
      mean += wex[0] * wex[0] * wpg * det;
    };
    };
  //printf("errl2=%f\n",sqrt(error) / (sqrt(mean)  + 1e-14));
  diag_vals[0]=sqrt(error);
  diag_vals[1]=sqrt(mean);
  diag_vals[2]=sqrt(error) / (sqrt(mean)  + 1e-14);
}

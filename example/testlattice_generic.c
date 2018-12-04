#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "global.h"
#include "lbm_generic.h"
#include "lbm_timeschemes.h"
#include "lbm_diagnostics.h"
//
typedef struct LbmSimuParams {
  int deg[3];
  int raf[3];
  schnaps_real cfl;
  schnaps_real dt;
  schnaps_real tmax;
  schnaps_real cref;
  schnaps_real tau;
  schnaps_real diag_2d_period;
} LbmSimuParams;
//! \brief parameters for the bump in a channel 2D flow
typedef struct BumpFlowParams {
  schnaps_real uref; //! flow velocity 
  schnaps_real ly; //! channel height (width)
  schnaps_real lx; //! channel length
  schnaps_real cx; //! x position  of center of circular bump
  schnaps_real cy; //! y position of center of circular bump
  schnaps_real r; //! radius of circular bump.
} BumpFlowParams;
typedef struct DoubleShearKHParams {
  schnaps_real kappa;
  schnaps_real delta;
  schnaps_real uref;
} DoubleShearKHParams;
typedef struct Linear2DWaveParams {
  int nkx;
  int nky;
  schnaps_real offset;
} Linear2DWaveParams;
//
int test_reverse_varindex(Simulation *simu);
//
int LBM_testmodels(void);
//
int TestLattice_BumpFlow(void);
void LBM_BumpFlow_InitData(schnaps_real x[3], schnaps_real w[]);
void LBM_BumpFlow_ImposedData(const schnaps_real x[3],
			      const schnaps_real t, schnaps_real w[]);
void LBM_BumpFlow_ImposedData_OneNode(const schnaps_real x[3],
				      const schnaps_real t,
				      schnaps_real w[]);
void LBM_BumpFlow_Plot_Fields(void *s, schnaps_real * w);
void LBM_BumpFlow_BoundaryFlux_OneNode(schnaps_real * x, schnaps_real t,
				       schnaps_real * wL,
				       schnaps_real * vnorm,
				       schnaps_real * flux);
//
int LBM_TestLattice_LinearWave2D(void);
void LBM_Linear2DWave_InitData(schnaps_real x[3], schnaps_real w[]);
void LBM_Linear2DWave_ImposedData(const schnaps_real x[3],
				  const schnaps_real t, schnaps_real w[]);
void LBM_Linear2DWave_ImposedData_OneNode(const schnaps_real x[3],
					  const schnaps_real t,
					  schnaps_real w[]);
void LBM_Linear2DWave_Periodic_BoundaryFlux_OneNode(schnaps_real * x,
						    schnaps_real t,
						    schnaps_real * wL,
						    schnaps_real * vnorm,
						    schnaps_real * flux);
void LBM_Linear2DWave_CollectDiags(void *s, schnaps_real * macro_diag_vals,
				   schnaps_real * micro_diag_vals);
void LBM_Linear2DWave_Plot_Fields(void *s, schnaps_real * w);
// test of local implicit newton scheme
int LBM_TestLattice_SteadyFlow_NLImplicit(schnaps_real alpha_deg);
void LBM_SteadyFlow_NLImplicit_MacInitData(schnaps_real x[3], schnaps_real w[]);
void LBM_SteadyFLow_NLImplicit_MacImposedData(const schnaps_real x[3],
			      const schnaps_real t, schnaps_real w[]);
void LBM_SteadyFlow_NLImplicit_MicInitData(schnaps_real x[3], schnaps_real w[]);
void LBM_SteadyFLow_NLImplicit_MicImposedData(const schnaps_real x[3],
			      const schnaps_real t, schnaps_real w[]);
void LBM_SteadyFlow_NLImplicit_BoundaryFlux(schnaps_real * x,
						    schnaps_real t,
						    schnaps_real * wL,
						    schnaps_real * vnorm,
						    schnaps_real * flux);
//
// global parameters with default values
LbmSimuParams SimParams = {
  .deg = {3, 3, 0},
  .raf = {3, 3, 1},
  .cfl = 1.0,.dt = 0.001,.tmax = 0.1,
  .tau = 0.00001,.cref = 1.0,
  .diag_2d_period = 1.0
};
DoubleShearKHParams DKHParams = {.kappa = 80.0,.delta = 0.05,.uref =
      0.05 };
Linear2DWaveParams LW2DParams = {.nkx = 1,.nky = 0,.offset = 0.0 };
BumpFlowParams BFParams = {.uref = 0.01,.lx = 6.0,.ly = 2.0,.cx = 3.0,.cy =
      -0.4,.r = 0.5 };
//
char simutag[4] = "TAG";
//
int main(void)
{
  //
  printf(" Lattice Boltzmann Model\n");
  SimParams = (LbmSimuParams) {
    .deg = {1, 1, 0},
    .raf = {2, 2, 1},
    .cfl = 1.0,
    .dt = 0.001,
    .tmax = 0.0005,
    .cref = 1.0,
    .tau = 0.01,
    .diag_2d_period = 0.5};
  // ***********************//
/*  int resu=LBM_testmodels();*/
/*  assert(1==2);*/
  // ***********************//
/*  BFParams.uref = 0.05;*/
/*  BFParams.lx = 4.0;*/
/*  BFParams.ly = 1.0;*/
/*  BFParams.cx = 2.0;*/
/*  BFParams.cy = -0.9;*/
/*  BFParams.r = 1.0;*/
/*  int resu = TestLattice_BumpFlow();*/
  //
  LW2DParams.nkx=1;
  LW2DParams.nky=1;
  LW2DParams.offset=0.0;
  int resu=LBM_TestLattice_LinearWave2D();
  //
/*  LW2DParams.nkx=1;*/
/*  LW2DParams.nky=1;*/
/*  LW2DParams.offset=0.0;*/
/*  int resu=LBM_TestLattice_LinearWave2D_NLImplicit();*/
  //
  //int resu = LBM_TestLattice_SteadyFlow_NLImplicit(0.0);
  if (resu)
    printf("lattice test OK !\n");
  else
    printf("lattice test failed !\n");
  return !resu;
}

//
int LBM_testmodels(void)
{
  int d=0;
  int nb_macro =0;
  int q = 0;
  LBModelDescriptor lbm = LBModelDescriptor_NULL;
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  d = 2;
  nb_macro = 2;
  q = 3;
  schnaps_real alpha_deg= 45.0;
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  printf(" D2Q3 TEST model for alpha=%f\n",alpha_deg);
  lsd->lb_model = &lbm;
  LBM_Set_D2Q3_TEST_model(&lbm, alpha_deg,SimParams.cref);
  DisplayLBModelDescriptorMomentPoly(&lbm);
  CheckLBModelDescriptorMacroConservation(&lbm, false);
  CheckLBMMomentMatrixInversion(&lbm,false);
  DisplayLBModelDescriptorMomentMatrix(&lbm);
  DestroyLBModelDescriptor(&lbm);
  //
  d = 2;
  nb_macro = 3;
  q = 9;
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  printf(" D2Q9 isothermal model\n");
  lsd->lb_model = &lbm;
  LBM_Set_D2Q9_ISOTH_model(&lbm, SimParams.cref);
  CheckLBModelDescriptorMacroConservation(&lbm, false);
  CheckLBMMomentMatrixInversion(&lbm,false);
  DisplayLBModelDescriptorMomentMatrix(&lbm);
  DestroyLBModelDescriptor(&lbm);
  //
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  printf(" D2Q9 isothermal INC model\n");
  lsd->lb_model = &lbm;
  LBM_Set_D2Q9_ISOTH_INC_model(&lbm, SimParams.cref);
  CheckLBModelDescriptorMacroConservation(&lbm, false);
  CheckLBMMomentMatrixInversion(&lbm,false);
  DisplayLBModelDescriptorMomentMatrix(&lbm);
  DestroyLBModelDescriptor(&lbm);
  //
  printf(" D2Q9 isothermal Linerarized model\n");
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  lsd->lb_model = &lbm;
  LBM_Set_D2Q9_ISOTH_LINEARIZED_model(&lbm, SimParams.cref);
  CheckLBModelDescriptorMacroConservation(&lbm, false);
  CheckLBMMomentMatrixInversion(&lbm,false);
  DisplayLBModelDescriptorMomentMatrix(&lbm);
  DestroyLBModelDescriptor(&lbm);
  //
  d=2;
  nb_macro=5;
  q=19;
  printf(" MDH D2Q9 ISOTH D2Q5 ISOTH model\n");
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  lsd->lb_model = &lbm;
  LBM_Set_MHD_D2Q9_2D2Q5_model(&lbm, SimParams.cref);
  //
  CheckLBModelDescriptorMacroConservation(&lbm, true);
  //DisplayLBModelDescriptorMomentMatrix(&lbm);
  DestroyLBModelDescriptor(&lbm);
  return 1;
}

//
int TestLattice_BumpFlow(void)
{
  int d = 2;
  int nb_macro = 3;
  int q = 9;
  LBModelDescriptor lbm = LBModelDescriptor_NULL;
  //
  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../geo/cbump2d.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  schnaps_real x0[3] = { 0, 0, 0 };
  AffineMapMacroMesh(&mesh, A, x0);
  // mesh preparation
  mesh.period[0] = 0.0;
  mesh.period[1] = 0.0;
  BuildConnectivity(&mesh);
  //
  CheckMacroMesh(&mesh, SimParams.deg, SimParams.raf);
  //
  // setup simulation paramaters in global shared LatticeBoltzmannSimData object
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  //LBM_Set_D2Q9_ISOTH_model(&lbm, SimParams.cref);
  LBM_Set_D2Q9_ISOTH_INC_model(&lbm, SimParams.cref);
  lsd->lb_model = &lbm;
  // setup LB Simulation object
  LBMSimulation lbsimu;
  lbsimu.macro_model.InitData = LBM_BumpFlow_InitData;
  lbsimu.macro_model.ImposedData = NULL;
  InitLBMSimulation(&lbsimu, lsd, &mesh, SimParams.deg, SimParams.raf);
  //
  lbsimu.vmax = lsd->lb_model->vmax;
  lbsimu.cfl = SimParams.cfl;
  lbsimu.micro_simu.cfl = lbsimu.cfl;
  lbsimu.dt = SimParams.dt;
  lbsimu.tmax = SimParams.tmax;
  lbsimu.micro_simu.tmax = lbsimu.tmax;
  lbsimu.micro_simu.cfl = lbsimu.cfl;
  lbsimu.macro_simu.cfl = lbsimu.cfl;
  lbsimu.micro_simu.tmax = lbsimu.tmax;
  lbsimu.macro_simu.tmax = lbsimu.tmax;
  //
  schnaps_real tau = SimParams.tau;
  lbm.s[0] = lbsimu.dt / (tau + 0.5 * lbsimu.dt);
  printf(" tau=%f s=%f\n", tau, lbm.s[0]);
  //
  lbsimu.macro_simu.nb_diags = 1;
  lbsimu.micro_simu.nb_diags = 0;
  lbsimu.pre_advec = LB_Relaxation_bgk_f;
  lbsimu.post_advec_one_node=NULL;
  lbsimu.post_advec = LB_ComputeMacroFromMicro;
  lbsimu.post_tstep = LBM_BumpFlow_Plot_Fields;
  //lbsimu.post_tstep=NULL;
  lbsimu.collect_diags = NULL;
  lbsimu.diag_2d_period = SimParams.diag_2d_period;
  //
  sprintf(simutag, "IMP");
  lbsimu.model_advec.m = 1;
  lbsimu.model_advec.InitData = LBM_Dummy_InitData_OneNode;
  lbsimu.model_advec.ImposedData = NULL;
  lbsimu.model_advec.NumFlux = LBM_OneNodeNumFlux;
  lbsimu.model_advec.BoundaryFlux = LBM_BumpFlow_BoundaryFlux_OneNode;
  lbsimu.model_advec.Source = NULL;
  //
  //
  schnaps_real hmin=lbsimu.micro_simu.hmin;
  schnaps_real vmax=lbsimu.vmax;
  printf(" hmin:%f \t vmax:%f \t dt:%f \t vmax*dt :%f\n",hmin,vmax,lbsimu.dt,vmax *lbsimu.dt);
  // Actual run
  LBMThetaTimeScheme(&lbsimu, 0.5, lbsimu.tmax, lbsimu.dt);
  //
  LBM_Dump_Lattice_Diagnostics(&lbsimu, "SIM");
  // cleanup
  FreeBMSimulation(&lbsimu);
  lsd->lb_model = NULL;
  DestroyLBModelDescriptor(&lbm);
  return 1;
}

void LBM_BumpFlow_InitData(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = 1.0;
  schnaps_real ux = 1.0;
  schnaps_real uy = 0.0;
  schnaps_real uref = BFParams.uref;
  //
  w[0] = rho;
  w[1] = ux * uref;
  w[2] = uy * uref;
}

void LBM_BumpFlow_ImposedData(const schnaps_real x[3],
			      const schnaps_real t, schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  schnaps_real rho = 1.0;
  schnaps_real ux = 0.0;
  schnaps_real uy = 0.0;
  schnaps_real uref = BFParams.uref;
  //
  if (x[0] < _SMALL) {
    ux = 1.0;
    rho = 1.0;
  }
/*  if (x[1] < _SMALL){*/
/*    ux=0.0;*/
/*  }*/
  w[0] = rho;
  w[1] = ux * uref;
  w[2] = uy * uref;
}

void LBM_BumpFlow_ImposedData_OneNode(const schnaps_real x[3],
				      const schnaps_real t,
				      schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = 1.0;
  schnaps_real ux = 1.0;
  schnaps_real uy = 0.0;
  //
/*  if (x[1] < _SMALL){*/
/*    ux=0.0;*/
/*  }*/
  schnaps_real wmac[lsd->lb_model->nb_macro];
  wmac[0] = rho;
  wmac[1] = ux * BFParams.uref;
  wmac[2] = uy * BFParams.uref;
  //
  int inode = lsd->current_node_index;
  w[0] = lsd->lb_model->feq_i(inode, lsd->lb_model->nb_macro, wmac);
}

void LBM_BumpFlow_BoundaryFlux_OneNode(schnaps_real * x, schnaps_real t,
				       schnaps_real * wL,
				       schnaps_real * vnorm,
				       schnaps_real * flux)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int i_node = lsd->current_node_index;
  //
  schnaps_real wR[1];
  // infinity state everywhere
/*    LBM_BumpFlow_ImposedData_OneNode(x, t, wR);*/
/*    LBM_OneNodeNumFlux(wL, wR, vnorm, flux);*/
/*    return;*/
  //
  // left flow inlet
  if (x[0] == 0) {
    LBM_BumpFlow_ImposedData_OneNode(x, t, wR);
    LBM_OneNodeNumFlux(wL, wR, vnorm, flux);
    return;
  }
/*  // right flow outlet infinity state*/
  if (x[0] == BFParams.lx) {
    LBM_BumpFlow_ImposedData_OneNode(x, t, wR);
    LBM_OneNodeNumFlux(wL, wR, vnorm, flux);
    return;
  }
  // right flow outlet=> continuity ?
  //if ((x[0] == BFParams.lx) && ((x[1] > 0.0) || (x[1] < BFParams.ly))){
/*  if (x[0] == BFParams.lx){*/
/*    wR[0] = wL[0];*/
/*    LBM_OneNodeNumFlux(wL, wR, vnorm, flux);*/
/*    return;*/
/*  }*/
  // corners ?
  //
  // uppper wall infinity state
  if (x[1] == BFParams.ly) {
    LBM_BumpFlow_ImposedData_OneNode(x, t, wR);
    LBM_OneNodeNumFlux(wL, wR, vnorm, flux);
    return;
  }
  schnaps_real ymax=BFParams.cy+BFParams.r;
  //if (x[1] <= ymax){
    // macro mirroring slip
/*  schnaps_real wmac[3];*/
/*  schnaps_real rho=lsd->current_lb_sim->wmac_buffer[0];*/
/*  schnaps_real ux= lsd->current_lb_sim->wmac_buffer[1];*/
/*  schnaps_real uy= lsd->current_lb_sim->wmac_buffer[2];*/
/*  schnaps_real udotn= ux *vnorm[0] + uy *vnorm[1];*/
/*  wmac[0]=rho;*/
/*  wmac[1]=ux - 2.0 * udotn *vnorm[0];*/
/*  wmac[2]=uy - 2.0 * udotn *vnorm[1];*/
/*  wR[0]=lsd->lb_model->feq(i_node, lsd->lb_model->nb_macro, wmac);*/
/*  LBM_OneNodeNumFlux(wL, wR, vnorm, flux);*/
  //
  // noslip on f
  int iopp = lsd->lb_model->iopposite[i_node];
  wR[0] = lsd->current_lb_sim->wmic_buffer[iopp];
  LBM_OneNodeNumFlux(wL, wR, vnorm, flux);
  //}
  //
  // test slip with full moments inversion
/*  schnaps_real Mom[lsd->lb_model->q];*/
/*  MatVect(lsd->lb_model->Msolv,lsd->current_lb_sim->wmic_buffer,Mom); // compute full Moments set*/
/*  schnaps_real jx = Mom[1];*/
/*  schnaps_real jy=  Mom[2];*/
/*  schnaps_real jdotn= jx *vnorm[0] + jy *vnorm[1];*/
/*  Mom[1] = jx - 2.0 * jdotn * vnorm[0];*/
/*  Mom[2] = jy - 2.0 * jdotn * vnorm[1];*/
/*  for (int i=0;i< lsd->lb_model->q;i++){*/
/*  lsd->lb_model->Msolv->rhs[i]=Mom[i];*/
/*  }*/
/*  SolveLinearSolver(lsd->lb_model->Msolv);*/
/*  wR[0] = lsd->lb_model->Msolv->sol[i_node];*/

/*  LBM_OneNodeNumFlux(wL, wR, vnorm, flux);*/
}

//
void LBM_BumpFlow_Plot_Fields(void *s, schnaps_real * w)
{
  LBMSimulation *lbsimu = s;
  Simulation *simu = &(lbsimu->macro_simu);
  Simulation *micsimu = &(lbsimu->micro_simu);
  schnaps_real period = lbsimu->diag_2d_period;
  schnaps_real dt = lbsimu->dt;
  int diagperiod = (int) (period / dt);
  int istep = lbsimu->micro_simu.iter_time_rk;
  schnaps_real t = lbsimu->micro_simu.tnow;
  int tmax = simu->tmax;
  int create_file = 0;
  LBM_Store_Lattice_diags(lbsimu);
  if (istep == 0) {
    create_file = 1;
  } else {
    create_file = 0;
  }
  if (diagperiod || create_file) {
    if ((istep % diagperiod == 0) || (t == tmax)) {
      istep = istep + 1;
      printf("Dumping fields at it=%i (period %i)\n", istep, diagperiod);
      int raf = simu->fd[0].raf[0];
      schnaps_real dt = simu->dt;
      char filename[sizeof("lbm_BumpFlow_TAG_raf000_cfl0.000.msh")];
      sprintf(filename, "lbm_BumpFlow_%s_raf%03d_dt%1.3f.msh", "MAC",
	      raf, dt);
      LBM_PlotFieldsBinSparseMultitime(0, false, simu, "rho", filename,
				       create_file, t, istep);
/*      LBM_PlotFieldsBinSparseMultitime(1, false, simu, "ux", filename, 0,*/
/*				       t, istep);*/
/*      LBM_PlotFieldsBinSparseMultitime(2, false, simu, "uy", filename, 0,*/
/*				       t, istep);*/
      LBM_PlotVecFieldsBinSparseMultitime((int[3]){1,2,-1},false,simu,"u",filename,0,t,istep);
      //
      char mic_filename[sizeof("lbm_BumpFlow_TAG_raf000_cfl0.000.msh")];
      sprintf(mic_filename, "lbm_BumpFlow_%s_raf%03d_dt%1.3f.msh", "MIC",
	      raf, dt);
      LBM_PlotFieldsBinSparseMultitime(0, false, micsimu, "f0",
				       mic_filename, create_file, t,
				       istep);
      for (int i = 1; i < lbsimu->q; i++) {
	char fieldname[32];
	sprintf(fieldname, "f%i", i);
	LBM_PlotFieldsBinSparseMultitime(i, false, micsimu, fieldname,
					 mic_filename, create_file, t,
					 istep);
      }
      //
    };
  }
};

/*******************************************************************************************************************************/
/************************* linearized D2Q9 test case ***************************************************************************/
int LBM_TestLattice_LinearWave2D(void)
{
  int d = 2;
  int nb_macro = 3;
  int q = 9;
  LBModelDescriptor lbm = LBModelDescriptor_NULL;
  //
  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  schnaps_real x0[3] = { 0, 0, 0 };
  AffineMapMacroMesh(&mesh, A, x0);

  // mesh preparation
  mesh.period[0] = 1.0;
  mesh.period[1] = 1.0;
  BuildConnectivity(&mesh);
  //
  CheckMacroMesh(&mesh, SimParams.deg, SimParams.raf);
  // setup simulation paramaters in global shared LatticeBoltzmannSimData object
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  lsd->lb_model = &lbm;
  LBM_Set_D2Q9_ISOTH_LINEARIZED_model(&lbm, SimParams.cref);
  // setup LB Simulation object
  LBMSimulation lbsimu;
  lbsimu.macro_model.InitData = LBM_Linear2DWave_InitData;
  lbsimu.macro_model.ImposedData = LBM_Linear2DWave_ImposedData;
  //
  lbsimu.micro_model.NumFlux = NULL;
  lbsimu.micro_model.InitData = LBM_Dummy_InitData;
  lbsimu.micro_model.ImposedData = NULL;
  lbsimu.micro_model.BoundaryFlux = NULL;
  lbsimu.micro_model.Source = NULL;
  //
  InitLBMSimulation(&lbsimu, lsd, &mesh, SimParams.deg, SimParams.raf);
  //
  lbsimu.vmax = lsd->lb_model->vmax;
  lbsimu.cfl = SimParams.cfl;
  lbsimu.micro_simu.cfl = lbsimu.cfl;
  lbsimu.tmax = SimParams.tmax;
  lbsimu.micro_simu.tmax = lbsimu.tmax;
  lbsimu.dt = SimParams.dt;
  lbsimu.micro_simu.cfl = lbsimu.cfl;
  lbsimu.macro_simu.cfl = lbsimu.cfl;
  lbsimu.micro_simu.tmax = lbsimu.tmax;
  lbsimu.macro_simu.tmax = lbsimu.tmax;
  //
  schnaps_real tau = SimParams.tau;
  // relaxation parameter for CN scheme
  lbm.s[0] = lbsimu.dt / (tau + 0.5 * lbsimu.dt); 
  printf(" tau=%f s=%f\n", tau, lbm.s[0]);
  //
  lbsimu.macro_simu.nb_diags = 4;
  lbsimu.micro_simu.nb_diags = 0;
  lbsimu.pre_advec = LB_Relaxation_bgk_f;
  lbsimu.post_advec_one_node=NULL;
  lbsimu.post_advec = LB_ComputeMacroFromMicro;
  //lbsimu.post_advec=NULL;
  lbsimu.post_tstep = LBM_Linear2DWave_Plot_Fields;
  lbsimu.collect_diags = LBM_Linear2DWave_CollectDiags;
  lbsimu.diag_2d_period = SimParams.diag_2d_period;
  //
  sprintf(simutag, "IMP");
  lbsimu.model_advec.m = 1;
  lbsimu.model_advec.InitData = LBM_Dummy_InitData_OneNode;
  lbsimu.model_advec.ImposedData = LBM_Linear2DWave_ImposedData_OneNode;
  lbsimu.model_advec.NumFlux = LBM_OneNodeNumFlux;
  lbsimu.model_advec.BoundaryFlux =
      LBM_Linear2DWave_Periodic_BoundaryFlux_OneNode;
  lbsimu.model_advec.Source = NULL;
  //
  // Actual run
  LBMThetaTimeScheme(&lbsimu, 0.5, lbsimu.tmax, lbsimu.dt);
  //LBMThetaTimeScheme_OneNode_Multistep(&lbsimu, 0.5, lbsimu.tmax, lbsimu.dt);
  //
  LBM_Dump_Lattice_Diagnostics(&lbsimu, "SIM");
  printf(" tau=%f s=%f\n", tau, lbm.s[0]);
  // cleanup
  FreeBMSimulation(&lbsimu);
  lsd->lb_model = NULL;
  DestroyLBModelDescriptor(&lbm);
  return 1;
}
/*int LBM_TestLattice_LinearWave2D_NLImplicit(void)*/
/*{*/
/*  int d = 2;*/
/*  int nb_macro = 3;*/
/*  int q = 9;*/
/*  LBModelDescriptor lbm = LBModelDescriptor_NULL;*/
/*  //*/
/*  MacroMesh mesh;*/
/*  ReadMacroMesh(&mesh, "../test/testcube.msh");*/
/*  Detect2DMacroMesh(&mesh);*/
/*  //*/
/*  schnaps_real A[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };*/
/*  schnaps_real x0[3] = { 0, 0, 0 };*/
/*  AffineMapMacroMesh(&mesh, A, x0);*/

/*  // mesh preparation*/
/*  mesh.period[0] = 1.0;*/
/*  mesh.period[1] = 1.0;*/
/*  BuildConnectivity(&mesh);*/
/*  //*/
/*  CheckMacroMesh(&mesh, SimParams.deg, SimParams.raf);*/
/*  // setup simulation paramaters in global shared LatticeBoltzmannSimData object*/
/*  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;*/
/*  NewLBModelDescriptor(&lbm, d, nb_macro, q);*/
/*  lsd->lb_model = &lbm;*/
/*  LBM_Set_D2Q9_ISOTH_LINEARIZED_model(&lbm, SimParams.cref);*/
/*  // setup LB Simulation object*/
/*  LBMSimulation lbsimu;*/
/*  lbsimu.macro_model.InitData = LBM_Linear2DWave_InitData;*/
/*  lbsimu.macro_model.ImposedData = LBM_Linear2DWave_ImposedData;*/
/*  //*/
/*  //*/
/*  InitLBMSimulation(&lbsimu, lsd, &mesh, SimParams.deg, SimParams.raf);*/
/*  lbsimu.micro_model.NumFlux=LBM_OneLatticeNumFlux_SelectVel;*/
/*  //*/
/*  lbsimu.vmax = lsd->lb_model->vmax;*/
/*  lbsimu.cfl = SimParams.cfl;*/
/*  lbsimu.micro_simu.cfl = lbsimu.cfl;*/
/*  lbsimu.tmax = SimParams.tmax;*/
/*  lbsimu.micro_simu.tmax = lbsimu.tmax;*/
/*  lbsimu.dt = SimParams.dt;*/
/*  lbsimu.micro_simu.cfl = lbsimu.cfl;*/
/*  lbsimu.macro_simu.cfl = lbsimu.cfl;*/
/*  lbsimu.micro_simu.tmax = lbsimu.tmax;*/
/*  lbsimu.macro_simu.tmax = lbsimu.tmax;*/
/*  //*/
/*  schnaps_real tau = SimParams.tau;*/
/*  lbm.s[0] = lbsimu.dt / (tau + 0.5 * lbsimu.dt);*/
/*  printf(" tau=%f s=%f\n", tau, lbm.s[0]);*/
/*  //*/
/*  lbsimu.macro_simu.nb_diags = 0;*/
/*  lbsimu.micro_simu.nb_diags = 0;*/
/*  lbsimu.pre_advec = NULL;*/
/*  lbsimu.post_advec_one_node=NULL;*/
/*  lbsimu.post_advec = NULL;*/
/*  //lbsimu.post_advec=NULL;*/
/*  lbsimu.post_tstep = LBM_Linear2DWave_Plot_Fields;*/
/*  lbsimu.collect_diags = LBM_Linear2DWave_CollectDiags;*/
/*  lbsimu.diag_2d_period = SimParams.diag_2d_period;*/
/*  //*/
/*  sprintf(simutag, "IMP");*/
/*  // Actual run*/
/*  LBMThetaTimeScheme_NL_LocalNewton(&lbsimu, 0.5, lbsimu.tmax, lbsimu.dt);*/
/*  //*/
/*  //LBM_Dump_Lattice_Diagnostics(&lbsimu, "SIM");*/
/*  //printf(" tau=%f s=%f\n", tau, lbm.s[0]);*/
/*  // cleanup*/
/*  FreeBMSimulation(&lbsimu);*/
/*  lsd->lb_model = NULL;*/
/*  DestroyLBModelDescriptor(&lbm);*/
/*  return 1;*/
/*}*/

//
void LBM_Linear2DWave_InitData(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real my_pi = 4.0 * atan(1.0);
  // wave mode numbers in half integer units
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset = LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky = 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix = kx * x[0];
  schnaps_real phiy = ky * x[1];
  //
  schnaps_real rho = offset + cos(phix) * cos(phiy);
  schnaps_real jx = 0.0;
  schnaps_real jy = 0.0;
  //
  w[0] = rho;
  w[1] = jx;
  w[2] = jy;
}

//*
void LBM_Linear2DWave_ImposedData(const schnaps_real x[3],
				  const schnaps_real t, schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real my_pi = 4.0 * atan(1.0);
  // wave mode numbers in half integer units
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset = LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky = 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix = kx * x[0];
  schnaps_real phiy = ky * x[1];
  //
  schnaps_real k = sqrt(kx * kx + ky * ky);
  schnaps_real c0 = lsd->lb_model->cref * sqrt(1 / 3.0);
  schnaps_real omega = k * c0;
  schnaps_real phit = omega * t;
  //
  schnaps_real rho = offset + cos(phix) * cos(phiy) * cos(phit);
  schnaps_real jx = c0 * kx * sin(phix) * cos(phiy) * sin(phit) / k;
  schnaps_real jy = c0 * ky * cos(phix) * sin(phiy) * sin(phit) / k;
  //
  w[0] = rho;
  w[1] = jx;
  w[2] = jy;
}

void LBM_Linear2DWave_ImposedData_OneNode(const schnaps_real x[3],
					  const schnaps_real t,
					  schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real my_pi = 4.0 * atan(1.0);
  // wave mode numbers in half integer units
  int nkx = LW2DParams.nkx;
  int nky = LW2DParams.nky;
  schnaps_real offset = LW2DParams.offset;
  // spatial frequencies
  schnaps_real kx = 2.0 * my_pi * (schnaps_real) nkx;
  schnaps_real ky = 2.0 * my_pi * (schnaps_real) nky;
  //
  schnaps_real phix = kx * x[0];
  schnaps_real phiy = ky * x[1];
  //
  schnaps_real k = sqrt(kx * kx + ky * ky);
  schnaps_real c0 = lsd->lb_model->cref * sqrt(1 / 3.0);
  schnaps_real omega = k * c0;
  schnaps_real phit = omega * t;
  //
  schnaps_real rho = offset + cos(phix) * cos(phiy) * cos(phit);
  schnaps_real jx = c0 * kx * sin(phix) * cos(phiy) * sin(phit) / k;
  schnaps_real jy = c0 * ky * cos(phix) * sin(phiy) * sin(phit) / k;
  //
  schnaps_real wmac[lsd->lb_model->nb_macro];
  wmac[0] = rho;
  wmac[1] = jx;
  wmac[2] = jy;
  //
  int inode = lsd->current_node_index;
  w[0] = lsd->lb_model->feq_i(inode, lsd->lb_model->nb_macro, wmac);
  //
}

void LBM_Linear2DWave_Periodic_BoundaryFlux_OneNode(schnaps_real * x,
						    schnaps_real t,
						    schnaps_real * wL,
						    schnaps_real * vnorm,
						    schnaps_real * flux)
{
  LatticeData *ld = &schnaps_lattice_data;
  int i_node = ld->current_node_index;
  schnaps_real wR[1];
  LBM_Linear2DWave_ImposedData_OneNode(x, t, wR);
  LBM_OneNodeNumFlux(wL, wR, vnorm, flux);
}

void LBM_Linear2DWave_CollectDiags(void *s, schnaps_real * macro_diag_vals,
				   schnaps_real * micro_diag_vals)
{
  LBMSimulation *lbsimu = s;
  Simulation *simu = &(lbsimu->macro_simu);
  schnaps_real error = 0;
  schnaps_real mean = 0;
  for (int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field *f = simu->fd + ie;
    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for (int ipg = 0; ipg < npg; ipg++) {
      schnaps_real w[f->model.m];
      for (int iv = 0; iv < f->model.m; iv++) {
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
      schnaps_ref2phy(f->physnode,	// phys. nodes
		      xpgref,	// xref
		      NULL, -1,	// dpsiref, ifa
		      xphy, dtau,	// xphy, dtau
		      codtau, NULL, NULL);	// codtau, dpsi, vnds
      det = dot_product(dtau[0], codtau[0]);
      // Get the exact value
      f->model.ImposedData(xphy, simu->tnow, wex);
      schnaps_real diff = w[0] - wex[0];
      error += diff * diff * wpg * det;
      mean += wex[0] * wex[0] * wpg * det;
    };
  };
  //printf("errl2=%f\n",sqrt(error) / (sqrt(mean)  + 1e-14));
  macro_diag_vals[0] = sqrt(error);
  macro_diag_vals[1] = sqrt(mean);
  macro_diag_vals[2] = sqrt(error) / (sqrt(mean) + 1e-14);
  //
  schnaps_real xphytest[3];
  schnaps_real wtest[3];
  xphytest[0] = 0.5;
  xphytest[1] = 0.5;
  xphytest[2] = 0.0;
  lbsimu->macro_model.ImposedData(xphytest, simu->tnow, wtest);
  macro_diag_vals[3] = wtest[0];
  //
}

void LBM_Linear2DWave_Plot_Fields(void *s, schnaps_real * w)
{
  LBMSimulation *lbsimu = s;
  Simulation *simu = &(lbsimu->macro_simu);
  schnaps_real period = lbsimu->diag_2d_period;
  schnaps_real dt = lbsimu->dt;
  int diagperiod = (int) (period / dt);
  int istep = lbsimu->micro_simu.iter_time_rk;
  schnaps_real t = lbsimu->micro_simu.tnow;
  int tmax = simu->tmax;
  int create_file = 0;
  LBM_Store_Lattice_diags(lbsimu);
  if (istep == 0) {
    create_file = 1;
  } else {
    create_file = 0;
  }
  if (diagperiod || create_file) {
    if ((istep % diagperiod == 0) || (t == tmax)) {
      istep = istep + 1;
      printf("Dumping fields at it=%i (period %i)\n", istep, diagperiod);
      int raf = simu->fd[0].raf[0];
      schnaps_real cfl = simu->cfl;
      char filename_rho[sizeof("lbm2DWave_rho_TAG_raf000_cfl0.000.msh")];
      sprintf(filename_rho, "lbm_2DWave_rho_%s_raf%03d_cfl%1.3f.msh",
	      simutag, raf, cfl);
      //char filename_rho_error[sizeof("lbm2DWave_rho_000.msh")];
      //sprintf(filename_rho_error,"lbm_2DWave_rho_error_%03d.msh",raf);
      //PlotFields(0,false,simu,"rho",filename_rho);
      LBM_PlotFieldsBinSparseMultitime(0, false, simu, "rho", filename_rho,
				       create_file, t, istep);
      LBM_PlotFieldsBinSparseMultitime(0, true, simu, "rho_error",
				       filename_rho, 0, t, istep);
      //LBM_PlotFieldsBinSparseMultitime(1,false,simu,"jx",filename_rho,0,t,istep);
      //LBM_PlotFieldsBinSparseMultitime(1,true,simu,"jx_error",filename_rho,0,t,istep);
    };
  }
};

//************************************************************************************************//
//************************************************************************************************//
//* basic test for NL solver
//************************************************************************************************//
//************************************************************************************************//
int LBM_TestLattice_SteadyFlow_NLImplicit(schnaps_real alpha_deg)
{
  int d = 2;
  int nb_macro = 2;
  int q = 3;
  LBModelDescriptor lbm = LBModelDescriptor_NULL;
  //
  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
  schnaps_real x0[3] = { 0, 0, 0 };
  AffineMapMacroMesh(&mesh, A, x0);

  // mesh preparation
  mesh.period[0] = 0.0;
  mesh.period[1] = 0.0;
  BuildConnectivity(&mesh);
  //
  CheckMacroMesh(&mesh, SimParams.deg,SimParams.raf);
  // setup simulation paramaters in global shared LatticeBoltzmannSimData object
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  NewLBModelDescriptor(&lbm, d, nb_macro, q);
  lsd->lb_model = &lbm;
  LBM_Set_D2Q3_TEST_model(&lbm, alpha_deg, SimParams.cref);
  // setup LB Simulation object
  LBMSimulation lbsimu;
  lbsimu.macro_model.InitData = LBM_SteadyFlow_NLImplicit_MacInitData;
  lbsimu.macro_model.ImposedData = LBM_SteadyFLow_NLImplicit_MacImposedData;
  //
  //
  lbsimu.micro_model.ImposedData = LBM_SteadyFLow_NLImplicit_MicImposedData;
  lbsimu.micro_model.InitData = LBM_SteadyFlow_NLImplicit_MacInitData;
  lbsimu.micro_model.NumFlux=LBM_OneLatticeNumFlux_SelectVel;
  lbsimu.micro_model.Source= LBM_BGK_Relaxation_Source_TestDecay;
  lbsimu.micro_model.SourceJac= LBM_BGK_Relaxation_Source_Jacobian_TestDecay;
  lbsimu.micro_model.BoundaryFlux = LBM_SteadyFlow_NLImplicit_BoundaryFlux;
  InitLBMSimulation(&lbsimu, lsd, &mesh, SimParams.deg, SimParams.raf);
  //
  int testvarindex= test_reverse_varindex(&(lbsimu.micro_simu));
  assert(testvarindex==1);
  if (testvarindex==1){
    printf(" ******* Reverse varindex ok ***********\n");
  }
  //
  int testmulti_index= CheckLBM_imem_multivarindex_constistency(&(lbsimu.micro_simu));
  assert(testmulti_index==1);
  if (testmulti_index==1){
    printf(" ******* Reverse Multindex ok ***********\n");
  }
  //
  lbsimu.vmax = lsd->lb_model->vmax;
  lbsimu.cfl = SimParams.cfl;
  lbsimu.micro_simu.cfl = lbsimu.cfl;
  lbsimu.tmax = SimParams.tmax;
  lbsimu.micro_simu.tmax = lbsimu.tmax;
  lbsimu.dt = SimParams.dt;
  lbsimu.micro_simu.cfl = lbsimu.cfl;
  lbsimu.macro_simu.cfl = lbsimu.cfl;
  lbsimu.micro_simu.tmax = lbsimu.tmax;
  lbsimu.macro_simu.tmax = lbsimu.tmax;
  //
  schnaps_real tau = SimParams.tau;
  //lbm.s[0] = lbsimu.dt / (tau + 0.5 * lbsimu.dt);
  lbm.s[0] = 1.0/tau;
  printf(" tau=%f s=%f\n", tau, lbm.s[0]);
  //
  lbsimu.macro_simu.nb_diags = 0;
  lbsimu.micro_simu.nb_diags = 0;
  lbsimu.pre_advec = NULL;
  lbsimu.post_advec_one_node=NULL;
  lbsimu.post_advec = NULL;
  //lbsimu.post_advec=NULL;
  lbsimu.post_tstep = NULL;
  lbsimu.collect_diags = NULL;
  lbsimu.diag_2d_period = SimParams.diag_2d_period;
  //
  sprintf(simutag, "IMP");
  // Actual run
  LBMThetaTimeScheme_NL_LocalNewton(&lbsimu, 0.5, lbsimu.tmax, lbsimu.dt);
  //
  // cleanup
  FreeBMSimulation(&lbsimu);
  lsd->lb_model = NULL;
  DestroyLBModelDescriptor(&lbm);
  return 1;
}
//
void LBM_SteadyFlow_NLImplicit_MacInitData(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  schnaps_real rho = 1.0;
  schnaps_real upar = 0.1;
  //
  w[0] = rho;
  w[1] = upar;
}
void LBM_SteadyFlow_NLImplicit_MicInitData(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  schnaps_real rho = 1.0;
  schnaps_real upar = 0.1;
  //
  schnaps_real wmac[2];
  wmac[0] = rho;
  wmac[1] = upar;
  int q = lsd->lb_model->q;
  lsd->lb_model->feq(q,2,wmac,w);
}
//*
void LBM_SteadyFLow_NLImplicit_MacImposedData(const schnaps_real x[3],
				  const schnaps_real t, schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  schnaps_real rho = 1.0;
  schnaps_real upar = 0.1;
  //
  w[0] = rho;
  w[1] = upar;
}
void LBM_SteadyFLow_NLImplicit_MicImposedData(const schnaps_real x[3],
				  const schnaps_real t, schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  schnaps_real wmac[2];
  schnaps_real rho = 1.0;
  schnaps_real upar = 0.1;
  //
  wmac[0] = rho;
  wmac[1] = upar;
  //
  int q = lsd->lb_model->q;
  lsd->lb_model->feq(q,2,wmac,w);
}
void LBM_SteadyFlow_NLImplicit_BoundaryFlux(schnaps_real * x,
						    schnaps_real t,
						    schnaps_real * wL,
						    schnaps_real * vnorm,
						    schnaps_real * flux)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real wR[lsd->lb_model->q];
  LBM_SteadyFLow_NLImplicit_MicImposedData(x,t,wR);
  LBM_OneLatticeNumFlux_SelectVel(wL,wR,vnorm,flux);
}
//
int test_reverse_varindex(Simulation *simu){
  for (int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field *f = simu->fd + ie;
    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for (int ipg = 0; ipg < npg; ipg++) {
      schnaps_real w[f->model.m];
      for (int iv = 0; iv < f->model.m; iv++) {
        int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
        int testiv=-1,testipg=-1;
        f->varindex_to_ipgiv(f->deg,f->raf,f->model.m,imem,&testipg,&testiv);
        //printf("(ipg,iv) in (%i,%i) \t imem %i \t (ipg,iv) out (%i,%i) \n",ipg,iv,imem,testipg,testiv);
        assert(testipg==ipg);
        assert(testiv==iv);
      } 
    }
  }
  return 1;
}
//

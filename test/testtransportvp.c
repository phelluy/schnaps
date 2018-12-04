#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"

int Test_TransportVP(void);
void UpdateVlasovPoisson(void *field);
void PlotVlasovPoisson(void *vf, schnaps_real *w);
void Test_TransportVP_ImposedData(const schnaps_real *x,
				  const schnaps_real t,
				  schnaps_real *w);
void Test_TransportVP_InitData(schnaps_real *x,
			       schnaps_real *w);
void Test_TransportVP_BoundaryFlux(schnaps_real *x,
				   schnaps_real t,
				   schnaps_real *wL,
				   schnaps_real *vnorm,
				   schnaps_real *flux);
schnaps_real TransportVP_ImposedKinetic_Data(const schnaps_real *x,
					     const schnaps_real t,
					     schnaps_real v);

int main(void) {

  // unit tests
  int resu = Test_TransportVP();

  if (resu) printf("poisson test OK !\n");
  else printf("poisson test failed !\n");

  return !resu;
}


int Test_TransportVP(void) {

  bool test = true;

#ifdef PARALUTION
  paralution_begin();
#endif

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect1DMacroMesh(&mesh);

  bool is1d = mesh.is1d;
  assert(is1d);

  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;
  schnaps_real degV=2;
  schnaps_real nbEV=24;
  KineticData * kd=&schnaps_kinetic_data;

  InitKineticData(kd,nbEV,degV);

  model.m=kd->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.InitData = Test_TransportVP_InitData;
  model.ImposedData = Test_TransportVP_ImposedData;
  model.BoundaryFlux = Test_TransportVP_BoundaryFlux;
  model.Source = VlasovP_Lagrangian_Source;

  int deg[]={2, 0, 0};
  int raf[]={16, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = kd->vmax; // maximal wave speed
  simu.cfl=0.2;
  simu.nb_diags = 4;
  simu.pre_dtfields = UpdateVlasovPoisson;
  simu.post_dtfields=NULL;
  simu.update_after_rk = PlotVlasovPoisson;

  schnaps_real tmax = 0.03;

  RK2(&simu, tmax);

   // save the results and the error
  int iel = 2 * kd->nb_elem_v / 3;
  int iloc = kd->deg_v;
  printf("Trace vi=%f\n", -kd->vmax + iel * kd->dv + kd->dv * glop(kd->deg_v, iloc));
  PlotFields(iloc + iel * kd->deg_v, false, &simu, "sol","dgvisu_kin.msh");
  PlotFields(iloc + iel * kd->deg_v, true, &simu, "error","dgerror_kin.msh");
  Plot_Energies(&simu, simu.dt);

  schnaps_real dd_Kinetic = L2_Kinetic_error(&simu);

  printf("erreur kinetic L2=%lf\n", dd_Kinetic);
  test= test && (dd_Kinetic < 1e-2);

  FreeMacroMesh(&mesh);

  return test;
}

void Test_TransportVP_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  KineticData * kd=&schnaps_kinetic_data;
  for(int i = 0; i <kd->index_max_kin + 1; i++) {
    int j = i % kd->deg_v; // local connectivity put in function
    int nel = i / kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax + nel * kd->dv + kd->dv * glop(kd->deg_v, j));
 
    w[i] = TransportVP_ImposedKinetic_Data(x, t, vi);
  }
  // exact value of the potential and electric field
  w[kd->index_phi] = -x[0];
  w[kd->index_ex] = 1.0;
  w[kd->index_rho] = 0.; //rho init
  w[kd->index_u] = 0; // u init
  w[kd->index_P] = 0; // p init
  w[kd->index_T] = 0; // e ou T init
}

void Test_TransportVP_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  Test_TransportVP_ImposedData(x, t, w);
}

schnaps_real TransportVP_ImposedKinetic_Data(const schnaps_real *x, const schnaps_real t, schnaps_real v) {
  schnaps_real f;
  schnaps_real pi = 4.0 * atan(1.0);
  schnaps_real xnew = 0, vnew = 0;
  f = exp(-(v - t) * (v - t)) *
    exp(-36 * ((x[0] - v * t + 0.5 * t * t) - 0.5)
	* ((x[0] - v * t + 0.5 * t * t) - 0.5));
  return f;
}

void Test_TransportVP_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  Test_TransportVP_ImposedData(x , t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
}

void UpdateVlasovPoisson(void *si) {
  Simulation *simu = si;
  KineticData * kd=&schnaps_kinetic_data;
  int type_bc = 1;
    
  // Computation_charge_density(simu,simu->w);
  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=kd->index_phi;
  
  InitContinuousSolver(&ps,simu,1,nb_var,listvar);

  ps.matrix_assembly=ContinuousOperator_Poisson1D;
  ps.rhs_assembly=RHSPoisson_Continuous;
  ps.bc_assembly= ExactDirichletContinuousMatrix;
  ps.postcomputation_assembly=Computation_ElectricField_Poisson;

#ifdef PARALUTION
  ps.lsol.solver_type = PAR_LU;
  ps.lsol.pc_type=NONE;
#else
  ps.lsol.solver_type = LU;
  ps.lsol.pc_type=NONE;
#endif

  SolveContinuous2D(&ps);
  //freeContinuousSolver(&ps);
}

void PlotVlasovPoisson(void *si, schnaps_real *w) {
  schnaps_real k_energy = 0, e_energy = 0, t_energy = 0, charge =0;
  Simulation *simu = si;

  Energies(simu, w, k_energy, e_energy, t_energy,1);
  Charge_total(simu, w, charge, 4);
  si = simu;
}

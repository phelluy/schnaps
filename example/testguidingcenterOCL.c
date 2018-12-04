#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "quantities_vp.h"
#include "gyro.h"
#include "solverpoisson.h"
#include "diagnostics_vp.h"

int TestGC_OCL(void);

void PlotVlasovPoisson(void* field, schnaps_real * w);

int main(void) {
  

  int resu=TestGC_OCL();
	 
  if (resu) printf("guiding center open CL test OK !\n");
  else printf("guiding center open CL test failed !\n");

  return !resu;
} 


int TestGC_OCL(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh"); 
  mesh.period[2]=1;
  BuildConnectivity(&mesh);
    
  int deg[]={2, 2, 1};
  int raf[]={8, 8, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  int nbelemv = 1;
  int deg_v = 0;
  
  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  kd->vmax = 0.5;
  kd->dv = 1;
  
  kd->solve_quasineutrality = true;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=GyroUpwindNumFlux;
  //model.NumFlux=GyroZeroNumFlux;
  model.BoundaryFlux=GC_OCLBoundaryFlux;
  model.InitData=GC_OCLInitData;
  model.ImposedData=GC_OCLImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;

  #if OCL
  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _NBELEMV=%d ",nbelemv);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _DEGV=%d ",deg_v);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _VMAX=%f ", kd->vmax);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D _PERIODZ=%f ", mesh.period[2]);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s","GC_OCLUpwindNumFlux");
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D BOUNDARYFLUX=%s","GC_OCLBoundaryFlux");
  strcat(cl_buildoptions, buf);

  schnaps_ocl_getcharge = true;
  #endif
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  //simu.pre_dtfields = UpdateGyroPoisson;
  simu.vmax = kd->vmax; // maximal wave speed
  simu.nb_diags = 6;
  /*simu.pre_dtfields =  UpdateGyroPoisson; XXXXXXX
  simu.post_dtfields = NULL;
  simu.update_after_rk = PlotVlasovPoisson; XXXXXXX
  */
 
  
  simu.cfl=0.7;
   schnaps_real dt = 0;
  schnaps_real tmax = 10;
  //UpdateGyroPoisson(&simu);
  GyroCFLVelocity(&simu);
#if OCL
 // OpenCL version
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return 0;
  }
  simu.update_after_rk = PlotVP;
  RK4_CL(&simu, tmax, dt,  0, 0, 0);
  CopyfieldtoCPU(&simu);
  printf("\nOpenCL Kernel time:\n");
  show_cl_timing(&simu);
  printf("\n");
#else
  simu.pre_dtfields =  UpdateGyroPoisson;
  simu.post_dtfields = NULL;
  simu.update_after_rk = PlotVP;
  UpdateGyroPoisson(&simu);
  RK2(&simu, tmax);
#endif
  
  Plot_Energies(&simu, simu.dt);

  PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","gc/rho.msh"); 
  PlotFields(kd->index_phi,(1==1),&simu,"sol_phi","gc/potential.msh"); 
  PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","gc/soluphi.msh");
  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  double dd_Kinetic=L2_Kinetic_error(&simu);
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd_Kinetic < 0.005);


  return test; 

};


  
void PlotVlasovPoisson(void *si, schnaps_real *w) {
  schnaps_real k_energy = 0, e_energy = 0, t_energy = 0, t_charge=0;
  schnaps_real taux_ins_m1 = 0;
  
  Simulation *simu = si;
  
  Energies(simu, w, k_energy, e_energy, t_energy, 1);
  Charge_total(simu, w, t_charge, 4);
  Taux_instability(simu, w, 7, taux_ins_m1, 5);
  //Taux_instability(simu, w, 4, taux_ins_m1, 6);
  si = simu; 
}
 

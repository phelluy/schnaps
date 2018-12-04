#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "quantities_vp.h"
#include "gyro.h"


int TestGyro(void);

int main(void) {
  
  // unit tests
    
  int resu=TestGyro();
	 
  if (resu) printf("gyro test OK !\n");
  else printf("gyro test failed !\n");

  return !resu;
} 


int TestGyro(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  //ReadMacroMesh(&mesh,"../test/unit-cube.msh");
  ReadMacroMesh(&mesh,"../geo/cylindre.msh");
  mesh.period[2]=1;
  BuildConnectivity(&mesh);

  
    
  int deg[]={2, 2, 2};
  int raf[]={16, 16, 16};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;
  KineticData *kd = &schnaps_kinetic_data;
  int deg_v = 3;
  int nbelemv = 6;

  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  kd->solve_quasineutrality = false;//true;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;


  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=GyroUpwindNumFlux;
  //model.NumFlux=GyroCenteredNumFlux;
  //model.NumFlux=GyroZeroNumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=GyroBoundaryFlux;
  model.InitData=GyroInitData;
  model.ImposedData=GyroImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;

  
  #if OCL
  bool cachememory = false;
  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _NBELEMV=%d ",nbelemv);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _DEGV=%d ",deg_v);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _VMAX=%f ", kd->vmax);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _CACHEMEMORY=%d ", cachememory);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D _PERIODZ=%f ", mesh.period[2]);
  strcat(cl_buildoptions, buf);
  //sprintf(buf," -D NUMFLUX=%s","GyroZeroNumFlux");
  sprintf(buf," -D NUMFLUX=%s","GyroUpwindNumFlux");
  strcat(cl_buildoptions, buf);
  sprintf(buf," -D BOUNDARYFLUX=%s", "GyroBoundaryFlux");
  strcat(cl_buildoptions, buf);
  
  schnaps_ocl_getcharge = false;
  #endif
  
  InitSimulation(&simu, &mesh, deg, raf, &model);
 /* CopyfieldtoCPU(&simu); */
 /* PlotFields(kd->index_phi,(1==0),&simu,"phi","init_phi.msh"); */
 /*  PlotFields(kd->index_rho,(1==0),&simu,"rho","init_rho.msh"); */
 /*  assert(1==2); */
  //PlotFields(0,(1==0),&simu,"init_f_p0","init_f.msh");
  //PlotFields(kd->index_phi,(1==0),&simu,"init_phi","init_potential.msh");
  //simu.pre_dtfields = UpdateGyroPoisson;
   simu.vmax = kd->vmax; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.3;
  schnaps_real dt = 0;
  schnaps_real tmax = 0.0001;
  //RK2(&simu,tmax);
  GyroCFLVelocity(&simu);
  simu.nb_diags = 6;
 #if OCL
  simu.update_after_rk = NULL;//PlotVP;
  RK4_CL(&simu,tmax, dt, 0, 0, 0);
  //Computation_charge_density(&simu);
  // save the results and the error
  //PlotFields(1,(1==0),&simu,"sol","dgvisu.msh");
  CopyfieldtoCPU(&simu);
  #else
  /*InitSimulation(&simu, &mesh, deg, raf, &model);
	simu.vmax = kd->vmax; // maximal wave speed 
  */
  simu.pre_dtfields =NULL;// UpdateGyroPoisson;
  simu.post_dtfields = NULL;
  simu.update_after_rk = NULL;//PlotVP;
  //simu.cfl=0.2;
  //schnaps_real tmax = 0.0005;
  //GyroCFLVelocity(&simu);
  //UpdateGyroPoisson(&simu);
  RK4(&simu,tmax);
  #endif
  /*Plot_Energies(&simu, simu.dt);
  PlotFields(0,(1==0),&simu,"init_f_p0","sol_f.msh");
  PlotFields(kd->index_phi,(1==0),&simu,"sol","sol_phi.msh");
  PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","sol_rho.msh");
  */
  PlotFields(0,(1==0),&simu,"f_p0","dis.msh");
  PlotFields(kd->index_phi,(1==0),&simu,"sol","sol_phi.msh");
  PlotFields(kd->index_rho,(1==0),&simu,"sol_rho","sol_rho.msh");
  double dd=L2error(&simu);
  //double dd_l2_vel =GyroL2VelError(&f)
  double dd_Kinetic=L2_Kinetic_error(&simu);
  
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd < 0.005);


  return test; 

};

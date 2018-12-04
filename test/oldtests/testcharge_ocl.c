#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "simulation.h"
#include "quantities_vp.h"
#include "gyro.h"

int TestChargeOCL(void);		      
		  			  			  
int main(void) {
  
  // unit tests
    
  int resu=TestChargeOCL();
	 
  if (resu) printf("chargeOCL test OK !\n");
  else printf("chargeOCL test failed !\n");

  return !resu;
} 


int TestChargeOCL(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh");
  mesh.period[2]=1;
  BuildConnectivity(&mesh);
    
  int deg[]={2, 2, 2};
  int raf[]={1, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  int nbelemv = 1;
  int deg_v = 2;
  InitKineticData(&schnaps_kinetic_data,nbelemv,deg_v);
  kd->solve_quasineutrality = true;
  kd->substract_mean_charge = false;
  kd->qn_damping = 0;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=GyroZeroNumFlux;
  model.BoundaryFlux=ChargeOCLBoundaryFlux;
  model.InitData=ChargeOCLInitData;
  model.ImposedData=ChargeOCLImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;

  #if OCL
  bool cachememory = true;
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
  sprintf(buf, " -D _CACHEMEMORY=%d ", cachememory);
  strcat(cl_buildoptions, buf);
  sprintf(buf," -D NUMFLUX=%s","GyroZeroNumFlux");
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D BOUNDARYFLUX=%s","ChargeOCLBoundaryFlux");
  strcat(cl_buildoptions, buf);

  schnaps_ocl_getcharge = true;
  #endif

  InitSimulation(&simu, &mesh, deg, raf, &model);
  //simu.pre_dtfields = UpdateGyroPoisson;
   simu.vmax = kd->vmax; // maximal wave speed 
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  simu.cfl=0.2;
  schnaps_real dt = 0;
  schnaps_real tmax = 0.0005;
 GyroCFLVelocity(&simu);
  simu.nb_diags = 6;
  
 #if OCL
 RK4_CL(&simu,tmax, dt, 0, 0, 0);
  CopyfieldtoCPU(&simu);
  #else
  simu.pre_dtfields =UpdateGyroPoisson;
  simu.post_dtfields = NULL;
  simu.update_after_rk = NULL;//PlotVP;
  //simu.cfl=0.2;
  //schnaps_real tmax = 0.0005;
  //GyroCFLVelocity(&simu);
  UpdateGyroPoisson(&simu);
  RK4(&simu,tmax);
  #endif


  PlotFields(kd->index_phi,(1==0),&simu,"sol_phi","sol_potential.msh");
  PlotFields(kd->index_phi,(1==1),&simu,"err_phi","err_potential.msh"); 

  PlotFields(kd->index_rho,(1==1),&simu,"sol_rho","sol_rho.msh"); 
  PlotFields(kd->index_ex,(1==1),&simu,"sol_ex","sol_ex.msh"); 
  PlotFields(kd->index_ey,(1==1),&simu,"sol_ey","sol_ey.msh"); 
  /* PlotFields(kd->index_rho,(1==1),&simu,"err_rho","err_rho.msh"); */
  //PlotFields(1,(1==1),&simu,"error","dgerror.msh");

  double dd=L2error(&simu);
  //double dd_l2_vel =ChargeOCLL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&f);
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",simu.tnow);
  Velocity_distribution_plot(simu.w);
  test= test && (dd < 0.3);


  return test; 

};

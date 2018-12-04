#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "acoustic.h"
#include "getopt.h"
#include "global.h"
#include "io.h"
#include "schnaps.h"
#include "../test/test.h"

#define HEADER_MAX_SIZE 9999
#define PATH_MAX_SIZE 9999
#define EXPORT_XDMF 0
#define EXPORT_CSV_RCP 1


void sn_rk2_spu(Simulation *simu, schnaps_real tmax){

  //Temporal infos
  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;
  simu->tmax = tmax;
  simu->itermax_rk = tmax / simu->dt + 1;
  simu->tnow = 0;
  
  //Setup diag freq
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10) ? 1 : simu->itermax_rk / 10;
  int iter = 0;
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

  assert(starpu_use);
 
  //Exporter
  schnaps_real rec_rho=0;
  schnaps_real rec_rho_ar[1];
  
  RegisterSimulation_SPU(simu);

  while(simu->tnow < tmax) {
      
    if (EXPORT_CSV_RCP) {
        UnregisterSimulation_SPU(simu);
        if(iter == 0){
            sn_compute_macro(simu, &rec_rho_ar[0]);
            printf("t=%f, Receptor value = %f \n",simu->tnow, rec_rho_ar[0]);
        } else {
            sn_compute_macro(simu, &rec_rho_ar[0]);
            printf("t=%f, Receptor value = %f \n",simu->tnow, rec_rho_ar[0]);
        }
        RegisterSimulation_SPU(simu);
    }
  
     
    
    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
#ifndef WORKER_ON_NODE
      ZeroBuffer_SPU(simu->fd[ie].dtwn_handle);
#else
      ZeroBuffer_SPU2(simu->fd[ie].dtwn_handle, ie);
#endif
    }

    DtFields_SPU(simu, NULL, NULL);

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      schnaps_real alpha = dt/2;
#ifndef WORKER_ON_NODE
      AddBuffer_SPU(alpha, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle);
#else
      AddBuffer_SPU2(alpha, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle, ie);
#endif
    }
    simu->tnow += 0.5 * dt;

    DtFields_SPU(simu, NULL, NULL);

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
#ifndef WORKER_ON_NODE
      AddBuffer_SPU(simu->dt, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle);
#else
      AddBuffer_SPU2(simu->dt, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle, ie);
#endif
    }

    simu->tnow += 0.5 * dt;

    iter++;
    simu->iter_time_rk = iter;
  }

  starpu_task_wait_for_all() ;

  // top chrono
  UnregisterSimulation_SPU(simu);
}


int test_sn_spec_gaussian_cube_spu(void){
    int test = true;

#ifdef _WITH_OPENCL
    if (!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
        printf("OpenCL device not acceptable.\n");
        return true;
    }
#endif 

    putenv("STARPU_PROFILING=1");
    putenv("STARPU_WORKER_STATS=1");
    putenv("STARPU_SCHED=dmda");

    /* schnaps STARPU boolean*/
    starpu_use = true;
    starpu_c_use = true;
    starpu_ocl_use = true;
    starpu_cuda_use = false;

    /*Read and init mesh struct*/
    MacroMesh mesh;
    ReadMacroMesh(&mesh, "../test/testcube.msh");
    BuildConnectivity(&mesh);

    /*DG parameters*/
    int thedeg  = 2;
    int therafx = 8;
    int therafy = 8;
    int therafz = 8;
    int deg[] = {thedeg, thedeg, thedeg};
    int raf[] = {therafx, therafy, therafz};
    CheckMacroMesh(&mesh, deg, raf);


    /*Init velocity space*/
    schnaps_acoustic_data = lebedev_50 ;
    acoustic_data *ad = &schnaps_acoustic_data ;

    /* Model parameters */
    Model model;

    model.m            = 50;
    model.InitData     = ocl_s50_init_gaussian_density;
    model.ImposedData  = ocl_s50_imposed_gaussian_density;
    model.NumFlux      = ocl_s50_upwind_numflux;
    model.BoundaryFlux = ocl_test_sn_mixed_boundary_flux;
    model.Source       = NULL;

    /* Init simulation */
    Simulation simu;
    EmptySimulation(&simu);
    InitSimulation(&simu, &mesh, deg, raf, &model);
    simu.cfl = 0.95;
    schnaps_real tmax = 1.0;
    
    /* OpenCL compilation options */
    char buf[1000];
#ifndef _DOUBLE_PRECISION
    sprintf(buf, " -D schnaps_real=float -D _M=%d -cl-single-precision-constant",
    model.m);
#else
    sprintf(buf, " -D schnaps_real=double -D _M=%d", model.m);
#endif
    strcat(cl_buildoptions, buf);
    sprintf(numflux_cl_name, "%s", "ocl_s50_upwind_numflux");
    sprintf(buf," -D NUMFLUX=");
    strcat(buf, numflux_cl_name);
    strcat(cl_buildoptions, buf);

    sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math",
    "ocl_test_sn_mixed_boundary_flux");
    strcat(cl_buildoptions, buf);

    //SmartPrefetch_SPU(&simu);

    /* Solve */
    sn_rk2_spu(&simu, tmax);
    CopyfieldtoCPU(&simu);

    schnaps_real rec_rho=0;
    schnaps_real rcp[1];
    printf("rcp = %f\n",rcp[0]);
    starpu_shutdown();
    return test;

}

int main(void) {
  int resu1 = test_sn_spec_gaussian_cube_spu();
  if (resu1) printf("Acoustic SN test OK !\n");
  return !resu1;
}



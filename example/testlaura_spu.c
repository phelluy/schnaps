#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test/test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "mhd.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <../test/test.h>
#define _XOPEN_SOURCE 700



// Headers declaration =========================================================
int TestLaura_SPU(int argc, char *argv[]) ;
// =============================================================================

#include <starpu.h>
#include <starpu_heteroprio.h>

void init_heteroprio(unsigned sched_ctx) {
    printf("[HETEROPRIO] Init\n");
   // Create queues for CPU tasks
   starpu_heteroprio_set_nb_prios(sched_ctx, STARPU_CPU_IDX, SSPU_NB_PRIO);

   // Set lookup order for CPU workers
   // 0 => 3
   // 1 => 2
   // ..
   // 3 => 0
   // Use simple mapping
   for (int bucketid=0; bucketid<SSPU_NB_PRIO; bucketid++) {
      starpu_heteroprio_set_mapping(sched_ctx, STARPU_CPU_IDX, bucketid, bucketid);
      starpu_heteroprio_set_faster_arch(sched_ctx, STARPU_CPU_IDX, bucketid);
   }

   // Create queues for CUDA tasks
   starpu_heteroprio_set_nb_prios(sched_ctx, STARPU_OPENCL_IDX, SSPU_NB_PRIO);
   for (int bucketid=0; bucketid<SSPU_NB_PRIO; bucketid++) {
      starpu_heteroprio_set_mapping(sched_ctx, STARPU_OPENCL_IDX, bucketid, SSPU_NB_PRIO-bucketid-1);
   }

   starpu_heteroprio_set_faster_arch(sched_ctx, STARPU_OPENCL_IDX, DGMass_SPU_PRIO);
   starpu_heteroprio_set_arch_slow_factor(sched_ctx, STARPU_OPENCL_IDX, DGMass_SPU_PRIO, 10.0f);

   starpu_heteroprio_set_faster_arch(sched_ctx, STARPU_OPENCL_IDX, DGSource_SPU_PRIO);
   starpu_heteroprio_set_arch_slow_factor(sched_ctx, STARPU_OPENCL_IDX, DGSource_SPU_PRIO, 10.0f);

   starpu_heteroprio_set_faster_arch(sched_ctx, STARPU_OPENCL_IDX, DGVolume_SPU_PRIO);
   starpu_heteroprio_set_arch_slow_factor(sched_ctx, STARPU_OPENCL_IDX, DGVolume_SPU_PRIO, 10.0f);
}


// Main function ===============================================================
int main(int argc, char *argv[]) {
  int resu = TestLaura_SPU(argc,argv);
  if (resu)
    printf("Laura test OK !\n");
  else
    printf("Laura test failed !\n");
  return !resu;
}
// =============================================================================



// Test function ===============================================================
int TestLaura_SPU(int argc, char *argv[]) {

  int test = true;

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif

  // Environmental variables................................
  // putenv("STARPU_PREFETCH=1");
  //putenv("STARPU_SCHED=dmdar");
  putenv("STARPU_PROFILING=1");
  putenv("STARPU_WORKER_STATS=1");
  //putenv("STARPU_NCPU=6");
  //putenv("STARPU_NOPENCL=0");
  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = true;
  //starpu_ocl_use = false;
  starpu_cuda_use = false;
  //........................................................


  // Mesh definition........................................
  MacroMesh mesh;
  // char *mshname =  "../geo/laura_cube.msh"; //TODO: ajouter le msh
  // char *mshname =  "../test/testcube.msh";
  char *mshname =  "../example/test_torus_400.msh";

  printf("\n\n---------------------------------------------------------\n");
  printf("1/ Reading Macromesh\n");
  ReadMacroMesh(&mesh, mshname);

  // Detect2DMacroMesh(&mesh);
  // bool is2d=mesh.is2d;
  // assert(!is2d);

  printf("\n\n---------------------------------------------------------\n");
  printf("2/ Building connectivity \n");
  BuildConnectivity(&mesh);
  int thedeg = 4;
  int theraf = 5;
  int deg[] = {thedeg, thedeg, thedeg};
  int raf[] = {theraf, theraf, theraf};

  printf("\n\n---------------------------------------------------------\n");
  printf("3/ Checking Macromesh\n");
  printf("      - Deg = {%d,%d,%d}\n", deg[0],deg[1],deg[2]);
  printf("      - Raf = {%d,%d,%d}\n", raf[0],raf[1],raf[2]);
  //CheckMacroMesh(&mesh, deg, raf);
  printf("---------------------------------------------------------\n");

  // schnaps_real periodsize = 6.2831853;
  schnaps_real periodsize = -1;
  // mesh.period[0]=periodsize;
  // mesh.period[1]=periodsize;
  // // mesh.period[2]=periodsize;
  //........................................................


  //........................................................
  schnaps_real cfl = 0.48*theraf/12;
  schnaps_real vcfl = 1.0;
  schnaps_real tmax = 0.01;
  //schnaps_real tmax = 0.02103;
  schnaps_real vmax = 1.0;
  //........................................................

  ////////////////////

  Model model;

  Simulation simu;
  EmptySimulation(&simu);

  simu.cfl = cfl;
  model.m = 8;

  strcpy(model.name,"Maxwell3DNumFluxClean_upwind");

  model.NumFlux      = Maxwell3DNumFluxClean_upwind;
  model.BoundaryFlux = Maxwell3DBoundaryFluxClean_upwind;
  model.InitData     = Maxwell3DInitData;
  model.ImposedData  = Maxwell3DImposedData;
  model.Source = NULL;

  /* model.m = 1; */
  /* model.NumFlux = TransNumFlux; */
  /* model.BoundaryFlux = TestTransBoundaryFlux; */
  /* model.InitData = TestTransInitData; */
  /* model.ImposedData = TestTransImposedData; */
  /* model.Source = NULL; */

  char buf[1000];
#ifndef _DOUBLE_PRECISION
  sprintf(buf, " -D schnaps_real=float -D _M=%d -cl-single-precision-constant",
	  model.m);
#else
  sprintf(buf, " -D schnaps_real=double -D _M=%d", model.m);
#endif
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "Maxwell3DNumFluxClean_upwind");
  //sprintf(numflux_cl_name, "%s", "TransNumFlux");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math",
	  "Maxwell3DBoundaryFluxClean_upwind");
  /* sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math", */
  /* 	  "TestTransBoundaryFlux"); */
  strcat(cl_buildoptions, buf);

  // Interesting environment variables
  printf("\n\n---------------------------------------------------------\n");
  printf("Environment variables...\n");
  printf("STARPU_SCHED               : %s\n", getenv("STARPU_SCHED"));
  printf("STARPU_WORKER_STATS        : %s\n", getenv("STARPU_WORKER_STATS"));
  printf("STARPU_PROFILING           : %s\n", getenv("STARPU_PROFILING"));
  printf("STARPU_NCPU                : %s\n", getenv("STARPU_NCPU"));
  printf("STARPU_NOPENCL             : %s\n", getenv("STARPU_NOPENCL"));
  printf("STARPU_NCUDA               : %s\n", getenv("STARPU_NCUDA"));
  printf("STARPU_NMIC                : %s\n", getenv("STARPU_NMIC"));
  printf("STARPU_OPENCL_ON_CPUS      : %s\n", getenv("STARPU_OPENCL_ON_CPUS"));
  printf("STARPU_OPENCL_ONLY_ON_CPUS : %s\n",
	 getenv("STARPU_OPENCL_ONLY_ON_CPUS"));

  if(starpu_is_init){
      destroy_global_arbiter();
      starpu_shutdown();
      starpu_is_init = false;
  }

  if (!starpu_is_init && starpu_use){
      struct starpu_conf conf;
      starpu_conf_init(&conf);

      conf.sched_policy_name = "heteroprio";
      conf.sched_policy_init = &init_heteroprio;

    int ret;
    ret = starpu_init(&conf);
    assert(ret != -ENODEV) ;
    starpu_is_init = true;
  }
  init_global_arbiter();

  
  // Workers
  const int nb_workers = starpu_worker_get_count();
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("Available workers...\n");
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  printf("OPENCL workers             : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());

  printf("\n\n---------------------------------------------------------\n");
  printf("4/ Initializing simulation\n");
  InitSimulation(&simu, &mesh, deg, raf, &model);



  //SmartPrefetch_SPU(&simu);


  simu.vmax = vmax;
  simu.vcfl = vcfl;

  //PlotFields(0, false, &simu, NULL, "init_dgvisu.msh");

  schnaps_real tps_deb = seconds();
  printf("\n\n---------------------------------------------------------\n");
  printf("5/ Solving\n");
  RK2_SPU(&simu, tmax);

  schnaps_real tps_fin = seconds();
  printf("temps RK2 total= %f\n", tps_fin-tps_deb);

  //CopyfieldtoCPU(&simu);

  // show_cl_timing(&simu);

  /* schnaps_real dd = L2error(&simu); */
  /* printf("erreur L2=%f\n", dd); */

  //PlotFields(0, false, &simu, NULL, "fin_dgvisu.msh");

  destroy_global_arbiter();
  starpu_shutdown(); /* Necessary to ensure perfmodels are written to disk */

  return test;
}
// =============================================================================

/*

The variable STARPU_SCHED can be set to one of the following strings:
modular-eager                 	-> eager modular policy
modular-eager-prefetching     	-> eager with prefetching modular policy
modular-prio                  	-> prio modular policy
modular-prio-prefetching      	-> prio prefetching modular policy
modular-random                	-> random modular policy
modular-random-prio           	-> random-prio modular policy
modular-random-prefetching    	-> random prefetching modular policy
modular-random-prio-prefetching	-> random-prio prefetching modular policy
modular-heft                  	-> heft modular policy
modular-heft-prio             	-> heft+prio modular policy
modular-heft2                 	-> heft modular2 policy
eager                         	-> eager policy with a central queue
prio                          	-> eager (with priorities)
random                        	-> weighted random based on worker overall performance
lws                           	-> locality work stealing
ws                            	-> work stealing
dm                            	-> performance model
dmda                          	-> data-aware performance model
dmdar                         	-> data-aware performance model (ready)
dmdas                         	-> data-aware performance model (sorted)
dmdasd                        	-> data-aware performance model (sorted decision)
pheft                         	-> parallel HEFT
peager                        	-> parallel eager policy
graph_test                    	-> test policy for using graphs in scheduling decisions
*/

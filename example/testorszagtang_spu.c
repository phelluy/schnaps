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
#define _XOPEN_SOURCE 700


static int
find_a_worker(enum starpu_worker_archtype type)
{
	int worker[STARPU_NMAXWORKERS];
	int ret = starpu_worker_get_ids_by_type(type, worker, STARPU_NMAXWORKERS);
	if (ret == 0)
		return -ENODEV;
	if (ret == -ERANGE)
		return worker[STARPU_NMAXWORKERS-1];
	return worker[ret-1];
}


int TestOrszagTang_SPU(int argc, char *argv[]) ;


int main(int argc, char *argv[]) {
  int resu = TestOrszagTang_SPU(argc,argv);
  if (resu)
    printf("OrszagTang test OK !\n");
  else 
    printf("OrszagTang test failed !\n");
  return !resu;
}

int TestOrszagTang_SPU(int argc, char *argv[]) {

  int test = true;

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif

  putenv("STARPU_SCHED=dmdar");
  putenv("STARPU_SCHED_BETA=10");
  putenv("STARPU_SCHED_ALPHA=1");

  // disable OPENCL for StarPU
  putenv("STARPU_NOPENCL=0");
  putenv("STARPU_NCPU=12");
  putenv("STARPU_PREFETCH=0");
    
  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = true;
  starpu_cuda_use = true;

  putenv("STARPU_NCUDA=0");
  putenv("STARPU_OPENCL_ON_CPU=0");
  putenv("STARPU_OPENCL_ONLY_ON_CPU=0");

  
  ////////////////////
  // init mesh
  MacroMesh mesh;
  char *mshname =  "../test/testOTgrid.msh";
  //char *mshname =  "../geo/grid.msh";
  
  ReadMacroMesh(&mesh, mshname);
  Detect2DMacroMesh(&mesh);
  bool is2d=mesh.is2d; 
  assert(is2d);  

  schnaps_real periodsize = 6.2831853;
  //schnaps_real periodsize = -1;
  mesh.period[0]=periodsize;
  mesh.period[1]=periodsize;

  BuildConnectivity(&mesh);

  int deg[]={1, 1, 0};
  int raf[]={50, 50, 1};
  CheckMacroMesh(&mesh, deg, raf);

  ////////////////////


  
  Model model;

  Simulation simu;
  EmptySimulation(&simu);

  schnaps_real cfl = 0.1;
  simu.cfl = cfl;
  model.m = 9;

  strcpy(model.name,"MHD");

  model.NumFlux=MHDNumFluxRusanov;
  model.BoundaryFlux=MHDBoundaryFluxOrszagTang;
  model.InitData=MHDInitDataOrszagTang;
  model.ImposedData=MHDImposedDataOrszagTang;
  model.Source = NULL;
  
  char buf[1000];
  sprintf(buf, " -D schnaps_real=double -D _M=%d -D _PERIODX=%f -D _PERIODY=%f",
          model.m,
          periodsize,
          periodsize);
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "MHDNumFluxRusanov");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math", "MHDBoundaryFluxOrszagTang");
  strcat(cl_buildoptions, buf);

  // Interesting environment variables
  printf("Environment variables...\n");
  printf("STARPU_NCPU                : %s\n", getenv("STARPU_NCPU"));
  printf("STARPU_NOPENCL             : %s\n", getenv("STARPU_NOPENCL"));
  printf("STARPU_NCUDA               : %s\n", getenv("STARPU_NCUDA"));
  printf("STARPU_NMIC                : %s\n", getenv("STARPU_NMIC"));
  printf("STARPU_OPENCL_ON_CPUS      : %s\n", getenv("STARPU_OPENCL_ON_CPUS"));
  printf("STARPU_OPENCL_ONLY_ON_CPUS : %s\n", getenv("STARPU_OPENCL_ONLY_ON_CPUS"));

  InitSimulation(&simu, &mesh, deg, raf, &model);

  // Workers
  const int nb_workers = starpu_worker_get_count();
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("Available workers...\n");
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  printf("OPENCL workers             : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());

  //SmartPrefetch_SPU(&simu);

  schnaps_real tmax = 10.;
  simu.vmax = 6.0;
  simu.vcfl = 6.0;
  schnaps_real dt = 0;

  schnaps_real tps_deb = seconds();
  RK2_SPU(&simu, tmax);

  schnaps_real tps_fin = seconds();
  printf("temps = %f\n", tps_fin-tps_deb);

  //CopyfieldtoCPU(&simu);

  //show_cl_timing(&simu);
  PlotFields(0, false, &simu, NULL, "dgvisu.msh");

  destroy_global_arbiter();
  starpu_shutdown(); /* Necessary to ensure perfmodels are written to disk */

  return test;
}

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




// schnaps_real seconds()
// {
//   struct timespec ts;
//   clock_gettime(CLOCK_MONOTONIC, &ts);
//   return (schnaps_real)ts.tv_sec + 1e-9 * (schnaps_real)ts.tv_nsec;
// }

int TestRiemann3D_SPU(int argc, char *argv[]) ;


int main(int argc, char *argv[]) {
  int resu = TestRiemann3D_SPU(argc,argv);
  if (resu)
    printf("Riemann3D test OK !\n");
  else 
    printf("Riemann3D test failed !\n");
  return !resu;
}

int TestRiemann3D_SPU(int argc, char *argv[]) {

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
  //putenv("STARPU_NCPU=0");
  putenv("STARPU_PREFETCH=1");
    
  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = false;
  starpu_cuda_use = false;

  //putenv("STARPU_NCUDA=0");
  //putenv("STARPU_OPENCL_ON_CPU=1");
  //putenv("STARPU_OPENCL_ONLY_ON_CPU=1");

  
  ////////////////////
  // init mesh
  MacroMesh mesh;
  //char *mshname =  "../test/testOTgrid.msh";
  //char *mshname =  "../test/testcartesiangrid1d.msh";
  char *mshname =  "../geo/grid.msh";
  
  ReadMacroMesh(&mesh, mshname);

//  Detect2DMacroMesh(&mesh);
//  bool is2d=mesh.is2d; 
//  assert(is2d);  

  
  schnaps_real periodsize = 10.;
  //schnaps_real periodsize = -1;
  //mesh.period[0]=periodsize;
  //mesh.period[1]=periodsize;
  //mesh.period[2]=periodsize;

  BuildConnectivity(&mesh);

  int deg[]={1, 1, 1};
  int raf[]={8, 8, 8};
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
  model.BoundaryFlux=MHDBoundaryFluxChocFort;
  model.InitData=MHDInitDataChocFort;
  model.ImposedData=MHDImposedDataChocFort;
  model.Source = NULL;
  
  char buf[1000];
  sprintf(buf, " -D schnaps_real=double -D _M=%d -D _PERIODY=%f -D _PERIODZ=%f",
          model.m,
          periodsize,
          periodsize);//,
          //          periodsize,
          //          periodsize);
  //  sprintf(buf, " -D schnaps_real=double -D _M=%d ", model.m);

  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "MHDNumFluxRusanov");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math", "MHDBoundaryFluxChocFort");
  strcat(cl_buildoptions, buf);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  //SmartPrefetch_SPU(&simu);
  
  schnaps_real tmax = 0.05;
  simu.vmax = 6.0;
  simu.vcfl = 6.0;
  schnaps_real dt = 0;
  
  // schnaps_real tps_deb = seconds();
  RK2_SPU(&simu, tmax);

  // schnaps_real tps_fin = seconds();
  // printf("temps = %f\n", tps_fin-tps_deb);
  
  //CopyfieldtoCPU(&simu);
 
  //show_cl_timing(&simu);
  PlotFields(0, false, &simu, NULL, "dgvisu.msh");

  starpu_shutdown(); /* Necessary to ensure perfmodels are written to disk */
  
  return test;
}

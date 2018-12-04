#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>



int test_user_model_spu(void){
  bool test = true;

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh


  //putenv("STARPU_NOPENCL=0");


  char *mshname =  "../test/disque2d.msh";

  MacroMesh mesh;
  //ReadMacroMesh(&mesh,"../test/testdisque.msh");
  ReadMacroMesh(&mesh,"../test/testcube2.msh");
  //Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  Model model;

  Simulation simu;
  EmptySimulation(&simu);

  char buf[1000];

  model.m = 1;
  schnaps_model_load("../test/model_trans.c", &model);

#ifdef _DOUBLE_PRECISION
  sprintf(buf, "-D schnaps_real=double -D _M=%d", model.m);
#else
  sprintf(buf, "-D schnaps_real=float -D _M=%d -cl-single-precision-constant", model.m);
#endif
  strcat(cl_buildoptions, buf);


  ////////////////// useless in the last version
  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TestTransBoundaryFlux");
  strcat(cl_buildoptions, buf);
  ////////////////////////////  useless in the last version

  int deg[]={2, 2, 2};
  int raf[]={3, 3, 3};
  //int raf[]={2, 2, 2};


  // 2015-01-19: the below parameters fail with testmacrocellinterface
  // but pass the test here (perhaps because the error is hidden by
  // the RK error?)
  /*
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement
  */

  //AffineMapMacroMesh(&(f.macromesh));

#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif

  CheckMacroMesh(&mesh, deg, raf);
  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = true;
  starpu_cuda_use = false;

  // Interesting environment variables
  printf("Environment variables...\n");
  printf("STARPU_NCPU                : %s\n", getenv("STARPU_NCPU"));
  printf("STARPU_NOPENCL             : %s\n", getenv("STARPU_NOPENCL"));
  printf("STARPU_NCUDA               : %s\n", getenv("STARPU_NCUDA"));
  printf("STARPU_NMIC                : %s\n", getenv("STARPU_NMIC"));
  printf("STARPU_OPENCL_ON_CPUS      : %s\n", getenv("STARPU_OPENCL_ON_CPUS"));
  printf("STARPU_OPENCL_ONLY_ON_CPUS : %s\n", getenv("STARPU_OPENCL_ONLY_ON_CPUS"));


  InitSimulation(&simu, &mesh, deg, raf, &model);
  // PlotFieldsBinSparse(0, false, &simu, NULL, "TESTTEST.msh");
  // Workers
  const int nb_workers = starpu_worker_get_count();
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("Available workers...\n");
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  printf("OPENCL workers             : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());


  schnaps_real tmax = .1;
  simu.cfl=0.1;
  simu.vmax=1;
  RK2_SPU(&simu,tmax);


  PlotFields(0, false, &simu, NULL, "dgvisu.msh");
  // PlotFields(0, true , &simu, "error", "dgerror.msh");

  schnaps_real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%.12f\n", dd);

//PlotFieldsBinSparse(0, false, &simu, NULL, "dgvisu.msh");
//PlotFieldsBinSparse(0, true , &simu, "error", "dgerror.msh");
  // pour tracer le graphe des tâches:
  // avant l'execution: export STARPU_FXT_PREFIX=~/schnaps/build/
  //  starpu_fxt_tool -i prof_file_*
  // dot -Tpdf dag.dot -o output.pdf
  // pour tracer le diagramme de Gant:
  // vite paje.trace
  // autres commandes utiles
  // export STARPU_GENERATE_TRACE=0 ou 1
  // voir aussi http://starpu.gforge.inria.fr/doc/html/ExecutionConfigurationThroughEnvironmentVariables.html
  // export STARPU_SCHED=dmda

schnaps_real tolerance = 0.02;

  test = dd < tolerance;

  FreeMacroMesh(&mesh);

  return test;
}

int main(void) {
  int resu = test_user_model_spu();
  if(resu)
    printf("test_user_model_spu OK !\n");
  else
    printf("test_user_model_spu failed !\n");
  return !resu;
}

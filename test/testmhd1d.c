#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "mhd.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#define _XOPEN_SOURCE 700

//real seconds()
//{
//  struct timespec ts;
//  clock_gettime(CLOCK_MONOTONIC, &ts);
//  return (real)ts.tv_sec + 1e-9 * (real)ts.tv_nsec;
//}

int main(int argc, char *argv[]) {
  int resu = TestMHD1D(argc,argv);
  if (resu)
    printf("MHD 1D test OK !\n");
  else 
    printf("MHD 1D test failed !\n");
  return !resu;
}

int TestMHD1D(int argc, char *argv[]) {
  real cfl = 0.2;
  real tmax = 1.0;
  bool writemsh = false;
  real vmax = 6.0;
  bool usegpu = false;
  real dt = 0.0;

  for (;;) {
    int cc = getopt(argc, argv, "c:t:w:D:P:g:s:");
    if (cc == -1) break;
    switch (cc) {
    case 0:
      break;
    case 'c':
      cfl = atof(optarg);
      break;
    case 'g':
      usegpu = atoi(optarg);
      break;
    case 't':
      tmax = atof(optarg);
      break;
    case 'w':
      writemsh = true;
      break;
    case 'D':
       ndevice_cl= atoi(optarg);
      break;
    case 'P':
      nplatform_cl = atoi(optarg);
      break;
    default:
      printf("Error: invalid option.\n");
      printf("Usage:\n");
      printf("./testmanyv -c <cfl> -d <deg> -n <nraf> -t <tmax> -C\n -P <cl platform number> -D <cl device number> FIXME");
      exit(1);
    }
  }

  bool test = true;
  field f;
  init_empty_field(&f);  

  f.varindex = GenericVarindex;
  f.model.m = 9;
  f.model.cfl = cfl;

  strcpy(f.model.name,"MHD");

  f.model.NumFlux=MHDNumFluxRusanov;
  f.model.BoundaryFlux=MHDBoundaryFlux;
  f.model.InitData=MHDInitData;
  f.model.ImposedData=MHDImposedData;
  
  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "MHDNumFluxRusanov");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "MHDBoundaryFlux");
  strcat(cl_buildoptions, buf);
  
  // Set the global parameters for the Vlasov equation
  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = 1; // x direction degree
  f.interp.interp_param[2] = 0; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 10; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement


  //set_vlasov_params(&(f.model));

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testcartesiangrid1d.msh");
  //ReadMacroMesh(&(f.macromesh), "test/testcube.msh");
  // Try to detect a 2d mesh
  Detect1DMacroMesh(&(f.macromesh));
  bool is1d=f.macromesh.is1d; 
  assert(is1d);  

  f.macromesh.period[0]=10;
  
  // Mesh preparation
  BuildConnectivity(&(f.macromesh));

  // Prepare the initial fields
  Initfield(&f);

  
  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  Plotfield(0, (1==0), &f, "Rho", "dginit.msh");

  f.vmax=vmax;

  if(dt <= 0)
    dt = set_dt(&f);

  real executiontime;
  if(usegpu) {
    printf("Using OpenCL:\n");
    //executiontime = seconds();
    assert(1==2);
    RK2(&f, tmax, dt);
    //executiontime = seconds() - executiontime;
  } else { 
    printf("Using C:\n");
    //executiontime = seconds();
    RK2(&f, tmax, dt);
    //executiontime = seconds() - executiontime;
  }

  Plotfield(0,false,&f, "Rho", "dgvisu.msh");
  Gnuplot(&f,0,0.0,"data1D.dat");

  printf("tmax: %f, cfl: %f\n", tmax, f.model.cfl);

  printf("deltax:\n");
  printf("%f\n", f.hmin);

  printf("deltat:\n");
  printf("%f\n", dt);

  printf("DOF:\n");
  printf("%d\n", f.wsize);

  printf("executiontime (s):\n");
  printf("%f\n", executiontime);

  printf("time per RK2 (s):\n");
  printf("%f\n", executiontime / (real)f.itermax);

  return test;
}

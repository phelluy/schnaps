#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "model.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#define _XOPEN_SOURCE 700

/* real seconds(struct timespec ta, struct timespec tb) { */
/*   return (real)(ta.tv_sec - tb.tv_sec)  */
/*     + 1e-9 * (real)(ta.tv_nsec - tb.tv_nsec); */
/* } */

int main(int argc, char *argv[]) {
  // Unit tests
  real cfl = 0.5;
  int deg = 3;
  int nx = 4;
  int ny = 4;
  real tmax = 1e-2;
  bool writemsh = false;
  real vmax = 1.0;
  int mx = 5;
  int my = 5;
  bool usegpu = false;
  real dt = 0.0;
  char *usage = "./testmanyv \
\n\t-c <float> set CFL\
\n\t-d <int> set interpolation degree\
\n\t-n <int> set number of subcells\
\n\t-g <1 or 0> GPU(1) instead of CPU(0) \
\n\t-t <float> set tmax\
\n\t-w write mesh output \
\n\t-P <int> set OpenCL platform number \
\n\t-D <int> set OpenCL device number \
\n\t-X <int> set number of velocities in x direction \
\n\t-Y <int> set number of velocities in y direction \
\n\t-x <int> set number of subcells in x direction \
\n\t-y <int> set number of subcells in y direction \
\n\t-s <float> set dt \
\n\t-h display usage information \n";

  for (;;) {
    int cc = getopt(argc, argv, "c:d:t:wD:P:X:Y:x:y:V:g:s:h");
    if (cc == -1) break;
    switch (cc) {
    case 0:
      break;
    case 'c':
      cfl = atof(optarg);
      break;
    case 'd':
      deg = atoi(optarg);
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
    case 'X':
      mx = atoi(optarg);
      break;
    case 'Y':
      my = atoi(optarg);
      break;
    case 'x':
      nx = atoi(optarg);
      break;
    case 'y':
      ny = atoi(optarg);
      break;
    case 'V':
      vmax = atof(optarg);
      break;
    case 's':
      dt = atof(optarg);
      break;
    case 'h':
      printf("Usage:\n%s", usage);
      exit(0);
      break;
    default:
      printf("Error: invalid option.\n");
      printf("Usage:\n%s", usage);
      exit(1);
    }
  }

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return 0;
  }

  bool test = true;
  field f;
  init_empty_field(&f);

  f.varindex = GenericVarindex;
  f.model.vlasov_mz = 1;
  f.model.cfl = cfl;
  f.model.NumFlux = vlaTransNumFlux2d;
  f.model.vlasov_mx = mx;
  f.model.vlasov_my = my;
  f.model.vlasov_vmax = vmax;
  f.model.m = f.model.vlasov_mx * f.model.vlasov_my * f.model.vlasov_mz;

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "vlaTransNumFlux2d");
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D vlasov_mx=%d", f.model.vlasov_mx);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D vlasov_my=%d", f.model.vlasov_my);
  strcat(cl_buildoptions, buf);
  sprintf(buf, " -D vlasov_vmax=%f", f.model.vlasov_vmax);
  strcat(cl_buildoptions, buf);

  f.model.BoundaryFlux = cemracs2014_TransBoundaryFlux;
  f.model.InitData = cemracs2014_TransInitData;
  f.model.ImposedData = cemcracs2014_imposed_data;
  sprintf(buf, " -D BOUNDARYFLUX=%s", "cemracs2014_TransBoundaryFlux");
  strcat(cl_buildoptions, buf);

  // Set the global parameters for the Vlasov equation
  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = deg; // x direction degree
  f.interp.interp_param[2] = deg; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = nx; // x direction refinement
  f.interp.interp_param[5] = ny; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  set_vlasov_params(&(f.model));

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "geo/square.msh");
  // Try to detect a 2d mesh
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);  

  // Mesh preparation
  BuildConnectivity(&(f.macromesh));
 
  // Prepare the initial fields

  Initfield(&f);
  f.vmax = f.model.vlasov_vmax;

  if(dt <= 0.0)
    dt = set_dt(&f);
  printf("dt: %f\n", dt);

  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  /* real executiontime; */
  /* struct timespec tstart, tend; */
  if(usegpu) {

    printf("Using OpenCL:\n");
    //clock_gettime(CLOCK_MONOTONIC, &tstart);
    RK2_CL(&f, tmax, dt,  0, NULL, NULL);
    //clock_gettime(CLOCK_MONOTONIC, &tend);

    CopyfieldtoCPU(&f);

    printf("\nOpenCL Kernel time:\n");
    show_cl_timing(&f);
    printf("\n");

  } else {
    printf("Using C:\n");
    //clock_gettime(CLOCK_MONOTONIC, &tstart);
    RK2(&f, tmax, dt);
    //clock_gettime(CLOCK_MONOTONIC, &tend);
  }
  /* executiontime = seconds(tend, tstart); */

  // Save the results and the error
  if(writemsh) {
    for(int ix = 0; ix < f.model.vlasov_mx; ++ix) {
      for(int iy = 0; iy < f.model.vlasov_my; ++iy) {
	int mplot = ix * f.model.vlasov_my + iy; 

	real vx = vlasov_vel(ix, f.model.vlasov_mx, f.model.vlasov_vmax);
	real vy = vlasov_vel(iy, f.model.vlasov_my, f.model.vlasov_vmax);
	char fieldname[100];
	sprintf(fieldname, "output field has v = (%f,%f)", vx, vy);
	//printf("%s\n", fieldname);
      
	char filename[100];
	sprintf(filename, "dgvisuix%diy%d.msh", ix, iy);
	//printf("ix: %d, iy: %d, fieldname: %s\n", ix, iy, fieldname);
	Plotfield(mplot, false, &f, fieldname, filename);
      }
    }
    /* Plotfield(mplot, true, &f, "dgerror.msh"); */
  }

  printf("tmax: %f, cfl: %f, deg: %d, nrafx: %d, nrafy: %d\n", 
	 tmax, f.model.cfl, deg, nx, ny);
  real dd = L2error(&f) / (f.model.vlasov_mx * f.model.vlasov_my);

  printf("deltax:\n");
  printf("%f\n", f.hmin);

  printf("deltat:\n");
  printf("%f\n", dt);

  printf("DOF:\n");
  printf("%d\n", f.wsize);

  printf("L2 error:\n");
  printf("%e\n", dd);

  /*
  printf("executiontime (s):\n");
  printf("%f\n", executiontime);
  */

  printf("itermax (s):\n");
  printf("%d\n", f.itermax);

  /*
  printf("perRK2time (s):\n");
  printf("%f\n", executiontime / (real)f.itermax);
  */

  real tolerance = 1e-2;
  test = test && (dd < tolerance);

  if(test) 
    printf("multiple velocity transport test OK !\n");
  else 
    printf("multiple velocity transport test failed !\n");
  return !test;
}

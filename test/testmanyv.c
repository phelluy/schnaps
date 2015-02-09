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

double seconds(struct timespec ta, struct timespec tb) {
  return (double)(ta.tv_sec - tb.tv_sec) 
    + 1e-9 * (double)(ta.tv_nsec - tb.tv_nsec);
}

int main(int argc, char *argv[]) {
  // Unit tests
  double cfl = 0.5;
  int deg = 3;
  int nraf = 4;
  double tmax = 1e-2;
  int cemracs = 0;
  bool writemsh = false;
  double vmax = 1.0;
  int mx = 5;
  int my = 5;
  bool usegpu = false;
  double dt = 0.0;
  for (;;) {
    int cc = getopt(argc, argv, "c:d:n:t:C:wD:P:X:Y:V:g:s:");
    if (cc == -1) break;
    switch (cc) {
    case 0:
      break;
    case 'C':
      cemracs = true;
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
    case 'n':
      nraf = atoi(optarg);
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
    case 'V':
      vmax = atof(optarg);
      break;
    case 's':
      dt = atof(optarg);
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

  if(cemracs > 0) {
    f.model.BoundaryFlux = cemracs2014_TransBoundaryFlux;
    if(cemracs == 1) {
      f.model.InitData = cemracs2014_TransInitData;
      f.model.ImposedData = cemcracs2014_imposed_data;
    }
    if(cemracs == 2) {
      f.model.InitData = cemracs2014a_TransInitData;
      f.model.ImposedData = cemcracs2014a_imposed_data;
    }

    sprintf(buf, " -D BOUNDARYFLUX=%s", "cemracs2014_TransBoundaryFlux");
    strcat(cl_buildoptions, buf);
  } else {
    // FIXME: set boundary flux.
    f.model.BoundaryFlux = vlaTransBoundaryFlux2d;
    f.model.InitData = vlaTransInitData2d;
    f.model.ImposedData = vlaTransImposedData2d;
  }

  // Set the global parameters for the Vlasov equation
  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = deg; // x direction degree
  f.interp.interp_param[2] = deg; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = nraf; // x direction refinement
  f.interp.interp_param[5] = nraf; // y direction refinement
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
  if(dt != 0.0)
    f.dt = dt;

  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  double executiontime;
  struct timespec tstart, tend;
  if(usegpu) {
    printf("Using OpenCL:\n");
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    RK2_CL(&f, tmax);
    clock_gettime(CLOCK_MONOTONIC, &tend);

  } else { 
    printf("Using C:\n");
    clock_gettime(CLOCK_MONOTONIC, &tstart);
    RK2(&f, tmax);
    clock_gettime(CLOCK_MONOTONIC, &tend);
  }
  executiontime = seconds(tend, tstart);

  // Save the results and the error
  if(writemsh) {
    for(int ix = 0; ix < f.model.vlasov_mx; ++ix) {
      for(int iy = 0; iy < f.model.vlasov_my; ++iy) {
	int mplot = ix * f.model.vlasov_my + iy; 

	double vx = vlasov_vel(ix, f.model.vlasov_mx, f.model.vlasov_vmax);
	double vy = vlasov_vel(iy, f.model.vlasov_my, f.model.vlasov_vmax);
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

  printf("tmax: %f, cfl: %f, deg: %d, nraf: %d\n", tmax, f.model.cfl, deg, nraf);
  double dd = L2error(&f) / (f.model.vlasov_mx * f.model.vlasov_my);

  printf("deltax:\n");
  printf("%f\n", f.hmin);

  printf("deltat:\n");
  printf("%f\n", f.dt);

  printf("DOF:\n");
  printf("%d\n", f.wsize);

  printf("L2 error:\n");
  printf("%e\n", dd);

  printf("executiontime (s):\n");
  printf("%f\n", executiontime);

  printf("itermax (s):\n");
  printf("%d\n", f.itermax);

  printf("perRK2time (s):\n");
  printf("%f\n", executiontime / (double)f.itermax);

  double tolerance = 1e-2;
  test = test && (dd < tolerance);

  if(test) 
    printf("multiple velocity transport test OK !\n");
  else 
    printf("multiple velocity transport test failed !\n");
  return !test;
}

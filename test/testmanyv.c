#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "model.h"
#include <math.h>
#include <string.h>


int main(int argc, char* argv[]) {
  // Unit tests
  double cfl = 0.05;
  int deg = 3;
  int nraf = 4;
  double tmax = 1e-2;
  bool cemracs = false;
  bool writemsh = false;
  double vmax = 1.0;
  int mx = 5;
  int my = 5;
  for (;;) {
    int cc = getopt(argc, argv, "c:d:n:t:CwD:P:X:Y:V:");
    if (cc == -1) break;
    switch (cc) {
    case 0:
      break;
    case 'C':
      cemracs = true;
      break;
    case 'c':
      cfl = atof(optarg);
      // set cfl
      break;
    case 'd':
      deg = atoi(optarg);
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
    default:
      printf("Error: invalid option.\n");
      printf("Usage:\n");
      printf("./testmanyv -c <cfl> -d <deg> -n <nraf> -t <tmax> -C\n -P <cl platform number> -D <cl device number>");
      exit(1);
    }
  }

  bool test = true;
  Field f;
  
  f.varindex = GenericVarindex;
  f.model.vlasov_mz = 1;
  f.model.cfl = 0.05;
  f.model.NumFlux = vlaTransNumFlux2d;
  f.model.vlasov_mx = mx;
  f.model.vlasov_my = my;
  f.model.vlasov_vmax = vmax;
  f.model.m = f.model.vlasov_mx * f.model.vlasov_my * f.model.vlasov_mz;

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);


  if(cemracs) {
    f.model.BoundaryFlux = cemracs2014_TransBoundaryFlux;
    f.model.InitData = cemracs2014_TransInitData;
    f.model.ImposedData = cemcracs2014_imposed_data;

    sprintf(numflux_cl_name, "%s", "vlaTransNumFlux2d");
    strcat(buf," -D NUMFLUX=");
    strcat(buf, numflux_cl_name);

    sprintf(buf, " -D vlasov_mx=%d",  f.model.vlasov_mx);
    strcat(cl_buildoptions, buf);
    sprintf(buf, " -D vlasov_my=%d",  f.model.vlasov_my);
    strcat(cl_buildoptions, buf);
    sprintf(buf, " -D vlasov_vmax=%f",  f.model.vlasov_vmax);
    strcat(cl_buildoptions, buf);
    
  } else {
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
  InitField(&f);

  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param: %f\n", f.hmin);

  RK2_CL(&f, tmax);
 
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
	PlotField(mplot, false, &f, fieldname, filename);
      }
    }
    /* PlotField(mplot, true, &f, "dgerror.msh"); */
  }

  printf("tmax: %f, cfl: %f, deg: %d, nraf: %d\n", tmax, cfl, deg, nraf);
  double dd = L2error(&f) / (f.model.vlasov_mx * f.model.vlasov_my);
  printf("L2 error:\n");
  printf("%e\n", dd);
  
  test = test && (dd < 1e-2);

  if(test) 
    printf("multiple velocity transport test OK !\n");
  else 
    printf("multiple velocity transport test failed !\n");
  return !test;
}

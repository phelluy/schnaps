#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */
#include "model.h"
#include <math.h>

void cemcracs2014_imposed_data(double x[3], double t, double *w)
{
  double PI = 4.0 * atan(1.0);
  double s2pi = sqrt(2.0 * PI);
  double xval = 1.0;

  double sigma = 1.0;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    double vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    double px = x[0] - vx * t;


    for(int iy = 0; iy < vlasov_my; ++iy) {
      double vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      double py = x[1] - vy * t;

      double r = sqrt(px * px + py * py);
      double pr = compact_bump(r);
      
      double vr = sqrt(vx * vx + vy * vy);
      double pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemracs2014_TransInitData(double x[3], double w[]) 
{
  double t = 0;
  cemcracs2014_imposed_data(x, t, w);
}

void cemracs2014_TransBoundaryFlux(double x[3], double t, 
			    double wL[], double *vnorm,
			    double* flux) 
{
  double wR[m];
  for(unsigned int i = 0; i < m; ++i) 
    wR[i] = 0;
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

int main(int argc, char* argv[]) {
  // Unit tests
  double cfl = 0.05;
  int deg = 3;
  int nraf = 4;
  double tmax = 1e-2;
  bool cemracs = false;
  bool writemsh = false;
  for (;;) {
    int cc = getopt(argc, argv, "c:d:n:t:CwD:P:");
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
  if(cemracs) {
    f.model.vlasov_mx = 64;
    f.model.vlasov_my = 64;
    f.model.vlasov_vmax = 1;

    f.model.BoundaryFlux = cemracs2014_TransBoundaryFlux;
    f.model.InitData = cemracs2014_TransInitData;
    f.model.ImposedData = cemcracs2014_imposed_data;
    
  } else {
    f.model.vlasov_mx = 5;
    f.model.vlasov_my = 5;
    f.model.vlasov_vmax = 0.5;

    f.model.BoundaryFlux = vlaTransBoundaryFlux2d;
    f.model.InitData = vlaTransInitData2d;
    f.model.ImposedData = vlaTransImposedData2d;
  }

  f.model.m = f.model.vlasov_mx * f.model.vlasov_my * f.model.vlasov_mz;
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

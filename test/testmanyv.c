#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include <stdlib.h>     /* atoi */

int test_manyv(int deg, int nraf, double cfl, double tmax) 
{
  bool test = true;
  Field f;
  f.model.cfl = 0.05;  
  f.model.vlasov_mx = 5;
  f.model.vlasov_my = 5;
  f.model.vlasov_mz = 1;
  f.model.m = f.model.vlasov_mx * f.model.vlasov_my * f.model.vlasov_mz;
  f.model.vlasov_vmax = 0.5;
  f.model.NumFlux = vlaTransNumFlux2d;
  f.model.BoundaryFlux = vlaTransBoundaryFlux2d;
  f.model.InitData = vlaTransInitData2d;
  f.model.ImposedData = vlaTransImposedData2d;
  f.varindex = GenericVarindex;

  // Set the global parameters for the Vlasov equation
  set_vlasov_params(&(f.model));
  

  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = deg; // x direction degree
  f.interp.interp_param[2] = deg; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = nraf; // x direction refinement
  f.interp.interp_param[5] = nraf; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  // Read the gmsh file
  ReadMacroMesh(&(f.macromesh), "geo/square.msh");
  // Try to detect a 2d mesh
  bool is2d = Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);

  // Mesh preparation
  BuildConnectivity(&(f.macromesh));
 
  // Prepare the initial fields
  InitField(&f);
  f.is2d = true;

  // Prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param: %f\n", f.hmin);

  RK2(&f, tmax);
 
  // Save the results and the error
  bool writemsh = false;
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
	printf("ix: %d, iy: %d, fieldname: %s\n", ix, iy, fieldname);
	PlotField(mplot, false, &f, fieldname, filename);
      }
    }
    /* PlotField(mplot, true, &f, "dgerror.msh"); */
  }

  printf("cfl: %f, deg: %d, nraf: %d\n", cfl, deg, nraf);
  double dd = L2error(&f);
  printf("L2 error\n");
  printf("%f\n", dd);
  
  test = test && (dd < 1e-3); // FIXME: reasonable precision?

  return test;
}

int main(int argc, char* argv[]) {
  // Unit tests
  double cfl = 0.05;
  int deg = 3;
  int nraf = 4;
  double tmax = 1e-2;
  for (;;) {
    int cc = getopt(argc, argv, "c:d:n:t:");
    if (cc == -1) break;
    switch (cc) {
    case 0:
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
    default:
      printf("Error: invalid option.\n");
      printf("Usage:\n");
      printf("./testmanyv -c <cfl> -d <deg> -n <nraf> -t <tmax>\n");
      exit(1);
    }
  }

  int resu = test_manyv(deg, nraf, cfl, tmax);
  if (resu) printf("multiple velocity transport test OK !\n");
  else printf("multiple velocity transport test failed !\n");
  return !resu;
}

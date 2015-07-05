#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "getopt.h"

int main(int argc, char *argv[]) 
{
  int dimension = 2; // Dimension of simulation
  int deg[3] = {3, 3, 3}; // Poynomial egree
  int raf[3] = {4, 4, 4}; // Number of subcells per macrocell
  real cfl = 0.05;
  real tmax = 0.1;
  real dt = 0;
  char *fluxdefault = "VecTransNumFlux2d";
  char *bfluxdefault = "TransBoundaryFlux2d";
  char *mshdefault = "disque.msh";

  int len;
  
  len = strlen(fluxdefault) + 1;
  char *fluxname = malloc(len);
  strncpy(fluxname, fluxdefault, len);

  len = strlen(bfluxdefault) + 1;
  char *bfluxname = malloc(len);
  strncpy(bfluxname, bfluxdefault, len);

  len = strlen(mshdefault) + 1;
  char *mshname = malloc(len);
  strncpy(mshname, mshdefault, len);

  bool usegpu = false;

  char *usage = "./schnaps\n \
\t-c <float> set CFL\n\
\t-n <int> Dimension of simulation (1, 2, or 3)\n\
\t-d <int> Interpolation degree\n\
\t-r <int> Number of subcells in each direction\n\
\t-f <string> Numerical flux dt\n\
\t-b <string> Boundary flux dt\n\
\t-m <string> gmsh filenam dt\n\
\t-s <float> dt\n\
\t-T <float> tmax\n";
  char *openclusage = "\t-g <0=false or 1=true> Use OpenCL\n\
\t-P <int> OpencL platform number\n \
\t-D <int> OpencL device number\n";

  for (;;) {
    int cc = getopt(argc, argv, "c:n:d:r:s:f:T:P:D:g:m:h");
    if (cc == -1) break;
    switch (cc) {
    case 0:
      break;
    case 'c':
      cfl = atof(optarg);
      break;
    case 'n':
      dimension = atoi(optarg);
      // TODO: check that argument is 1, 2, or 3
      break;
    case 'd':
      deg[0] = deg[1] = deg[2] = atoi(optarg);
      // TODO: check that arg is >= 0
      break;
    case 'r':
      raf[0] = raf[1] = raf[2] = atoi(optarg);
      // TODO: check that arg is >= 0
      break;
    case 's':
      dt = atof(optarg);
      // TODO: check that arg is >= 0
      break;
    case 'T':
      tmax = atof(optarg);
      // TODO: check that arg is >= 0
      break;
    case 'f':
      {
	int len = strlen(optarg) + 1;
	free(fluxname);
	fluxname = malloc(len);
	strncpy(fluxname, optarg, len);
      }
      break;
    case 'b':
      {
	int len = strlen(optarg) + 1;
	free(bfluxname);
	bfluxname = malloc(len);
	strncpy(bfluxname, optarg, len);
      }
      break;
    case 'm':
      {
	int len = strlen(optarg) + 1;
	free(mshname);
	mshname = malloc(len);
	strncpy(mshname, optarg, len);
      }
      break;
#ifdef _WITH_OPENCL
    case 'g':
      usegpu = atoi(optarg);
      break;
    case 'D':
      ndevice_cl = atoi(optarg);
      break;
    case 'P':
      nplatform_cl = atoi(optarg);
      break;
#endif
    case 'h':
      printf("Usage:\n%s", usage);
#ifdef _WITH_OPENCL 
      printf("%s", openclusage);
#endif
      exit(0);
      break;
    default:
      printf("Error: invalid option.\n");
      printf("Usage:\n%s", usage);
#ifdef _WITH_OPENCL 
      printf("%s", openclusage);
#endif
      exit(1);
    }
  }

  if(dimension < 3) {
    deg[2] = 0;
    raf[2] = 1;
  }
  if(dimension < 2) {
    deg[1] = 0;
    raf[1] = 1;
  }

  field f;
  init_empty_field(&f);

  f.model.cfl = cfl;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = numflux(fluxname);
  f.model.BoundaryFlux = bflux(bfluxname);
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = deg[0]; // x direction degree
  f.interp.interp_param[2] = deg[1]; // y direction degree
  f.interp.interp_param[3] = deg[2]; // z direction degree
  f.interp.interp_param[4] = raf[0]; // x direction refinement
  f.interp.interp_param[5] = raf[1]; // y direction refinement
  f.interp.interp_param[6] = raf[2]; // z direction refinement

#ifdef _WITH_OPENCL
  char buf[1000];

  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", fluxname);
  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", bfluxname);
  strcat(cl_buildoptions, buf);
#endif

  // Read the gmsh file
  ReadMacroMesh(&f.macromesh, mshname);
  //ReadMacroMesh(&(f.macromesh), "geo/cube.msh");

  if(dimension == 2) {
    Detect2DMacroMesh(&f.macromesh);
    assert(f.macromesh.is2d);
  }

  if(dimension == 1) {
    Detect2DMacroMesh(&f.macromesh);
    assert(f.macromesh.is1d);
  }

  // Mesh preparation
  BuildConnectivity(&f.macromesh);
  //PrintMacroMesh(&f.macromesh);

  // Prepare the initial fields
  Initfield(&f);

  // AffineMapMacroMesh(&f.macromesh);
  CheckMacroMesh(&f.macromesh, f.interp.interp_param + 1);

  f.vmax = 0.1;
  if(dt <= 0.0)
    dt = set_dt(&f);

  printf("\n\n");

  printf("Working in dimension %d\n", dimension);
  printf("Polynomial degree: %d, %d, %d\n", deg[0], deg[1], deg[2]);
  printf("Number of subcells: %d, %d, %d\n", raf[0], raf[1], raf[2]);
  printf("cfl param: %f\n", f.hmin);
  printf("dt: %f\n", dt);
  printf("tmax: %f\n", tmax);
  printf("Numerical flux: %s\n", fluxname);
  printf("Boundary flux: %s\n", bfluxname);  
  printf("gmsh file: %s\n", mshname);

  printf("\n\n");

  if(usegpu) {
#ifdef _WITH_OPENCL
    RK2_CL(&f, tmax, dt,  0, NULL, NULL);
    CopyfieldtoCPU(&f);
    printf("\nOpenCL Kernel time:\n");
    show_cl_timing(&f);
    printf("\n");
#else
    printf("OpenCL not enabled!\n");
    exit(1);
#endif
  } else {
    RK2(&f, tmax, dt);
  }

  // Save the results and the error
  Plotfield(0, false, &f, NULL, "dgvisu.msh");
  Plotfield(0, true, &f, "Error", "dgerror.msh");

  real dd = L2error(&f);
 
  printf("\n");
  printf("L2 error: %f\n", dd);
  return 0;
}

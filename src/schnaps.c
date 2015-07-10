#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "getopt.h"
#ifdef _WITH_OPENCL
#include "clutils.h"
#include "clinfo.h"
#endif
int main(int argc, char *argv[]) 
{
  int dimension = 2; // Dimension of simulation
#ifdef _WITH_OPENCL
  bool usegpu = true;
#else
  bool usegpu = false;
#endif
  int deg[3] = {3, 3, 3}; // Poynomial degree
  int raf[3] = {2, 2, 2}; // Number of subcells per macrocell
  real cfl = 0.05;
  real tmax = 0.1;
  real dt = 0;
  char *fluxdefault = "TransNumFlux2d";
  char *bfluxdefault = "TransBoundaryFlux2d";
  char *initdatadefault = "TransInitData2d";
  char *imposeddatadefault = "TransImposedData2d";
  char *mshdefault = "disque.msh";
  int m = 1;
  bool writeout = true;

  int len;
  
  len = strlen(fluxdefault) + 1;
  char *fluxname = malloc(len);
  strncpy(fluxname, fluxdefault, len);

  len = strlen(bfluxdefault) + 1;
  char *bfluxname = malloc(len);
  strncpy(bfluxname, bfluxdefault, len);

  len = strlen(bfluxdefault) + 1;
  char *initdataname = malloc(len);
  strncpy(initdataname, initdatadefault, len);

  len = strlen(bfluxdefault) + 1;
  char *imposeddataname = malloc(len);
  strncpy(imposeddataname, imposeddatadefault, len);

  len = strlen(mshdefault) + 1;
  char *mshname = malloc(len);
  strncpy(mshname, mshdefault, len);

  char *usage = "./schnaps\n \
\t-c <float> set CFL\n\
\t-n <int> Dimension of simulation (1, 2, or 3)\n\
\t-m <int> Number of hyperbolicly conserved variables\n\
\t-d <int> Interpolation degree\n\
\t-r <int> Number of subcells in each direction\n\
\t-f <string> Numerical flux\n\
\t-b <string> Boundary flux\n\
\t-i <string> Initialization function\n\
\t-I <string> Imposed data function\n\
\t-G <string> gmsh filename\n\
\t-s <float> dt\n\
\t-T <float> tmax\n\
\t-w <0=false or 1=true> Write output to disk.\n";
  char *openclusage = "\t-g <0=false or 1=true> Use OpenCL\n\
\t-P <int> OpencL platform number\n \
\t-D <int> OpencL device number\n";

  for (;;) {
    int cc = getopt(argc, argv, "c:n:m:d:r:s:f:b:w:i:I:T:P:D:g:G:h");
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
    case 'm':
      m = atoi(optarg);
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
    case 'i':
      {
	int len = strlen(optarg) + 1;
	free(initdataname);
	initdataname = malloc(len);
	strncpy(initdataname, optarg, len);
      }
      break;
    case 'I':
      {
	int len = strlen(optarg) + 1;
	free(imposeddataname);
	imposeddataname = malloc(len);
	strncpy(imposeddataname, optarg, len);
      }
      break;
    case 'G':
      {
	int len = strlen(optarg) + 1;
	free(mshname);
	mshname = malloc(len);
	strncpy(mshname, optarg, len);
      }
      break;
    case 'w':
      writeout = atoi(optarg);
      break;
#ifdef _WITH_OPENCL
    case 'g':
      usegpu = atoi(optarg);
      break;
    case 'P':
      nplatform_cl = atoi(optarg);
      break;
    case 'D':
      ndevice_cl = atoi(optarg);
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
  f.model.m = m; // only one conservative variable
  if(!usegpu) {
    f.model.NumFlux = numflux(fluxname);
    f.model.BoundaryFlux = bflux(bfluxname);
  }
  f.model.InitData = initdata(initdataname);
  f.model.ImposedData = imposeddata(imposeddataname);
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
    Detect1DMacroMesh(&f.macromesh);
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
  
  if(!usegpu) {
    printf("C version\n");
  } else {
    printf("OpenCL version\n");
    printf("OpenCL platform: %d\n", nplatform_cl);
    printf("OpenCL device: %d\n", ndevice_cl);
  }
  printf("Working in dimension %d\n", dimension);
  printf("m: %d\n", m);
  printf("Numerical flux: %s\n", fluxname);
  printf("Boundary flux: %s\n", bfluxname);  
  printf("Init function: %s\n", initdataname);
  printf("Imposed data function: %s\n", imposeddataname);

  printf("gmsh file: %s\n", mshname);
  printf("Polynomial degree: %d, %d, %d\n", deg[0], deg[1], deg[2]);
  printf("Number of subcells: %d, %d, %d\n", raf[0], raf[1], raf[2]);
  printf("cfl param: %f\n", f.hmin);
  printf("dt: %f\n", dt);
  printf("tmax: %f\n", tmax);
  printf("Buffer size (GB): %f\n", f.wsize * sizeof(real) * 1e-9);

  printf("\n\n");

  if(usegpu) {
#ifdef _WITH_OPENCL
    RK2_CL(&f, tmax, dt, 0, NULL, NULL);

    cl_int status = clFinish(f.cli.commandqueue);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    
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
  if(writeout) {
    Plotfield(0, false, &f, NULL, "dgvisu.msh");
    Plotfield(0, true, &f, "Error", "dgerror.msh");
  }

  real dd = L2error(&f);
 
  printf("\n");
  printf("L2 error: %f\n", dd);
  return 0;
}

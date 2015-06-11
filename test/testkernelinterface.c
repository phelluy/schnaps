#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestKernelInterface(void){
  bool test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  field f;
  init_empty_field(&f);
  
  // Original:
  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 3; // x direction refinement
  f.interp.interp_param[5] = 3; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  BuildConnectivity(&(f.macromesh));
  PrintMacroMesh(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  Initfield(&f);

  MacroFace mface[f.macromesh.nbfaces];
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ++ifa){
    mface[ifa].first = ifa;
    mface[ifa].last_p1 = ifa + 1;
  }

  void* chkptr;
  cl_int status;
  chkptr = clEnqueueMapBuffer(f.cli.commandqueue,
			      f.dtwn_cl,
			      CL_TRUE,
			      CL_MAP_WRITE,
			      0, // offset
			      sizeof(real) * (f.wsize),
			      0, NULL, NULL, // events management
			      &status);
  assert(status == CL_SUCCESS);
  assert(chkptr == f.dtwn);

  for(int i = 0; i < f.wsize; i++)
    f.dtwn[i] = 0.0;

  status=clEnqueueUnmapMemObject(f.cli.commandqueue,
				 f.dtwn_cl,
				 f.dtwn,
				 0, NULL, NULL);

  assert(status == CL_SUCCESS);
  status = clFinish(f.cli.commandqueue);
  assert(status == CL_SUCCESS);

  // OpenCL version
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    printf("ifa: %d\n", ifa);
    DGMacroCellInterface_CL((void*) (mface + ifa), &f, &(f.wn_cl), 
			    0, NULL, NULL);
    //clFinish(f.cli.commandqueue);
  }
  clFinish(f.cli.commandqueue);
  CopyfieldtoCPU(&f);
  //Displayfield(&f);
  real *fdtwn_opencl = f.dtwn;

  // OpenMP version
  f.dtwn = calloc(f.wsize, sizeof(real));
  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++)
    DGMacroCellInterface((void*) (mface + ifa), &f, f.wn, f.dtwn);
  //Displayfield(&f);
  real *fdtwn_openmp = f.dtwn;

  real maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++) {
    real error = fabs(fdtwn_openmp[i] - fdtwn_opencl[i]);
    printf("error: %f \t%f \t%f\n", error, fdtwn_openmp[i], fdtwn_opencl[i]);
    maxerr = fmax(error, maxerr);
  }
 
  printf("error: %f\n", maxerr);

  test = (maxerr < 1e-8);

  return test;
}

int main(void) {
  int resu = TestKernelInterface();
  if(resu) 
    printf("Interface Kernel test OK !\n");
  else 
    printf("Interface Kernel test failed !\n");
  return !resu;
}

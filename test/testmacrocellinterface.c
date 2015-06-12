#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "clutils.h"

int TestMacroFace(void){
  bool test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  field f;
  init_empty_field(&f);
  
  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  char *mshname =  "test/disque2d.msh";
  
  ReadMacroMesh(&(f.macromesh), mshname);

  Detect2DMacroMesh(&(f.macromesh));
  BuildConnectivity(&(f.macromesh));

#if 1
  // 2D version
  f.model.cfl = 0.05;
  f.model.m = 1;

  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  assert(f.macromesh.is2d);
#else
  // 3D version
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 2; // z direction degree
  f.interp.interp_param[4] = 3; // x direction refinement
  f.interp.interp_param[5] = 3; // y direction refinement
  f.interp.interp_param[6] = 3; // z direction refinement
  
  printf("is2d: %d\n", f.macromesh.is2d);
#endif

  // From testfieldrk2:
  
  //char *mshname =  "test/testdisque.msh";
  /* f.model.cfl = 0.05; */
  /* f.model.m = 1; */
  /* f.model.NumFlux = TransNumFlux; */
  /* f.model.BoundaryFlux = TestTransBoundaryFlux; */
  /* f.model.InitData = TestTransInitData; */
  /* f.model.ImposedData = TestTransImposedData; */
  /* f.varindex = GenericVarindex; */

  /* f.interp.interp_param[0] = f.model.m; */
  /* f.interp.interp_param[1] = 3; // x direction degree */
  /* f.interp.interp_param[2] = 3; // y direction degree */
  /* f.interp.interp_param[3] = 3; // z direction degree */
  /* f.interp.interp_param[4] = 1; // x direction refinement */
  /* f.interp.interp_param[5] = 1; // y direction refinement */
  /* f.interp.interp_param[6] = 1; // z direction refinement */
  
  Initfield(&f);

  MacroFace mface[f.macromesh.nbfaces];
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    mface[ifa].first = ifa;
    mface[ifa].last_p1 = ifa + 1;
  }
 
  //f.is2d = true;

  // OpenCL method
  // NB: Initfield expects a certain address for dtwn, so the OpenCL
  // version must come before the other versions.
  cl_int status;
  void* chkptr = clEnqueueMapBuffer(f.cli.commandqueue,
  				    f.dtwn_cl,
  				    CL_TRUE,
  				    CL_MAP_WRITE,
  				    0, // offset
  				    sizeof(real) * 60,
  				    0, NULL, NULL,
  				    &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  assert(chkptr == f.dtwn);

  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;

  status = clEnqueueUnmapMemObject(f.cli.commandqueue,
  				   f.dtwn_cl,
  				   f.dtwn,
  				   0, NULL, NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  clFinish(f.cli.commandqueue);

  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++)
    DGMacroCellInterface_CL((void*) (mface + ifa), &f, &(f.wn_cl),
			    0, NULL, NULL);

  CopyfieldtoCPU(&f);
  real *fdtwn_opencl = f.dtwn;

  // OpenMP, new method
  f.dtwn = calloc(f.wsize, sizeof(real));
  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++)
    DGMacroCellInterface((void*) (mface + ifa), &f, f.wn, f.dtwn);
  real *fdtwn_openmp = f.dtwn;

  // Check that the results are the same
  test = true;
  real tolerance = 1e-8;

  real maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++) {
    real error = fabs(fdtwn_openmp[i] - fdtwn_opencl[i]);
    //printf("error: %f \t%f \t%f\n", error, fdtwn_openmp[i], fdtwn_opencl[i]);
    maxerr = fmax(error, maxerr);
  }
  printf("Max difference between OpenCL and OpenMP: %f\n", maxerr);
  test = test && (maxerr < tolerance);

  // OpenMP, slow method
  /* MacroCell mcell[f.macromesh.nbelems]; */
  /* for(int ie = 0; ie < f.macromesh.nbelems; ie++) { */
  /*   mcell[ie].field = &f; */
  /*   mcell[ie].first = ie; */
  /*   mcell[ie].last_p1 = ie + 1; */
  /* } */
  /* f.dtwn = calloc(f.wsize, sizeof(real)); */
  /* for(int iw = 0; iw < f.wsize; iw++) */
  /*   f.dtwn[iw] = 0; */
  /* for(int ie = 0; ie < f.macromesh.nbelems; ++ie) { */
  /*   MacroCell *mcelli = mcell + ie; */
  /*   DGMacroCellInterfaceSlow(mcelli); */
  /* } */
  /* real *fdtwn_slow = f.dtwn; */

  /* maxerr = 0.0; */
  /* for(int i = 0; i < f.wsize; i++) { */
  /*   maxerr = fmax(fabs(fdtwn_openmp[i] - fdtwn_slow[i]), maxerr); */
  /* } */
  /* printf("Max difference between OpenMP and OpenMP-slow: %f\n", maxerr); */
  /* test = test && (maxerr < tolerance); */

  /* maxerr = 0.0; */
  /* for(int i = 0; i < f.wsize; i++) { */
  /*   maxerr = fmax(fabs(fdtwn_opencl[i] - fdtwn_slow[i]), maxerr); */
  /* } */
  /* printf("Max difference between OpenCL and OpenMP-slow: %f\n", maxerr); */
  /* test = test && (maxerr < tolerance); */

  return test;
}

int main(void) {
  int resu = TestMacroFace();
  if(resu) 
    printf("MacroFace test OK\n");
  else 
    printf("MacroFace test FAILED\n");
  return !resu;
} 

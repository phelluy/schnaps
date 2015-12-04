#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "clutils.h"

int TestMacroFace()
{
  int retval = 0;

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

  char *mshname =  "../test/disque2d.msh";
  
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

  //f.is2d = true;

  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;

  real tnow = 0.0;
  
  CopyfieldtoGPU(&f);
  
  clFinish(f.cli.commandqueue);

  const int ninterfaces = f.macromesh.nmacrointerfaces;
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f.macromesh.macrointerface[i];
    DGMacroCellInterface_CL(f.mface + ifa, &f, f.wn_cl,
			    0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }
  
  const int nboundaryfaces = f.macromesh.nboundaryfaces;
  for(int i = 0; i < nboundaryfaces; ++i) {
    int ifa = f.macromesh.boundaryface[i];
    DGBoundary_CL(f.mface + ifa, &f, f.wn_cl, tnow, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  CopyfieldtoCPU(&f);
  real *fdtwn_opencl = f.dtwn;

  // OpenMP, new method
  f.dtwn = calloc(f.wsize, sizeof(real));

  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;
  
  {
    // Macrocell interfaces
    const int ninterfaces = f.macromesh.nmacrointerfaces;
    for(int i = 0; i < ninterfaces; ++i) {
      int ifa = f.macromesh.macrointerface[i];
      MacroFace *mface = f.mface + ifa;

      int ieL = mface->ieL;
      MacroCell *mcellL = f.mcell + ieL;
      real *wL = f.wn + mcellL->woffset;
      real *dtwL = f.dtwn + mcellL->woffset;

      int ieR = mface->ieR;
      MacroCell *mcellR = f.mcell + ieR;
      real *wR = f.wn + mcellR->woffset;
      real *dtwR = f.dtwn + mcellR->woffset;

      DGMacroCellInterface(f.mface + ifa, &f, wL, wR, dtwL, dtwR);
    }
  
    // Macrocell boundaries
    const int nboundaryfaces = f.macromesh.nboundaryfaces;
    for(int i = 0; i < nboundaryfaces; ++i) {
      int ifa = f.macromesh.boundaryface[i];
      MacroFace *mface = f.mface + ifa;
      int ie = mface->ieL;
      MacroCell *mcell = f.mcell + ie;
      real *wmc = f.wn + mcell->woffset;
      real *dtwmc = f.dtwn + mcell->woffset;
      DGMacroCellBoundary(mface, &f, wmc, dtwmc);
    }
  }

  real *fdtwn_openmp = f.dtwn;

  // Check that the results are the same
  real tolerance;
  if(sizeof(real) == sizeof(double))
    tolerance = 1e-8;
  else
    tolerance = 1e-4;

  real maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++) {
    real error = fabs(fdtwn_openmp[i] - fdtwn_opencl[i]);
    //printf("error: %f \t%f \t%f\n", error, fdtwn_openmp[i], fdtwn_opencl[i]);
    maxerr = fmax(error, maxerr);
  }
  printf("Max difference between OpenCL and OpenMP: %f\n", maxerr);
  if(maxerr > tolerance)
    retval += 1;

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

  return retval;
}

int main()
{
  int retval = TestMacroFace();
  if(retval == 0) 
    printf("MacroFace test OK\n");
  else 
    printf("MacroFace test FAILED\n");
  return retval;
} 

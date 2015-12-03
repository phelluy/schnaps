#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestKernelInterface()
{
  int retval = 0;

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

  ReadMacroMesh(&f.macromesh, "../test/testmacromesh.msh");
  //ReadMacroMesh(&f.macromesh, "test/testcube.msh");
  Detect2DMacroMesh(&f.macromesh);
  assert(f.macromesh.is2d);

  BuildConnectivity(&f.macromesh);
  PrintMacroMesh(&f.macromesh);

  //AffineMapMacroMesh(&f.macromesh);
 
  Initfield(&f);
  free(f.dtwn);

  real* dtwn_extract = malloc(f.wsize * sizeof(real));
  f.dtwn = dtwn_extract;
  
  for(int i = 0; i < f.wsize; i++)
    f.dtwn[i] = 0.0;
  
  const int ninterfaces = f.macromesh.nmacrointerfaces;

  // Test interface extraction
  printf("OpenCL extraction:\n");
  for(int i = 0; i < f.wsize; i++)
    f.dtwn[i] = 0.0;

  CopyfieldtoGPU(&f);

  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f.macromesh.macrointerface[i];
    MacroFace *mface = f.mface + ifa;
    //ExtractInterface_CL(mface, &f, f.wn_cl, 0, NULL, NULL);
    // FIXME
    
    clFinish(f.cli.commandqueue);
  }

  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f.macromesh.macrointerface[i];
    MacroFace *mface = f.mface + ifa;
    compute_extracted_DGInterface_CL(mface, &f, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }
  
  
  // FIXME: add actual computation, addition to dtwn_cl, and comparisons
  
  CopyfieldtoCPU(&f);

   
  // OpenCL version
  printf("OpenCL version:\n");

  real* dtwn_cl = malloc(f.wsize * sizeof(real));
  f.dtwn = dtwn_cl;
  
  for(int i = 0; i < f.wsize; i++)
    f.dtwn[i] = 0.0;

  CopyfieldtoGPU(&f);

  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f.macromesh.macrointerface[i];
    MacroFace *mface = f.mface + ifa;
    DGMacroCellInterface_CL(mface, &f, f.wn_cl, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }
  
  const int nboundaryfaces = f.macromesh.nboundaryfaces;
  for(int i = 0; i < nboundaryfaces; ++i) {
    int ifa = f.macromesh.boundaryface[i];
    MacroFace *mface = f.mface + ifa;
    DGBoundary_CL(mface, &f, f.wn_cl, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  CopyfieldtoCPU(&f);
  //Displayfield(&f);


  // OpenMP version
  printf("OpenMP version:\n");

  real* dtwn_omp = malloc(f.wsize * sizeof(real));
  f.dtwn = dtwn_omp;
  
  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0.0;

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
      int ie = f.macromesh.face2elem[4 * ifa + 0]; // FIXME: put in ifa
      MacroCell *mcell = f.mcell + ie;
      real *wmc = f.wn + mcell->woffset;
      real *dtwmc = f.dtwn + mcell->woffset;
      DGMacroCellBoundary(f.mface + ifa, &f, wmc, dtwmc);
    }
  }

  real maxerr = 0.0;
  for(int i = 0; i < f.wsize; i++) {
    real error = fabs(dtwn_omp[i] - dtwn_cl[i]);
    //printf("error: %f \t%f \t%f\n", error, dtwn_omp[i], dtwn_cl[i]);
    if(error > maxerr) 
      maxerr = error;
  }
 
  printf("maxerr: %f\n", maxerr);

  real tolerance;
  if(sizeof(real) == sizeof(double))
    tolerance = 1e-8;
  else
    tolerance = 1e-6;
  if(maxerr > tolerance) {
    retval += 1;
  }

  return retval;
}

int main()
{
  int retval = TestKernelInterface();
  if(retval == 0) 
    printf("Interface Kernel test OK !\n");
  else 
    printf("Interface Kernel test failed !\n");
  return retval;
}

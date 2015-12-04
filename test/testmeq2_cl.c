#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "getopt.h"
#include "model.h"
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

real maxerr(real *a, real *b, int n) 
{
  real err = 0.0;
  for(int i = 0; i < n; ++i) {
    err = fmax(fabs(a[i] - b[i]), err);
  }
  return err;
}

int TestmEq2()
{
  int retval = 0;

  field f;
  init_empty_field(&f);
  
  f.model.cfl = 0.05;  
  f.model.m = 2; 
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 1; // x direction degree
  f.interp.interp_param[2] = 1; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  f.model.NumFlux = VecTransNumFlux2d;
  f.model.BoundaryFlux = VecTransBoundaryFlux2d;
  f.model.InitData = VecTransInitData2d;
  f.model.ImposedData = VecTransImposedData2d;
  f.varindex = GenericVarindex;

  char buf[1000];
  sprintf(buf, "-D _M=%d", f.model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "VecTransNumFlux2d");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "VecTransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);

  ReadMacroMesh(&f.macromesh, "../test/testcube.msh");

  // Try to detect a 2d mesh
  Detect2DMacroMesh(&f.macromesh);
  assert(f.macromesh.is2d);
  BuildConnectivity(&f.macromesh);

  //AffineMapMacroMesh(&f.macromesh);
 
  Initfield(&f);
  
  CheckMacroMesh(&f.macromesh, f.interp.interp_param + 1);

  real *dtwn_cl = f.dtwn;
  real *dtwn = calloc(f.wsize, sizeof(real));

  real tnow = 0.0;
  
  real err;
    real tolerance;
  if(sizeof(real) == sizeof(double))
    tolerance = 1e-8;
  else
    tolerance = 1e-6;

  printf("Test volume terms\n");
  
  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }
  
  f.dtwn = dtwn_cl;

  const int nmacro = f.macromesh.nbelems;

  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f.mcell + ie;
    set_buf_to_zero_cl(f.dtwn_cl + ie, mcell, &f,
		       0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    DGVolume_CL(f.mcell + ie, &f, f.wn_cl, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);
  
  f.dtwn = dtwn;
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    DGVolume(f.mcell + ie, &f, f.wn, f.dtwn);
  }

  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  if(err > tolerance) {
    retval += 1;
    printf("Error!\n");
  }

  printf("Test flux terms\n");
  
  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }
  f.dtwn = dtwn_cl;
  
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f.mcell + ie;
    set_buf_to_zero_cl(f.dtwn_cl + ie, mcell, &f,
		       0, NULL, NULL);
  }
  clFinish(f.cli.commandqueue);

  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    int *param = f.interp_param;
    int nraf[3] = {param[4], param[5], param[6]};

    MacroCell *mcell = f.mcell + ie;
    
    if(nraf[0] > 1) {
      DGFlux_CL(mcell, &f, 0, f.wn_cl, 0, NULL, NULL);
      clFinish(f.cli.commandqueue);
    }
      
    if(!f.macromesh.is1d && nraf[1] > 1) {
      DGFlux_CL(mcell, &f, 1, f.wn_cl, 0, NULL, NULL);
      clFinish(f.cli.commandqueue);
    }
      
    if(!f.macromesh.is2d && nraf[2] > 1) {
      DGFlux_CL(mcell, &f, 2, f.wn_cl, 0, NULL, NULL);
      clFinish(f.cli.commandqueue);
    }
  }
  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);

  f.dtwn = dtwn;
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    DGSubCellInterface(f.mcell + ie, &f, f.wn, f.dtwn);
  }

  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  if(err > tolerance) {
    retval += 1;
    printf("Error!\n");
  }
  

  printf("Test macrocell interfaces\n");

  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }

  MacroFace mface[f.macromesh.nbfaces];
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    mface[ifa].ifa = ifa;
  }
  
  f.dtwn = dtwn_cl;
  
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f.mcell + ie;
    set_buf_to_zero_cl(f.dtwn_cl, mcell, &f,
		       0, NULL, NULL);
  }
  clFinish(f.cli.commandqueue);

  const int ninterfaces = f.macromesh.nmacrointerfaces;
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f.macromesh.macrointerface[i];
    DGMacroCellInterface_CL(mface + ifa, &f, f.wn_cl,
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
  clFinish(f.cli.commandqueue);
  
  f.dtwn = dtwn;

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
  
  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  if(err > tolerance) {
    retval += 1;
    printf("Error!\n");
  }

  printf("Test all terms\n");
  for(int i = 0; i < f.wsize; ++i) {
    dtwn_cl[i] = 0.0;
    dtwn[i] = 0.0;
  }

 
  f.dtwn = dtwn_cl;
  dtfield_CL(&f, tnow, f.wn_cl, 0, NULL, NULL);
  clFinish(f.cli.commandqueue);
  
  CopyfieldtoCPU(&f);
  clFinish(f.cli.commandqueue);
  
  f.dtwn = dtwn;
  dtfield(&f, tnow, f.wn, f.dtwn);

  err = maxerr(dtwn, dtwn_cl, f.wsize);
  printf("\tmax error: %f\n", err);
  if(err > tolerance) {
    retval += 1;
    printf("Error!\n");
  }

  return retval;
}

int main()
{
  int retval = TestmEq2();
  if (retval == 0) 
    printf("OpenCL m greater than 1 test OK!\n");
  else 
    printf("OpenCL m greater than 1 test failed !\n");
  return retval;
}

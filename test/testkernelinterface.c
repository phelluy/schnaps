#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

int TestKernelInterface(void){
  bool test = true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  /* field f; */
  /* init_empty_field(&f); */
  
  /* // Original: */
  /* f.model.cfl = 0.05; */
  /* f.model.m = 1; // only one conservative variable */
  /* f.model.NumFlux = TransNumFlux2d; */
  /* f.model.BoundaryFlux = TransBoundaryFlux2d; */
  /* f.model.InitData = TransInitData2d; */
  /* f.model.ImposedData = TransImposedData2d; */
  /* f.varindex = GenericVarindex; */

  /* f.interp.interp_param[0] = f.model.m; */
  /* f.interp.interp_param[1] = 2; // x direction degree */
  /* f.interp.interp_param[2] = 2; // y direction degree */
  /* f.interp.interp_param[3] = 0; // z direction degree */
  /* f.interp.interp_param[4] = 3; // x direction refinement */
  /* f.interp.interp_param[5] = 3; // y direction refinement */
  /* f.interp.interp_param[6] = 1; // z direction refinement */

  /* ReadMacroMesh(&(f.macromesh),"../test/testmacromesh.msh"); */
  /* //ReadMacroMesh(&(f.macromesh),"../test/testcube.msh"); */
  /* Detect2DMacroMesh(&(f.macromesh)); */
  /* assert(f.macromesh.is2d); */

  /* BuildConnectivity(&(f.macromesh)); */
  /* PrintMacroMesh(&(f.macromesh)); */

  /* //AffineMapMacroMesh(&(f.macromesh)); */
 
  /* Initfield(&f); */
  Model model;
  
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  model.Source = NULL;

  char buf[1000];
  sprintf(buf, "-D _M=%d", model.m);
  strcat(cl_buildoptions, buf);

  sprintf(buf," -D NUMFLUX=%s", "TransNumFlux2d");
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s", "TransBoundaryFlux2d");
  strcat(cl_buildoptions, buf);

  int deg[]={3, 3, 0};
  int raf[]={2, 2, 2};

  MacroMesh mesh;

  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");
  //ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  assert(mesh.is2d);
  BuildConnectivity(&mesh);
  
  PrintMacroMesh(&mesh);


  Simulation simu;
  //AffineMapMacroMesh(&(f.macromesh));
  InitSimulation(&simu, &mesh, deg, raf, &model);

  init_field_cl(&simu);

  
  void* chkptr;
  cl_int status;
  chkptr = clEnqueueMapBuffer(simu.cli.commandqueue,
			      simu.dtw_cl,
			      CL_TRUE,
			      CL_MAP_WRITE,
			      0, // offset
			      sizeof(schnaps_real) * (simu.wsize),
			      0, NULL, NULL, // events management
			      &status);
  assert(status == CL_SUCCESS);
  assert(chkptr == simu.dtw);

  for(int i = 0; i < simu.wsize; i++)
    simu.dtw[i] = 0;

  status = clEnqueueUnmapMemObject(simu.cli.commandqueue,
				   simu.dtw_cl,
				   simu.dtw,
				   0, NULL, NULL);
  assert(status == CL_SUCCESS);
  status = clFinish(simu.cli.commandqueue);
  assert(status == CL_SUCCESS);

  // OpenCL version
  
  const int ninterfaces = simu.macromesh.nmacrointerfaces;
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = simu.macromesh.macrointerface[i];
    DGMacroCellInterface_CL(ifa, &simu, &(simu.w_cl), 
    			    0, NULL, NULL);
    clFinish(simu.cli.commandqueue);
  }
  
  const int nboundaryfaces = simu.macromesh.nboundaryfaces;
  for(int i = 0; i < nboundaryfaces; ++i) {
    int ifa = simu.macromesh.boundaryface[i];
    printf("ocl ifa=%d\n",ifa);
    DGBoundary_CL(ifa, &simu, &(simu.w_cl),
			    0, NULL, NULL);
    clFinish(simu.cli.commandqueue);
  }
  clFinish(simu.cli.commandqueue);
  CopyfieldtoCPU(&simu);
  //Displayfield(&f);
  schnaps_real *fdtw_opencl = simu.dtw;

  // OpenMP version
  simu.dtw = calloc(simu.wsize, sizeof(schnaps_real));
  for(int iw = 0; iw < simu.wsize; iw++) simu.dtw[iw] = 0;
  
  for(int ifa = 0; ifa < simu.macromesh.nbfaces; ifa++){
    int fsize =  simu.wsize / simu.macromesh.nbelems;
    int ieL = simu.macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu.macromesh.face2elem[4 * ifa + 1];
    int ieR = simu.macromesh.face2elem[4 * ifa + 2];
    field *fL = simu.fd + ieL;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu.fd + ieR;
      offsetR = fsize * ieR;
    }
    
    printf("cpu ifa=%d\n",ifa);
    DGMacroCellInterface(locfaL, fL, offsetL,
			 fR, offsetR, simu.w, simu.dtw);
  }
  //Displayfield(&f);
  schnaps_real *fdtw_openmp = simu.dtw;

  schnaps_real maxerr = 0.0;
  for(int i = 0; i < simu.wsize; i++) {
    schnaps_real error = fabs(fdtw_openmp[i] - fdtw_opencl[i]);
    printf("error: %f \t%f \t%f\n", error, fdtw_openmp[i], fdtw_opencl[i]);
    maxerr = fmax(error, maxerr);
  }
 
  printf("error: %f\n", maxerr);

  test = (maxerr < _SMALL);

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

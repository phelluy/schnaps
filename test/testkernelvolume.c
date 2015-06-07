#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestKernelVolume(void){
  bool test=true;

  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }

  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1;  // _M
  f.interp.interp_param[1] = 1;  // x direction degree
  f.interp.interp_param[2] = 1;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 3;  // x direction refinement
  f.interp.interp_param[5] = 3;  // y direction refinement
  f.interp.interp_param[6] = 1;  // z direction refinement

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  Initfield(&f);

  /* // set dtwn to 1 for testing */
  
  /* void* chkptr; */
  /* cl_int status; */
  /* chkptr=clEnqueueMapBuffer(f.cli.commandqueue, */
  /*       		    f.dtwn_cl,  // buffer to copy from */
  /*       		    CL_TRUE,  // block until the buffer is available */
  /*       		     CL_MAP_WRITE,  */
  /*       		    0, // offset */
  /*       		    sizeof(real)*(f.wsize),  // buffersize */
  /*       		    0,NULL,NULL, // events management */
  /*       		    &status); */
  /*   assert(status == CL_SUCCESS); */
  /*   assert(chkptr == f.dtwn); */

  /* for(int i=0;i<f.wsize;i++){ */
  /*   f.dtwn[i]=1; */
  /* } */

  /* status=clEnqueueUnmapMemObject (f.cli.commandqueue, */
  /*       			  f.dtwn_cl, */
  /*       			  f.dtwn, */
  /*   			     0,NULL,NULL); */

  /* assert(status == CL_SUCCESS); */
  /* status=clFinish(f.cli.commandqueue); */
  /* assert(status == CL_SUCCESS); */

  clFinish(f.cli.commandqueue);
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    update_physnode_cl(&f, ie, f.physnode_cl, f.physnode, NULL,
		       0, NULL, NULL);
    clFinish(f.cli.commandqueue);

    DGVolume_CL((void*) &(f.mcell[ie]), &f, &(f.wn_cl), 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }
  CopyfieldtoCPU(&f);

  Displayfield(&f);

  // save the dtwn pointer
  real *dtwn_cl = f.dtwn;

  // malloc a new dtwn.
  f.dtwn = calloc(f.wsize, sizeof(real));
 
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    //DGSubCellInterface((void*) &(f.mcell[ie]), &f, f.wn, f.dtwn);
    DGVolume((void*) &(f.mcell[ie]), &f, f.wn, f.dtwn);
  }

  Displayfield(&f);

  //check that the results are the same
  real maxerr = 0.0;
  printf("\nDifference\tC\t\tOpenCL\n");
  for(int i = 0; i < f.wsize; ++i) {
    printf("%f\t%f\t%f\n", f.dtwn[i] - dtwn_cl[i], f.dtwn[i], dtwn_cl[i]);
    maxerr = fmax(fabs(f.dtwn[i] - dtwn_cl[i]), maxerr);
  }
  printf("max error: %f\n",maxerr);

  real tolerance = 1e-8;

  test = (maxerr < tolerance);

  return test;
}

int main(void) {
  // Unit tests
  int resu = TestKernelVolume();
  if (resu) printf("Volume Kernel test OK !\n");
  else printf("Volume Kernel test failed !\n");
  return !resu;
} 

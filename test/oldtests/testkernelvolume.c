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

  //field f;
  //init_empty_field(&f);
  Model model;
  
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;
  model.Source = NULL;
  
  int deg[]={2, 2, 0};
  int raf[]={3, 3, 1};

  MacroMesh mesh;

  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");
  //ReadMacroMesh(&(f.macromesh),"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  assert(mesh.is2d);
  BuildConnectivity(&mesh);
  
  PrintMacroMesh(&mesh);


  Simulation simu;
  //AffineMapMacroMesh(&(f.macromesh));
  InitSimulation(&simu, &mesh, deg, raf, &model);


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

  clFinish(simu.cli.commandqueue);
  for(int ie = 0; ie < simu.macromesh.nbelems; ++ie) {
    /* update_physnode_cl(&f, ie, f.physnode_cl, f.physnode, NULL, */
    /* 		       0, NULL, NULL); */
    /* clFinish(f.cli.commandqueue); */

    printf("macrocell %d\n",ie);
    DGVolume_CL(ie, &simu, &(simu.w_cl), 0, NULL, NULL);
    clFinish(simu.cli.commandqueue);
  }
  CopyfieldtoCPU(&simu);

  DisplaySimulation(&simu);

  // save the dtwn pointer
  schnaps_real *dtwn_cl = simu.dtw;

  // malloc a new dtwn.
  simu.dtw = calloc(simu.wsize, sizeof(schnaps_real));
 
  int fsize =  simu.wsize / simu.macromesh.nbelems;

  for(int ie = 0; ie < simu.macromesh.nbelems; ++ie) {
     //DGSubCellInterface((void*) &(simu.mcell[ie]), &f, simu.wn, simu.dtwn);
    DGVolume(simu.fd + ie, simu.w + ie * fsize, simu.dtw + ie * fsize);
  }

  DisplaySimulation(&simu);

  //check that the results are the same
  schnaps_real maxerr = 0.0;
  printf("\nDifference\tC\t\tOpenCL\n");
  for(int i = 0; i < simu.wsize; ++i) {
    printf("%f\t%f\t%f\n", simu.dtw[i] - dtwn_cl[i], simu.dtw[i], dtwn_cl[i]);
    maxerr = fmax(fabs(simu.dtw[i] - dtwn_cl[i]), maxerr);
  }
  printf("max error: %f\n",maxerr);

  schnaps_real tolerance = _SMALL;

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

#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int main(void) {
  // unit tests
  int resu=TestKernelInterface();
  if (resu) printf("Interface Kernel test OK !\n");
  else printf("Interface Kernel test failed !\n");
  return !resu;
} 

int TestKernelInterface(void){
  bool test = true;

  Field f;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;


  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=3;  // x direction refinement
  f.interp.interp_param[5]=3;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement


  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  bool is2d=Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  InitField(&f);
  f.is2d=true;


  //
  MacroFace mface[f.macromesh.nbfaces];

  for(int ifa=0;ifa<f.macromesh.nbfaces;ifa++){
    mface[ifa].field = &f;
    mface[ifa].first = ifa;
    mface[ifa].last_p1 = ifa + 1;
  }

  /* // set dtwn to 1 for testing */
  
  /* void* chkptr; */
  /* cl_int status; */
  /* chkptr=clEnqueueMapBuffer(f.cli.commandqueue, */
  /*       		    f.dtwn_cl,  // buffer to copy from */
  /*       		    CL_TRUE,  // block until the buffer is available */
  /*       		     CL_MAP_WRITE,  */
  /*       		    0, // offset */
  /*       		    sizeof(double)*(f.wsize),  // buffersize */
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


 
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
     DGMacroCellInterface_CL((void*) (mface + ifa));
  }

  CopyFieldtoCPU(&f);

  DisplayField(&f);

  // save the dtwn pointer
  double* saveptr=f.dtwn;

  // malloc a new dtwn.
  f.dtwn=calloc(f.wsize,sizeof(double));
 
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    DGMacroCellInterface2((void*) (mface + ifa));
  }

  DisplayField(&f);

  //check that the results are the same
  double maxerr=0;
  for(int i=0;i<f.wsize;i++){
    printf("error=%f %f %f\n",f.dtwn[i]-saveptr[i],f.dtwn[i],saveptr[i]);
    maxerr=fmax(fabs(f.dtwn[i]-saveptr[i]),maxerr);
  }
  printf("max error=%f\n",maxerr);

  test=(maxerr < 1e-8);

  return test;

}

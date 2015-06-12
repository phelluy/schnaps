#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestKernel(void)
{
  bool test = true;

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

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 1; // x direction degree
  f.interp.interp_param[2] = 1; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testmacromesh.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);  
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  Initfield(&f);

  printf("&f=%p\n",&f);

  // Set dtwn to 1 for testing
  { 
    void *chkptr;
    cl_int status;
    chkptr = clEnqueueMapBuffer(f.cli.commandqueue,
				f.dtwn_cl,  // buffer to copy from
				CL_TRUE,  // block until the buffer is available
				CL_MAP_WRITE, 
				0, // offset
				sizeof(real) * (f.wsize), // buffersize
				0,
				NULL,
				NULL, // events management
				&status);
    assert(status == CL_SUCCESS);
    assert(chkptr == f.dtwn);

    for(int i = 0; i < f.wsize; i++){
      f.dtwn[i] = 1;
    }

    status = clEnqueueUnmapMemObject(f.cli.commandqueue,
				     f.dtwn_cl,
				     f.dtwn,
				     0,
				     NULL,
				     NULL);
    assert(status == CL_SUCCESS);

    status = clFinish(f.cli.commandqueue);
    assert(status == CL_SUCCESS);
  }
 
  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    printf("ie: %d\n", ie);
    update_physnode_cl(&f, ie, f.physnode_cl, f.physnode, NULL,
		       0, NULL, NULL);
    clFinish(f.cli.commandqueue);

    DGMass_CL((void*) &(f.mcell[ie]), &f, 0, NULL, NULL);
    clFinish(f.cli.commandqueue);
  }

  CopyfieldtoCPU(&f);

  Displayfield(&f);

  // save the dtwn pointer
  real *saveptr = f.dtwn;

  // malloc a new dtwn.
  f.dtwn = calloc(f.wsize,sizeof(real));
  for(int i = 0; i < f.wsize; i++){
    f.dtwn[i] = 1;
  }

  for(int ie = 0; ie < f.macromesh.nbelems; ++ie)
    DGMass((void*) &(f.mcell[ie]), &f, f.dtwn);

  Displayfield(&f);

  //check that the results are the same
  real maxerr = 0;
  for(int i = 0; i < f.wsize; i++){
    printf("error=%f %f %f\n", f.dtwn[i]-saveptr[i], f.dtwn[i],saveptr[i]);
    maxerr=fmax(fabs(f.dtwn[i] - saveptr[i]), maxerr);
  }
  printf("max error=%f\n",maxerr);

  test = (maxerr < 1e-8);

  return test;
}

int main(void) {
  // Unit tests
  int resu = TestKernel();
  if (resu) 
    printf("Kernel test OK !\n");
  else 
    printf("Kernel test failed !\n");
  return !resu;
} 

#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  
  // unit tests
    
  int resu=TestKernel();
	 
  if (resu) printf("Kernel test OK !\n");
  else printf("Kernel test failed !\n");

  return !resu;
} 




int TestKernel(void){

  bool test=true;

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux2d;
  f.model.BoundaryFlux=TransportBoundaryFlux2d;
  f.model.InitData=TransportInitData2d;
  f.model.ImposedData=TransportImposedData2d;
  f.varindex=GenericVarindex;


  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=4;  // x direction refinement
  f.interp.interp_param[5]=4;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement


  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  bool is2d=Detect2DMacroMesh(&(f.macromesh));
  assert(is2d);
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  InitField(&f);
  f.is2d=true;

  MacroCell mcell[f.macromesh.nbelems];

  for(int ie=0;ie<f.macromesh.nbelems;ie++){
    mcell[ie].field=&f;
    mcell[ie].first_cell=ie;
    mcell[ie].last_cell_p1=ie+1;
  }

  for(int ie=0; ie < f.macromesh.nbelems; ++ie) {
    DGMass_CL((void*) (mcell+ie));
  }

  return test;

}

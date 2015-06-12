#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>


int TestfieldDG(void){

  int test = true;

  field f;
  init_empty_field(&f);
  
  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.model.Source = NULL;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 2; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 2; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testcube2.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));
  //AffineMapMacroMesh(&(f.macromesh));
  PrintMacroMesh(&(f.macromesh));
  
  Initfield(&f);
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  dtfield(&f, f.wn, f.dtwn);
  
  Displayfield(&f);

  Plotfield(0, false, &f, NULL, "visu.msh");
  Plotfield(0, true, &f, "error", "error.msh");

  // Test the time derivative with the exact solution
  for(int i = 0; 
      i < f.model.m * f.macromesh.nbelems * NPG(f.interp.interp_param+1); 
      i++){
    test = test && fabs(4 * f.wn[i] - pow(f.dtwn[i], 2)) < 1e-2;
    printf("i=%d err=%f \n",i,4 * f.wn[i] - pow(f.dtwn[i], 2));
    assert(test);
  }
  
  return test;
};

int main(void) {
  // Unit tests
  int resu = TestfieldDG();
  if (resu) printf("field DG test OK !\n");
  else printf("field DG test failed !\n");
  return !resu;
} 

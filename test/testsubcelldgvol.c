#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldSubCellDGVol(void){
  int test = true;

  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 2; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testdisque.msh");
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));
  //AffineMapMacroMesh(&(f.macromesh));
  PrintMacroMesh(&(f.macromesh));
  
  Initfield(&f);
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  for(int ie = 0;ie < f.macromesh.nbelems; ie++)
    DGMacroCellInterfaceSlow((void*) (f.mcell+ie), &f, f.wn, f.dtwn);
  for(int ie = 0; ie < f.macromesh.nbelems; ie++) {
    DGSubCellInterface((void*) (f.mcell+ie), &f, f.wn, f.dtwn);
    DGVolume((void*) (f.mcell+ie), &f, f.wn, f.dtwn);
    DGMass((void*) (f.mcell+ie), &f, f.dtwn);
    DGSource((void*) (f.mcell+ie), &f, f.wn, f.dtwn);
  }

  /* DGMacroCellInterfaceSlow(&f); */
  /* DGSubCellInterface(&f); */
  /* DGVolume(&f); */
  /* DGMass(&f); */
  
  Displayfield(&f);  

  Plotfield(0, false, &f, NULL, "visu.msh");
  Plotfield(0, true, &f, "error", "error.msh");

  // test the time derivative with the exact solution
  for(int i=0;
      i < f.model.m * f.macromesh.nbelems * NPG(f.interp.interp_param+1);
      i++) {
    test = test && fabs(4 * f.wn[i] - pow(f.dtwn[i] , 2)) < 1e-2;
    assert(test);
  }
  
  return test;
}

int main(void) {
  int resu = TestfieldSubCellDGVol();
  if(resu) 
    printf("field DG Subcell Vol test OK !\n");
  else 
    printf("field DG Subcell Vol test failed !\n");
  return !resu;
} 

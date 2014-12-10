#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  // unit tests
  int resu=TestFieldSubCellDGVol();
  if (resu) printf("Field DG Subcell Vol test OK !\n");
  else printf("Field DG Subcell Vol test failed !\n");
  return !resu;
} 

int TestFieldSubCellDGVol(void){
  int test = true;

  Field f;
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

  InitField(&f);
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  MacroCell mcell[f.macromesh.nbelems];

  for(int ie = 0; ie < f.macromesh.nbelems; ie++){
    mcell[ie].field = &f;
    mcell[ie].first = ie;
    mcell[ie].last_p1 = ie + 1;
  }

  for(int ie = 0;ie < f.macromesh.nbelems; ie++){
    DGMacroCellInterface((void*) (mcell+ie));
  }
  for(int ie = 0; ie < f.macromesh.nbelems; ie++){
    DGSubCellInterface((void*) (mcell+ie));
    DGVolume((void*) (mcell+ie));
    DGMass((void*) (mcell+ie));
  }

  /* DGMacroCellInterface(&f); */
  /* DGSubCellInterface(&f); */
  /* DGVolume(&f); */
  /* DGMass(&f); */
  
  DisplayField(&f);  

  PlotField(0, false, &f, NULL, "visu.msh");
  PlotField(0, true, &f, "error", "error.msh");

  // test the time derivative with the exact solution
  for(int i=0;
      i < f.model.m * f.macromesh.nbelems * NPG(f.interp.interp_param+1);
      i++) {
    test = test && fabs(4 * f.wn[i] - pow(f.dtwn[i] , 2)) < 1e-2;
    assert(test);
  }
  
  return test;
};

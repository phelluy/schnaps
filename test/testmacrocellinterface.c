#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

int TestMacroFace(void);

int main(void) {
  // unit tests
  int resu=TestMacroFace();
  if (resu) printf("MacroFace test OK !\n");
  else printf("MacroFace test failed !\n");
  return !resu;
} 

int TestMacroFace(void){
  bool test = true;

  Field f;
  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testdisque.msh");
  BuildConnectivity(&(f.macromesh));
  InitField(&f);

  MacroFace mface[f.macromesh.nbfaces];
  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    mface[ifa].field = &f;
    mface[ifa].first = ifa;
    mface[ifa].last_p1 = ifa + 1;
  }
  
  printf("2");
  MacroCell mcell[f.macromesh.nbelems];
  for(int ie = 0; ie < f.macromesh.nbelems; ie++) {
    mcell[ie].field = &f;
    mcell[ie].first = ie;
    mcell[ie].last_p1 = ie + 1;
  }

  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;

  for(int ifa = 0; ifa < f.macromesh.nbfaces; ifa++) {
    DGMacroCellInterface((void*) (mface + ifa));
  }

  // save the dtwn pointer
  double *saveptr = f.dtwn;
  // malloc a new dtwn.
  f.dtwn = calloc(f.wsize, sizeof(double));

  for(int iw = 0; iw < f.wsize; iw++)
    f.dtwn[iw] = 0;

  for(int ie = 0; ie < f.macromesh.nbelems; ++ie) {
    MacroCell *mcelli = mcell + ie;
    DGMacroCellInterfaceSlow(mcelli);
  }
 
  // Check that the results are the same
  double maxerr = 0;
  for(int i = 0; i < f.wsize; i++){
    printf("error=%f %f %f\n", f.dtwn[i] - saveptr[i], f.dtwn[i], saveptr[i]);
    maxerr = fmax(fabs(f.dtwn[i] - saveptr[i]), maxerr);
  }
  printf("max error=%f\n", maxerr);

  test=(maxerr < 1e-8);

  return test;
}

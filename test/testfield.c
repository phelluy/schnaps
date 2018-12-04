#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int Testfield(void){
  int test = true;

  Model model;
  
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;

  int deg[]={2, 2, 2};
  int raf[]={2, 2, 2};
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  BuildConnectivity(&mesh);

  CheckMacroMesh(&mesh, deg, raf);

  field f;

  schnaps_real physnode[20][3];

  for(int inoloc = 0; inoloc < 20; inoloc++) {
    int ino = mesh.elem2node[20 * 0 + inoloc];
    physnode[inoloc][0] = mesh.node[3 * ino + 0];
    physnode[inoloc][1] = mesh.node[3 * ino + 1];
    physnode[inoloc][2] = mesh.node[3 * ino + 2];
  }
  
  init_empty_field(&f);
  Initfield(&f, model, physnode, deg, raf, NULL, NULL);


  int ipg = 4;

  schnaps_real xref[3],xphy[3],wtest[1];

  ref_pg_vol(f.deg, f.raf, ipg, xref, NULL, NULL);
  schnaps_ref2phy(f.physnode,
	  xref,
	  NULL, -1, // dphiref, ifa
	  xphy, NULL,
	  NULL, NULL, NULL); // codtau, dphi, vnds

  int imem = f.varindex(f.deg, f.raf, f.model.m,ipg,0);
  f.model.InitData(xphy,wtest);
  test = fabs(f.wn[imem] - wtest[0]) < _SMALL;

  FreeMacroMesh(&mesh);
  
  return test;
}

int main(void) {
  // Unit tests
  int resu = Testfield();
  if (resu)
    printf("field test OK !\n");
  else 
    printf("field test failed !\n");
  return !resu;
} 

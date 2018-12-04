#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int TestGeometry(void);

int main(void) {
  // unit tests
  int resu=TestGeometry();
  if (resu) printf("Geometry test OK !\n");
  else printf("Geometry test failed !\n");
  return !resu;
} 

int TestGeometry(void){
  int test = true;
  MacroMesh mc;

  ReadMacroMesh(&mc,"../test/testgeometry.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  schnaps_real xref[3],v[3];
  schnaps_real xref0[3] = {0.1, 0.3, 0.7};
  schnaps_real v0[3] = {0.1, 0.3, 0.7};
  schnaps_real physnode[20][3];

  schnaps_real xphy[3], dtau[3][3];

  for(int inoloc = 0; inoloc <20; inoloc++){
    int ino = mc.elem2node[0 * 20 + inoloc];
    physnode[inoloc][0] = mc.node[3 * ino + 0]; //x
    physnode[inoloc][1] = mc.node[3 * ino + 1]; //y
    physnode[inoloc][2] = mc.node[3 * ino + 2]; //z
  }

  schnaps_ref2phy(physnode, xref0, 0, -1, xphy, dtau, 0, 0, 0);

  printf("xphy= %f %f %f \n", xphy[0], xphy[1], xphy[2]);

  schnaps_phy2ref(physnode, xphy, xref);

  printf("xref= %f %f %f \n", xref[0], xref[1], xref[2]);

  v[0] = xref[0]-xref0[0];
  v[1] = xref[1]-xref0[1];
  v[2] = xref[2]-xref0[2];

  schnaps_real d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  test = (d < _SMALL);


  // again with the other function
  schnaps_ref2phy(physnode, xref0, 0, -1, xphy, dtau, 0, 0, 0);
  RobustPhy2Ref(physnode, xphy, xref);
  v[0] = xref[0]-xref0[0];
  v[1] = xref[1]-xref0[1];
  v[2] = xref[2]-xref0[2];

  d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  test = test && (d < _SMALL);

  FreeMacroMesh(&mc);


  return test;
}

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

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

  ReadMacroMesh(&mc,"test/testgeometry.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  real xref[3],v[3];
  real xref0[3] = {0.1, 0.3, 0.7};
  real v0[3] = {0.1, 0.3, 0.7};
  real physnode[20][3];

  real xphy[3], dtau[3][3];

  for(int inoloc = 0; inoloc <20; inoloc++){
    int ino = mc.elem2node[0 * 20 + inoloc];
    physnode[inoloc][0] = mc.node[3 * ino + 0]; //x
    physnode[inoloc][1] = mc.node[3 * ino + 1]; //y
    physnode[inoloc][2] = mc.node[3 * ino + 2]; //z
  }

  Ref2Phy(physnode, xref0, 0, -1, xphy, dtau, 0, 0, 0);

  printf("xphy= %f %f %f \n", xphy[0], xphy[1], xphy[2]);

  Phy2Ref(physnode, xphy, xref);

  printf("xref= %f %f %f \n", xref[0], xref[1], xref[2]);

  v[0] = xref[0]-xref0[0];
  v[1] = xref[1]-xref0[1];
  v[2] = xref[2]-xref0[2];

  real d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  test = (d < 1e-12);


  // again with the other function
  Ref2Phy(physnode, xref0, 0, -1, xphy, dtau, 0, 0, 0);
  RobustPhy2Ref(physnode, xphy, xref);
  v[0] = xref[0]-xref0[0];
  v[1] = xref[1]-xref0[1];
  v[2] = xref[2]-xref0[2];

  d = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  test = test && (d < 1e-12);




  return test;
}

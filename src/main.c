#include "macromesh.h"
#include <stdio.h>
#include <assert.h>
#include <h20.h>

int main(void) {

  MacroMesh m;

  ReadMacroMesh(&m,"../geo/disque.msh");
  BuildConnectivity(&m);
  PrintMacroMesh(&m);

  MacroMesh mc;

  ReadMacroMesh(&mc,"../geo/cube.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  double xref[3]={0.5_F,0.5_F,0.5_F};

  return 0;

} 

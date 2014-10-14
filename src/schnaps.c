#include "schnaps.h"
#include <stdio.h>
#include <assert.h>


int main(void) {
 
  MacroMesh mm1;
  ReadMacroMesh(&mm1,"test/permutcube.msh");
  bool is2d=Detect2DMacroMesh(&mm1);
  assert(is2d);
  BuildConnectivity(&mm1);
  CheckMacroMesh(&mm1);
 
  MacroMesh mm3;
  ReadMacroMesh(&mm3,"test/unit-cube.msh");
  is2d=Detect2DMacroMesh(&mm3);
  assert(is2d);
  BuildConnectivity(&mm3);
  CheckMacroMesh(&mm3);


  MacroMesh mm2;
  ReadMacroMesh(&mm2,"test/disque2d.msh");
  is2d=Detect2DMacroMesh(&mm2);
  assert(is2d);
  BuildConnectivity(&mm2);
  CheckMacroMesh(&mm2);
  
  return 0;



};



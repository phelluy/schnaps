#include "schnaps.h"
#include <stdio.h>
#include <assert.h>


int main(void) {

  //printf("%d\n",12%1);
  //assert(1==2);
 
  MacroMesh mm1;
  ReadMacroMesh(&mm1,"test/permutcube.msh");
  bool is2d=Detect2DMacroMesh(&mm1);
  assert(is2d);
  BuildConnectivity(&mm1);
  int param[]={4,4,4,1,1,1,0};
  CheckMacroMesh(&mm1,param);
 
  MacroMesh mm3;
  ReadMacroMesh(&mm3,"test/unit-cube.msh");
  is2d=Detect2DMacroMesh(&mm3);
  assert(is2d);
  BuildConnectivity(&mm3);
  CheckMacroMesh(&mm3,param);


  MacroMesh mm2;
  ReadMacroMesh(&mm2,"test/disque2d.msh");
  is2d=Detect2DMacroMesh(&mm2);
  assert(is2d);
  BuildConnectivity(&mm2);
  CheckMacroMesh(&mm2,param);
  
  return 0;



};



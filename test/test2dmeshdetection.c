#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  int resu = Test2DMeshDetection();
  
  if(resu) 
    printf("2D detect test OK !\n");
  else 
    printf("2D detect failed !\n");

  return !resu;
} 

int Test2DMeshDetection(void) {
  int test = true;
 
  int param[] = {4, // deg x
		 4, // deg y
		 4, // deg z
		 1, // ref x
		 1, // ref y
		 1, // ref z
		 0}; // nraf z

  // 2D mesh tests:

  MacroMesh mm1;
  ReadMacroMesh(&mm1, "test/permutcube.msh");
  Detect2DMacroMesh(&mm1);
  BuildConnectivity(&mm1);
  CheckMacroMesh(&mm1, param);

  test = test && mm1.is2d;

 
  MacroMesh mm2;
  ReadMacroMesh(&mm2, "test/disque2d.msh");
  Detect2DMacroMesh(&mm2);
  BuildConnectivity(&mm2);
  CheckMacroMesh(&mm2, param);

  test = test && mm2.is2d;

  MacroMesh mm3;
  ReadMacroMesh(&mm3, "test/unit-cube.msh");
  Detect2DMacroMesh(&mm3);
  BuildConnectivity(&mm3);
  CheckMacroMesh(&mm3, param);
  
  test = test && mm3.is2d;

  // 3D mesh tests:
  param[2] = 4;

  MacroMesh mm4;
  ReadMacroMesh(&mm4, "test/testdisque.msh");
  Detect2DMacroMesh(&mm4);
  BuildConnectivity(&mm4);
  CheckMacroMesh(&mm4, param);
    
  test = test && !mm4.is2d;

  return test;
};

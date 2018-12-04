#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int Test2DMeshDetection(void);


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

  int deg[] = {4, // deg x
	       4, // deg y
	       4}; // deg z
  
  int raf[] = {1, // raf x
	       1, // raf y
	       1}; // raf z
  MacroMesh mm1;

  // 2D mesh tests:
  ReadMacroMesh(&mm1, "../test/permutcube.msh");
  Detect2DMacroMesh(&mm1);
  BuildConnectivity(&mm1);
	       CheckMacroMesh(&mm1, deg,raf);

  test = test && mm1.is2d;

  FreeMacroMesh(&mm1);

 
  MacroMesh mm2;
  ReadMacroMesh(&mm2, "../test/disque2d.msh");
  Detect2DMacroMesh(&mm2);
  BuildConnectivity(&mm2);
  CheckMacroMesh(&mm2, deg,raf);

  test = test && mm2.is2d;

  FreeMacroMesh(&mm2);

  MacroMesh mm3;
  ReadMacroMesh(&mm3, "../test/unit-cube.msh");
  Detect2DMacroMesh(&mm3);
  BuildConnectivity(&mm3);
  CheckMacroMesh(&mm3, deg,raf);
  
  test = test && mm3.is2d;

  FreeMacroMesh(&mm3);


  MacroMesh mm4;
  ReadMacroMesh(&mm4, "../test/testdisque.msh");
  Detect2DMacroMesh(&mm4);
  BuildConnectivity(&mm4);
  CheckMacroMesh(&mm4, deg,raf);
    
  test = test && !mm4.is2d;

  FreeMacroMesh(&mm4);

  

  return test;
};

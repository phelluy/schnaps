#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>



int main(void) {
  
  // unit tests
    
  int resu=Test2DMeshDetection();
	 
  if (resu) printf("2D detect test OK !\n");
  else printf("2D detect failed !\n");

  return !resu;
} 


int Test2DMeshDetection(void) {

  int test = (1==1);
 
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
  
  return test;



};



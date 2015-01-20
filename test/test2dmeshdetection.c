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
 
  int param[]={4,4,4,1,1,1,0};

  MacroMesh mm1;
  ReadMacroMesh(&mm1, "test/permutcube.msh");
  Detect2DMacroMesh(&mm1);
  assert(mm1.is2d);
  BuildConnectivity(&mm1);
  CheckMacroMesh(&mm1,param);
 
  MacroMesh mm3;
  ReadMacroMesh(&mm3, "test/unit-cube.msh");
  Detect2DMacroMesh(&mm3);
  assert(mm3.is2d);
  BuildConnectivity(&mm3);
  CheckMacroMesh(&mm3,param);

  MacroMesh mm2;
  ReadMacroMesh(&mm2, "test/disque2d.msh");
  Detect2DMacroMesh(&mm2);
  assert(mm2.is2d);
  BuildConnectivity(&mm2);
  CheckMacroMesh(&mm2,param);
  
  return test;
};

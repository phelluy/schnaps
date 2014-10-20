#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  
  // unit tests
    
  int resu=TestMacroMesh();
	 

  if (resu) printf("Mesh test OK !\n");
  else printf("Mesh test failed !\n");

  return !resu;
} 



// some unit tests of the macromesh code
int TestMacroMesh(void){

  int test=1;
  MacroMesh m;

  int param[]={4,4,4,1,1,1,0};
  
  ReadMacroMesh(&m,"test/testmacromesh.msh");
  BuildConnectivity(&m);
  CheckMacroMesh(&m,param);
  PrintMacroMesh(&m);

  test = (m.nbelems == 5);
  test =  (test && m.nbnodes == 50);

  return test;

}

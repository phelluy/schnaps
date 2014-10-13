#include "macromesh.h"
//#include "geometry.h"
//#include "interpolation.h"
#include "test.h"
//#include "model.h"
//#include "field.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

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
  
  ReadMacroMesh(&m,"test/testmacromesh.msh");
  BuildConnectivity(&m);
  PrintMacroMesh(&m);

  test = (m.nbelems == 5);
  test =  (test && m.nbnodes == 50);

  return test;

}

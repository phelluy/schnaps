#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  // Unit tests
  int resu=TestMacroMesh();
  if (resu) printf("Mesh test OK !\n");
  else printf("Mesh test failed !\n");
  return !resu;
} 

// some unit tests of the macromesh code
int TestMacroMesh(void)
{
  MacroMesh m;

  int param[]={4, 4, 4, 1, 1, 1, 0};
  
  // test gmsh file reading
  ReadMacroMesh(&m, "test/testmacromesh.msh");
  BuildConnectivity(&m);
  CheckMacroMesh(&m, param);
  PrintMacroMesh(&m);

  int test = (m.nbelems == 5);
  test = (test && m.nbnodes == 50);



  // test search methods
  real xphy[3]={1,1.1,0.5};
  real xref[3];

  test= test && IsInElem(&m,0,xphy,xref);

  printf("xphy=%f %f %f xref=%f %f %f \n",xphy[0],xphy[1],xphy[2],
	 xref[0],xref[1],xref[2]);

  xphy[2]=-0.5;

  test= test && !IsInElem(&m,0,xphy,xref);

  int num=NumElemFromPoint(&m,xphy,NULL);
  printf("xphy=%f %f %f is in elem=%d\n",xphy[0],xphy[1],xphy[2],num);
  test=test && (num == -1);

  xphy[2]=0.5;
  num=NumElemFromPoint(&m,xphy,NULL);
  printf("xphy=%f %f %f is in elem=%d\n",xphy[0],xphy[1],xphy[2],num);
  test=test && (num == 0);

  real xphy2[3]={1,0,0.33};
  num=NumElemFromPoint(&m,xphy2,NULL);
  printf("xphy=%f %f %f is in elem=%d\n",xphy2[0],xphy2[1],xphy2[2],num);
  test=test && (num == 3);





  return test;
}

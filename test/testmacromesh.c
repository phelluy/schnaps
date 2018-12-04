#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int TestMacroMesh(void);

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

  int deg[]={4, 4, 4};
  int raf[]={2, 2, 2};
  
  // test gmsh file reading
  ReadMacroMesh(&m, "../test/testmacromesh.msh");
  BuildConnectivity(&m);
  CheckMacroMesh(&m, deg, raf);
  PrintMacroMesh(&m);

  int test = (m.nbelems == 5);
  test = (test && m.nbnodes == 50);



  // test search methods
  schnaps_real xphy[3]={1,1.1,0.5};
  schnaps_real xref[3];

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

  schnaps_real xphy2[3]={1,0,0.33};
  num=NumElemFromPoint(&m,xphy2,NULL);
  printf("xphy=%f %f %f is in elem=%d\n",xphy2[0],xphy2[1],xphy2[2],num);
  test=test && (num == 3);


  FreeMacroMesh(&m);


  return test;
}

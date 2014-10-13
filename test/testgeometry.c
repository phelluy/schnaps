#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  
  // unit tests
    
  int resu=TestGeometry();
	 

  if (resu) printf("Geometry test OK !\n");
  else printf("Geometry test failed !\n");

  return !resu;
} 




int TestGeometry(void){

  int test=1;
  MacroMesh mc;

  ReadMacroMesh(&mc,"test/testgeometry.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  double xref[3]={0.1,0.3,0.7};
  double v[3]={0.1,0.3,0.7};
  double physnode[20][3];

  double xphy[3],dtau[3][3];

  for(int inoloc=0;inoloc<20;inoloc++){
    int ino=mc.elem2node[0*20+inoloc];
    physnode[inoloc][0]=mc.node[3*ino+0]; //x
    physnode[inoloc][1]=mc.node[3*ino+1]; //y
    physnode[inoloc][2]=mc.node[3*ino+2]; //z
  }

  Ref2Phy(physnode,xref,0,-1,xphy,dtau,0,0,0);

  printf("xphy= %f %f %f \n",xphy[0],xphy[1],xphy[2]);

  Phy2Ref(physnode,xphy,xref);

  printf("xref= %f %f %f \n",xref[0],xref[1],xref[2]);

  v[0]-=xref[0];
  v[1]-=xref[1];
  v[2]-=xref[2];

  double d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  test = (d < 1e-12);

  return test;

}

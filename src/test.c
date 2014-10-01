#include "macromesh.h"
#include "geometry.h"

#include <stdio.h>
#include <math.h>

// some unit tests of the macromesh code
int TestMacroMesh(void){

  int test=1;
  MacroMesh m;
  
  ReadMacroMesh(&m,"../geo/disque.msh");
  BuildConnectivity(&m);
  PrintMacroMesh(&m);

  test = (m.nbelems == 5);
  test = (m.nbnodes == 50);

  return test;

}



int TestGeometry(void){

  int test=1;
  MacroMesh mc;

  ReadMacroMesh(&mc,"../geo/cube.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  double xref[3]={0.1,0.3,0.7};
  double v[3]={0.1,0.3,0.7};
  double physnode[20*3];

  double xphy[3],dtau[9];

  for(int inoloc=0;inoloc<20;inoloc++){
    int ino=mc.elem2node[0*20+inoloc];
    physnode[3*inoloc+0]=mc.node[3*ino+0]; //x
    physnode[3*inoloc+1]=mc.node[3*ino+1]; //y
    physnode[3*inoloc+2]=mc.node[3*ino+2]; //z
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

int TestInterpolation(void){
  
  int test= (1==1);

  return test;

}

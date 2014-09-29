#include "macromesh.h"
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "global.h"

int main(void) {

  MacroMesh m;

  ReadMacroMesh(&m,"../geo/disque.msh");
  BuildConnectivity(&m);
  PrintMacroMesh(&m);

  MacroMesh mc;

  ReadMacroMesh(&mc,"../geo/cube.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  double xref[3]={0.1,0.3,0.7};
  double physnode[20*3];

  double xphy[3],dtau[9];

  for(int inoloc=0;inoloc<20;inoloc++){
    int ino=m.elem2node[0*20+inoloc];
    physnode[3*inoloc+0]=m.node[3*ino+0]; //x
    physnode[3*inoloc+1]=m.node[3*ino+1]; //y
    physnode[3*inoloc+2]=m.node[3*ino+2]; //z
  }

  Ref2Phy(physnode,xref,0,-1,xphy,dtau,0,0,0);

  printf("xphy= %f %f %f \n",xphy[0],xphy[1],xphy[2]);

  Phy2Ref(physnode,xphy,xref);

  printf("xref= %f %f %f \n",xref[0],xref[1],xref[2]);


  return 0;

} 

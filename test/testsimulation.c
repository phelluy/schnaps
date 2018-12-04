#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestSimulation(void){
  int test = true;

  Model model;
  
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;

  int deg[]={2, 2, 2};
  int raf[]={4, 4, 4};
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  BuildConnectivity(&mesh);

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);


  for(int ie = 0; ie < simu.macromesh.nbelems; ++ie){

    int ipg = 4;

    schnaps_real xref[3],xphy[3],wtest[1];

    field f = simu.fd[ie];

    ref_pg_vol(f.deg, f.raf, ipg, xref, NULL, NULL);

    schnaps_ref2phy(f.physnode,
	    xref,
	    NULL, -1, // dphiref, ifa
	    xphy, NULL,
	    NULL, NULL, NULL); // codtau, dphi, vnds


    int imem = f.varindex(f.deg, f.raf, f.model.m,ipg,0);
    f.model.InitData(xphy,wtest);
    schnaps_real val = f.wn[imem];
    test = test && fabs(val - wtest[0]) < _SMALL;
    printf("x= %f %f %f val = %f valex = %f\n",xphy[0],xphy[1],xphy[2],val,wtest[0]);
  }
  

  //PlotFields(0, false, &simu, "trans", "visu.msh");
  //PlotFields(0, true, &simu, "error", "error.msh");

  freeSimulation(&simu);
  FreeMacroMesh(&mesh);
  
  return test;
}

int main(void) {
  // Unit tests
  int resu = TestSimulation();
  if (resu)
    printf("Simulation test OK !\n");
  else 
    printf("Simulation test failed !\n");
  return !resu;
} 

#include "field.h"
#include <stdio.h>
#include <assert.h>


int main(void) {
 
  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TestTransportBoundaryFlux;
  f.model.InitData=TestTransportInitData;
  f.model.ImposedData=TestTransportImposedData;
  f.varindex=GenericVarindex;

  ReadMacroMesh(&(f.macromesh),"test/testdisque.msh");
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
  CheckMacroMesh(&(f.macromesh));
 
  InitField(&f);

  printf("cfl param =%f\n",f.hmin);


  RK2(&f,.5);
 
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);

  
  return 0;



};



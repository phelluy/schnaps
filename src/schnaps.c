#include "schnaps.h"
#include <stdio.h>
#include <assert.h>


int main(void) {

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TransportBoundaryFlux;
  f.model.InitData=TransportInitData;
  f.model.ImposedData=TransportImposedData;
  f.varindex=GenericVarindex;


  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=2;  // z direction degree
  f.interp.interp_param[4]=1;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement


  ReadMacroMesh(&(f.macromesh),"test/testdisque.msh");
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
  //CheckMacroMesh(&(f.macromesh));
 
  InitField(&f);

  printf("cfl param =%f\n",f.hmin);


  RK2(&f,1);
 
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);

  return 0;

};



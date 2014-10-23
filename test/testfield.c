#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  
  // unit tests
    
  int resu=TestField();
	 
  if (resu) printf("Field test OK !\n");
  else printf("Field test failed !\n");

  return !resu;
} 




int TestField(void){

  int test = (1==1);

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TestTransportBoundaryFlux;
  f.model.InitData=TestTransportInitData;
  f.model.ImposedData=TestTransportImposedData;
  f.varindex=GenericVarindex;

  f.interp.interp_param[0]=1;  // _M
  f.interp.interp_param[1]=2;  // x direction degree
  f.interp.interp_param[2]=2;  // y direction degree
  f.interp.interp_param[3]=2;  // z direction degree
  f.interp.interp_param[4]=2;  // x direction refinement
  f.interp.interp_param[5]=2;  // y direction refinement
  f.interp.interp_param[6]=2;  // z direction refinement

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));

  InitField(&f);
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  PlotField(0,(1==0),&f,"testvisufield.msh");
  
  return test;


};

#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
int main(void) {
  
  // unit tests
    
  int resu=TestFieldRK2();
	 
  if (resu) printf("Field RK2 test OK !\n");
  else printf("Field RK2 test failed !\n");

  return !resu;
} 



int TestFieldRK2(void){

  int test = (1==1);

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


  RK2(&f,.1);
 
  PlotField(0,(1==0),&f,"dgvisu.msh");
  PlotField(0,(1==1),&f,"dgerror.msh");

  double dd=L2error(&f);

  printf("erreur L2=%f\n",dd);

  test=test && (dd < 0.0009);
  
  return test;



};




#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>
int main(void) {
  
  // unit tests
    
  int resu=TestFieldDG();
	 
  if (resu) printf("Field DG test OK !\n");
  else printf("Field DG test failed !\n");

  return !resu;
} 




int TestFieldDG(void){

  int test = (1==1);

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TestTransportBoundaryFlux;
  f.model.InitData=TestTransportInitData;
  f.model.ImposedData=TestTransportImposedData;
  f.varindex=GenericVarindex;

  ReadMacroMesh(&(f.macromesh),"test/testcube.msh");
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));
  //CheckMacroMesh(&(f.macromesh));
  //AffineMapMacroMesh(&(f.macromesh));
  PrintMacroMesh(&(f.macromesh));


  InitField(&f);

  dtField(&f);
  
  DisplayField(&f);  

  int yes_compare = 1;
  int no_compare = 0;

  PlotField(0,no_compare,&f,"visu.msh");
  PlotField(0,yes_compare,&f,"error.msh");

  // test the time derivative that has to be -1
  for(int i=0;i<f.model.m * f.macromesh.nbelems * 
	(_DEGX+1)*(_DEGY+1)*(_DEGZ+1);i++){
    test = test && fabs(4*f.wn[i]-pow(f.dtwn[i],2))<1e-2;
  }
  
  return test;



};

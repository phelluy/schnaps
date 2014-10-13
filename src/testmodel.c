//#include "macromesh.h"
//#include "geometry.h"
//#include "interpolation.h"
#include "test.h"
#include "model.h"
//#include "field.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(void) {
  
  // unit tests
    
  int resu=TestModel();
	 

  if (resu) printf("Model test OK !\n");
  else printf("Model test failed !\n");

  return !resu;
} 

int TestModel(void){

  int test= (1==1);
  // creation of a simple transport model
  Model tr;
  tr.m=1; // only one conservative variable
  tr.NumFlux=TransportNumFlux;
  tr.BoundaryFlux=TestTransportBoundaryFlux;
  tr.InitData=TestTransportInitData;
  tr.ImposedData=TestTransportImposedData;

  double wL[tr.m];
  double wR[tr.m];
  double flux1[tr.m],flux2[tr.m];

  double x[3]={1,1,2};
  double t=0;
  double vn[3]={1/sqrt(3),1/sqrt(3),-1/sqrt(3)};

  tr.InitData(x,wR);
  tr.NumFlux(wL,wR,vn,flux1);
  printf("NumFlux %f \n",flux1[0]);
  tr.BoundaryFlux(x,t,wL,vn,flux2);
  printf("BoundaryFlux %f \n",flux2[0]);

  double err=fabs(flux2[0]-flux1[0]);

  test=(err < 1e-8);

  return test;

};

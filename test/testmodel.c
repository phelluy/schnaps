//#include "macromesh.h"
//#include "geometry.h"
//#include "interpolation.h"
#include "test.h"
#include "model.h"
//#include "field.h"
#include <stdbool.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(void) {
  // Unit tests
  int resu = TestModel();
  if (resu) printf("Model test OK !\n");
  else printf("Model test failed !\n");
  return !resu;
} 

int TestModel(void){
  int test = true;
  // Creation of a simple transport model
  Model tr;
  tr.cfl = 0.05;
  tr.m = 1; // only one conservative variable
  tr.NumFlux = TransNumFlux;
  tr.BoundaryFlux = TestTransBoundaryFlux;
  tr.InitData = TestTransInitData;
  tr.ImposedData = TestTransImposedData;

  real wL[tr.m];
  real wR[tr.m];
  real flux1[tr.m], flux2[tr.m];

  real x[3] = {1, 1, 2};
  real t = 0;
  real vn[3] = {sqrt(1.0 / 3.0), sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};

  real vmax = 1.0;
  tr.InitData(x, wR);
  tr.NumFlux(wL, wR, vn, flux1);
  printf("NumFlux %f \n", flux1[0]);
  tr.BoundaryFlux(x, t, wL, vn, flux2);
  printf("BoundaryFlux %f \n", flux2[0]);

  real err = fabs(flux2[0] - flux1[0]);
  test = (err < 1e-8);
  return test;
};

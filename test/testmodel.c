#include "test.h"
#include "model.h"
#include <stdbool.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>

int main()
{
  // Unit tests
  int retval = TestModel();
  if(retval == 0)
    printf("Model test OK !\n");
  else
    printf("Model test failed !\n");
  return retval;
} 

int TestModel()
{
  int retval = 0;
  
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
  real nflux[tr.m];
  real bflux[tr.m];

  for(int iv = 0; iv < tr.m; ++iv) {
    wL[iv] = 0.0;
    wR[iv] = 0.0;
    nflux[iv] = 0.0;
    bflux[iv] = 0.0;
  }

  real x[3] = {1.0, 1.0, 2.0};
  real t = 0.0;
  real vn[3] = {sqrt(1.0 / 3.0), sqrt(1.0 / 3.0), sqrt(1.0 / 3.0)};
  
  tr.InitData(x, wR);
  tr.NumFlux(wL, wR, vn, nflux);
  printf("NumFlux %f \n", nflux[0]);
  tr.BoundaryFlux(x, t, wL, vn, bflux);
  printf("BoundaryFlux %f \n", bflux[0]);

  real err = fabs(bflux[0] - nflux[0]);
  if(err > 1e-8)
    retval++;
  
  return retval;
}

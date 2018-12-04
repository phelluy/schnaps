//#include "../src/user_model.h"
#define schnaps_real float   // !!!!!!!!!!!!!!  to be improved !!!!!!!!!
#include <math.h>
#include <stdio.h>

void NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real* vnorm, schnaps_real* flux)
{
  // printf("NumFlux\n");
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
  schnaps_real vn
    = transport_v[0] * vnorm[0]
    + transport_v[1] * vnorm[1]
    + transport_v[2] * vnorm[2];
  schnaps_real vnp = vn > 0 ? vn : 0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
}

void ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
  //const real transport_v[] = {1,0,0};
  // printf("ImposedData\n");
  schnaps_real vx
    = transport_v[0] * x[0]
    + transport_v[1] * x[1]
    + transport_v[2] * x[2];
  schnaps_real xx = vx - t;
  //w[0] = 1;
  w[0] = sin(3*xx);
  w[0] = xx * xx;
  //w[0]=0;
}

void BoundaryFlux(schnaps_real *x, schnaps_real t,
			       schnaps_real *wL, schnaps_real* vnorm,
			       schnaps_real* flux) {
  // printf("BoundaryFlux\n");
  schnaps_real wR[1];
  ImposedData(x, t, wR);
  NumFlux(wL, wR, vnorm, flux);
}

void InitData(schnaps_real *x, schnaps_real *w)
{
  // printf("InitData\n");
  schnaps_real t = 0;
  ImposedData(x, t, w);
}

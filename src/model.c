#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

const double transport_v[] = {ONE_OVER_SQRT_3,
			      ONE_OVER_SQRT_3,
			      ONE_OVER_SQRT_3};

//const double transport_v[] = {1, 0, 0};

const double transport_v2d[] = {ONE_OVER_SQRT_2,
				ONE_OVER_SQRT_2,
				0};

void TransNumFlux(double wL[], double wR[], double* vnorm, double* flux) {
  double vn 
    = transport_v[0] * vnorm[0] 
    + transport_v[1] * vnorm[1]
    + transport_v[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
};

void TransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux) {
  double vn 
    = transport_v2d[0] * vnorm[0]
    + transport_v2d[1] * vnorm[1]
    + transport_v2d[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  /* if (fabs(vnorm[2])>1e-6) { */
  /*   printf("vnds %lf %lf %lf \n", vnorm[0], vnorm[1], vnorm[2]); */
  /* } */
  // verify that 2d computations are actually
  // activated
  assert(fabs(vnorm[2]) < 1e-8);
};

void VecTransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux) {
  double vn 
    = transport_v2d[0] * vnorm[0]
    + transport_v2d[1] * vnorm[1]
    + transport_v2d[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  flux[1] = vnp * wL[1] + vnm * wR[1];
  /* if (fabs(vnorm[2])>1e-6) { */
  /*   printf("vnds %lf %lf %lf \n", vnorm[0], vnorm[1], vnorm[2]); */
  /* } */
  // Verify that 2d computations are actually activated
  assert(fabs(vnorm[2]) < 1e-8);
};

void TransBoundaryFlux(double x[3], double t, double wL[], double* vnorm,
		       double* flux) {
  double wR[1];
  TransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
};

void TransBoundaryFlux2d(double x[3], double t, double wL[], double *vnorm,
			 double *flux) {
  double wR[1];
  TransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
};

void VecTransBoundaryFlux2d(double x[3], double t, 
			    double wL[], double *vnorm,
			    double* flux) {
  double wR[2];
  VecTransImposedData2d(x, t, wR);
  VecTransNumFlux2d(wL, wR, vnorm, flux);
};

void TransInitData(double x[3], double w[]) {
  double t = 0;
  TransImposedData(x, t, w);
};

void TransInitData2d(double x[3], double w[]) {
  double t = 0;
  TransImposedData2d(x, t, w);
};

void VecTransInitData2d(double x[3], double w[]) {
  double t = 0;
  VecTransImposedData2d(x, t, w);
};

void TransImposedData(double x[3], double t, double w[]) {
  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];
  double xx = vx - t;
  w[0]=cos(xx);
};

void TransImposedData2d(double x[3], double t, double w[]) {
  double vx 
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = cos(xx);
};

// m = 2 test-case
void VecTransImposedData2d(double x[3], double t, double* w) {
  double vx 
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
  w[1] = xx * xx;
};


void TestTransBoundaryFlux(double x[3], double t, 
			       double wL[], double* vnorm,
			       double* flux) {
  double wR[1];
  TestTransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
};

void TestTransBoundaryFlux2d(double x[3], double t, double wL[],
				 double* vnorm, double* flux) {
  double wR[1];
  TestTransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
};

void TestTransInitData(double x[3], double w[]) {
  double t = 0;
  TestTransImposedData(x, t, w);
};

void TestTransInitData2d(double x[3], double w[]) {
  double t = 0;
  TestTransImposedData2d(x, t, w);
};

void TestTransImposedData(double x[3], double t, double w[]) {
  double vx 
    = transport_v[0] * x[0] 
    + transport_v[1] * x[1] 
    + transport_v[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

void TestTransImposedData2d(double x[3], double t, double w[]) {
  double vx 
    = transport_v2d[0] * x[0] 
    + transport_v2d[1] * x[1] 
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

// Vlasov 2D transport equation functions

// Set the parameters for the Vlasov equation that are stored in the
// global space of model.h
void set_vlasov_params(int mx0, int my0, int mz0, double vmax0) 
{
  mx = mx0;
  my = my0;
  mz = mz0;
  vmax = vmax0;
}

void vlaTransInitData2d(double x[3], double w[]) 
{
  double t = 0;
  vlaTransImposedData2d(x, t, w);
};

void vlaTransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux) 
{
  int im = 0;
  for(int ix = 0; ix < mx; ++ix) {
    double vx = vmax * (ix - (mx / 2));

    for(int iy = 0; iy < my; ++iy) {
      double vy = vmax * (iy - (my / 2));
      
      double vn = vx * vnorm[0]	+ vy * vnorm[1];
      double vnp = vn > 0 ? vn : 0;
      double vnm = vn - vnp;
      
      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      flux[im] = vnp * wL[im] + vnm * wR[im];

      ++im;
    }
  }
  // Verify that 2d computations are actually activated
  assert(fabs(vnorm[2]) < 1e-8);
};

void vlaTransBoundaryFlux2d(double x[3], double t, 
			    double wL[], double* vnorm,
			    double* flux) 
{
  double wR[2]; // FIXME implicit assumption m = 2
  vlaTransImposedData2d(x, t, wR);
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
};

// 6th-degree polynomial with compact support
double compact_poly6(double r)
{
  if (fabs(2 * r) > 1)
    return 0;
  double rrm1 = 2 * r - 1;
  double rrp1 = 2 * r + 1;
  return -35.0 / 16.0 * rrm1 * rrm1 * rrm1 * rrp1 * rrp1 * rrp1;
}

// Impose a compact-supporrt 6th degree polynomial in 2+2D
void vlaTransImposedData2d(double x[3], double t, double* w) 
{
  double PI = 4.0 * atan(1.0);
  double s2pi = sqrt(2.0 * PI);
  double xval = 1.0;

  printf("mx: %d\n", mx);
  printf("my: %d\n", my);

  double r = sqrt(x[0] * x[0] + x[1] * x[1]);
  double pr = compact_poly6(r);

  int im = 0;
  for(int ix = 0; ix < mx; ++ix) {
    double vx = vmax * (ix - (mx / 2));

    for(int iy = 0; iy < my; ++iy) {
      double vy = vmax * (iy - (my / 2));

      double vr = sqrt(vx * vx + vy * vy);
      double pvr = compact_poly6(vr);

      printf("im: %d", im);
      w[im++] = pr * pvr;
    }
  }
};

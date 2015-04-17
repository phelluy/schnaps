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
}

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
}

#pragma start_opencl
void VecTransNumFlux2d(double *wL, double *wR, double *vnorm, double *flux) 
{
  double vn = sqrt(0.5) * vnorm[0] + sqrt(0.5) * vnorm[1];
  double vnp = vn > 0.0 ? vn : 0.0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  flux[1] = vnp * wL[1] + vnm * wR[1];
}
#pragma end_opencl

void TransBoundaryFlux(double x[3], double t, double wL[], double *vnorm,
		       double *flux)
{
  double wR[1];
  TransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
}

void TransBoundaryFlux2d(double x[3], double t, double wL[], double *vnorm,
			 double *flux) {
  double wR[1];
  TransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
}

// m = 2 test-case
#pragma start_opencl
void VecTransImposedData2d(double *x, double t, double* w) 
{
  double vx = sqrt(0.5) * x[0] + sqrt(0.5) * x[1];
  double xx = vx - t;
  w[0] = xx * xx;
  w[1] = 2 * xx * xx;
}
#pragma end_opencl

#pragma start_opencl
void VecTransBoundaryFlux2d(double *x, double t, 
			    double *wL, double *vnorm,
			    double *flux) 
{
  double wR[2];
  VecTransImposedData2d(x, t, wR);
  VecTransNumFlux2d(wL, wR, vnorm, flux);
}
#pragma end_opencl

void TransInitData(double x[3], double w[]) {
  double t = 0;
  TransImposedData(x, t, w);
}

void TransInitData2d(double x[3], double w[]) {
  double t = 0;
  TransImposedData2d(x, t, w);
}

void VecTransInitData2d(double x[3], double w[]) {
  double t = 0;
  VecTransImposedData2d(x, t, w);
}

void TransImposedData(double x[3], double t, double w[]) {
  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];
  double xx = vx - t;
  w[0]=cos(xx);
}

void TransImposedData2d(double *x, double t, double *w) 
{
  double vx  = sqrt(0.5) * x[0] + sqrt(0.5) * x[1];
  double xx = vx - t;
  w[0] = cos(xx);
}

void TestTransBoundaryFlux(double x[3], double t, 
			       double wL[], double* vnorm,
			       double* flux) {
  double wR[1];
  TestTransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
}

void TestTransBoundaryFlux2d(double x[3], double t, double wL[],
				 double* vnorm, double* flux) {
  double wR[1];
  TestTransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
}

void TestTransInitData(double x[3], double w[]) {
  double t = 0;
  TestTransImposedData(x, t, w);
}

void TestTransInitData2d(double x[3], double w[]) {
  double t = 0;
  TestTransImposedData2d(x, t, w);
}

void TestTransImposedData(double x[3], double t, double w[]) {
  double vx 
    = transport_v[0] * x[0] 
    + transport_v[1] * x[1] 
    + transport_v[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
}

void TestTransImposedData2d(double x[3], double t, double w[]) {
  double vx 
    = transport_v2d[0] * x[0] 
    + transport_v2d[1] * x[1] 
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
}

// Vlasov 2D transport equation functions

// Set the parameters for the Vlasov equation that are stored in the
// global space of model.h
void set_vlasov_params(Model *mod) 
{
  m = mod->m;
  assert(m > 0);
  vlasov_mx = mod->vlasov_mx;
  vlasov_my = mod->vlasov_my;
  vlasov_mz = mod->vlasov_mz;
  assert(m == vlasov_mx * vlasov_my * vlasov_mz);
  vlasov_vmax = mod->vlasov_vmax;
}

void vlaTransInitData2d(double x[3], double w[]) 
{
  double t = 0;
  vlaTransImposedData2d(x, t, w);
}

// Return the component of the vlasov velocity with index id.
double vlasov_vel(const int id, const int md, double vlasov_vmax)
{
  int mid = md / 2;
  double dv = vlasov_vmax / mid;
  return (id - mid) * dv;
}

void vlaTransNumFlux2d(double wL[], double wR[], double *vnorm, double *flux) 
{
  for(int ix = 0; ix < vlasov_mx; ++ix) {
    double vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);

    for(int iy = 0; iy < vlasov_my; ++iy) {
      double vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      
      double vn = vx * vnorm[0]	+ vy * vnorm[1];
      double vnp = vn > 0 ? vn : 0;
      double vnm = vn - vnp;
      
      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      flux[im] = vnp * wL[im] + vnm * wR[im];
    }
  }
  // Verify that 2d computations are actually activated
  assert(fabs(vnorm[2]) < 1e-8);
}

void vlaTransBoundaryFlux2d(double x[3], double t, 
			    double wL[], double *vnorm,
			    double* flux) 
{
  double wR[m];
  vlaTransImposedData2d(x, t, wR);
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

// compact support bump (C-infinity, but not analytic):
double compact_bump(double r)
{
  if(fabs(r) >= 0.5)
    return 0;
  return exp(-1.0 / (1.0 - 4.0 * r * r));
}

double icgaussian(double r, double sigma)
{
  double PI = 4.0 * atan(1.0);
  return exp(-r * r / sigma) / (sigma * sqrt(2.0 * PI));
}

// 6th-degree polynomial with compact support
double compact_poly6(double r)
{
  if (fabs(2 * r) > 1)
    return 0;
  double rrm1 = 2 * r - 1;
  double rrp1 = 2 * r + 1;
  return -35.0 / 16.0 * rrm1 * rrm1 * rrm1 * rrp1 * rrp1 * rrp1;
}

void vlaTransImposedData2d(double x[3], double t, double *w) 
{
  double PI = 4.0 * atan(1.0);
  double s2pi = sqrt(2.0 * PI);
  double xval = 1.0;
  double sigma = 0.1;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    double vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    double px = x[0] - vx * t;

    for(int iy = 0; iy < vlasov_my; ++iy) {
      double vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      double py = x[1] - vy * t;

      double r = sqrt(px * px + py * py);
      double pi = 4.0 * atan(1.0);
      double pr = icgaussian(r, sigma);

      double vr = sqrt(vx * vx + vy * vy);
      double pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemcracs2014_imposed_data(double x[3], double t, double *w)
{
  double sigma = 1.0;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    double vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    double px = x[0] - vx * t;

    for(int iy = 0; iy < vlasov_my; ++iy) {
      double vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      double py = x[1] - vy * t;

      double r = sqrt(px * px + py * py);
      double pr = compact_bump(r);
      //double pr = compact_poly6(r);
      
      double vr = sqrt(vx * vx + vy * vy);
      double pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemcracs2014a_imposed_data(double x[3], double t, double *w)
{
  double sigma = 1.0;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    double vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    double px = x[0] - vx * t;

    for(int iy = 0; iy < vlasov_my; ++iy) {
      double vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      double py = x[1] - vy * t;

      double r = sqrt(px * px + py * py);
      //double pr = compact_bump(r);
      double pr = compact_poly6(r);
      
      double vr = sqrt(vx * vx + vy * vy);
      double pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemracs2014_TransInitData(double x[3], double w[]) 
{
  double t = 0;
  cemcracs2014_imposed_data(x, t, w);
}

void cemracs2014a_TransInitData(double x[3], double w[]) 
{
  double t = 0;
  cemcracs2014a_imposed_data(x, t, w);
}

void cemracs2014_TransBoundaryFlux(double x[3], double t, 
				   double wL[], double *vnorm,
				   double *flux) 
{
  double wR[m];
  for(unsigned int i = 0; i < m; ++i)
    wR[i] = 0;
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

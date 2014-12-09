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

void TransNumFlux(double wL[], double wR[], 
		      double* vnorm, double* flux, 
		      const int mx, const int my, const int mz, 
		      const double vmax) {
  double vn 
    = transport_v[0] * vnorm[0] 
    + transport_v[1] * vnorm[1]
    + transport_v[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
};

void TransNumFlux2d(double wL[], double wR[], 
			double* vnorm, double* flux, 
			const int mx, const int my, const int mz, 
			const double vmax) {
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

void vTransNumFlux2d(double wL[], double wR[], 
			 double* vnorm, double* flux,
			 const int mx, const int my, const int mz, 
			 const double vmax) 
{
  int im = 0;
  for(int ix = 0; ix < mx; ++ix) {
    double vx = vmax * (ix - (mx / 2));
    for(int iy = 0; iy < my; ++iy) {
      double vy = vmax * (iy - (my / 2));
      // The 2D velocity is (vx, vy)
      
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

void VecTransNumFlux2d(double wL[], double wR[], 
		       double* vnorm, double* flux,
		       const int mx, const int my, const int mz, 
		       const double vmax) {
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

void TransBoundaryFlux(double x[3], double t, 
			   double wL[], double* vnorm,
			   double* flux,
			   int mx, int my, int mz, double vmax) {
  double wR[1];
  TransImposedData(x, t, wR, mx, my, mz, vmax);
  TransNumFlux(wL, wR, vnorm, flux, 0, 0, 0, 0.0);
};

void TransBoundaryFlux2d(double x[3], double t, 
			     double wL[], double* vnorm,
			     double* flux,
			     int mx, int my, int mz, double vmax) {
  double wR[1];
  TransImposedData2d(x, t, wR, mx, my, mz, vmax);
  TransNumFlux2d(wL, wR, vnorm, flux, mx, my, mz, vmax);
};

void VecTransBoundaryFlux2d(double x[3], double t, 
			    double wL[], double* vnorm,
			    double* flux,
			    int mx, int my, int mz, double vmax) {
  double wR[2];
  VecTransImposedData2d(x, t, wR, mx, my, mz, vmax);
  VecTransNumFlux2d(wL, wR, vnorm, flux, mx, my, mz, vmax);
};

void TransInitData(double x[3], double w[], 
		       int mx, int my, int mz, double vmax) {
  double t = 0;
  TransImposedData(x, t, w, mx, my, mz, vmax);
};

void TransInitData2d(double x[3], double w[], 
			 int mx, int my, int mz, double vmax) {
  double t = 0;
  TransImposedData2d(x, t, w, mx, my, mz, vmax);
};

void vTransInitData2d(double x[3], double w[], 
			  int mx, int my, int mz, double vmax) {
  double t = 0;
  vTransImposedData2d(x, t, w, mx, my, mz, vmax);
};

void VecTransInitData2d(double x[3], double w[], 
			int mx, int my, int mz, double vmax) {
  double t = 0;
  VecTransImposedData2d(x, t, w, mx, my, mz, vmax);
};

void TransImposedData(double x[3], double t, double w[], 
			  int mx, int my, int mz, double vmax) {
  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];
  double xx = vx - t;
  w[0]=cos(xx);
};

void TransImposedData2d(double x[3], double t, double w[], 
			    int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = cos(xx);
};

// m = 2 test-case
void VecTransImposedData2d(double x[3], double t, double* w, 
			  int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
  w[1] = xx * xx;
};

// zero-centered Gaussian distribution in both x and v.
void vTransImposedData2d(double x[3], double t, double* w, 
			 int mx, int my, int mz, double vmax) {
  double PI = 4.0 * atan(1.0);
  double s2pi = sqrt(2.0 * PI);
  double xval = 1.0;


  int im = 0;
  for(int ix = 0; ix < mx; ++ix) {
    double vx = vmax * (ix - (mx / 2));
        
    for(int iy = 0; iy < my; ++iy) {
      double vy = vmax * (iy - (my / 2));
      // The 2D velocity is (vx, vy)
      double v = vx * x[0] + vy * x[1];
      double xx = v - t;
      w[im++] = xx * xx; // FIXME: choose a function.
    }
  }
  //v = vx - y
  //f = (1 + (- 3 + (3 - v * v ) * v * v) * v * v) * hh
};

void TestTransBoundaryFlux(double x[3], double t, 
			       double wL[], double* vnorm,
			       double* flux, 
			       int mx, int my, int mz, double vmax) {
  double wR[1];
  TestTransImposedData(x, t, wR, mx, my, mz, vmax);
  TransNumFlux(wL, wR, vnorm, flux, 0 ,0, 0, 0.0);
};

void TestTransBoundaryFlux2d(double x[3], double t, double wL[],
				 double* vnorm, double* flux, 
				 int mx, int my, int mz, double vmax) {
  double wR[1];
  TestTransImposedData2d(x, t, wR, mx, my, mz, vmax);
  TransNumFlux2d(wL, wR, vnorm, flux, 0 ,0, 0, 0.0);
};

void TestTransInitData(double x[3], double w[], 
			   int mx, int my, int mz, double vmax) {
  double t = 0;
  TestTransImposedData(x, t, w, mx, my, mz, vmax);
};

void TestTransInitData2d(double x[3], double w[], 
			     int mx, int my, int mz, double vmax) {
  double t = 0;
  TestTransImposedData2d(x, t, w, mx, my, mz, vmax);
};

void TestTransImposedData(double x[3], double t, double w[], 
			      int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v[0] * x[0] 
    + transport_v[1] * x[1] 
    + transport_v[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

void TestTransImposedData2d(double x[3], double t, double w[], 
				int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v2d[0] * x[0] 
    + transport_v2d[1] * x[1] 
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

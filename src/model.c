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

void TransportNumFlux(double wL[], double wR[], 
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

void TransportNumFlux2d(double wL[], double wR[], 
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

void vTransportNumFlux2d(double wL[], double wR[], 
			 double* vnorm, double* flux,
			 const int mx, const int my, const int mz, 
			 const double vmax) 
{
  // FIXME: how to get mx, my, and vmax from the model to here?
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

void TransportBoundaryFlux(double x[3], double t, 
			   double wL[], double* vnorm,
			   double* flux,
			   int mx, int my, int mz, double vmax) {
  double wR[1];
  TransportImposedData(x, t, wR, mx, my, mz, vmax);
  TransportNumFlux(wL, wR, vnorm, flux, 0, 0, 0, 0.0);
};

void TransportBoundaryFlux2d(double x[3], double t, 
			     double wL[], double* vnorm,
			     double* flux,
			     int mx, int my, int mz, double vmax) {
  double wR[1];
  TransportImposedData2d(x, t, wR, mx, my, mz, vmax);
  TransportNumFlux2d(wL, wR, vnorm, flux, mx, my, mz, vmax);
};

void VecTransBoundaryFlux2d(double x[3], double t, 
			    double wL[], double* vnorm,
			    double* flux,
			    int mx, int my, int mz, double vmax) {
  double wR[2];
  VecTransImposedData2d(x, t, wR, mx, my, mz, vmax);
  VecTransNumFlux2d(wL, wR, vnorm, flux, mx, my, mz, vmax);
};

void TransportInitData(double x[3], double w[], 
		       int mx, int my, int mz, double vmax) {
  double t = 0;
  TransportImposedData(x, t, w, mx, my, mz, vmax);
};

void TransportInitData2d(double x[3], double w[], 
			 int mx, int my, int mz, double vmax) {
  double t = 0;
  TransportImposedData2d(x, t, w, mx, my, mz, vmax);
};

void vTransportInitData2d(double x[3], double w[], 
			  int mx, int my, int mz, double vmax) {
  double t = 0;
  // TransportImposedData2d(x, t, w);
  // FIXME
};

void VecTransInitData2d(double x[3], double w[], 
			int mx, int my, int mz, double vmax) {
  double t = 0;
  VecTransImposedData2d(x, t, w, mx, my, mz, vmax);
};

void TransportImposedData(double x[3], double t, double w[], 
			  int mx, int my, int mz, double vmax) {
  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];
  double xx = vx - t;
  w[0]=cos(xx);
};

void TransportImposedData2d(double x[3], double t, double w[], 
			    int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = cos(xx);
};

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

void TestTransportBoundaryFlux(double x[3], double t, 
			       double wL[], double* vnorm,
			       double* flux, 
			       int mx, int my, int mz, double vmax) {
  double wR[1];
  TestTransportImposedData(x, t, wR, mx, my, mz, vmax);
  TransportNumFlux(wL, wR, vnorm, flux, 0 ,0, 0, 0.0);
};

void TestTransportBoundaryFlux2d(double x[3], double t, double wL[],
				 double* vnorm, double* flux, 
				 int mx, int my, int mz, double vmax) {
  double wR[1];
  TestTransportImposedData2d(x, t, wR, mx, my, mz, vmax);
  TransportNumFlux2d(wL, wR, vnorm, flux, 0 ,0, 0, 0.0);
};

void TestTransportInitData(double x[3], double w[], 
			   int mx, int my, int mz, double vmax) {
  double t = 0;
  TestTransportImposedData(x, t, w, mx, my, mz, vmax);
};

void TestTransportInitData2d(double x[3], double w[], 
			     int mx, int my, int mz, double vmax) {
  double t = 0;
  TestTransportImposedData2d(x, t, w, mx, my, mz, vmax);
};

void TestTransportImposedData(double x[3], double t, double w[], 
			      int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v[0] * x[0] 
    + transport_v[1] * x[1] 
    + transport_v[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

void TestTransportImposedData2d(double x[3], double t, double w[], 
				int mx, int my, int mz, double vmax) {
  double vx 
    = transport_v2d[0] * x[0] 
    + transport_v2d[1] * x[1] 
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

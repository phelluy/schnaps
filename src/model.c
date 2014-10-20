#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)
/* const double transport_v[] = { */
/*   ONE_OVER_SQRT_3, */
/*   ONE_OVER_SQRT_3, */
/*   ONE_OVER_SQRT_3}; */


const double transport_v[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};


void TransportNumFlux(double wL[],double wR[],double* vnorm,double* flux){
  
  double vn =
    transport_v[0] * vnorm[0] +
    transport_v[1] * vnorm[1] +
    transport_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];
};

void TransportBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TransportImposedData(x,t,wR);
  TransportNumFlux(wL,wR,vnorm,flux);
};


void TransportInitData(double x[3],double w[]){

  double t=0;
  TransportImposedData(x,t,w);

};

void TransportImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};


void TestTransportBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestTransportImposedData(x,t,wR);
  TransportNumFlux(wL,wR,vnorm,flux);
};


void TestTransportInitData(double x[3],double w[]){

  double t=0;
  TestTransportImposedData(x,t,w);

};

void TestTransportImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};

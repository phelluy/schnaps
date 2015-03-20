#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>



// Centered flux if eps=0, uncentered flux if eps=1
void Maxwell2DNumFlux(double wL[], double wR[], double* vnorm, double* flux) {

  double r = sqrt(vnorm[0]*vnorm[0]+vnorm[1]*vnorm[1]);
  double overr = 1.0 / (r+1e-14);
  double nx=vnorm[0];
  double ny=vnorm[1];
  double eps=1;
#define _KHI 1

  flux[0] = - ny * (wR[2]+wL[2]) + _KHI * nx * (wR[3]+wL[3])
            - eps * (ny * ny + _KHI * nx * nx) * overr * (wR[0]-wL[0])
            - eps * (ny * nx * (_KHI - 1)) * overr * (wR[1]-wL[1]);

  flux[1] =   nx * (wR[2]+wL[2]) + _KHI * ny * (wR[3]+wL[3])
            - eps * (ny * nx * (_KHI - 1)) * overr * (wR[0]-wL[0])
            - eps * (nx * nx + _KHI * ny * ny)* overr  * (wR[1]-wL[1]);

  flux[2] = - ny * (wR[0]+wL[0]) + nx * (wR[1]+wL[1]) - eps * r * (wR[2]-wL[2]);

  flux[3] = _KHI * nx * (wR[0]+wL[0]) + _KHI * ny * (wR[1]+wL[1]) - eps * _KHI * r * (wR[3]-wL[3]);

  flux[0]/=2;
  flux[1]/=2;
  flux[2]/=2;
  flux[3]/=2;




}


void Maxwell2DImposedData(double x[3], double t, double w[]) {

    double c,u,v,theta,pi;
    double un=1;
    pi=4*atan(un);
    //pi=_PI;
    double k=2*pi;
    k=0;
    double r=1;
    k=2*pi;
    theta=pi/4;

    u=cos(theta);
    v=sin(theta); 

    k=2*pi/v;
  
    c=-cos(k * (u * x[0] + v * x[1] - t));
    
    w[0]= -v * c/r;
    w[1]= u * c/r;
    w[2]= c/r;
    w[3]= 0;

}


void Maxwell2DBoundaryFlux(double x[3], double t, double wL[], double *vnorm,
		       double *flux) {
  double wR[4];
  Maxwell2DImposedData(x, t, wR);
  Maxwell2DNumFlux(wL, wR, vnorm, flux);
}


void Maxwell2DInitData(double x[3], double w[]) {
  double t = 0;
  Maxwell2DImposedData(x, t, w);
}



#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define _KHI 1


// Centered flux if eps=0, uncentered flux if eps=1
void Maxwell2DNumFlux(real wL[], real wR[], real* vnorm, real* flux) {

  real r = sqrt(vnorm[0]*vnorm[0]+vnorm[1]*vnorm[1]);
  real overr = 1.0 / (r+1e-14);
  real nx=vnorm[0];
  real ny=vnorm[1];
  real eps=1;

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
  flux[4]=0;
  flux[5]=0;
  flux[6]=0;




}


void Maxwell2DImposedData(real x[3], real t, real w[]) {

    real c,u,v,theta,pi;
    real un=1;
    pi=4*atan(un);
    //pi=_PI;
    real k=2*pi;
    k=0;
    real r=1;
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
    w[4]=0;
    w[5]=0;
    w[6]=0;
}


void Maxwell2DBoundaryFlux(real x[3], real t, real wL[], real *vnorm,
		       real *flux) {
  real wR[7];
  Maxwell2DImposedData(x, t, wR);
  Maxwell2DNumFlux(wL, wR, vnorm, flux);
}


void Maxwell2DInitData(real x[3], real w[]) {
  real t = 0;
  Maxwell2DImposedData(x, t, w);
}

void Maxwell2DSource(real x[3],real t,real w[],real source[]){

  source[0]=w[4];
  source[1]=w[5];
  source[2]=0;
  source[3]= _KHI * w[6]; 
  source[4]=0;
  source[5]=0;
  source[6]=0;


}



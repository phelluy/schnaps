#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

// Centered flux if eps=0, uncentered flux if eps=1
#pragma start_opencl
void Maxwell2DNumFlux(real *wL, real *wR, real *vnorm, 
			real *flux) 
{
  real r = sqrt(vnorm[0] *vnorm[0] + vnorm[1] * vnorm[1]);
  real overr = 1.0 / (r + 1e-14);
  real nx = vnorm[0];
  real ny = vnorm[1];
  real eps = 1;
  real khi = 1.0;

  flux[0] = 
    - ny * (wR[2] + wL[2]) + khi * nx * (wR[3] + wL[3])
    - eps * (ny * ny + khi * nx * nx) * overr * (wR[0] - wL[0])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[1] - wL[1]);

  flux[1] =   
    nx * (wR[2]+wL[2]) + khi * ny * (wR[3]+wL[3])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[0]-wL[0])
    - eps * (nx * nx + khi * ny * ny)* overr  * (wR[1]-wL[1]);

  flux[2] = - ny * (wR[0] + wL[0]) 
    + nx * (wR[1] + wL[1]) 
    - eps * r * (wR[2] - wL[2]);

  flux[3] = 
    khi * nx * (wR[0]+wL[0]) 
    + khi * ny * (wR[1]+wL[1]) 
    - eps * khi * r * (wR[3]-wL[3]);

  flux[0] *= 0.5;
  flux[1] *= 0.5;
  flux[2] *= 0.5;
  flux[3] *= 0.5;
  flux[4] = 0;
  flux[5] = 0;
  flux[6] = 0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DImposedData(const real *x, const real t, real *w) 
{
  real pi = 4.0 * atan(1.0);
  real r = 1.0;
  real theta = pi / 4.0;
  real u = cos(theta);
  real v = sin(theta); 
  real k = 2.0 * pi / v;
  real c = -cos(k * (u * x[0] + v * x[1] - t));
  
  w[0] = -v * c / r;
  w[1] = u * c / r;
  w[2] = c / r;
  w[3] = 0;
  w[4] = 0;
  w[5] = 0;
  w[6] = 0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DBoundaryFlux(real *x, real t, 
			   real *wL, real *vnorm, real *flux)
{
  real wR[7];
  Maxwell2DImposedData(x, t, wR);
  Maxwell2DNumFlux(wL, wR, vnorm, flux);
}
#pragma end_opencl

void Maxwell2DInitData(real *x, real *w) 
{
  real t = 0;
  Maxwell2DImposedData(x, t, w);
}

#pragma start_opencl
void Maxwell2DSource(const real *x, const real t, const real *w, real *source)
{
  real khi = 1.0;
  source[0] = w[4];
  source[1] = w[5];
  source[2] = 0;
  source[3] = khi * w[6]; 
  source[4] = 0;
  source[5] = 0;
  source[6] = 0;
}
#pragma end_opencl


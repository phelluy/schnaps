#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#pragma start_opencl
void Maxwell2DNumFlux_centered(real *wL, real *wR, real *vnorm, 
			       real *flux) 
{
  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real khi = 1.0;

  const real s0 = wR[0] + wL[0];
  const real s1 = wR[1] + wL[1];
  const real s2 = wR[2] + wL[2];
  const real s3 = wR[3] + wL[3];

  flux[0] = 0.5 * (-ny * s2 + khi * nx * s3);
  flux[1] = 0.5 * ( nx * s2 + khi * ny * s3);
  flux[2] = 0.5 * (-ny * s0 + nx * s1);
  flux[3] = 0.5 * khi * (nx * s0 + ny * s1);
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DNumFlux_uncentered(real *wL, real *wR, real *vnorm, 
			real *flux) 
{
  const real r = sqrt(vnorm[0] * vnorm[0] + vnorm[1] * vnorm[1]);
  const real overr = 1.0 / (r + 1e-16);
  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real khi = 1.0;

  const real s0 = wR[0] + wL[0];
  const real s1 = wR[1] + wL[1];
  const real s2 = wR[2] + wL[2];
  const real s3 = wR[3] + wL[3];

  const real d0 = wR[0] - wL[0];
  const real d1 = wR[1] - wL[1];
  const real d2 = wR[2] - wL[2];
  const real d3 = wR[3] - wL[3];

  flux[0] = 0.5 * (-ny * s2 + khi * nx * s3
		   -overr * ( d0 * (ny * ny + khi * nx * nx)
			      + d1 * nx * ny * (khi - 1)
			      )
		   );
  flux[1] = 0.5 * (nx * s2 + khi * ny * s3
		   -overr * ((nx * ny * (khi - 1)) * d0
			     +(nx * nx + khi * ny * ny) * d1)
		   );
  flux[2] = 0.5 * (-ny * s0 + nx * s1 - r * d2);
  flux[3] = 0.5 * khi * (nx * s0 + ny * s1 -r * d3);
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DImposedData(const real *x, const real t, real *w) 
{
  const real pi = 4.0 * atan(1.0);
  const real r = 1.0;
  const real theta = pi / 4.0;
  const real u = cos(theta);
  const real v = sin(theta); 
  const real k = 2.0 * pi / v;
  const real c = -cos(k * (u * x[0] + v * x[1] - t));
  
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
  Maxwell2DNumFlux_centered(wL, wR, vnorm, flux);
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
  const real khi = 1.0;
  source[0] = w[4];
  source[1] = w[5];
  source[2] = 0;
  source[3] = khi * w[6]; 
  source[4] = 0;
  source[5] = 0;
  source[6] = 0;
}
#pragma end_opencl

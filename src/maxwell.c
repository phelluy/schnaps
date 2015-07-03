#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#pragma start_opencl
void Maxwell2DNumFlux_centered(real *wL, real *wR, real *vnorm, real *flux) 
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation
  
  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real khi = 1.0;

  // FIXME: improve names for these (refer to E, H, etc).
  const real s0 = 0.5 * ( wR[0] + wL[0] );
  const real s1 = 0.5 * ( wR[1] + wL[1] );
  const real s2 = 0.5 * ( wR[2] + wL[2] );
  const real s3 = 0.5 * ( wR[3] + wL[3] );

  flux[0] = -ny * s2 + khi * nx * s3;
  flux[1] =  nx * s2 + khi * ny * s3;
  flux[2] = -ny * s0 + nx * s1;
  flux[3] = khi * (nx * s0 + ny * s1);
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DNumFlux_uncentered(real *wL, real *wR, real *vnorm, real *flux) 
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real r = sqrt(nx * nx + ny * ny);
  const real overr = 1.0 / (r + 1e-16);
  const real khi = 1.0;

  const real s0 = 0.5 * ( wR[0] + wL[0] );
  const real s1 = 0.5 * ( wR[1] + wL[1] );
  const real s2 = 0.5 * ( wR[2] + wL[2] );
  const real s3 = 0.5 * ( wR[3] + wL[3] );

  const real d0 = 0.5 * ( wR[0] - wL[0] );
  const real d1 = 0.5 * ( wR[1] - wL[1] );
  const real d2 = 0.5 * ( wR[2] - wL[2] );
  const real d3 = 0.5 * ( wR[3] - wL[3] );

  flux[0] = 
    -ny * s2 + khi * nx * s3
    -overr * ( d0 * (ny * ny + khi * nx * nx)
	       + d1 * nx * ny * (khi - 1) );
  flux[1] = 
    nx * s2 + khi * ny * s3
    -overr * ((nx * ny * (khi - 1)) * d0
	      +(nx * nx + khi * ny * ny) * d1);
  flux[2] = -ny * s0 + nx * s1 - r * d2;
  flux[3] = khi * (nx * s0 + ny * s1 -r * d3);
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DNumFlux_unoptimised(real *wL, real *wR, real *vnorm, real *flux) 
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real r = sqrt(nx * nx + ny * ny);
  const real overr = 1.0 / (r + 1e-16);
  // Centered flux if eps=0, uncentered flux if eps=1
  const real eps = 1;
  const real khi = 1.0;

  flux[0] = 
    - ny * (wR[2] + wL[2]) + khi * nx * (wR[3] + wL[3])
    - eps * (ny * ny + khi * nx * nx) * overr * (wR[0] - wL[0])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[1] - wL[1]);

  flux[1] =   
    nx * (wR[2] + wL[2]) + khi * ny * (wR[3] + wL[3])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[0] - wL[0])
    - eps * (nx * nx + khi * ny * ny) * overr  * (wR[1] - wL[1]);

  flux[2] = - ny * (wR[0] + wL[0]) 
    + nx * (wR[1] + wL[1]) 
    - eps * r * (wR[2] - wL[2]);

  flux[3] = 
    khi * nx * (wR[0] + wL[0]) 
    + khi * ny * (wR[1] + wL[1]) 
    - eps * khi * r * (wR[3] - wL[3]);

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
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation
  
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
void Maxwell2DBoundaryFlux_uncentered(real *x, real t, 
				      real *wL, real *vnorm, real *flux)
{
  real wR[7];
  Maxwell2DImposedData(x, t, wR);
  Maxwell2DNumFlux_uncentered(wL, wR, vnorm, flux);
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
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)
  
  // FIXME add documentation
  
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

#pragma start_opencl
void Maxwell3DNumFlux_uncentered(real *wL, real *wR, real *vnorm, real *flux) 
{
  // Uncentered flux (upwind) for Maxwell's equations

  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz}

  // Let {{E}} = ( ER + EL ) / 2, [[E]] = ( ER - EL ) / 2 
  // The first three components of the flux are
  // - n x {{H}} + n x n x [[E]] / r 
  // and the last three are
  //   n x {{E}} + n x n x [[H]] / r 
  
  const real n0 = vnorm[0];
  const real n1 = vnorm[1];
  const real n2 = vnorm[2];

  const real overr = 1.0 / ( sqrt( n0 * n0 + n1 * n1 + n2 * n2 ) + 1e-16 );
  const real n01 = overr * n0 * n1;
  const real n02 = overr * n0 * n2;
  const real n12 = overr * n1 * n2;
  const real n00 = overr * n0 * n0;
  const real n11 = overr * n1 * n1;
  const real n22 = overr * n2 * n2;
    
  const real Es0 = 0.5 * ( wR[0] + wL[0] );
  const real Es1 = 0.5 * ( wR[1] + wL[1] );
  const real Es2 = 0.5 * ( wR[2] + wL[2] );

  const real Hs0 = 0.5 * ( wR[3] + wL[3] );
  const real Hs1 = 0.5 * ( wR[4] + wL[4] );
  const real Hs2 = 0.5 * ( wR[5] + wL[5] );

  const real Ed0 = 0.5 * ( wR[0] - wL[0] );
  const real Ed1 = 0.5 * ( wR[1] - wL[1] );
  const real Ed2 = 0.5 * ( wR[2] - wL[2] );

  const real Hd0 = 0.5 * ( wR[3] - wL[3] );
  const real Hd1 = 0.5 * ( wR[4] - wL[4] );
  const real Hd2 = 0.5 * ( wR[5] - wL[5] );

  // E flux
  flux[0] = n01 * Ed2 - n22 * Ed0 - n22 * Ed0 + n02 * Ed2 
    - n1 * Hs2 + n2 * Hs1;
  flux[1] = n12 * Ed2 - n22 * Ed1 - n11 * Ed2 + n01 * Ed0 
    - n2 * Hs0 + n0 * Hs2;
  flux[2] = n02 * Ed0 - n00 * Ed2 - n11 * Ed2 + n12 * Ed1 
    - n0 * Hs1 + n1 * Hs0;

  // H flux
  flux[3] = n01 * Hd2 - n22 * Hd0 - n22 * Hd0 + n02 * Hd2 
    + n1 * Es2 - n2 * Es1;
  flux[4] = n12 * Hd2 - n22 * Hd1 - n11 * Hd2 + n01 * Hd0 
    + n2 * Es0 - n0 * Es2;
  flux[5] = n02 * Hd0 - n00 * Hd2 - n11 * Hd2 + n12 * Hd1 
    + n0 * Es1 - n1 * Es0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DNumFluxClean_uncentered(real *wL, real *wR, real *vnorm, 
				      real *flux) 
{
  // Uncentered flux (upwind) for Maxwell's equations

  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz, lambda_E, lambda_H}

  // Let {{E}} = ( ER + EL ) / 2, [[E]] = ( ER - EL ) / 2 
  // The E flux is:
  // - n x {{H}} + n x n x [[E]] / r 
  //   + c1 n {{lambda_E}} + c1 n ( n . [[E]] ) / r 
  // The H flux is:
  //   n x {{E}} + n x n x [[H]] / r 
  //   + c2 n {{lambda_H}} + c2 n ( n . [[H]] ) / r 
  // The lambda_E flux is
  // c1 * ( n . {{E}} + r [[lambda_e]] )
  // The lambda_H flux is
  // c2 * ( n . {{H}} + r [[lambda_H]] )

  // We first compute the uncleaned flux, and then add the cleaning
  // (which is cleaner, though perhaps uses a few extra operations).

  // Add the normal uncentered flux.
  Maxwell3DNumFlux_uncentered(wL, wR, vnorm, flux);

  // FIXME: how do we set these?  What are good values?
  const real c1 = 1.0; // E-cleaning parameter
  const real c2 = 1.0; // H-cleaning parameter

  // Consts based on vnorm
  const real n0 = vnorm[0];
  const real n1 = vnorm[1];
  const real n2 = vnorm[2];
  const real r = sqrt( n0 * n0 + n1 * n1 + n2 * n2 );
  
  // Consts based on the mean
  const real Es0 = 0.5 * ( wR[0] + wL[0] );
  const real Es1 = 0.5 * ( wR[1] + wL[1] );
  const real Es2 = 0.5 * ( wR[2] + wL[2] );

  const real Hs0 = 0.5 * ( wR[3] + wL[3] );
  const real Hs1 = 0.5 * ( wR[4] + wL[4] );
  const real Hs2 = 0.5 * ( wR[5] + wL[5] );

  const real lEs = 0.5 * ( wR[6] + wL[6] );
  const real lHs = 0.5 * ( wR[7] + wL[7] );

  // Consts based on the jump
  const real Ed0 = 0.5 * ( wR[0] - wL[0] );
  const real Ed1 = 0.5 * ( wR[1] - wL[1] );
  const real Ed2 = 0.5 * ( wR[2] - wL[2] );

  const real Hd0 = 0.5 * ( wR[3] - wL[3] );
  const real Hd1 = 0.5 * ( wR[4] - wL[4] );
  const real Hd2 = 0.5 * ( wR[5] - wL[5] );

  const real lEd = 0.5 * ( wR[6] - wL[6] );
  const real lHd = 0.5 * ( wR[7] - wL[7] );

  // Add correction term to E flux
  // c_1 * (n \cdot \jump{E} + \mean{\lambda_E} )
  const real Ec = c1 * (n0 * Ed0 + n1 * Ed1 + n2 * Ed2 + lEs );
  flux[0] += n0 * Ec;
  flux[1] += n1 * Ec;
  flux[2] += n2 * Ec;

  // Add correction term to H flux
  // c_2 * (n \cdot \jump{E} + \mean{\lambda_E})
  const real Hc = c2 * (n0 * Hd0 + n1 * Hd1 + n2 * Hd2 + lHs);
  flux[3] += n0 * Hc;
  flux[4] += n1 * Hc;
  flux[5] += n2 * Hc;

  // Flux for lambda_E
  // c_1 * ( n \cdot \mean{E} + r \jump{\lambda_E} )
  flux[6] = c1 * (n0 * Es0 + n1 * Es1 + n2 * Es2 + r * lEd );

  // Flux for lambda_H
  // c_2 * ( n \cdot \mean{H} + r \jump{\lambda_H} )
  flux[7] = c2 * (n0 * Hs0 + n1 * Hs1 + n2 * Hs2 + r * lHd );
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DImposedData(const real *x, const real t, real *w) 
{
  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz, lambda_E, lambda_H}

  // FIXME add documentation
  
  const real pi = 4.0 * atan(1.0);
  const real theta = pi / 4.0;
  const real r = 1.0;

  const real u = cos(theta);
  const real v = sin(theta); 
  const real k = 2.0 * pi / v;
  const real c = -cos(k * (u * x[0] + v * x[1] - t));

  // set E
  w[0] = -v * c / r;
  w[1] = u * c / r;
  w[2] = 0.0;
  // set H
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = c / r;
  // set cleaners
  w[6] = 0.0;
  w[7] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DInitData(real *x, real *w) 
{
  real t = 0;
  Maxwell3DImposedData(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DBoundaryFlux_uncentered(real *x, real t, 
				      real *wL, real *vnorm, real *flux)
{
  real wR[8];
  Maxwell3DImposedData(x, t, wR);
  Maxwell3DNumFluxClean_uncentered(wL, wR, vnorm, flux);
}
#pragma end_opencl

// TODO: add 3D clean source.

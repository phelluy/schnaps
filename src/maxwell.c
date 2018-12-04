#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#pragma start_opencl
void
Maxwell2DNumFlux_centered (schnaps_real * wL, schnaps_real * wR,
			   schnaps_real * vnorm, schnaps_real * flux)
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const schnaps_real nx = vnorm[0];
  const schnaps_real ny = vnorm[1];
  const schnaps_real khi = 1.0;

  // FIXME: improve names for these (refer to E, H, etc).
  const schnaps_real s0 = 0.5 * (wR[0] + wL[0]);
  const schnaps_real s1 = 0.5 * (wR[1] + wL[1]);
  const schnaps_real s2 = 0.5 * (wR[2] + wL[2]);
  const schnaps_real s3 = 0.5 * (wR[3] + wL[3]);

  flux[0] = -ny * s2 + khi * nx * s3;
  flux[1] = nx * s2 + khi * ny * s3;
  flux[2] = -ny * s0 + nx * s1;
  flux[3] = khi * (nx * s0 + ny * s1);
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell2DNumFlux_upwind (schnaps_real * wL, schnaps_real * wR,
			 schnaps_real * vnorm, schnaps_real * flux)
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const schnaps_real nx = vnorm[0];
  const schnaps_real ny = vnorm[1];
  const schnaps_real r = sqrt (nx * nx + ny * ny);
  const schnaps_real overr = 1.0 / (r + 1e-16);
  const schnaps_real khi = 1.0;

  const schnaps_real s0 = 0.5 * (wR[0] + wL[0]);
  const schnaps_real s1 = 0.5 * (wR[1] + wL[1]);
  const schnaps_real s2 = 0.5 * (wR[2] + wL[2]);
  const schnaps_real s3 = 0.5 * (wR[3] + wL[3]);

  const schnaps_real d0 = 0.5 * (wR[0] - wL[0]);
  const schnaps_real d1 = 0.5 * (wR[1] - wL[1]);
  const schnaps_real d2 = 0.5 * (wR[2] - wL[2]);
  const schnaps_real d3 = 0.5 * (wR[3] - wL[3]);

  flux[0] =
    -ny * s2 + khi * nx * s3
    - overr * (d0 * (ny * ny + khi * nx * nx) + d1 * nx * ny * (khi - 1));
  flux[1] =
    nx * s2 + khi * ny * s3
    - overr * ((nx * ny * (khi - 1)) * d0 + (nx * nx + khi * ny * ny) * d1);
  flux[2] = -ny * s0 + nx * s1 - r * d2;
  flux[3] = khi * (nx * s0 + ny * s1 - r * d3);
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell2DNumFlux_unoptimised (schnaps_real * wL, schnaps_real * wR,
			      schnaps_real * vnorm, schnaps_real * flux)
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const schnaps_real nx = vnorm[0];
  const schnaps_real ny = vnorm[1];
  const schnaps_real r = sqrt (nx * nx + ny * ny);
  const schnaps_real overr = 1.0 / (r + 1e-16);
  // Centered flux if eps=0, upwind flux if eps=1
  const schnaps_real eps = 1;
  const schnaps_real khi = 1.0;

  flux[0] =
    -ny * (wR[2] + wL[2]) + khi * nx * (wR[3] + wL[3])
    - eps * (ny * ny + khi * nx * nx) * overr * (wR[0] - wL[0])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[1] - wL[1]);

  flux[1] =
    nx * (wR[2] + wL[2]) + khi * ny * (wR[3] + wL[3])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[0] - wL[0])
    - eps * (nx * nx + khi * ny * ny) * overr * (wR[1] - wL[1]);

  flux[2] = -ny * (wR[0] + wL[0])
    + nx * (wR[1] + wL[1]) - eps * r * (wR[2] - wL[2]);

  flux[3] =
    khi * nx * (wR[0] + wL[0])
    + khi * ny * (wR[1] + wL[1]) - eps * khi * r * (wR[3] - wL[3]);

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
void
Maxwell2DImposedData (const schnaps_real * x, const schnaps_real t,
		      schnaps_real * w)
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const schnaps_real pi = 4.0 * atan (1.0);
  const schnaps_real r = 1.0;
  const schnaps_real theta = pi / 4.0;
  const schnaps_real u = cos (theta);
  const schnaps_real v = sin (theta);
  const schnaps_real k = pi / 2;
  const schnaps_real c = -cos (k * (u * x[0] + v * x[1] - t));

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
void
Maxwell2DBoundaryFlux_upwind (schnaps_real * x, schnaps_real t,
			      schnaps_real * wL, schnaps_real * vnorm,
			      schnaps_real * flux)
{
  schnaps_real wR[7];
  Maxwell2DImposedData (x, t, wR);
  Maxwell2DNumFlux_upwind (wL, wR, vnorm, flux);
}

#pragma end_opencl

void
Maxwell2DInitData (schnaps_real * x, schnaps_real * w)
{
  schnaps_real t = 0;
  Maxwell2DImposedData (x, t, w);
}

#pragma start_opencl
void
Maxwell2DSource (const schnaps_real * x, const schnaps_real t,
		 const schnaps_real * w, schnaps_real * source)
{
  // w: (Ex, Ey, Hz, \lambda,  Jx, Jy, rho)

  // FIXME add documentation

  const schnaps_real khi = 1.0;
  source[0] = -w[4];
  source[1] = -w[5];
  source[2] = 0;
  source[3] = khi * w[6];
  source[4] = 0;
  source[5] = 0;
  source[6] = 0;
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell3DNumFlux_upwind (schnaps_real * wL, schnaps_real * wR,
			 schnaps_real * vnorm, schnaps_real * flux)
{
  // Upwind flux (upwind) for Maxwell's equations

  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz}

  // Let {{E}} = ( ER + EL ) / 2, [[E]] = ( ER - EL ) / 2
  // The first three components of the flux are
  // - n x {{H}} + n x n x [[E]] / r
  // and the last three are
  //   n x {{E}} + n x n x [[H]] / r

  const schnaps_real nx = vnorm[0];
  const schnaps_real ny = vnorm[1];
  const schnaps_real nz = vnorm[2];

  const schnaps_real overr =
    1.0 / (sqrt (nx * nx + ny * ny + nz * nz) + 1e-16);
  const schnaps_real nxy = overr * nx * ny;
  const schnaps_real nxz = overr * nx * nz;
  const schnaps_real nyz = overr * ny * nz;
  const schnaps_real nxx = overr * nx * nx;
  const schnaps_real nyy = overr * ny * ny;
  const schnaps_real nzz = overr * nz * nz;

  const schnaps_real Esx = 0.5 * (wR[0] + wL[0]);
  const schnaps_real Esy = 0.5 * (wR[1] + wL[1]);
  const schnaps_real Esz = 0.5 * (wR[2] + wL[2]);

  const schnaps_real Hsx = 0.5 * (wR[3] + wL[3]);
  const schnaps_real Hsy = 0.5 * (wR[4] + wL[4]);
  const schnaps_real Hsz = 0.5 * (wR[5] + wL[5]);

  const schnaps_real Edx = 0.5 * (wR[0] - wL[0]);
  const schnaps_real Edy = 0.5 * (wR[1] - wL[1]);
  const schnaps_real Edz = 0.5 * (wR[2] - wL[2]);

  const schnaps_real Hdx = 0.5 * (wR[3] - wL[3]);
  const schnaps_real Hdy = 0.5 * (wR[4] - wL[4]);
  const schnaps_real Hdz = 0.5 * (wR[5] - wL[5]);

  // E flux
  flux[0] = nz * Hsy - ny * Hsz - (nyy + nzz) * Edx + nxy * Edy + nxz * Edz;
  flux[1] = -nz * Hsx + nx * Hsz + nxy * Edx - (nxx + nzz) * Edy + nyz * Edz;
  flux[2] = -nx * Hsy + ny * Hsx + nxz * Edx + nyz * Edy - (nxx + nyy) * Edz;

  // H flux
  flux[3] = -nz * Esy + ny * Esz - (nyy + nzz) * Hdx + nxy * Hdy + nxz * Hdz;
  flux[4] = nz * Esx - nx * Esz + nxy * Hdx - (nxx + nzz) * Hdy + nyz * Hdz;
  flux[5] = -ny * Esx + nx * Esy + nxz * Hdx + nyz * Hdy - (nxx + nyy) * Hdz;

}

#pragma end_opencl


void Maxwell3DNumFluxClean_Maple(schnaps_real * wL, schnaps_real * wR,
				 schnaps_real * vnorm, schnaps_real * Fnum){


  const schnaps_real nx = vnorm[0];
  const schnaps_real ny = vnorm[1];
  const schnaps_real nz = vnorm[2];

  const schnaps_real cg = sqrt (nx * nx + ny * ny + nz * nz) + 1e-16;
  const schnaps_real c = 1.1;
  
  Fnum[0] = -(c * nx * nx + ny * ny + nz * nz) / cg * wR[0] / 0.2e1 - (c - 0.1e1) * nx * ny / cg * wR[1] / 0.2e1 - (c - 0.1e1) * nx * nz / cg * wR[2] / 0.2e1 + nz * wR[4] / 0.2e1 - ny * wR[5] / 0.2e1 + c * nx * wR[6] / 0.2e1 + (c * nx * nx + ny * ny + nz * nz) / cg * wL[0] / 0.2e1 + (c - 0.1e1) * nx * ny / cg * wL[1] / 0.2e1 + (c - 0.1e1) * nx * nz / cg * wL[2] / 0.2e1 + nz * wL[4] / 0.2e1 - ny * wL[5] / 0.2e1 + c * nx * wL[6] / 0.2e1;
  Fnum[1] = -(c - 0.1e1) * nx * ny / cg * wR[0] / 0.2e1 - (c * ny * ny + nx * nx + nz * nz) / cg * wR[1] / 0.2e1 - (c - 0.1e1) * ny * nz / cg * wR[2] / 0.2e1 - nz * wR[3] / 0.2e1 + nx * wR[5] / 0.2e1 + c * ny * wR[6] / 0.2e1 + (c - 0.1e1) * nx * ny / cg * wL[0] / 0.2e1 + (c * ny * ny + nx * nx + nz * nz) / cg * wL[1] / 0.2e1 + (c - 0.1e1) * ny * nz / cg * wL[2] / 0.2e1 - nz * wL[3] / 0.2e1 + nx * wL[5] / 0.2e1 + c * ny * wL[6] / 0.2e1;
  Fnum[2] = -(c - 0.1e1) * nx * nz / cg * wR[0] / 0.2e1 - (c - 0.1e1) * ny * nz / cg * wR[1] / 0.2e1 - (c * nz * nz + nx * nx + ny * ny) / cg * wR[2] / 0.2e1 + ny * wR[3] / 0.2e1 - nx * wR[4] / 0.2e1 + c * nz * wR[6] / 0.2e1 + (c - 0.1e1) * nx * nz / cg * wL[0] / 0.2e1 + (c - 0.1e1) * ny * nz / cg * wL[1] / 0.2e1 + (c * nz * nz + nx * nx + ny * ny) / cg * wL[2] / 0.2e1 + ny * wL[3] / 0.2e1 - nx * wL[4] / 0.2e1 + c * nz * wL[6] / 0.2e1;
  Fnum[3] = -nz * wR[1] / 0.2e1 + ny * wR[2] / 0.2e1 - (c * nx * nx + ny * ny + nz * nz) / cg * wR[3] / 0.2e1 - (c - 0.1e1) * nx * ny / cg * wR[4] / 0.2e1 - (c - 0.1e1) * nx * nz / cg * wR[5] / 0.2e1 + c * nx * wR[7] / 0.2e1 - nz * wL[1] / 0.2e1 + ny * wL[2] / 0.2e1 + (c * nx * nx + ny * ny + nz * nz) / cg * wL[3] / 0.2e1 + (c - 0.1e1) * nx * ny / cg * wL[4] / 0.2e1 + (c - 0.1e1) * nx * nz / cg * wL[5] / 0.2e1 + c * nx * wL[7] / 0.2e1;
  Fnum[4] = nz * wR[0] / 0.2e1 - nx * wR[2] / 0.2e1 - (c - 0.1e1) * nx * ny / cg * wR[3] / 0.2e1 - (c * ny * ny + nx * nx + nz * nz) / cg * wR[4] / 0.2e1 - (c - 0.1e1) * ny * nz / cg * wR[5] / 0.2e1 + c * ny * wR[7] / 0.2e1 + nz * wL[0] / 0.2e1 - nx * wL[2] / 0.2e1 + (c - 0.1e1) * nx * ny / cg * wL[3] / 0.2e1 + (c * ny * ny + nx * nx + nz * nz) / cg * wL[4] / 0.2e1 + (c - 0.1e1) * ny * nz / cg * wL[5] / 0.2e1 + c * ny * wL[7] / 0.2e1;
  Fnum[5] = -ny * wR[0] / 0.2e1 + nx * wR[1] / 0.2e1 - (c - 0.1e1) * nx * nz / cg * wR[3] / 0.2e1 - (c - 0.1e1) * ny * nz / cg * wR[4] / 0.2e1 - (c * nz * nz + nx * nx + ny * ny) / cg * wR[5] / 0.2e1 + c * nz * wR[7] / 0.2e1 - ny * wL[0] / 0.2e1 + nx * wL[1] / 0.2e1 + (c - 0.1e1) * nx * nz / cg * wL[3] / 0.2e1 + (c - 0.1e1) * ny * nz / cg * wL[4] / 0.2e1 + (c * nz * nz + nx * nx + ny * ny) / cg * wL[5] / 0.2e1 + c * nz * wL[7] / 0.2e1;
  Fnum[6] = c * nx * wR[0] / 0.2e1 + c * ny * wR[1] / 0.2e1 + c * nz * wR[2] / 0.2e1 - cg * c * wR[6] / 0.2e1 + c * nx * wL[0] / 0.2e1 + c * ny * wL[1] / 0.2e1 + c * nz * wL[2] / 0.2e1 + cg * c * wL[6] / 0.2e1;
  Fnum[7] = c * nx * wR[3] / 0.2e1 + c * ny * wR[4] / 0.2e1 + c * nz * wR[5] / 0.2e1 - cg * c * wR[7] / 0.2e1 + c * nx * wL[3] / 0.2e1 + c * ny * wL[4] / 0.2e1 + c * nz * wL[5] / 0.2e1 + cg * c * wL[7] / 0.2e1;

  
}


#pragma start_opencl
void
Maxwell3DNumFluxClean_upwind (schnaps_real * wL, schnaps_real * wR,
			      schnaps_real * vnorm, schnaps_real * flux)
{
  // Upwind flux (upwind) for Maxwell's equations

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
  Maxwell3DNumFlux_upwind (wL, wR, vnorm, flux);

  // FIXME: temp
  /* flux[6] = 0.0; */
  /* flux[7] = 0.0; */
  /* return;  */

  // FIXME: how do we set these?  What are good values?
  const schnaps_real c1 = 0.1;	// E-cleaning parameter
  const schnaps_real c2 = 0.1;	// H-cleaning parameter

  // FIXME: temp
  /* const real c1 = 0.0; */
  /* const real c2 = 0.0; */

  // Consts based on vnorm
  const schnaps_real nx = vnorm[0];
  const schnaps_real ny = vnorm[1];
  const schnaps_real nz = vnorm[2];
  const schnaps_real r = sqrt (nx * nx + ny * ny + nz * nz);

  // Consts based on the mean
  const schnaps_real Esx = 0.5 * (wR[0] + wL[0]);
  const schnaps_real Esy = 0.5 * (wR[1] + wL[1]);
  const schnaps_real Esz = 0.5 * (wR[2] + wL[2]);

  const schnaps_real Hsx = 0.5 * (wR[3] + wL[3]);
  const schnaps_real Hsy = 0.5 * (wR[4] + wL[4]);
  const schnaps_real Hsz = 0.5 * (wR[5] + wL[5]);

  const schnaps_real lEs = 0.5 * (wR[6] + wL[6]);
  const schnaps_real lHs = 0.5 * (wR[7] + wL[7]);

  // Consts based on the jump
  const schnaps_real Edx = 0.5 * (wR[0] - wL[0]);
  const schnaps_real Edy = 0.5 * (wR[1] - wL[1]);
  const schnaps_real Edz = 0.5 * (wR[2] - wL[2]);

  const schnaps_real Hdx = 0.5 * (wR[3] - wL[3]);
  const schnaps_real Hdy = 0.5 * (wR[4] - wL[4]);
  const schnaps_real Hdz = 0.5 * (wR[5] - wL[5]);

  const schnaps_real lEd = 0.5 * (wR[6] - wL[6]);
  const schnaps_real lHd = 0.5 * (wR[7] - wL[7]);

  // Add correction term to E flux
  // c_1 * ( n \cdot \jump{E} + \mean{\lambda_E} )
  const schnaps_real Ec = c1 * (nx * Edx + ny * Edy + nz * Edz + lEs);
  flux[0] += nx * Ec;
  flux[1] += ny * Ec;
  flux[2] += nz * Ec;

  // Add correction term to H flux
  // c_2 * ( n \cdot \jump{E} + \mean{\lambda_E} )
  const schnaps_real Hc = c2 * (nx * Hdx + ny * Hdy + nz * Hdz + lHs);
  flux[3] += nx * Hc;
  flux[4] += ny * Hc;
  flux[5] += nz * Hc;

  // Flux for lambda_E
  // c_1 * ( n \cdot \mean{E} + r \jump{\lambda_E} )
  flux[6] = c1 * (nx * Esx + ny * Esy + nz * Esz + r * lEd);

  // Flux for lambda_H
  // c_2 * ( n \cdot \mean{H} + r \jump{\lambda_H} )
  flux[7] = c2 * (nx * Hsx + ny * Hsy + nz * Hsz + r * lHd);
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell3DImposedData (const schnaps_real * x, const schnaps_real t,
		      schnaps_real * w)
{
  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz, lambda_E, lambda_H}

#if 0

#ifndef M_PI
#define M_PI (3.141592653589793)
#endif

#define MAXWELL_THETA (M_PI / 4)
#define MAXWELL_PHI   (M_PI / 4)

  schnaps_real Vn[] = {
    sin (MAXWELL_THETA) * cos (MAXWELL_PHI),
    sin (MAXWELL_THETA) * sin (MAXWELL_PHI),
    cos (MAXWELL_THETA)
  };
  schnaps_real Ws[] = {
    cos (MAXWELL_THETA) * cos (MAXWELL_PHI),
    cos (MAXWELL_THETA) * sin (MAXWELL_PHI),
    -sin (MAXWELL_THETA),
    -sin (MAXWELL_PHI),
    cos (MAXWELL_PHI),
    0
  };
  schnaps_real test_cos_frequency = 1;

  schnaps_real xdotVn = Vn[0] * x[0] + Vn[1] * x[1] + Vn[2] * x[2];
  schnaps_real magnitude =
    cos (M_PI * 2.0 * test_cos_frequency * (xdotVn - t));

  for (int ii = 0; ii < 6; ii++)
    {
      w[ii] = Ws[ii] * magnitude;
    }
#else

  const schnaps_real pi = 4.0 * atan (1.0);
  const schnaps_real theta = pi / 4.0;
  const schnaps_real r = 1.0;

  const schnaps_real u = cos (theta);
  const schnaps_real v = sin (theta);
  const schnaps_real k = pi/2;
  const schnaps_real c = -cos (k * (u * x[0] + v * x[1] - t));

  // set E
  w[0] = -v * c / r;
  w[1] = u * c / r;
  w[2] = 0.0;
  // set H
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = c / r;

#endif
  // set cleaners
  w[6] = 0.0;
  w[7] = 0.0;
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell3DInitData (schnaps_real * x, schnaps_real * w)
{
  schnaps_real t = 0;
  Maxwell3DImposedData (x, t, w);
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell3DBoundaryFlux_upwind (schnaps_real * x, schnaps_real t,
			      schnaps_real * wL, schnaps_real * vnorm,
			      schnaps_real * flux)
{
  schnaps_real wR[8];
  Maxwell3DImposedData (x, t, wR);
  Maxwell3DNumFlux_upwind (wL, wR, vnorm, flux);
}

#pragma end_opencl

#pragma start_opencl
void
Maxwell3DBoundaryFluxClean_upwind (schnaps_real * x, schnaps_real t,
				   schnaps_real * wL, schnaps_real * vnorm,
				   schnaps_real * flux)
{
  schnaps_real wR[8];
  Maxwell3DImposedData (x, t, wR);
  Maxwell3DNumFluxClean_upwind (wL, wR, vnorm, flux);
}

#pragma end_opencl

// TODO: add 3D clean source.

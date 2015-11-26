#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "global.h"

#include "interpolation.h"

#pragma start_opencl
//! Gauss LObatto Points (GLOP) up to order 4
__constant real gauss_lob_point[] = {
  0.5,
  0.0,
  1.0,
  0.0,
  0.5,
  1.0,
  0.0,
  0.276393202250021030359082633127,
  0.723606797749978969640917366873,
  1.0,
  0.0,
  0.172673164646011428100853771877,
  0.5,
  0.827326835353988571899146228123,
  1.0
};

//! GLOP weights up to order 4
__constant real gauss_lob_weight[] = {
  1.0,
  0.5,
  0.5,
  0.166666666666666666666666666667,
  0.666666666666666666666666666668,
  0.166666666666666666666666666667,
  0.0833333333333333333333333333333,
  0.416666666666666666666666666666,
  0.416666666666666666666666666666,
  0.0833333333333333333333333333333,
  0.05,
  0.272222222222222222222222222223,
  0.355555555555555555555555555556,
  0.272222222222222222222222222219,
  0.05
};

//! indirection for finding the GLOP
//! data for a given degree in the previous arrays
__constant int gauss_lob_offset[] = {0, 1, 3, 6, 10};

__constant real gauss_lob_dpsi[] = {
  0.0,
  -1.0,
  -1.0,
  1.0,
  1.0,
  -3.0,
  -1.0,
  1.0,
  4.0,
  0.0,
  -4.0,
  -1.0,
  1.0,
  3.0,
  -6,
  -1.61803398874989484820458683436,
  .618033988749894848204586834362,
  -1,
  8.09016994374947424102293417177,
  0,
  -2.23606797749978969640917366872,
  3.09016994374947424102293417184,
  -3.09016994374947424102293417182,
  2.23606797749978969640917366872,
  0,
  -8.09016994374947424102293417177,
  1,
  -.618033988749894848204586834362,
  1.61803398874989484820458683436,
  6,
  -10,
  -2.48198050606196571569743868436,
  .75,
  -.518019493938034284302561315632,
  1,
  13.5130049774484800076860550594,
  0,
  -2.67316915539090667050969419631,
  1.52752523165194666886268239794,
  -2.82032835588485332564727827404,
  -5.33333333333333333333333333336,
  3.49148624377587810025755976667,
  0,
  -3.49148624377587810025755976662,
  5.33333333333333333333333333336,
  2.82032835588485332564727827399,
  -1.52752523165194666886268239791,
  2.67316915539090667050969419635,
  0,
  -13.5130049774484800076860550594,
  -1,
  .518019493938034284302561315631,
  -.75,
  2.48198050606196571569743868437,
  10
};

//! indirection for finding the GLOP
//! data for a given degree in the previous arrays
__constant int gauss_lob_dpsi_offset[] = {0, 1, 5, 14, 30};

#pragma end_opencl

//! \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
real wglop(int deg, int i)
{
  return gauss_lob_weight[gauss_lob_offset[deg] + i];
}

real glop(int deg, int i)
{
  return gauss_lob_point[gauss_lob_offset[deg] + i];
}

void lagrange_polynomial(real* p, const real* subdiv,
			 int deg, int ii, real x)
{
  *p = 1;
  const int npg = deg + 1;
  for(int j = 0; j < npg; j++) {
    if (j != ii) {
      *p *= (x - subdiv[j]) / (subdiv[ii] - subdiv[j]);
    }
  }
}

void dlagrange_polynomial(real* dp, const real* subdiv,
			  int deg, int i, real x)
{
  *dp = 0;
  const int npg = deg + 1;
  for(int k = 0; k < npg; k++) {
    if (k != i) {
      real xk = subdiv[k];
      real dploc = 1.0 / (subdiv[i] -  xk);
      for(int j = 0; j < (deg) + 1; j++) {
	if (j != i && j != k) {
	  real xj = subdiv[j];
	  dploc *= (x - xj) / (subdiv[i] - xj);
	}
      }
      *dp += dploc;
    }
  }
}

// Number of Gauss Lobatto Points (GLOPS) in a macrocell
int NPG(int *raf, int *deg)
{
  return raf[0] * raf[1] * raf[2] * (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1);
}

// Number of interpolation points for each face of a subcell
int NPGF(int *raf, int *deg, int ifa)
{
  // For each face, give the dimension index i
  int permut[6][4] = { {0, 2, 1, 0},
		       {1, 2, 0, 1},
		       {2, 0, 1, 1},
		       {2, 1, 0, 0},
		       {0, 1, 2, 1},
		       {1, 0, 2, 0} };
  int i0 = permut[ifa][0];
  int i1 = permut[ifa][1];
  return raf[i0] * raf[i1] * (deg[i0] + 1) * (deg[i1] + 1);
}

#pragma start_opencl
int xyz_to_ipg(const int *raf, const int *deg, const int *ic, const int *ix) 
{
  const int nc = ic[0] + raf[0] * (ic[1] + raf[1] * ic[2]);
  const int offset = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1) * nc;

  return ix[0] + (deg[0] + 1) * (ix[1] + (deg[1] + 1) * ix[2]) + offset;
}
#pragma end_opencl

#pragma start_opencl
void ipg_to_xyz(const int *raf, const int *deg, int *ic, int *ix, 
		const int *pipg)
{
  int ipg = *pipg;

  ix[0] = ipg % (deg[0] + 1);
  ipg /= (deg[0] + 1);

  ix[1] = ipg % (deg[1] + 1);
  ipg /= (deg[1] + 1);

  ix[2] = ipg % (deg[2] + 1);
  ipg /= (deg[2] + 1);

  ic[0] = ipg % raf[0];
  ipg /= raf[0];

  ic[1] = ipg % raf[1];
  ipg /= raf[1];

  ic[2] = ipg;
}
#pragma end_opencl

// From a reference point find the nearest gauss point
// Warning: works only  degree 1, 2, or 3
int ref_ipg(int *raf, int *deg, real *xref)
{
  real h[3] = {1.0 / raf[0], 1.0 / raf[1], 1.0 / raf[2]};
  
  // get the subcell id
  int ic[3] = {floor(xref[0] * raf[0]),
	       floor(xref[1] * raf[1]),
	       floor(xref[2] * raf[2]) };

  assert(ic[0] >=0 && ic[0] < raf[0]);
  assert(ic[1] >=0 && ic[1] < raf[1]);
  assert(ic[2] >=0 && ic[2] < raf[2]);

  // round to the nearest integer
  int ix[3];
  for(int i = 0; i < 3; ++i) {
    ix[i] = floor((xref[i] - ic[i] * h[i]) / h[i] * deg[i] + 0.5);
  }
  
  int ipg = xyz_to_ipg(raf, deg, ic, ix);
  return ipg;
}

// Return the reference coordinates xpg[3] and weight wpg of the GLOP ipg
void ref_pg_vol(int* raf, int* deg,
		int ipg, real *xpg, real *wpg, real *xpg_in)
{
  int ic[3];
  int ix[3];
  ipg_to_xyz(raf, deg, ic, ix, &ipg);

  real h[3] = {1.0 / (real) raf[0],
	       1.0 / (real) raf[1],
	       1.0 / (real) raf[2] };

  int offset[3] = {gauss_lob_offset[deg[0]] + ix[0],
		   gauss_lob_offset[deg[1]] + ix[1],
		   gauss_lob_offset[deg[2]] + ix[2] };

  if (xpg != NULL) {
    xpg[0] = h[0] * (ic[0] + gauss_lob_point[offset[0]]);
    xpg[1] = h[1] * (ic[1] + gauss_lob_point[offset[1]]);
    xpg[2] = h[2] * (ic[2] + gauss_lob_point[offset[2]]);
  }
  
  if (wpg != NULL) {
    *wpg = h[0] * h[1] * h[2] *
      gauss_lob_weight[offset[0]]*
      gauss_lob_weight[offset[1]]*
      gauss_lob_weight[offset[2]];
  }
    
  if (xpg_in != NULL) {
    real small = 1e-5; //1e-3;

    for(int i = 0; i < 3; ++i) {
      xpg_in[i] = xpg[i];
      if (ix[i] == 0)
	xpg_in[i] += h[i] * small;
      if (ix[i] == deg[i])
	xpg_in[i] -= h[i] * small;
    }
  }
}

// Return the reference coordinates xpg[3] and weight wpg of the GLOP
// ipg on the face ifa.
int ref_pg_face(int *raf, int *deg, int ifa, int ipgf, 
		 real *xpg, real *wpg, real *xpgin)
{
  // For each face, give the dimension index i
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };

  // number of subcells in each permuted direction
  int praf[3] = {raf[axis_permut[ifa][0]],
		 raf[axis_permut[ifa][1]],
		 raf[axis_permut[ifa][2]] };
  
  // approximation degree in each permuted direction
  int pdeg[3] = {deg[axis_permut[ifa][0]],
		 deg[axis_permut[ifa][1]],
		 deg[axis_permut[ifa][2]] };

  // Compute permuted indices
  int pix[3];
  pix[0] = ipgf % (pdeg[0] + 1);
  ipgf /= (pdeg[0] + 1);
  pix[1] = ipgf % (pdeg[1] + 1);
  ipgf /= (pdeg[1] + 1);
  // pix[2] is 0 or d depending on the face
  pix[2] = axis_permut[ifa][3] * pdeg[2];

  // Compute permuted indices of the subface
  int pic[3];
  pic[0] = ipgf % praf[0];
  ipgf /= praf[0];
  pic[1] = ipgf;
  // pic[2] is 0 or praf-1 depending on the face
  pic[2] = axis_permut[ifa][3] * (praf[2] - 1);

  real h[3] = {1.0 / (real) praf[0],
	       1.0 / (real) praf[1],
	       1.0 / (real) praf[2] };
  
  // Compute non permuted indices for points and subfaces
  int ic[3];
  ic[axis_permut[ifa][0]] = pic[0];
  ic[axis_permut[ifa][1]] = pic[1];
  ic[axis_permut[ifa][2]] = pic[2];

  int ix[3];
  ix[axis_permut[ifa][0]] = pix[0];
  ix[axis_permut[ifa][1]] = pix[1];
  ix[axis_permut[ifa][2]] = pix[2];

  // Compute the global index of the Gauss-Lobatto point in the volume
  int ipgv = ix[0] + (deg[0] + 1) *
    (ix[1] + (deg[1] + 1) *
     (ix[2] + (deg[2] + 1) *
      (ic[0] + raf[0] *
       (ic[1] + raf[1] *
	ic[2]))));

  // Compute the reference coordinates of the Gauss-Lobatto point in
  // the volume
  int offset[2] = {gauss_lob_offset[pdeg[0]] + pix[0],
		   gauss_lob_offset[pdeg[1]] + pix[1] };

  if(xpg != NULL) {
    xpg[axis_permut[ifa][0]] = h[0] * (pic[0] + gauss_lob_point[offset[0]]);
    xpg[axis_permut[ifa][1]] = h[1] * (pic[1] + gauss_lob_point[offset[1]]);
    xpg[axis_permut[ifa][2]] = axis_permut[ifa][3];
  }

  if(wpg != NULL) {
    *wpg
      = h[0] * gauss_lob_weight[offset[0]]
      * h[1] * gauss_lob_weight[offset[1]];
  }
    
  // If xpgin exists, compute a point slightly INSIDE the opposite
  // subcell along the face.
  if(xpgin != NULL) {
    real small = 1e-3;  //0.001
    real vsmall = 1e-5; //0.000001;

    xpgin[axis_permut[ifa][0]] = h[0] * (pic[0] + gauss_lob_point[offset[0]]);
    xpgin[axis_permut[ifa][1]] = h[1] * (pic[1] + gauss_lob_point[offset[1]]);

    if(axis_permut[ifa][3] == 0)
      xpgin[axis_permut[ifa][2]] = -vsmall;
    if(axis_permut[ifa][3] == 1)
      xpgin[axis_permut[ifa][2]] = 1 + vsmall;

    if(pix[0] == 0)
      xpgin[axis_permut[ifa][0]]
	= h[0] * (pic[0] + gauss_lob_point[offset[0]] + small);
    if(pix[0] == pdeg[0])
      xpgin[axis_permut[ifa][0]]
	= h[0] * (pic[0] + gauss_lob_point[offset[0]] - small);

    if(pix[1] == 0)
      xpgin[axis_permut[ifa][1]]
	= h[1] * (pic[1] + gauss_lob_point[offset[1]] + small);
    if(pix[1] == pdeg[1])
      xpgin[axis_permut[ifa][1]]
	= h[1] * (pic[1] + gauss_lob_point[offset[1]] - small);
  }
  
  return ipgv;
}

// return the 1d derivative of lagrange polynomial ib at glop ipg
real dlag(int deg, int ib, int ipg) 
{
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg] + ib * (deg + 1) + ipg];
}

// Return the value psi and the gradient dpsi[3] of the basis function
// ib at point xref[3].
// Warning: the value of the gradient is not reliable if xref is on
// the boundary of a subcell (because the gradient is discontinuous)
void psi_ref(int *param, int ib, real *xref, real *psi, real *dpsi)
{
  // number of subcells in each direction
  int raf[3] = {param[3], param[4], param[5]};

  // approximation degree in each direction
  int deg[3] = {param[0], param[1], param[2]};

  int ic[3];
  int ix[3];
  ipg_to_xyz(raf, deg, ic, ix, &ib); 

  real h[3] = {1.0 / (real) raf[0], 1.0 / (real) raf[1], 1.0 / (real) raf[2]};

  // Starting Gauss-Lobatto point in each direction
  int offset[3] = {gauss_lob_offset[deg[0]],
		   gauss_lob_offset[deg[1]],
		   gauss_lob_offset[deg[2]] };

  real psib[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i) {
    lagrange_polynomial(&psib[i], gauss_lob_point + offset[i],
			deg[i], ix[i], xref[i] / h[i] - ic[i]);
  }
  
  int is[3] = {xref[0] * raf[0], xref[1] * raf[1], xref[2] * raf[2]};
  for(int ii = 0; ii < 3; ii++) {
    assert(is[ii] < raf[ii] && is[ii]>= 0);
  }
  
  int is_in_subcell= (ic[0] == is[0]) && (ic[1] == is[1]) && (ic[2] == is[2]);

  *psi = psib[0] * psib[1] * psib[2] * is_in_subcell;

  if (dpsi != NULL) {
    real dpsib[3];
    for(int i = 0; i < 3; ++i) {
      dlagrange_polynomial(&dpsib[i], gauss_lob_point + offset[i],
                         deg[i], ix[i], xref[i]);
    }
    
    dpsi[0] = dpsib[0] *  psib[1] *  psib[2] * is_in_subcell;
    dpsi[1] =  psib[0] * dpsib[1] *  psib[2] * is_in_subcell;
    dpsi[2] =  psib[0] *  psib[1] * dpsib[2] * is_in_subcell;
  }
}

// Return the value psi and the gradient dpsi[3] of the basis function
// ib at point xref[3] given the subcell indices is[3].
// The computation is reliable.
void psi_ref_subcell(int *param, int *is, int ib,
		     real *xref, real *psi, real *dpsi)
{
  // Number of subcells in each direction
  int raf[3] = {param[3], param[4], param[5]};

  // Approximation degree in each direction
  int deg[3] = {param[0], param[1], param[2]};

  // Starting Gauss-Lobatto point in each direction
  int offset[3] = {gauss_lob_offset[deg[0]],
		   gauss_lob_offset[deg[1]],
		   gauss_lob_offset[deg[2]] };

  int ic[3];
  int ix[3];
  ipg_to_xyz(raf, deg, ic, ix, &ib);

  real h[3] = {1.0 / (real) raf[0], 1.0 / (real) raf[1], 1.0 / (real) raf[2]};

  int is_in_subcell = (ic[0] == is[0]) && (ic[1] == is[1]) && (ic[2] == is[2]);

  real psib[3] = {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; ++i) {
    lagrange_polynomial(&psib[i], gauss_lob_point + offset[i],
			deg[i], ix[i], xref[i] / h[i] - ic[i]);
  }
  
  // might be useful for the future subcell case
  /* psibx *= (xref[0] <= (ncbx + 1) * hx)&&(xref[0] > ncbx * hx); */
  /* psiby *= (xref[1] <= (ncby + 1) * hy)&&(xref[1] > ncby * hy); */
  /* psibz *= (xref[2] <= (ncbz + 1) * hz)&&(xref[2] > ncbz * hz); */

  *psi = psib[0] * psib[1] * psib[2] * is_in_subcell ;

  if (dpsi != NULL) {
    real dpsib[3];
    for(int i = 0; i < 3; ++i) {
      dlagrange_polynomial(&dpsib[i], gauss_lob_point + offset[i],
                         deg[i], ix[i], xref[i]);
    }
    
    dpsi[0] = dpsib[0] *  psib[1] *  psib[2] * is_in_subcell;
    dpsi[1] =  psib[0] * dpsib[1] *  psib[2] * is_in_subcell;
    dpsi[2] =  psib[0] *  psib[1] * dpsib[2] * is_in_subcell;
  }
}

// Return the gradient dpsi[0..2] of the basis function ib at GLOP
// ipg.
void grad_psi_pg(int *param, int ib, int ipg, real *dpsi)
{
  // number of subcells in each direction
  int nraf[3] = {param[3], param[4], param[5]};
  
  // approximation degree in each direction
  int deg[3] = {param[0], param[1], param[2]};

  // glop 3d indices
  int ic[3];
  int ix[3];
  ipg_to_xyz(nraf, deg, ic, ix, &ipg);

  // basis function 3d indices
  int ibc[3];
  int ibx[3];
  ipg_to_xyz(nraf, deg, ibc, ibx, &ib);

  real h[3] = {1.0 / (real) nraf[0],
	       1.0 / (real) nraf[1],
	       1.0 / (real) nraf[2] };

  // Number of Gauss-Lobatto points in each direction
  int offset[3] = {gauss_lob_dpsi_offset[deg[0]],
		   gauss_lob_dpsi_offset[deg[1]],
		   gauss_lob_dpsi_offset[deg[2]] };

  // Computation of the value of the interpolation polynomial gradient
  real psib[3] = {(ix[0] == ibx[0]) * (ic[0] == ibc[0]),
		  (ix[1] == ibx[1]) * (ic[1] == ibc[1]),
		  (ix[2] == ibx[2]) * (ic[2] == ibc[2]) };
  real dpsib[3];
  for(int i = 0; i < 3; ++i) {
    dpsib[i] = (ic[i] == ibc[i])
      * gauss_lob_dpsi[offset[i] + ibx[i] * (deg[i] + 1) + ix[i]] / h[i];
  }
  
  dpsi[0] = dpsib[0] *  psib[1] *  psib[2];
  dpsi[1] =  psib[0] * dpsib[1] *  psib[2];
  dpsi[2] =  psib[0] *  psib[1] * dpsib[2];
}

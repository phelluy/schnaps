#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "global.h"

#include "interpolation.h"

#pragma start_opencl
//! Gauss LObatto Points (GLOP) up to order 4
__constant schnaps_real gauss_lob_point[] = {
  0.5,
  0,
  1,
  0,
  0.5,
  1,
  0,
  0.276393202250021030359082633127,
  0.723606797749978969640917366873,
  1,
  0,
  0.172673164646011428100853771877,
  0.5,
  0.827326835353988571899146228123,
  1
};

//! GLOP weights up to order 4
__constant schnaps_real gauss_lob_weight[] = {
  1,
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

__constant schnaps_real gauss_lob_dpsi[] = {
  0.0,
  -1.,
  -1.,
  1.,
  1.,
  -3.,
  -1.,
  1.,
  4.,
  0.,
  -4.,
  -1.,
  1.,
  3.,
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

#pragma start_opencl
//! \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
schnaps_real wglop(int deg, int i) {
  return gauss_lob_weight[gauss_lob_offset[deg] + i];
}

schnaps_real glop(int deg, int i){
  return gauss_lob_point[gauss_lob_offset[deg] + i];
}
#pragma end_opencl

void lagrange_polynomial(schnaps_real* p, const schnaps_real* subdiv,
			 int deg, int ii, schnaps_real x) {
  *p = 1;
  const int npg = deg + 1;
  for(int j = 0; j < npg; j++) {
    if (j != ii) {
      *p *= (x - subdiv[j]) / (subdiv[ii] - subdiv[j]);
    }
  }
}

void dlagrange_polynomial(schnaps_real* dp, const schnaps_real* subdiv,
			  int deg, int i, schnaps_real x) {
  schnaps_real xj;
  *dp = 0;
  const int npg = deg + 1;
  for(int k = 0; k < npg; k++) {
    if (k != (i)) {
      schnaps_real xk = (subdiv)[k];
      schnaps_real dploc = ((schnaps_real)1) / ((subdiv)[i] -  xk);
      for(int j=0;j<(deg)+1;j++) {
	if (j != (i) && j != k) {
	  xj = (subdiv)[j];
	  dploc *= ((x) - xj) / ((subdiv)[i] - xj);
	}
      }
      *dp += dploc;
    }
  }
}

// Number of Gauss Lobatto Points (GLOPS) in a macrocell
#pragma start_opencl
int NPG(const int deg[], const int raf[]) {
  return
    (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1)
    * raf[0] * raf[1] * raf[2];
}
#pragma end_opencl

// Number of Gauss Lobatto Points (GLOPS) in a continous macrocell
int NPG_CG(int deg[], int raf[]) {
  return  (deg[0] * raf[0] + 1) * (deg[1] * raf[1] + 1)
    * (deg[2] * raf[2] + 1);
}

// Number of interpolation points for each face of a subcell
#pragma start_opencl
int NPGF(const int deg[], const int raf[], int ifa) {
  // For each face, give the dimension index i
  int permut[6][4] = {
    {0, 2, 1, 0},
    {1, 2, 0, 1},
    {2, 0, 1, 1},
    {2, 1, 0, 0},
    {0, 1, 2, 1},
    {1, 0, 2, 0}
  };
  int i0 = permut[ifa][0];
  int i1 = permut[ifa][1];
  /* printf("locfa=%d deg= %d %d %d raf=%d %d %d\n",ifa,deg[0],deg[1],deg[2], */
  /* 	raf[0],raf[1],raf[2]); */
  return (deg[i0] + 1) * (deg[i1] + 1) * raf[i0] * raf[i1];
}
#pragma end_opencl

// Number of interpolation points for each continuous face of a subcell
int NPGF_CG(int deg[], int raf[], int ifa) {
  // For each face, give the dimension index i
  int permut[6][4] = {
    {0, 2, 1, 0},
    {1, 2, 0, 1},
    {2, 0, 1, 1},
    {2, 1, 0, 0},
    {0, 1, 2, 1},
    {1, 0, 2, 0}
  };
  int i0 = permut[ifa][0];
  int i1 = permut[ifa][1];
  return (deg[i0] * raf[i0] + 1) * (deg[i1] * raf[i1] + 1);
}



#pragma start_opencl
void xyz_to_ipg(const int *raf, const int *deg, const int *ic, const int *ix,
		int *ipg)
{
  const int nc = ic[0] + raf[0] * (ic[1] + raf[1] * ic[2]);
  const int offset = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1)*nc;

  *ipg= ix[0] + (deg[0] + 1) * (ix[1] + (deg[1] + 1) * ix[2]) + offset;
}
#pragma end_opencl

#pragma start_opencl
void ipg_to_xyz(const int *raf, const int *deg, int *ic, int *ix,
		const int *pipg) {
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
// Warning: works only for  degrees 1,2 or 3
int ref_ipg(int *deg, int *nraf, schnaps_real *xref) {

  schnaps_real hh[3] = {1./nraf[0],1./nraf[1],1./nraf[2]};

  int ic[3],ix[3];

  // get the subcell id
  ic[0] = floor(xref[0] * nraf[0]);
  ic[1] = floor(xref[1] * nraf[1]);
  ic[2] = floor(xref[2] * nraf[2]);

  if (!(ic[0] >=0 && ic[0]<nraf[0]) ||
      !(ic[1] >=0 && ic[1]<nraf[1]) ||
      !(ic[2] >=0 && ic[2]<nraf[2]) ){
    printf("x=%.10e ic[0]=%d nrafx=%d\n",xref[0], ic[0],nraf[0]);
    printf("y=%.10e ic[1]=%d nrafy=%d\n",xref[1], ic[1],nraf[1]);
    printf("z=%.10e ic[2]=%d nrafz=%d\n",xref[2], ic[2],nraf[2]);
    printf("correction...\n");
    for(int ii = 0; ii < 2; ii++){
      if (ic[ii] < 0) ic[ii] = 0;
      if (ic[ii] >= nraf[ii]) ic[ii] = nraf[ii] - 1;
    }
  }
  
  assert(ic[0] >=0 && ic[0]<nraf[0]);
  assert(ic[1] >=0 && ic[1]<nraf[1]);
  assert(ic[2] >=0 && ic[2]<nraf[2]);

  // subcell index in the macrocell
  //int nc = ic[0] + nraf[0] * (ic[1] + nraf[1] * ic[2]);
  //int offset = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1)*nc;

  // round to the nearest integer
  ix[0] = floor((xref[0] - ic[0] * hh[0]) / hh[0] * deg[0] + 0.5);
  ix[1] = floor((xref[1] - ic[1] * hh[1]) / hh[1] * deg[1] + 0.5);
  ix[2] = floor((xref[2] - ic[2] * hh[2]) / hh[2] * deg[2] + 0.5);
  //int ix[2]=floor(xref[2]*deg[2]+0.5);

  //printf("xref %f %f %f ix[0]=%d ix[1]=%d ix[2]=%d\n",
  //	 xref[0],xref[1],xref[2],ix[0],ix[1],ix[2]);

  int ipg;
  xyz_to_ipg(nraf,deg,ic,ix,&ipg);

  //return ix[0] + (deg[0] + 1) * (ix[1] + (deg[1] + 1) * ix[2]) + offset;
  return ipg;
} // ref_ipg


int ref_ipg_CG(int *deg, int *nraf, schnaps_real *xref) {

  schnaps_real hh[3] = {1./nraf[0],1./nraf[1],1./nraf[2]};

  int ic[3],ix[3];

  // get the subcell id
  ic[0] = floor(xref[0] * nraf[0]);
  ic[1] = floor(xref[1] * nraf[1]);
  ic[2] = floor(xref[2] * nraf[2]);

  //printf("x=%f ic[0]=%d nrafx=%d\n",xref[0], ic[0],nraf[0]);
  //printf("y=%f ic[1]=%d nrafy=%d\n",xref[1], ic[1],nraf[1]);
  //printf("z=%f ic[2]=%d nrafz=%d\n",xref[2], ic[2],nraf[2]);
  assert(ic[0] >=0 && ic[0]<nraf[0]);
  assert(ic[1] >=0 && ic[1]<nraf[1]);
  assert(ic[2] >=0 && ic[2]<nraf[2]);

  // subcell index in the macrocell
  //int nc = ic[0] + nraf[0] * (ic[1] + nraf[1] * ic[2]);
  //int offset = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1)*nc;

  // round to the nearest integer
  ix[0] = floor((xref[0] - ic[0] * hh[0]) / hh[0] * deg[0] + 0.5);
  ix[1] = floor((xref[1] - ic[1] * hh[1]) / hh[1] * deg[1] + 0.5);
  ix[2] = floor((xref[2] - ic[2] * hh[2]) / hh[2] * deg[2] + 0.5);
  //int ix[2]=floor(xref[2]*deg[2]+0.5);

  //printf("xref %f %f %f ix[0]=%d ix[1]=%d ix[2]=%d\n",
  //	 xref[0],xref[1],xref[2],ix[0],ix[1],ix[2]);

  int nx[3] = {
    deg[0] * nraf[0] + 1,
    deg[1] * nraf[1] + 1,
    deg[2] * nraf[2] + 1};

  int irefx[3]= {
    deg[0] * ic[0] + ix[0],
    deg[1] * ic[1] + ix[1],
    deg[2] * ic[2] + ix[2]};

  int ipg = irefx[0] + nx[0] * (irefx[1] + nx[1] * irefx[2]);
  //xyz_to_ipg(nraf,deg,ic,ix,&ipg);

  //return ix[0] + (deg[0] + 1) * (ix[1] + (deg[1] + 1) * ix[2]) + offset;
  return ipg;
} // ref_ipg_CG

void ref_pg_vol(int *deg, int *nraf, int ipg, schnaps_real *xpg, schnaps_real *wpg, schnaps_real *xpg_in) {

  int offset[3];


  int ix[3], ic[3];

  ipg_to_xyz(nraf,deg,ic,ix,&ipg);

  schnaps_real hx = 1 / (schnaps_real) nraf[0];
  schnaps_real hy =1 / (schnaps_real) nraf[1];
  schnaps_real hz = 1 / (schnaps_real) nraf[2];

  //printf("h=%f %f %f\n",hx,hy,hz);

  offset[0] = gauss_lob_offset[deg[0]] + ix[0];
  offset[1] = gauss_lob_offset[deg[1]] + ix[1];
  offset[2] = gauss_lob_offset[deg[2]] + ix[2];

  if (xpg != NULL){
    xpg[0] = hx * (ic[0] + gauss_lob_point[offset[0]]);
    xpg[1] = hy * (ic[1] + gauss_lob_point[offset[1]]);
    xpg[2] = hz * (ic[2] + gauss_lob_point[offset[2]]);
  }

  if (wpg != NULL) *wpg = hx * hy * hz *
		     gauss_lob_weight[offset[0]]*
		     gauss_lob_weight[offset[1]]*
		     gauss_lob_weight[offset[2]];

  if (xpg_in !=0) {
    schnaps_real small = _SMALL;
    xpg_in[0] = xpg[0];
    xpg_in[1] = xpg[1];
    xpg_in[2] = xpg[2];

    if (ix[0] == 0) xpg_in[0] += hx * _SMALL;
    if (ix[0] == deg[0]) xpg_in[0] -= hx * _SMALL;
    if (ix[1] == 0) xpg_in[1] += hy * _SMALL;
    if (ix[1] == deg[1]) xpg_in[1] -= hy * _SMALL;
    if (ix[2] == 0) xpg_in[2] += hz * _SMALL;
    if (ix[2] == deg[2]) xpg_in[2] -= hz * _SMALL;

    /* printf("xpg %f %f %f\n",xpg[0],xpg[1],xpg[2]); */
    /*  printf("xpg_in %f %f %f %d %d %d\n",xpg_in[0],xpg_in[1],xpg_in[2], */
    /* 	   ix,iy,iz); */
  }
}

void ref_pg_vol_CG(int *deg, int *nraf, int ipg, schnaps_real *xpg, schnaps_real *wpg, schnaps_real *xpg_in) {

  int offset[3], irefx[3], ic[3], ix[3];
  schnaps_real hx[3];

  int nx[3] = {
    deg[0] * nraf[0] + 1,
    deg[1] * nraf[1] + 1,
    deg[2] * nraf[2] + 1};


  hx[0] = 1 / (schnaps_real) nraf[0];
  hx[1] = 1 / (schnaps_real) nraf[1];
  hx[2] = 1 / (schnaps_real) nraf[2];

  //printf("h=%f %f %f\n",hx,hy,hz);

  irefx[0] = ipg % nx[0];
  ipg /= nx[0];
  irefx[1] = ipg % nx[1];

  ic[0] = irefx[0] / deg[0];
  ic[1] = irefx[1] / deg[1];

  ix[0] = irefx[0] % deg[0];
  ix[1] = irefx[1] % deg[1];

  if (deg[2] !=0) {
    ipg /= nx[1];
    irefx[2] = ipg % nx[2];
    ic[2] = irefx[2] / deg[2];
    ix[2] = irefx[2] % deg[2];
  } else{
    assert(nraf[2] == 1);
    ic[2] = 0;
    ix[2] = 0;
  }

  offset[0] = gauss_lob_offset[deg[0]] + ix[0];
  offset[1] = gauss_lob_offset[deg[1]] + ix[1];
  offset[2] = gauss_lob_offset[deg[2]] + ix[2];

  if (xpg != NULL){
    xpg[0] = hx[0] * (ic[0] + gauss_lob_point[offset[0]]);
    xpg[1] = hx[1] * (ic[1] + gauss_lob_point[offset[1]]);
    xpg[2] = hx[2] * (ic[2] + gauss_lob_point[offset[2]]);
  }

  if (wpg != NULL) *wpg = hx[0] * hx[1] * hx[2] *
		     gauss_lob_weight[offset[0]]*
		     gauss_lob_weight[offset[1]]*
		     gauss_lob_weight[offset[2]];

  if (xpg_in != NULL) {
    schnaps_real small = 1e-5;
    xpg_in[0] = xpg[0];
    xpg_in[1] = xpg[1];
    xpg_in[2] = xpg[2];

    for(int ii=0; ii < 3; ii++){
      if (fabs(xpg[ii]) < small)  xpg_in[ii] += hx[ii] * small;
      if (fabs(xpg[ii] - 1) < small)  xpg_in[ii] -= hx[ii] * small;
    }
    /* if (ix[0] == 0) xpg_in[0] += hx[0] * _SMALL; */
    /* if (ix[0] == deg[0]) xpg_in[0] -= hx[0] * _SMALL; */
    /* if (ix[1] == 0) xpg_in[1] += hx[1] * _SMALL; */
    /* if (ix[1] == deg[1]) xpg_in[1] -= hx[1] * _SMALL; */
    /* if (ix[2] == 0) xpg_in[2] += hx[2] * _SMALL; */
    /* if (ix[2] == deg[2]) xpg_in[2] -= hx[2] * _SMALL; */

    /* printf("xpg %f %f %f\n",xpg[0],xpg[1],xpg[2]); */
    /*  printf("xpg_in %f %f %f %d %d %d\n",xpg_in[0],xpg_in[1],xpg_in[2], */
    /* 	   ix,iy,iz); */
  }
} // ref_pg_vol_CG



// compute the reference coordinates xpg[3] and weight wpg of the GLOP
// ipg on the face ifa.
// and returns the index of the volume gauss point
int ref_pg_face(int deg3d[], int nraf3d[], int ifa, int ipg,
		 schnaps_real *xpg, schnaps_real *wpg, schnaps_real *xpgin) {
  // For each face, give the dimension index i
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };

  int deg[3], offset[2],nraf[3];
  schnaps_real h[3];
  int ipgxyz[3], ncpgxyz[3];
  //int ipgf=ipg;

  // approximation degree in each direction
  deg[0] = deg3d[axis_permut[ifa][0]];
  deg[1] = deg3d[axis_permut[ifa][1]];
  deg[2] = deg3d[axis_permut[ifa][2]];

  // number of subcells in each direction
  nraf[0] = nraf3d[axis_permut[ifa][0]];
  nraf[1] = nraf3d[axis_permut[ifa][1]];
  nraf[2] = nraf3d[axis_permut[ifa][2]];

  // Compute permuted indices
  int ix = ipg % (deg[0] + 1);
  ipg /= (deg[0] + 1);

  int iy = ipg % (deg[1] + 1);
  ipg /= (deg[1] + 1);

  // Equals 0 or d depending on the face
  int iz = axis_permut[ifa][3] * deg[2];

  // Compute permuted indices of the subface
  int ncx = ipg % nraf[0];
  h[0] = 1.0 / (schnaps_real) nraf[0];
  ipg /= nraf[0];

  int ncy = ipg;
  h[1] = 1.0 / (schnaps_real) nraf[1];

  // Equals 0 or nraf-1 depending on the face
  int ncz = axis_permut[ifa][3] * (nraf[2] - 1);
  h[2] = 1.0 / (schnaps_real) nraf[2];

  // Compute non permuted indices for points and subfaces
  ipgxyz[axis_permut[ifa][0]] = ix;
  ipgxyz[axis_permut[ifa][1]] = iy;
  ipgxyz[axis_permut[ifa][2]] = iz;

  ncpgxyz[axis_permut[ifa][0]] = ncx;
  ncpgxyz[axis_permut[ifa][1]] = ncy;
  ncpgxyz[axis_permut[ifa][2]] = ncz;

  // Compute the global index of the
  // Gauss-Lobatto point in the volume
  // TODO: call rather xyz_to_ipg !!!!!!!!!!!!!!!!!!!!!
  int ipg3d = ipgxyz[0] + (deg3d[0] + 1) *
    (ipgxyz[1] + (deg3d[1] + 1) *
     (ipgxyz[2] + (deg3d[2] + 1) *
      (ncpgxyz[0] + nraf3d[0] *
       (ncpgxyz[1] + nraf3d[1] *
	ncpgxyz[2]))));

  // Compute the reference coordinates of the Gauss-Lobatto point in
  // the volume
  offset[0] = gauss_lob_offset[deg[0]] + ix;
  offset[1] = gauss_lob_offset[deg[1]] + iy;

  if (xpg != NULL) {
    xpg[axis_permut[ifa][0]] = h[0] * (ncx + gauss_lob_point[offset[0]]);
    xpg[axis_permut[ifa][1]] = h[1] * (ncy + gauss_lob_point[offset[1]]);
    xpg[axis_permut[ifa][2]] = axis_permut[ifa][3];
  }

  if (wpg != NULL) *wpg = h[0] * h[1] *
		     gauss_lob_weight[offset[0]] *
		     gauss_lob_weight[offset[1]];

  // If xpgin exists, compute a point slightly INSIDE the opposite
  // subcell along the face.
  if(xpgin != NULL) {
    schnaps_real small = 1e-3;//0.001
    schnaps_real vsmall = 1e-5;//0.000001;

    xpgin[axis_permut[ifa][0]] = h[0] * (ncx + gauss_lob_point[offset[0]]);
    xpgin[axis_permut[ifa][1]] = h[1] * (ncy + gauss_lob_point[offset[1]]);

    if(axis_permut[ifa][3] == 0)
      xpgin[axis_permut[ifa][2]] = -vsmall;
    if(axis_permut[ifa][3] == 1)
      xpgin[axis_permut[ifa][2]] = 1 + vsmall;

    if(ix == 0)
      xpgin[axis_permut[ifa][0]]
	= h[0] * (ncx + gauss_lob_point[offset[0]] + small);
    if(ix == deg[0])
      xpgin[axis_permut[ifa][0]]
	= h[0] * (ncx + gauss_lob_point[offset[0]] - small);

    if(iy == 0)
      xpgin[axis_permut[ifa][1]]
	= h[1] * (ncy + gauss_lob_point[offset[1]] + small);
    if(iy == deg[1])
      xpgin[axis_permut[ifa][1]]
	= h[1] * (ncy + gauss_lob_point[offset[1]] - small);
  }

  return ipg3d;
} // ref_pg_face

// compute the reference coordinates xpg[3] and weight wpg of the GLOP
// ipg on the face ifa.
// and returns the index of the volume gauss point
int ref_pg_face_CG(int deg3d[], int nraf3d[], int ifa, int ipg,
		 schnaps_real *xpg, schnaps_real *wpg, schnaps_real *xpgin) {
  // For each face, give the dimension index i
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };

  int deg[3], offset[2],nraf[3];

  int ipgxyz[3], ncpgxyz[3];
  //int ipgf=ipg;

  // approximation degree in each direction
  deg[0] = deg3d[axis_permut[ifa][0]];
  deg[1] = deg3d[axis_permut[ifa][1]];
  deg[2] = deg3d[axis_permut[ifa][2]];


  // number of subcells in each direction
  nraf[0] = nraf3d[axis_permut[ifa][0]];
  nraf[1] = nraf3d[axis_permut[ifa][1]];
  nraf[2] = nraf3d[axis_permut[ifa][2]];

  schnaps_real h[3];
  h[0] = 1.0 / (schnaps_real) nraf[0];
  h[1] = 1.0 / (schnaps_real) nraf[1];
  h[2] = 1.0 / (schnaps_real) nraf[2];

 int nx[3] = {
    deg[0] * nraf[0] + 1,
    deg[1] * nraf[1] + 1,
    deg[2] * nraf[2] + 1};

  // permuted point indices in each direction in [0..deg * raf + 1[
  int ix = ipg % nx[0];
  ipg /= nx[0];

  int iy = ipg % nx[1];
  ipg /= nx[1];

  // Equals 0 or nx (a side of the cube)
  int iz = axis_permut[ifa][3] * (nx[2] - 1);


  // point index in each subcell in [0..deg]
  int igx = ix % deg[0];
  int igy = iy % deg[1];
  int igz = 0;

  // Compute permuted indices of the subface
  int ncx = ix / deg[0];
  int ncy = iy / deg[1];

  // Equals 0 or raf depending on the side
  int ncz = axis_permut[ifa][3] * nraf[2];


  // Compute non permuted indices for points and subfaces
  ipgxyz[axis_permut[ifa][0]] = igx;
  ipgxyz[axis_permut[ifa][1]] = igy;
  ipgxyz[axis_permut[ifa][2]] = igz;

  int iref[3];
  iref[axis_permut[ifa][0]] = ix;
  iref[axis_permut[ifa][1]] = iy;
  iref[axis_permut[ifa][2]] = iz;

  printf("ifa=%d iref=%d %d %d \n",ifa,iref[0],iref[1],iref[2]);
  printf("ifa=%d ix=%d %d %d \n",ifa,ix,iy,iz);
  printf("ifa=%d nx=%d %d %d \n",ifa,nx[0],nx[1],nx[2]);


  ncpgxyz[axis_permut[ifa][0]] = ncx;
  ncpgxyz[axis_permut[ifa][1]] = ncy;
  ncpgxyz[axis_permut[ifa][2]] = ncz;

  // Compute the global index of the
  // Gauss-Lobatto point in the volume
  /* int ipg3d = ipgxyz[0] + (deg3d[0] + 1) * */
  /*   (ipgxyz[1] + (deg3d[1] + 1) * */
  /*    (ipgxyz[2] + (deg3d[2] + 1) * */
  /*     (ncpgxyz[0] + nraf3d[0] * */
  /*      (ncpgxyz[1] + nraf3d[1] * */
  /* 	ncpgxyz[2])))); */

  int nxp[3];
  nxp[axis_permut[ifa][0]] = nx[0];
  nxp[axis_permut[ifa][1]] = nx[1];
  nxp[axis_permut[ifa][2]] = nx[2];


  int ipg3d = iref[0] + nxp[0] * (iref[1] + nxp[1] * iref[2]);

  // Compute the reference coordinates of the Gauss-Lobatto point in
  // the volume
  offset[0] = gauss_lob_offset[deg[0]] + ix;
  offset[1] = gauss_lob_offset[deg[1]] + iy;

  if (xpg != NULL) {
    xpg[axis_permut[ifa][0]] = h[0] * (ncx + gauss_lob_point[offset[0]]);
    xpg[axis_permut[ifa][1]] = h[1] * (ncy + gauss_lob_point[offset[1]]);
    xpg[axis_permut[ifa][2]] = axis_permut[ifa][3];
  }

  if (wpg != NULL) *wpg = h[0] * h[1] *
		     gauss_lob_weight[offset[0]] *
		     gauss_lob_weight[offset[1]];

  // If xpgin exists, compute a point slightly INSIDE the opposite
  // subcell along the face.
  if(xpgin != NULL) {
    schnaps_real small = 1e-4;//0.001
    schnaps_real vsmall = 1e-6;//0.000001;

    xpgin[axis_permut[ifa][0]] = h[0] * (ncx + gauss_lob_point[offset[0]]);
    xpgin[axis_permut[ifa][1]] = h[1] * (ncy + gauss_lob_point[offset[1]]);

    if(axis_permut[ifa][3] == 0)
      xpgin[axis_permut[ifa][2]] = -vsmall;
    if(axis_permut[ifa][3] == 1)
      xpgin[axis_permut[ifa][2]] = 1 + vsmall;

    if(ix == 0)
      xpgin[axis_permut[ifa][0]]
	= h[0] * (ncx + gauss_lob_point[offset[0]] + small);
    if(ix == deg[0])
      xpgin[axis_permut[ifa][0]]
	= h[0] * (ncx + gauss_lob_point[offset[0]] - small);

    if(iy == 0)
      xpgin[axis_permut[ifa][1]]
	= h[1] * (ncy + gauss_lob_point[offset[1]] + small);
    if(iy == deg[1])
      xpgin[axis_permut[ifa][1]]
	= h[1] * (ncy + gauss_lob_point[offset[1]] - small);
  }

  return ipg3d;
} // ref_pg_face_CG


// return the 1d derivative of lagrange polynomial ib at glop ipg
schnaps_real dlag(int deg, int ib, int ipg)
{
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg] + ib * (deg + 1) + ipg];
}

// Return the value psi and the gradient dpsi[3] of the basis function
// ib at point xref[3].
// Warning: the value of the gradient is not reliable if xref is on
// the boundary of a subcell (because the gradient is discontinuous)
void psi_ref(int *deg, int *nraf, int ib, schnaps_real *xref, schnaps_real *psi, schnaps_real *dpsi)
{
  schnaps_real dpsibx;
  schnaps_real dpsiby;
  schnaps_real dpsibz;

  int offset[3];

  // Starting Gauss-Lobatto point in each direction
  offset[0] = gauss_lob_offset[deg[0]];
  offset[1] = gauss_lob_offset[deg[1]];
  offset[2] = gauss_lob_offset[deg[2]];

  int ix[3],ic[3];

  ipg_to_xyz(nraf,deg,ic,ix,&ib);

  /* int ibx = ib % (deg[0] + 1); */
  /* ib /= (deg[0] + 1); */

  /* int iby = ib % (deg[1] + 1); */
  /* ib /= (deg[1] + 1); */

  /* int ibz = ib % (deg[2] + 1); */
  /* ib /= (deg[2] + 1); */

  /* int ncbx= ib % nraf[0]; */
  schnaps_real hx=1 / (schnaps_real) nraf[0];
  //ib /= nraf[0];

  //int ncby= ib % nraf[1];
  schnaps_real hy=1 / (schnaps_real) nraf[1];
  //ib /= nraf[1];

  //int ncbz= ib;
  schnaps_real hz=1 / (schnaps_real) nraf[2];

  schnaps_real psibx = 0;
  schnaps_real psiby = 0;
  schnaps_real psibz = 0;

  lagrange_polynomial(&psibx, gauss_lob_point + offset[0],
                      deg[0], ix[0], xref[0]/hx-ic[0]);
  lagrange_polynomial(&psiby, gauss_lob_point + offset[1],
                      deg[1], ix[1], xref[1]/hy-ic[1]);
  lagrange_polynomial(&psibz, gauss_lob_point + offset[2],
                      deg[2], ix[2], xref[2]/hz-ic[2]);

  int is[3];
  for(int ii=0;ii<3;ii++){
    is[ii]=xref[ii]*nraf[ii];
    assert(is[ii] < nraf[ii] && is[ii]>= 0);
  }

  int is_in_subcell= (ic[0] == is[0]) && (ic[1] == is[1])
    && (ic[2] == is[2]);


  *psi = psibx * psiby * psibz * is_in_subcell;

  if (dpsi != NULL) {

    dlagrange_polynomial(&dpsibx, gauss_lob_point + offset[0],
                         deg[0], ix[0], xref[0]);
    dlagrange_polynomial(&dpsiby, gauss_lob_point + offset[1],
                         deg[1], ix[1], xref[1]);
    dlagrange_polynomial(&dpsibz, gauss_lob_point + offset[2],
                         deg[2], ix[2], xref[2]);

    dpsi[0] = dpsibx * psiby * psibz * is_in_subcell;
    dpsi[1] = psibx * dpsiby * psibz * is_in_subcell ;
    dpsi[2] = psibx * psiby * dpsibz * is_in_subcell;
  }
}

// ibx iby ibz ncbx ncby ncbz

// Return the value psi and the gradient dpsi[3] of the basis function
// ib at point xref[3] given the subcell indices is[3].
// The computation is reliable.
void psi_ref_subcell(int *deg, int *nraf, int *is, int ib,
		     schnaps_real *xref, schnaps_real *psi, schnaps_real *dpsi) {
  schnaps_real dpsibx;
  schnaps_real dpsiby;
  schnaps_real dpsibz;


  int offset[3];

  // Starting Gauss-Lobatto point in each direction
  offset[0] = gauss_lob_offset[deg[0]];
  offset[1] = gauss_lob_offset[deg[1]];
  offset[2] = gauss_lob_offset[deg[2]];


  int ix[3],ic[3];

  ipg_to_xyz(nraf,deg,ic,ix,&ib);

  schnaps_real hx=1 / (schnaps_real) nraf[0];
  schnaps_real hy=1 / (schnaps_real) nraf[1];
  schnaps_real hz=1 / (schnaps_real) nraf[2];

  int is_in_subcell= (ic[0] == is[0]) && (ic[1] == is[1])
    && (ic[2] == is[2]);

  schnaps_real psibx = 0;
  schnaps_real psiby = 0;
  schnaps_real psibz = 0;

  lagrange_polynomial(&psibx, gauss_lob_point + offset[0],
                      deg[0], ix[0], xref[0]/hx-ic[0]);
  lagrange_polynomial(&psiby, gauss_lob_point + offset[1],
                      deg[1], ix[1], xref[1]/hy-ic[1]);
  lagrange_polynomial(&psibz, gauss_lob_point + offset[2],
                      deg[2], ix[2], xref[2]/hz-ic[2]);

  // might be useful for the future subcell case
  /* psibx *= (xref[0] <= (ncbx + 1) * hx)&&(xref[0] > ncbx * hx); */
  /* psiby *= (xref[1] <= (ncby + 1) * hy)&&(xref[1] > ncby * hy); */
  /* psibz *= (xref[2] <= (ncbz + 1) * hz)&&(xref[2] > ncbz * hz); */

  *psi = psibx * psiby * psibz * is_in_subcell ;

  if (dpsi != NULL) {

    dlagrange_polynomial(&dpsibx, gauss_lob_point + offset[0],
                         deg[0], ix[0], xref[0]);
    dlagrange_polynomial(&dpsiby, gauss_lob_point + offset[1],
                         deg[1], ix[1], xref[1]);
    dlagrange_polynomial(&dpsibz, gauss_lob_point + offset[2],
                         deg[2], ix[2], xref[2]);

    dpsi[0] = dpsibx *  psiby *  psibz * is_in_subcell;
    dpsi[1] =  psibx * dpsiby *  psibz * is_in_subcell;
    dpsi[2] =  psibx *  psiby * dpsibz * is_in_subcell;
  }
}

// Return the gradient dpsi[0..2] of the basis function ib at GLOP
// ipg.
void grad_psi_pg(int *deg, int *nraf, int ib, int ipg, schnaps_real *dpsi) {
  int offset[3];

  // glop 3d indices
  int ix[3],ic[3];
  ipg_to_xyz(nraf,deg,ic,ix,&ipg);

  // basis function 3d indices
  int ibx[3],ibc[3];
  ipg_to_xyz(nraf,deg,ibc,ibx,&ib);

  // // indices of Gauss-Lobatto points in each subcell
  // int ipgx = ipg % (deg[0] + 1);
  // ipg /= (deg[0] + 1);

  // int ipgy = ipg % (deg[1] + 1);
  // ipg /= (deg[1] + 1);

  // int ipgz = ipg % (deg[2] + 1);
  // ipg /= (deg[2] + 1);

  // // indices of each subcell and space step in each direction
  // int ncpgx= ipg % nraf[0];
  schnaps_real hx=1 / (schnaps_real) nraf[0];
  // ipg /= nraf[0];

  // int ncpgy= ipg % nraf[1];
  schnaps_real hy=1 / (schnaps_real) nraf[1];
  // ipg /= nraf[1];

  // int ncpgz= ipg;
  schnaps_real hz=1 / (schnaps_real) nraf[2];

  // basis functions indices
  // int ibx = ib % (deg[0] + 1);
  // ib /= (deg[0] + 1);

  // int iby = ib % (deg[1] + 1);
  // ib /= (deg[1] + 1);

  // int ibz = ib % (deg[2] + 1);
  // ib /= (deg[2] + 1);

  // int ncbx= ib % nraf[0];
  // ib /= nraf[0];

  // int ncby= ib % nraf[1];
  // ib /= nraf[1];

  // int ncbz= ib;

  // Number of Gauss-Lobatto points in each direction
  offset[0] = gauss_lob_dpsi_offset[deg[0]];
  offset[1] = gauss_lob_dpsi_offset[deg[1]];
  offset[2] = gauss_lob_dpsi_offset[deg[2]];

  // Computation of the value of the interpolation polynomial gradient

  schnaps_real psibx,psiby,psibz,dpsibx,dpsiby,dpsibz;

  psibx = (ix[0] == ibx[0]) * (ic[0] == ibc[0]);
  dpsibx = (ic[0] == ibc[0]) * gauss_lob_dpsi[offset[0]+ibx[0]*(deg[0]+1)+ix[0]] / hx;

  psiby = (ix[1] == ibx[1]) * (ic[1] == ibc[1]);
  dpsiby = (ic[1] == ibc[1]) * gauss_lob_dpsi[offset[1]+ibx[1]*(deg[1]+1)+ix[1]] / hy;

  psibz = (ix[2] == ibx[2]) * (ic[2] == ibc[2]);
  dpsibz = (ic[2] == ibc[2]) * gauss_lob_dpsi[offset[2]+ibx[2]*(deg[2]+1)+ix[2]] / hz;

  dpsi[0] = dpsibx*psiby*psibz;
  dpsi[1] = psibx*dpsiby*psibz;
  dpsi[2] = psibx*psiby*dpsibz;
}

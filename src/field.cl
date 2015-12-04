// A kernel used solely for managing events (cf: SOCL)
__kernel
void empty_kernel()
{
}

__kernel
void set_buffer_to_zero(__global real *w)
{
  w[get_global_id(0)] = 0.0;
}

// First and second columns: 2D loop dimensions
// Third column: face dimension
// Fourth column: 0 -> negative orientation, 1 -> positive orientation
__constant int axis_permut[6][4] = { {0, 2, 1, 0},
				     {1, 2, 0, 1},
				     {2, 0, 1, 1},
				     {2, 1, 0, 0},
				     {0, 1, 2, 1},
				     {1, 0, 2, 0} };

real dlag(int deg, int ib, int ipg)
{
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg] + ib * (deg + 1) + ipg];
}

#ifndef VARINDEX
#define VARINDEX GenericVarindex
#endif

// From a reference point find the nearest gauss point
// Warning: works only  degree 1, 2, or 3 (FIXME: why?)
int ref_ipg(int *raf, int *deg, real *xref) 
{
  real h[3] = {1.0 / raf[0], 1.0 / raf[1], 1.0 / raf[2]};

  // Get the subcell indices
  int ic[3] = {floor(xref[0] * raf[0]),
	       floor(xref[1] * raf[1]), 
	       floor(xref[2] * raf[2]) };

  // Point indices in the subcell
  int ix[3] = {floor((xref[0] - ic[0] * h[0]) / h[0] * deg[0] + 0.5),
	       floor((xref[1] - ic[1] * h[1]) / h[1] * deg[1] + 0.5),
	       floor((xref[2] - ic[2] * h[2]) / h[2] * deg[2] + 0.5) };

  return xyz_to_ipg(raf, deg, ic, ix); 
}

void compute_gradphi(const real xref[3], real gradphi[20][4]) 
{
  real x = xref[0];
  real y = xref[1];
  real z = xref[2];

  real t1 = -1.0 + z;
  real t2 = -1.0 + y;
  real t3 = t1 * t2;
  real t4 = 2 * y;
  real t5 = 2 * z;
  real t6 = 4 * x;
  real t9 = -1.0 + x;
  real t10 = t1 * t9;
  real t11 = 2 * x;
  real t12 = 4 * y;
  real t15 = t2 * t9;
  real t16 = 4 * z;
  real t24 = x * t1;
  real t27 = x * t2;
  real t33 = y * t1;
  real t38 = x * y;
  real t48 = y * t9;
  real t54 = z * t2;
  real t57 = z * t9;
  real t67 = x * z;
  real t75 = y * z;
  real t94 = t11 - 1.0;
  real t98 = 4 * t24 * t9;
  real t100 = 4 * t27 * t9;
  real t104 = 4 * t33 * t2;
  real t105 = t4 - 1.0;
  real t111 = 4 * y * t2 * t9;
  real t114 = z * t1;
  real t116 = 4 * t114 * t2;
  real t118 = 4 * t114 * t9;
  real t119 = t5 - 1.0;
  real t128 = 4 * t38 * t2;
  real t132 = 4 * t67 * t1;
  real t141 = 4 * t38 * t9;
  real t145 = 4 * t75 * t1;
  real t158 = 4 * t67 * t9;
  real t162 = 4 * t75 * t2;

  gradphi[0][0] = t3 * (t4 + t5 - 3 + t6);
  gradphi[0][1] = t10 * (t11 + t5 - 3 + t12);
  gradphi[0][2] = t15 * (t11 + t4 - 3 + t16);
  gradphi[0][3] = t3 * t9 * (t11 + t4 + t5 - 1);
  gradphi[1][0] = t3 * (-t4 - t5 - 1 + t6);
  gradphi[1][1] = t24 * (-t5 + 1 + t11 - t12);
  gradphi[1][2] = t27 * (-t4 + 1 + t11 - t16);
  gradphi[1][3] = t24 * t2 * (-t4 - t5 + t11 - 1);
  gradphi[2][0] = -t33 * (t4 - t5 - 3 + t6);
  gradphi[2][1] = -t24 * (-t5 - 3 + t11 + t12);
  gradphi[2][2] = -t38 * (t4 - 1 + t11 - t16);
  gradphi[2][3] = -t38 * t1 * (t4 - t5 - 3 + t11);
  gradphi[3][0] = -t33 * (-t4 + t5 - 1 + t6);
  gradphi[3][1] = -t10 * (t11 + t5 + 1 - t12);
  gradphi[3][2] = -t48 * (t11 - 1 - t4 + t16);
  gradphi[3][3] = -t33 * t9 * (t11 + t5 - t4 + 1);
  gradphi[4][0] = -t54 * (t4 - t5 - 1 + t6);
  gradphi[4][1] = -t57 * (t11 - t5 - 1 + t12);
  gradphi[4][2] = -t15 * (t11 + 1 + t4 - t16);
  gradphi[4][3] = -t54 * t9 * (t11 + t4 - t5 + 1);
  gradphi[5][0] = -t54 * (-t4 + t5 - 3 + t6);
  gradphi[5][1] = -t67 * (t5 - 1 + t11 - t12);
  gradphi[5][2] = -t27 * (-t4 - 3 + t11 + t16);
  gradphi[5][3] = -t67 * t2 * (-t4 + t5 - 3 + t11);
  gradphi[6][0] = t75 * (t4 + t5 - 5 + t6);
  gradphi[6][1] = t67 * (t5 - 5 + t11 + t12);
  gradphi[6][2] = t38 * (t4 - 5 + t11 + t16);
  gradphi[6][3] = t38 * z * (t4 + t5 - 5 + t11);
  gradphi[7][0] = t75 * (-t4 - t5 + 1 + t6);
  gradphi[7][1] = t57 * (t11 - t5 + 3 - t12);
  gradphi[7][2] = t48 * (t11 - t4 + 3 - t16);
  gradphi[7][3] = t75 * t9 * (t11 - t5 - t4 + 3);
  gradphi[8][0] = -4 * t3 * t94;
  gradphi[8][1] = -t98;
  gradphi[8][2] = -t100;
  gradphi[8][3] = -4 * t24 * t15;
  gradphi[9][0] = -t104;
  gradphi[9][1] = -4 * t1 * t105 * t9;
  gradphi[9][2] = -t111;
  gradphi[9][3] = -4 * t33 * t15;
  gradphi[10][0] = -t116;
  gradphi[10][1] = -t118;
  gradphi[10][2] = -4 * t119 * t2 * t9;
  gradphi[10][3] = -4 * t114 * t15;
  gradphi[11][0] = t104;
  gradphi[11][1] = 4 * t24 * t105;
  gradphi[11][2] = t128;
  gradphi[11][3] = 4 * t38 * t3;
  gradphi[12][0] = t116;
  gradphi[12][1] = t132;
  gradphi[12][2] = 4 * x * t119 * t2;
  gradphi[12][3] = 4 * t67 * t3;
  gradphi[13][0] = 4 * t33 * t94;
  gradphi[13][1] = t98;
  gradphi[13][2] = t141;
  gradphi[13][3] = 4 * t38 * t10;
  gradphi[14][0] = -t145;
  gradphi[14][1] = -t132;
  gradphi[14][2] = -4 * t38 * t119;
  gradphi[14][3] = -4 * t38 * t114;
  gradphi[15][0] = t145;
  gradphi[15][1] = t118;
  gradphi[15][2] = 4 * y * t119 * t9;
  gradphi[15][3] = 4 * t75 * t10;
  gradphi[16][0] = 4 * t54 * t94;
  gradphi[16][1] = t158;
  gradphi[16][2] = t100;
  gradphi[16][3] = 4 * t67 * t15;
  gradphi[17][0] = t162;
  gradphi[17][1] = 4 * z * t105 * t9;
  gradphi[17][2] = t111;
  gradphi[17][3] = 4 * t75 * t15;
  gradphi[18][0] = -t162;
  gradphi[18][1] = -4 * t67 * t105;
  gradphi[18][2] = -t128;
  gradphi[18][3] = -4 * t38 * t54;
  gradphi[19][0] = -4 * t75 * t94;
  gradphi[19][1] = -t158;
  gradphi[19][2] = -t141;
  gradphi[19][3] = -4 * t38 * t57;
}

inline void compute_xphy(__constant real *physnode,
		  real gradphi[20][4],
		  real xphy[3])
{
  // FIXME: FP_FAST_FMA checks for double, we should also look at
  // FP_FAST_FMAF for single-precision comparison.
#ifdef FP_FAST_FMA
  xphy[0] = 0;
  xphy[1] = 0;
  xphy[2] = 0;
  for(int i = 0; i < 20; ++i) {
    real gp = gradphi[i][3];
    int i3 = 3 * i;
    xphy[0] = fma(physnode[i3], gp, xphy[0]);
    xphy[1] = fma(physnode[i3 + 1], gp, xphy[1]);
    xphy[2] = fma(physnode[i3 + 2], gp, xphy[2]);
  }
#else
  for(int ii = 0; ii < 3; ++ii) {
    xphy[ii] = 0;
    for(int i = 0; i < 20; ++i) {
      xphy[ii] += physnode[3 * i + ii] * gradphi[i][3];
    }
  }
#endif
}

inline void compute_dtau(__constant real *physnode,
			 real gradphi[20][4],
			 real dtau[3][3])
{
  for(int ii = 0; ii < 3; ii++) {
    dtau[ii][0] = 0;
    dtau[ii][1] = 0;
    dtau[ii][2] = 0;

    for(int i = 0; i < 20; ++i) {
      real pn = physnode[3 * i + ii];
#ifdef FP_FAST_FMA
      dtau[ii][0] = fma(pn, gradphi[i][0], dtau[ii][0]);
      dtau[ii][1] = fma(pn, gradphi[i][1], dtau[ii][1]);
      dtau[ii][2] = fma(pn, gradphi[i][2], dtau[ii][2]);
#else
      for(int jj = 0; jj < 3; jj++) {
	dtau[ii][jj] += pn * gradphi[i][jj];
      }
#endif
    }
  }
}

inline void compute_codtau(real dtau[3][3], real codtau[3][3])
{
#ifdef FP_FAST_FMA
  codtau[0][0] = fma( dtau[1][1], dtau[2][2], -dtau[1][2] * dtau[2][1]);
  codtau[0][1] = fma(-dtau[1][0], dtau[2][2],  dtau[1][2] * dtau[2][0]);
  codtau[0][2] = fma( dtau[1][0], dtau[2][1], -dtau[1][1] * dtau[2][0]);
  codtau[1][0] = fma(-dtau[0][1], dtau[2][2],  dtau[0][2] * dtau[2][1]);
  codtau[1][1] = fma( dtau[0][0], dtau[2][2], -dtau[0][2] * dtau[2][0]);
  codtau[1][2] = fma(-dtau[0][0], dtau[2][1],  dtau[0][1] * dtau[2][0]);
  codtau[2][0] = fma( dtau[0][1], dtau[1][2], -dtau[0][2] * dtau[1][1]);
  codtau[2][1] = fma(-dtau[0][0], dtau[1][2],  dtau[0][2] * dtau[1][0]);
  codtau[2][2] = fma( dtau[0][0], dtau[1][1], -dtau[0][1] * dtau[1][0]);
#else
  codtau[0][0] =  dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
  codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
  codtau[0][2] =  dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
  codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
  codtau[1][1] =  dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
  codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
  codtau[2][0] =  dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
  codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
  codtau[2][2] =  dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
#endif
}

inline void compute_dphi(real dphiref[3], real codtau[3][3], real dphi[3])
{

#ifdef FP_FAST_FMA
  dphi[0] = fma(codtau[0][0], dphiref[0], 
		fma(codtau[0][1], dphiref[1],
		    codtau[0][2] * dphiref[2]) );
    
  dphi[1] = fma(codtau[1][0], dphiref[0], 
		fma(codtau[1][1], dphiref[1],
		    codtau[1][2] * dphiref[2]) );
    
  dphi[2] = fma(codtau[2][0], dphiref[0], 
		fma(codtau[2][1], dphiref[1],
		    codtau[2][2] * dphiref[2]) );
#else
  for(int ii = 0; ii < 3; ii++) {
    dphi[ii]=0;
    for(int jj = 0; jj < 3; jj++) {
      dphi[ii] += codtau[ii][jj] * dphiref[jj];
    }
  }
#endif
}

inline void ComputeNormal(real codtau[3][3], int ifa, real vnds[3])
{
  int h20_refnormal[6][3]={ {0, -1,  0},
			    {1,  0,  0},
			    {0,  1,  0},
			    {-1, 0,  0},
			    {0,  0,  1},
			    {0,  0, -1}};
  // FIXME: replace some of this code with select to avoid unnecessary
  // multiplies with zero.
#ifdef FP_FAST_FMA
  vnds[0] = fma(codtau[0][0], h20_refnormal[ifa][0],
		fma(codtau[0][1], h20_refnormal[ifa][1],
		    codtau[0][2] * h20_refnormal[ifa][2]) );
    
  vnds[1] = fma(codtau[1][0], h20_refnormal[ifa][0],
		fma(codtau[1][1], h20_refnormal[ifa][1],
		    codtau[1][2] * h20_refnormal[ifa][2]) );

  vnds[2] = fma(codtau[2][0], h20_refnormal[ifa][0],
		fma(codtau[2][1], h20_refnormal[ifa][1],
		    codtau[2][2] * h20_refnormal[ifa][2]) );
#else
  for(int ii = 0; ii < 3; ii++) {
    vnds[ii] = 0;
    for(int jj = 0; jj < 3; jj++) {
      vnds[ii] += codtau[ii][jj] * h20_refnormal[ifa][jj];
    }
  }
#endif
}

inline void Phy2Ref(__constant real *physnode, real *xphy, real *xref) 
{
#define ITERNEWTON 10
  real dxref[3], dxphy[3];
  xref[0] = 0.5;
  xref[1] = 0.5;
  xref[2] = 0.5;
  for(int iter = 0; iter < ITERNEWTON; ++iter ) {
    real dtau[3][3];
    real codtau[3][3];
    
    real gradphi[20][4];
    compute_gradphi(xref, gradphi);
    compute_xphy(physnode, gradphi, dxphy);
    compute_dtau(physnode, gradphi, dtau);
    compute_codtau(dtau, codtau);
    
    dxphy[0] -= xphy[0];
    dxphy[1] -= xphy[1];
    dxphy[2] -= xphy[2];

#ifdef FP_FAST_FMA
    real det = fma(dtau[0][0], codtau[0][0],
		   fma(dtau[0][1], codtau[0][1],
		       dtau[0][2] * codtau[0][2]) );
    real overdet = 1.0 / det;
    
    dxref[0] = fma(codtau[0][0], dxphy[0],
		   fma(codtau[1][0], dxphy[1],
		       codtau[2][0] * dxphy[2]) );
    xref[0] = fma(-dxref[0], overdet, xref[0]);
    
    dxref[1] = fma(codtau[0][1], dxphy[0],
		   fma(codtau[1][1], dxphy[1],
		       codtau[2][1] * dxphy[2]) );
    xref[1] = fma(-dxref[1], overdet, xref[1]);
    
    dxref[2] = fma(codtau[0][2], dxphy[0],
		   fma(codtau[1][2], dxphy[1],
		       codtau[2][2] * dxphy[2]) );
    xref[2] = fma(-dxref[2], overdet, xref[2]);
#else
    real overdet = 1.0 / (  dtau[0][0] * codtau[0][0]
			    + dtau[0][1] * codtau[0][1]
			    + dtau[0][2] * codtau[0][2] );
 
    for(int ii = 0; ii < 3; ++ii) {
      dxref[ii]
	= codtau[0][ii] * dxphy[0] 
	+ codtau[1][ii] * dxphy[1] 
	+ codtau[2][ii] * dxphy[2];
      xref[ii] -= dxref[ii] * overdet;
    }
#endif
  }
}

// Given parameters deg and nraf and input ipg, compute the reference
// coordinages (xpg) and the weght of the Gauss piont (wpg).
inline void ref_pg_vol(const int *deg, const int *nraf, 
		       const int ipg, real *xpg, real *wpg)
{
  int ix[3], ic[3];
  ipg_to_xyz(nraf, deg, ic, ix, &ipg);

  real h[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  int offset[3] = {gauss_lob_offset[deg[0]] + ix[0],
		   gauss_lob_offset[deg[1]] + ix[1],
		   gauss_lob_offset[deg[2]] + ix[2] };

  xpg[0] = h[0] * (ic[0] + gauss_lob_point[offset[0]]);
  xpg[1] = h[1] * (ic[1] + gauss_lob_point[offset[1]]);
  xpg[2] = h[2] * (ic[2] + gauss_lob_point[offset[2]]);
  
  *wpg = h[0] * h[1] * h[2] *
    gauss_lob_weight[offset[0]]*
    gauss_lob_weight[offset[1]]*
    gauss_lob_weight[offset[2]];
}

inline void permute_indices(const int *i, int *pi, const int ifa)
{
  pi[0] = i[axis_permut[ifa][0]];
  pi[1] = i[axis_permut[ifa][1]];
  pi[2] = i[axis_permut[ifa][2]];
}

inline void unpermute_indices(int *i, const int *pi, const int ifa)
{
  i[axis_permut[ifa][0]] = pi[0];
  i[axis_permut[ifa][1]] = pi[1];
  i[axis_permut[ifa][2]] = pi[2];
}

void compute_xpgin(const int *raf, const int *deg,
		   const int* ic, const int* ix,
		   const int ifa, real *xpgin)
{
  const int paxis[4] = {axis_permut[ifa][0],
			axis_permut[ifa][1],
			axis_permut[ifa][2],
			axis_permut[ifa][3]};
  int praf[3];
  permute_indices(raf, praf, ifa);
  
  int pdeg[3];
  permute_indices(deg, pdeg, ifa);

  int pic[3];
  permute_indices(ic, pic, ifa);
  
  int pix[3];
  permute_indices(ix, pix, ifa);
  
  const int poffset[2] = {gauss_lob_offset[pdeg[0]] + pix[0],
			  gauss_lob_offset[pdeg[1]] + pix[1]};
  
  const real h[3] = {1.0 / praf[0], 1.0 / praf[1], 1.0 / praf[2] };
  
  const real  small = 1e-3;
  const real vsmall = 1e-6;

  xpgin[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[poffset[0]]);
  if(pix[0] == 0)
    xpgin[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[poffset[0]] + small);
  if(pix[0] == pdeg[0])
    xpgin[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[poffset[0]] - small);
  
  xpgin[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[poffset[1]]);
  if(pix[1] == 0)
    xpgin[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[poffset[1]] + small);
  if(pix[1] == pdeg[1])
    xpgin[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[poffset[1]] - small);

  xpgin[paxis[2]] = paxis[3] == 0 ? -vsmall: 1.0 + vsmall;
}

// Compute the permuted volumic indices based corresponding to the
// facial index.
void ipgf_to_pxyz(const int *deg,
		  const int *raf,
		  const int ifa,
		  int ipgf,
		  int *pic,
		  int *pix)
{
  const int paxis[4] = {axis_permut[ifa][0],
			axis_permut[ifa][1],
			axis_permut[ifa][2],
			axis_permut[ifa][3]};

  // number of subcells in each permuted direction
  int praf[3];
  permute_indices(raf, praf, ifa);
  
  // Polynomial degree in each permuted direction
  int pdeg[3];
  permute_indices(deg, pdeg, ifa);
  
  // permuted point index in subcell
  pix[0] = ipgf % (pdeg[0] + 1);
  ipgf /= (pdeg[0] + 1);
  pix[1] = ipgf % (pdeg[1] + 1);
  ipgf /= (pdeg[1] + 1);
  // pix[2] equals 0 or d depending on the face
  pix[2] = select(0, pdeg[2], paxis[3]); 

  // Compute permuted subcell  indices of the subface
  pic[0] = ipgf % praf[0];
  pic[1] = ipgf / praf[0];
  pic[2] = select(0, praf[2] - 1, paxis[3]); // Equals 0 or raf-1
}

int ref_pg_face(const int *deg,
		const int *raf,
		const int ifa,  // face index
		int ipgf,       // facial index of the point
                real *xpg,      // reference-space coordinates
		real *wpg)      // weight for the point
{
  const int paxis[4] = {axis_permut[ifa][0],
			axis_permut[ifa][1],
			axis_permut[ifa][2],
			axis_permut[ifa][3]};
  
  int pic[3];
  int pix[3];
  ipgf_to_pxyz(deg, raf, ifa, ipgf, pic, pix);

  // number of subcells in each permuted direction
  int praf[3];
  permute_indices(raf, praf, ifa);
  
  // Polynomial degree in each permuted direction
  int pdeg[3];
  permute_indices(deg, pdeg, ifa);
  
  // non-permuted subcell index
  int ic[3];
  unpermute_indices(ic, pic, ifa);

  // non-permuted point index
  int ix[3];
  unpermute_indices(ix, pix, ifa);

  const real ph[3] = {1.0 / praf[0], 1.0 / praf[1], 1.0 / praf[2] };
  
  // Compute the global index of the Gauss-Lobatto point in the volume
  int ipgv = xyz_to_ipg(raf, deg, ic, ix); 

  // Compute the reference coordinates of the Gauss-Lobatto point in
  // the volume
  const int poffset[2] = {gauss_lob_offset[pdeg[0]] + pix[0],
			  gauss_lob_offset[pdeg[1]] + pix[1]};

  xpg[paxis[0]] = ph[0] * (pic[0] + gauss_lob_point[poffset[0]]);
  xpg[paxis[1]] = ph[1] * (pic[1] + gauss_lob_point[poffset[1]]);
  xpg[paxis[2]] = paxis[3];

  *wpg = ph[0] * ph[1] *
    gauss_lob_weight[poffset[0]] * gauss_lob_weight[poffset[1]];

  return ipgv;
}

#ifndef NUMFLUX
#define NUMFLUX NumFlux
#endif

void NumFlux(real wL[], real wR[], real *vnorm, real *flux) {
  real vn = sqrt(0.5) * (vnorm[0] + vnorm[1]);

  real vnp = vn > 0 ? vn : 0;
  real vnm = vn - vnp;

  flux[0] = vnp * wL[0] + vnm * wR[0];
}

// Return the component of the vlasov velocity with index id.
inline real vlasov_vel(const int id, const int md)
{
#ifndef vlasov_vmax
#define vlasov_vmax 0.5
#endif
  int mid = md / 2;
  real dv = vlasov_vmax / mid;
  return (id - mid) * dv;
}

// Sample flux for 2D Vlasov equation
inline void vlaTransNumFlux2d(real wL[], real wR[], real *vnorm, real *flux)
{
#ifndef vlasov_mx
#define vlasov_mx 1
#endif

#ifndef vlasov_my
#define vlasov_my 1
#endif
  for(int ix = 0; ix < vlasov_mx; ++ix) {
    real vx = vlasov_vel(ix, vlasov_mx);

    for(int iy = 0; iy < vlasov_my; ++iy) {
      real vy = vlasov_vel(iy, vlasov_my);
      
      real vn = vx * vnorm[0]	+ vy * vnorm[1];
      real vnp = vn > 0 ? vn : 0;
      real vnm = vn - vnp;
      
      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      flux[im] = vnp * wL[im] + vnm * wR[im];
    }
  }
}

#ifndef _M
#define _M 1
#endif

void cemracs2014_TransBoundaryFlux(real x[3], real t, 
				   real wL[], real *vnorm,
				   real *flux)
{
  real wR[_M];
  int m = vlasov_mx * vlasov_my;
  for(int i = 0; i < m; ++i)
    wR[i] = 0;
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

// Sample boundary flux
inline void BoundaryFlux(real x[3], real t, real *wL, real *vnorm,
                  real *flux) 
{
  real wR[_M];
  real vx = sqrt(0.5) * (x[0] + x[1]);
  wR[0] = cos(vx - t);

  NUMFLUX(wL, wR, vnorm, flux);
}

//! \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
inline real wglop(int deg, int i) 
{
  return gauss_lob_weight[gauss_lob_offset[deg] + i];
}

void get_dtau(real x, real y, real z,
	      __constant real *p, real dtau[][3]) 
{
  // Gradient of the shape functions and value (4th component) of the
  // shape functions

  /* real gradphi[20][3]; */
  /* //real x,y,z; */
  /* // this fills the values of gradphi */
  /* gradphi[0][0] = (-1 + z) * (-1 + y) * (2 * y + 2 * z + 4 * x - 3); */
  /* gradphi[0][1] = (-1 + z) * (-1 + x) * (2 * x + 2 * z - 3 + 4 * y); */
  /* gradphi[0][2] = (-1 + y) * (-1 + x) * (2 * x + 2 * y - 3 + 4 * z); */
  /* gradphi[1][0] = (-1 + z) * (-1 + y) * (-2 * y - 2 * z - 1 + 4 * x); */
  /* gradphi[1][1] = x * (-1 + z) * (-2 * z + 1 - 4 * y + 2 * x); */
  /* gradphi[1][2] = x * (-1 + y) * (-2 * y + 1 - 4 * z + 2 * x); */
  /* gradphi[2][0] = -y * (-1 + z) * (2 * y - 2 * z - 3 + 4 * x); */
  /* gradphi[2][1] = -x * (-1 + z) * (-2 * z - 3 + 4 * y + 2 * x); */
  /* gradphi[2][2] = -x * y * (2 * y - 4 * z - 1 + 2 * x); */
  /* gradphi[3][0] = -y * (-1 + z) * (-2 * y + 2 * z - 1 + 4 * x); */
  /* gradphi[3][1] = -(-1 + z) * (-1 + x) * (2 * x + 2 * z + 1 - 4 * y); */
  /* gradphi[3][2] = -y * (-1 + x) * (2 * x - 1 + 4 * z - 2 * y); */
  /* gradphi[4][0] = -z * (-1 + y) * (2 * y - 2 * z - 1 + 4 * x); */
  /* gradphi[4][1] = -z * (-1 + x) * (2 * x - 2 * z - 1 + 4 * y); */
  /* gradphi[4][2] = -(-1 + y) * (-1 + x) * (2 * x - 4 * z + 2 * y + 1); */
  /* gradphi[5][0] = -z * (-1 + y) * (-2 * y + 2 * z + 4 * x - 3); */
  /* gradphi[5][1] = -x * z * (2 * z - 4 * y - 1 + 2 * x); */
  /* gradphi[5][2] = -x * (-1 + y) * (-2 * y - 3 + 4 * z + 2 * x); */
  /* gradphi[6][0] = y * z * (2 * y + 2 * z + 4 * x - 5); */
  /* gradphi[6][1] = x * z * (2 * z + 4 * y - 5 + 2 * x); */
  /* gradphi[6][2] = x * y * (2 * y + 4 * z - 5 + 2 * x); */
  /* gradphi[7][0] = y * z * (-2 * y - 2 * z + 1 + 4 * x); */
  /* gradphi[7][1] = z * (-1 + x) * (2 * x - 2 * z + 3 - 4 * y); */
  /* gradphi[7][2] = y * (-1 + x) * (2 * x - 2 * y + 3 - 4 * z); */
  /* gradphi[8][0] = -4 * (-1 + z) * (-1 + y) * (2 * x - 1); */
  /* gradphi[8][1] = -4 * x * (-1 + z) * (-1 + x); */
  /* gradphi[8][2] = -4 * x * (-1 + y) * (-1 + x); */
  /* gradphi[9][0] = -4 * y * (-1 + z) * (-1 + y); */
  /* gradphi[9][1] = -4 * (-1 + z) * (2 * y - 1) * (-1 + x); */
  /* gradphi[9][2] = -4 * y * (-1 + y) * (-1 + x); */
  /* gradphi[10][0] = -4 * z * (-1 + z) * (-1 + y); */
  /* gradphi[10][1] = -4 * z * (-1 + z) * (-1 + x); */
  /* gradphi[10][2] = -4 * (2 * z - 1) * (-1 + y) * (-1 + x); */
  /* gradphi[11][0] = 4 * y * (-1 + z) * (-1 + y); */
  /* gradphi[11][1] = 4 * x * (-1 + z) * (2 * y - 1); */
  /* gradphi[11][2] = 4 * x * y * (-1 + y); */
  /* gradphi[12][0] = 4 * z * (-1 + z) * (-1 + y); */
  /* gradphi[12][1] = 4 * x * z * (-1 + z); */
  /* gradphi[12][2] = 4 * x * (2 * z - 1) * (-1 + y); */
  /* gradphi[13][0] = 4 * y * (-1 + z) * (2 * x - 1); */
  /* gradphi[13][1] = 4 * x * (-1 + z) * (-1 + x); */
  /* gradphi[13][2] = 4 * x * y * (-1 + x); */
  /* gradphi[14][0] = -4 * y * z * (-1 + z); */
  /* gradphi[14][1] = -4 * x * z * (-1 + z); */
  /* gradphi[14][2] = -4 * x * y * (2 * z - 1); */
  /* gradphi[15][0] = 4 * y * z * (-1 + z); */
  /* gradphi[15][1] = 4 * z * (-1 + z) * (-1 + x); */
  /* gradphi[15][2] = 4 * y * (2 * z - 1) * (-1 + x); */
  /* gradphi[16][0] = 4 * z * (-1 + y) * (2 * x - 1); */
  /* gradphi[16][1] = 4 * x * z * (-1 + x); */
  /* gradphi[16][2] = 4 * x * (-1 + y) * (-1 + x); */
  /* gradphi[17][0] = 4 * y * z * (-1 + y); */
  /* gradphi[17][1] = 4 * z * (2 * y - 1) * (-1 + x); */
  /* gradphi[17][2] = 4 * y * (-1 + y) * (-1 + x); */
  /* gradphi[18][0] = -4 * y * z * (-1 + y); */
  /* gradphi[18][1] = -4 * x * z * (2 * y - 1); */
  /* gradphi[18][2] = -4 * x * y * (-1 + y); */
  /* gradphi[19][0] = -4 * y * z * (2 * x - 1); */
  /* gradphi[19][1] = -4 * x * z * (-1 + x); */
  /* gradphi[19][2] = -4 * x * y * (-1 + x); */
  /* for(int ii=0;ii<3;ii++){ */
  /*   for(int jj=0;jj<3;jj++){ */
  /*     dtau[ii][jj]=0; */
  /*   } */
  /*   for(int i=0;i<20;i++){ */
  /*     //printf("xyzphy=%f %f %f \n",physnode[3*i+0],physnode[3*i+1],physnode[3*i+2]); */
  /*     for(int jj=0;jj<3;jj++){ */
  /*       dtau[ii][jj]+=physnode[3*i+ii]*gradphi[i][jj];; */
  /*     } */
  /*   } */
  /* } */

  dtau[0][0]
    = 2 * (-1 + z) * (-1 + y) * (-1 + x) * p[0]
    + 2 * x * (-1 + z) * (-1+ y ) * p[3]
    - 2 * x * y * (-1 + z) * p[6]
    - 2 * y * (-1 + z) * (-1 + x) * p[9]
    - 2 * z * (-1 + y) * (-1 + x) * p[12]
    - 2 * x * z * (-1 + y) * p[15]
    + 2  * x * y * z * p[18]
    +2 * y * z * (-1+x) * p[21]
    +(-1+z) * (-1+y) * (2 * x + 2 * y + 2 * z - 1) * p[0]
    +(-1+z) * (-1+y) * (-2 * y-2 * z+2 * x-1) * p[3]
    -y * (-1+z) * (2 * y-2 * z-3 + 2 * x) * p[6]
    -y * (-1+z) * (2 * x+2 * z+1 - 2 * y) * p[9]
    -z * (-1+y) * (2 * x+2 * y-2 * z+1) * p[12]
    -z * (-1+y) * (-2 * y+2 * z-3 + 2 * x) * p[15]
    +y * z * (2 * y+2 * z-5+2 * x) * p[18]
    +y * z * (2 * x-2 * z+3-2 * y) * p[21]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[24]
    -4 * x * (-1+z) * (-1+y) * p[24]
    -4 * y * (-1+z) * (-1+y) * p[27]
    -4 * z * (-1+z) * (-1+y) * p[30]
    +4 * y * (-1+z) * (-1+y) * p[33]
    +4 * z * (-1+z) * (-1+y) * p[36]
    +4 * y * (-1+z) * (-1+x) * p[39]
    +4 * x * y * (-1+z) * p[39]
    -4 * y * z * (-1+z) * p[42]
    +4 * y * z * (-1+z) * p[45]
    +4 * z * (-1+y) * (-1+x) * p[48]
    +4 * x * z * (-1+y) * p[48]
    +4 * y * z * (-1+y) * p[51]
    -4 * y * z * (-1+y) * p[54]
    -4 * y * z * (-1+x) * p[57]
    -4 * x * y * z * p[57];

  dtau[0][1]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[0]
    -2 * x * (-1+z) * (-1+y) * p[3]
    -2 * x * y * (-1+z) * p[6]
    +2 * y * (-1+z) * (-1+x) * p[9]
    -2 * z * (-1+y) * (-1+x) * p[12]
    +2 * x * z * (-1+y) * p[15]
    +2 * x * y * z * p[18]
    -2 * y * z * (-1+x) * p[21]
    +(-1+z) * (-1+x) * (2 * x+2 * y+2 * z-1) * p[0]
    +x * (-1+z) * (-2 * y-2 * z+2 * x-1) * p[3]
    -x * (-1+z) * (2 * y-2 * z-3+2 * x) * p[6]
    -(-1+z) * (-1+x) * (2 * x+2 * z+1-2 * y) * p[9]
    -z * (-1+x) * (2 * x+2 * y-2 * z+1) * p[12]
    -x * z * (-2 * y+2 * z-3+2 * x) * p[15]
    +x * z * (2 * y+2 * z-5+2 * x) * p[18]
    +z * (-1+x) * (2 * x-2 * z+3-2 * y) * p[21]
    -4 * x * (-1+z) * (-1+x) * p[24]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[27]
    -4 * y * (-1+z) * (-1+x) * p[27]
    -4 * z * (-1+z) * (-1+x) * p[30]
    +4 * x * (-1+z) * (-1+y) * p[33]
    +4 * x * y * (-1+z) * p[33]
    +4 * x * z * (-1+z) * p[36]
    +4 * x * (-1+z) * (-1+x) * p[39]
    -4 * x * z * (-1+z) * p[42]
    +4 * z * (-1+z) * (-1+x) * p[45]
    +4 * x * z * (-1+x) * p[48]
    +4 * z * (-1+y) * (-1+x) * p[51]
    +4 * y * z * (-1+x) * p[51]
    -4 * x * z * (-1+y) * p[54]
    -4 * x * y * z * p[54]-4 * x * z * (-1+x) * p[57];

  dtau[0][2]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[0]
    -2 * x * (-1+z) * (-1+y) * p[3]
    +2 * x * y * (-1+z) * p[6]
    -2 * y * (-1+z) * (-1+x) * p[9]
    +2 * z * (-1+y) * (-1+x) * p[12]
    -2 * x * z * (-1+y) * p[15]
    +2 * x * y * z * p[18]
    -2 * y * z * (-1+x) * p[21]
    +(-1+y) * (-1+x) * (2 * x+2 * y+2 * z-1) * p[0]
    +x * (-1+y) * (-2 * y-2 * z+2 * x-1) * p[3]
    -x * y * (2 * y-2 * z-3+2 * x) * p[6]
    -y * (-1+x) * (2 * x+2 * z+1-2 * y) * p[9]
    -(-1+y) * (-1+x) * (2 * x+2 * y-2 * z+1) * p[12]
    -x * (-1+y) * (-2 * y+2 * z-3+2 * x) * p[15]
    +x * y * (2 * y+2 * z-5+2 * x) * p[18]
    +y * (-1+x) * (2 * x-2 * z+3-2 * y) * p[21]
    -4 * x * (-1+y) * (-1+x) * p[24]
    -4 * y * (-1+y) * (-1+x) * p[27]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[30]
    -4 * z * (-1+y) * (-1+x) * p[30]
    +4 * x * y * (-1+y) * p[33]
    +4 * x * (-1+z) * (-1+y) * p[36]
    +4 * x * z * (-1+y) * p[36]
    +4 * x * y * (-1+x) * p[39]
    -4 * x * y * (-1+z) * p[42]
    -4 * x * y * z * p[42]
    +4 * y * (-1+z) * (-1+x) * p[45]
    +4 * y * z * (-1+x) * p[45]
    +4 * x * (-1+y) * (-1+x) * p[48]
    +4 * y * (-1+y) * (-1+x) * p[51]
    -4 * x * y * (-1+y) * p[54]
    -4 * x * y * (-1+x) * p[57];

  dtau[1][0]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[1]
    +2 * x * (-1+z) * (-1+y) * p[4]
    -2 * x * y * (-1+z) * p[7]
    -2 * y * (-1+z) * (-1+x) * p[10]
    -2 * z * (-1+y) * (-1+x) * p[13]
    -2 * x * z * (-1+y) * p[16]
    +2 * x * y * z * p[19]
    +2 * y * z * (-1+x) * p[22]
    +(-1+z) * (-1+y) * (2 * x+2 * y+2 * z-1) * p[1]
    +(-1+z) * (-1+y) * (-2 * y-2 * z+2 * x-1) * p[4]
    -y * (-1+z) * (2 * y-2 * z-3+2 * x) * p[7]
    -y * (-1+z) * (2 * x+2 * z+1-2 * y) * p[10]
    -z * (-1+y) * (2 * x+2 * y-2 * z+1) * p[13]
    -z * (-1+y) * (-2 * y+2 * z-3+2 * x) * p[16]
    +y * z * (2 * y+2 * z-5+2 * x) * p[19]
    +y * z * (2 * x-2 * z+3-2 * y) * p[22]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[25]
    -4 * x * (-1+z) * (-1+y) * p[25]
    -4 * y * (-1+z) * (-1+y) * p[28]
    -4 * z * (-1+z) * (-1+y) * p[31]
    +4 * y * (-1+z) * (-1+y) * p[34]
    +4 * z * (-1+z) * (-1+y) * p[37]
    +4 * y * (-1+z) * (-1+x) * p[40]
    +4 * x * y * (-1+z) * p[40]
    -4 * y * z * (-1+z) * p[43]
    +4 * y * z * (-1+z) * p[46]
    +4 * z * (-1+y) * (-1+x) * p[49]
    +4 * x * z * (-1+y) * p[49]
    +4 * y * z * (-1+y) * p[52]
    -4 * y * z * (-1+y) * p[55]
    -4 * y * z * (-1+x) * p[58]
    -4 * x * y * z * p[58];

  dtau[1][1]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[1]
    -2 * x * (-1+z) * (-1+y) * p[4]
    -2 * x * y * (-1+z) * p[7]
    +2 * y * (-1+z) * (-1+x) * p[10]
    -2 * z * (-1+y) * (-1+x) * p[13]
    +2 * x * z * (-1+y) * p[16]
    +2 * x * y * z * p[19]
    -2 * y * z * (-1+x) * p[22]
    +(-1+z) * (-1+x) * (2 * x+2 * y+2 * z-1) * p[1]
    +x * (-1+z) * (-2 * y-2 * z+2 * x-1) * p[4]
    -x * (-1+z) * (2 * y-2 * z-3+2 * x) * p[7]
    -(-1+z) * (-1+x) * (2 * x+2 * z+1-2 * y) * p[10]
    -z * (-1+x) * (2 * x+2 * y-2 * z+1) * p[13]
    -x * z * (-2 * y+2 * z-3+2 * x) * p[16]
    +x * z * (2 * y+2 * z-5+2 * x) * p[19]
    +z * (-1+x) * (2 * x-2 * z+3-2 * y) * p[22]
    -4 * x * (-1+z) * (-1+x) * p[25]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[28]
    -4 * y * (-1+z) * (-1+x) * p[28]
    -4 * z * (-1+z) * (-1+x) * p[31]
    +4 * x * (-1+z) * (-1+y) * p[34]
    +4 * x * y * (-1+z) * p[34]
    +4 * x * z * (-1+z) * p[37]
    +4 * x * (-1+z) * (-1+x) * p[40]
    -4 * x * z * (-1+z) * p[43]
    +4 * z * (-1+z) * (-1+x) * p[46]
    +4 * x * z * (-1+x) * p[49]
    +4 * z * (-1+y) * (-1+x) * p[52]
    +4 * y * z * (-1+x) * p[52]
    -4 * x * z * (-1+y) * p[55]
    -4 * x * y * z * p[55]
    -4 * x * z * (-1+x) * p[58];

  dtau[1][2]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[1]
    -2 * x * (-1+z) * (-1+y) * p[4]
    +2 * x * y * (-1+z) * p[7]
    -2 * y * (-1+z) * (-1+x) * p[10]
    +2 * z * (-1+y) * (-1+x) * p[13]
    -2 * x * z * (-1+y) * p[16]
    +2 * x * y * z * p[19]
    -2 * y * z * (-1+x) * p[22]
    +(-1+y) * (-1+x) * (2 * x+2 * y+2 * z-1) * p[1]
    +x * (-1+y) * (-2 * y-2 * z+2 * x-1) * p[4]
    -x * y * (2 * y-2 * z-3+2 * x) * p[7]
    -y * (-1+x) * (2 * x+2 * z+1-2 * y) * p[10]
    -(-1+y) * (-1+x) * (2 * x+2 * y-2 * z+1) * p[13]
    -x * (-1+y) * (-2 * y+2 * z-3+2 * x) * p[16]
    +x * y * (2 * y+2 * z-5+2 * x) * p[19]
    +y * (-1+x) * (2 * x-2 * z+3-2 * y) * p[22]
    -4 * x * (-1+y) * (-1+x) * p[25]
    -4 * y * (-1+y) * (-1+x) * p[28]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[31]
    -4 * z * (-1+y) * (-1+x) * p[31]
    +4 * x * y * (-1+y) * p[34]
    +4 * x * (-1+z) * (-1+y) * p[37]
    +4 * x * z * (-1+y) * p[37]
    +4 * x * y * (-1+x) * p[40]
    -4 * x * y * (-1+z) * p[43]
    -4 * x * y * z * p[43]
    +4 * y * (-1+z) * (-1+x) * p[46]
    +4 * y * z * (-1+x) * p[46]
    +4 * x * (-1+y) * (-1+x) * p[49]
    +4 * y * (-1+y) * (-1+x) * p[52]
    -4 * x * y * (-1+y) * p[55]-4 * x * y * (-1+x) * p[58];

  dtau[2][0]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[2]
    +2 * x * (-1+z) * (-1+y) * p[5]
    -2 * x * y * (-1+z) * p[8]
    -2 * y * (-1+z) * (-1+x) * p[11]
    -4 * y * (-1+z) * (-1+y) * p[29]
    +(-1+z) * (-1+y) * (2 * x+2 * y+2 * z-1) * p[2]
    +(-1+z) * (-1+y) * (-2 * y-2 * z+2 * x-1) * p[5]
    -y * (-1+z) * (2 * y-2 * z-3+2 * x) * p[8]
    -y * (-1+z) * (2 * x+2 * z+1-2 * y) * p[11]
    -z * (-1+y) * (2 * x+2 * y-2 * z+1) * p[14]
    -z * (-1+y) * (-2 * y+2 * z-3+2 * x) * p[17]
    +y * z * (2 * y+2 * z-5+2 * x) * p[20]
    +y * z * (2 * x-2 * z+3-2 * y) * p[23]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[26]
    -4 * x * (-1+z) * (-1+y) * p[26]
    -2 * z * (-1+y) * (-1+x) * p[14]
    -2 * x * z * (-1+y) * p[17]
    +2 * x * y * z * p[20]
    +2 * y * z * (-1+x) * p[23]
    -4 * z * (-1+z) * (-1+y) * p[32]
    +4 * y * (-1+z) * (-1+y) * p[35]
    +4 * z * (-1+z) * (-1+y) * p[38]
    +4 * y * (-1+z) * (-1+x) * p[41]
    +4 * x * y * (-1+z) * p[41]
    -4 * y * z * (-1+z) * p[44]
    +4 * y * z * (-1+z) * p[47]
    +4 * z * (-1+y) * (-1+x) * p[50]
    +4 * x * z * (-1+y) * p[50]
    +4 * y * z * (-1+y) * p[53]
    -4 * y * z * (-1+y) * p[56]
    -4 * y * z * (-1+x) * p[59]
    -4 * x * y * z * p[59];

  dtau[2][1]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[2]
    -2 * x * (-1+z) * (-1+y) * p[5]
    -2 * x * y * (-1+z) * p[8]
    +2 * y * (-1+z) * (-1+x) * p[11]
    -2 * z * (-1+y) * (-1+x) * p[14]
    +2 * x * z * (-1+y) * p[17]
    +2 * x * y * z * p[20]
    -2 * y * z * (-1+x) * p[23]
    +(-1+z) * (-1+x) * (2 * x+2 * y+2 * z-1) * p[2]
    +x * (-1+z) * (-2 * y-2 * z+2 * x-1) * p[5]
    -x * (-1+z) * (2 * y-2 * z-3+2 * x) * p[8]
    -(-1+z) * (-1+x) * (2 * x+2 * z+1-2 * y) * p[11]
    -z * (-1+x) * (2 * x+2 * y-2 * z+1) * p[14]
    -x * z * (-2 * y+2 * z-3+2 * x) * p[17]
    +x * z * (2 * y+2 * z-5+2 * x) * p[20]
    +z * (-1+x) * (2 * x-2 * z+3-2 * y) * p[23]
    -4 * x * (-1+z) * (-1+x) * p[26]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[29]
    -4 * y * (-1+z) * (-1+x) * p[29]
    -4 * z * (-1+z) * (-1+x) * p[32]
    +4 * x * (-1+z) * (-1+y) * p[35]
    +4 * x * y * (-1+z) * p[35]
    +4 * x * z * (-1+z) * p[38]
    +4 * x * (-1+z) * (-1+x) * p[41]
    -4 * x * z * (-1+z) * p[44]
    +4 * z * (-1+z) * (-1+x) * p[47]
    +4 * x * z * (-1+x) * p[50]
    +4 * z * (-1+y) * (-1+x) * p[53]
    +4 * y * z * (-1+x) * p[53]
    -4 * x * z * (-1+y) * p[56]
    -4 * x * y * z * p[56]
    -4 * x * z * (-1+x) * p[59];

  dtau[2][2]
    =2 * (-1+z) * (-1+y) * (-1+x) * p[2]
    -2 * x * (-1+z) * (-1+y) * p[5]
    +2 * x * y * (-1+z) * p[8]
    -2 * y * (-1+z) * (-1+x) * p[11]
    +2 * z * (-1+y) * (-1+x) * p[14]
    -2 * x * z * (-1+y) * p[17]
    +2 * x * y * z * p[20]
    -2 * y * z * (-1+x) * p[23]
    +(-1+y) * (-1+x) * (2 * x+2 * y+2 * z-1) * p[2]
    +x * (-1+y) * (-2 * y-2 * z+2 * x-1) * p[5]
    -x * y * (2 * y-2 * z-3+2 * x) * p[8]
    -y * (-1+x) * (2 * x+2 * z+1-2 * y) * p[11]
    -(-1+y) * (-1+x) * (2 * x+2 * y-2 * z+1) * p[14]
    -x * (-1+y) * (-2 * y+2 * z-3+2 * x) * p[17]
    +x * y * (2 * y+2 * z-5+2 * x) * p[20]
    +y * (-1+x) * (2 * x-2 * z+3-2 * y) * p[23]
    -4 * x * (-1+y) * (-1+x) * p[26]
    -4 * y * (-1+y) * (-1+x) * p[29]
    -4 * (-1+z) * (-1+y) * (-1+x) * p[32]
    -4 * z * (-1+y) * (-1+x) * p[32]
    +4 * x * y * (-1+y) * p[35]
    +4 * x * (-1+z) * (-1+y) * p[38]
    +4 * x * z * (-1+y) * p[38]
    +4 * x * y * (-1+x) * p[41]
    -4 * x * y * (-1+z) * p[44]
    -4 * x * y * z * p[44]
    +4 * y * (-1+z) * (-1+x) * p[47]
    +4 * y * z * (-1+x) * p[47]
    +4 * x * (-1+y) * (-1+x) * p[50]
    +4 * y * (-1+y) * (-1+x) * p[53]
    -4 * x * y * (-1+y) * p[56]
    -4 * x * y * (-1+x) * p[59];
}

// Get the logical index of the gaussian point given the coordinate
// p[] of the point in the subcell and the index of the subcell icell.
inline int ipg(const int npg[], const int p[], const int icell) 
{
  return npg[0] * npg[1] * npg[2] * icell
    + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
}

// Compute the surface terms inside one macrocell
__kernel
void DGFlux(__constant int *param,     // interp param
	    int dim0,                  // face direction
	    __constant real *physnode, // macrocell nodes
	    __global   real *wn,       // field values
	    __global   real *dtwn,     // time derivative
	    __local    real *wnloc     // wn and dtwn in local memory
	    )
{
  // Use __local memory in DGFlux kernel?
#define DGFLUX_LOCAL 1

  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};

  int dim1 = (dim0 + 1) % 3;
  int dim2 = (dim1 + 1) % 3;

  // Subcell id
  int icL[3];
  int icR[3];
  int icell = get_group_id(0);

  icL[dim0] = icell % (nraf[dim0] - 1);
  icL[dim1] = (icell / (nraf[dim0] - 1)) % nraf[dim1];
  icL[dim2] = icell / (nraf[dim0]-1) / nraf[dim1];

  icR[dim0] = icL[dim0] + 1;
  icR[dim1] = icL[dim1];
  icR[dim2] = icL[dim2];

  __local real *wnlocL = wnloc;
  __local real *wnlocR = wnloc + get_local_size(0) * m;
  __local real *dtwnlocL = wnloc + 2 * get_local_size(0) * m;
  __local real *dtwnlocR = wnloc + 3 * get_local_size(0) * m;

  int pL[3];
  
#if DGFLUX_LOCAL
  for(int i = 0; i < m; i++) {
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipg = iread / m;
    
    int p[3];
    p[dim1] = ipg % npg[dim1];
    p[dim2] = ipg / npg[dim1];

    // Left point
    p[dim0] = deg[dim0];
    int ipgL = xyz_to_ipg(nraf, deg, icL, p);
    int imemL = VARINDEX(param, ipgL, iv);
    // wnlocL[iread] = wn[imemL];
    wnlocL[ipg * m + iv] = wn[imemL];

    // Right point
    p[dim0] = 0;
    int ipgR =  xyz_to_ipg(nraf, deg, icR, p);
    int imemR = VARINDEX(param, ipgR, iv);
    // wnlocR[iread] = wn[imemR];
    wnlocR[ipg * m + iv] = wn[imemR];
  }

  //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  barrier(CLK_LOCAL_MEM_FENCE);
#else
  // Gauss point id where we compute the jacobian

  int pR[3];
  
  //ipg_to_xyz(get_local_id(0), p, npg
  {
    int ipg = get_local_id(0);
    pL[dim0] = deg[dim0];
    pL[dim1] = ipg % npg[dim1];
    pL[dim2] = (ipg / npg[dim1]);

    pR[dim0] = 0;
    pR[dim1] = pL[dim1];
    pR[dim2] = pL[dim2];
  }
  
  real wL[_M], wR[_M];
  int ipgL = xyz_to_ipg(nraf, deg, icL, pL);
  int ipgR = xyz_to_ipg(nraf, deg, icR, pR);
  for(int iv = 0; iv < m; iv++) {
    int imemL = VARINDEX(param, ipgL, iv);
    wL[iv] = wn[imemL];
    int imemR = VARINDEX(param, ipgR, iv);
    wR[iv] = wn[imemR];
  }
#endif
 
  // Gauss point id where we compute the Jacobian
  //int pL[3], pR[3];
  //ipg_to_xyz(get_local_id(0), p, npg
  {
    int ipg = get_local_id(0);
    pL[dim0] = deg[dim0];
    pL[dim1] = ipg % npg[dim1];
    pL[dim2] = ipg / npg[dim1];
  }

  real h[3] = {1.0 / (real) nraf[0],
	       1.0 / (real) nraf[1],
	       1.0 / (real) nraf[2] };

  int offset[3] = {gauss_lob_offset[deg[0]] + pL[0],
  		   gauss_lob_offset[deg[1]] + pL[1],
  		   gauss_lob_offset[deg[2]] + pL[2]};

  real x = h[0] * (icL[0] + gauss_lob_point[offset[0]]);
  real y = h[1] * (icL[1] + gauss_lob_point[offset[1]]);
  real z = h[2] * (icL[2] + gauss_lob_point[offset[2]]);

  real codtau[3][3];
  {
    real dtau[3][3];
    get_dtau(x, y, z, physnode, dtau);
    compute_codtau(dtau, codtau);
  }

#if DGFLUX_LOCAL
  real wL[_M], wR[_M]; // TODO: remove?
  __local real *wnL = wnlocL + get_local_id(0) * m;
  __local real *wnR = wnlocR + get_local_id(0) * m;
  for(int iv = 0; iv < m; iv++) {
    wL[iv] = wnL[iv];
    wR[iv] = wnR[iv];
  }
#else

  ipgL = xyz_to_ipg(nraf, deg, icL, pL);
  ipgR = xyz_to_ipg(nraf, deg, icR, pR);

  for(int iv = 0; iv < m; iv++) {
    int imemL = VARINDEX(param, ipgL, iv);
    wL[iv] = wn[imemL];
    int imemR = VARINDEX(param, ipgR, iv);
    wR[iv] = wn[imemR];
  }
#endif
  
  const real h1h2 = 1.0 / nraf[dim1] / nraf[dim2];
  real vnds[3] = {codtau[0][dim0] * h1h2,
		  codtau[1][dim0] * h1h2,
		  codtau[2][dim0] * h1h2 };

  // TODO: wL and wR could be passed without a copy to __private.  (ie
  // we can just pass *wnL and *wnR).
  real flux[_M];
  NUMFLUX(wL, wR, vnds, flux);

  real wpgs = wglop(deg[dim1], pL[dim1]) * wglop(deg[dim2], pL[dim2]);

#if DGFLUX_LOCAL
    // write flux to local memory
  __local real *dtwnL = dtwnlocL + get_local_id(0) * m;
  __local real *dtwnR = dtwnlocR + get_local_id(0) * m;
  for(int iv = 0; iv < m; ++iv) {
    real fluxivwpgs = flux[iv] * wpgs; 
    dtwnL[iv] = -fluxivwpgs;
    dtwnR[iv] =  fluxivwpgs;
  }

  //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  barrier(CLK_LOCAL_MEM_FENCE);

  for(int i = 0; i < m; i++) {
    int p[3];
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipg = iread / m;
    
    p[dim1] = ipg % npg[dim1];
    p[dim2] = ipg / npg[dim1];

    // Left point
    p[dim0] = deg[dim0];
    int ipgL = xyz_to_ipg(nraf, deg, icL, p);
    int imemL = VARINDEX(param, ipgL, iv);
    // wnlocL[iread] = wn[imemL];
    dtwn[imemL] += dtwnlocL[ipg * m + iv];
    
    // Right point
    p[dim0] = 0;
    int ipgR =  xyz_to_ipg(nraf, deg, icR, p);
    int imemR = VARINDEX(param, ipgR, iv);
    // wnlocR[iread] = wn[imemR];
    dtwn[imemR] += dtwnlocR[ipg * m + iv];
  }
#else
  for(int iv = 0; iv < m; ++iv) {
    //int ipgL = ipg(npg, p, icell);
    //int imemL = VARINDEX(param, ie, ipgL, iv);
    real fluxivwpgs = flux[iv] * wpgs; 
    
    int imemL = VARINDEX(param, ipgL, iv);
    dtwn[imemL] -= fluxivwpgs;

    int imemR = VARINDEX(param, ipgR, iv);
    dtwn[imemR] += fluxivwpgs;
  }
#endif
}

inline void prefetch_macrocell(const __global real *in,
			__local real *out,
			__constant int *param)
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m ; ++i) {
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m;
    int ipgL = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipgL, iv);
    int imemloc = iv + ipgloc * m;
    
    out[imemloc] = in[imem];
  }
}

inline void zero_macrocell_buffer(__local real *out,
			   __constant int *param
			   )
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m ; ++i) {
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m;
    int ipgL = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipgL, iv);
    int imemloc = iv + ipgloc * m;
    
    out[imemloc] = 0.0;
  }
}

inline void postfetch_macrocell(const __local real *in,
			 __global real *out,
			 __constant int *param
			 )
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m; ++i){
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m ;
    int ipg = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipg, iv);
    int imemloc = ipgloc * m + iv;
    out[imem] += in[imemloc];
  }
}

inline void postfetch_macrocell_overwrite(const __local real *in,
			 __global real *out,
			 __constant int *param
			 )
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m; ++i){
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m ;
    int ipg = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipg, iv);
    int imemloc = ipgloc * m + iv;
    out[imem] = in[imemloc];
  }
}

inline void compute_volume(__constant int *param,     // interp param
		    __constant real *physnode, // macrocell nodes
		    __local real *wnloc,       // cache for wn
		    __local real *dtwnloc      // cache for dtwn
		    )
{
  const int m = param[0];
  const int deg[3] = {param[1],param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int raf[3] = {param[4], param[5], param[6]};
  
  // subcell id
  int icell = get_group_id(0);
  int ic[3] = {icell % raf[0],
		(icell / raf[0]) % raf[1],
		icell / raf[0] / raf[1] };

  // gauss point id where we compute the jacobian
  int localid = get_local_id(0);
  int ix[3] = {localid % npg[0],
	      (localid / npg[0]) % npg[1],
	      localid / npg[0] / npg[1] };

  real h[3] = {1.0 / raf[0], 1.0 / raf[1], 1.0 / raf[2]};

  int offset[3] = {gauss_lob_offset[deg[0]] + ix[0],
		   gauss_lob_offset[deg[1]] + ix[1],
		   gauss_lob_offset[deg[2]] + ix[2]};

  real x[3] = {h[0] * (ic[0] + gauss_lob_point[offset[0]]),
	       h[1] * (ic[1] + gauss_lob_point[offset[1]]),
	       h[2] * (ic[2] + gauss_lob_point[offset[2]]) };

  real wpg
    = h[0] * gauss_lob_weight[offset[0]]
    * h[1] * gauss_lob_weight[offset[1]]
    * h[2] * gauss_lob_weight[offset[2]];

  real codtau[3][3];
  {
    real dtau[3][3];
    get_dtau(x[0], x[1], x[2], physnode, dtau);
    compute_codtau(dtau, codtau);
  }

  real wL[_M];
  int ipgL = ipg(npg, ix, 0);
  __local real *wnloc0 = wnloc + ipgL * m;

  for(int iv = 0; iv < m; iv++) {
    wL[iv] = wnloc0[iv];
  }

  real flux[_M];
  for(int dim0 = 0; dim0 < 3; dim0++) {
    int iy[3] = {ix[0], ix[1], ix[2]};

    // Loop on the "cross" points
    for(int iq = 0; iq < npg[dim0]; iq++) {
      iy[dim0] = (ix[dim0] + iq) % npg[dim0];
      real dphiref[3] = {0, 0, 0};
      dphiref[dim0] = dlag(deg[dim0], iy[dim0], ix[dim0]) * raf[dim0];
      real dphi[3];
      compute_dphi(dphiref, codtau, dphi);
      
      NUMFLUX(wL, wL, dphi, flux);

      int ipgR = ipg(npg, iy, 0);

      int imemR0loc = ipgR * m;
      __local real *dtwnloc0 =  dtwnloc + imemR0loc;
      for(int iv = 0; iv < m; iv++) {
	// TODO: is it worth fma when it's += ?
#ifdef FP_FAST_FMA
	dtwnloc0[iv] = fma(flux[iv], wpg, dtwnloc0[iv]);
#else
	dtwnloc0[iv] += flux[iv] * wpg;
#endif
      }
    }

  } // dim0 loop
}

inline void compute_volume_global(__constant int *param,     // interp param
				  __constant real *physnode, // macrocell nodes
				  __global real *wn,         // field values
				  __global real *dtwn       // time derivative
				  )
{
  
  const int m = param[0];
  const int deg[3] = {param[1],param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};
  
  // subcell id
  int icell = get_group_id(0);
  int ic[3] = {icell % nraf[0],
		(icell / nraf[0]) % nraf[1],
		icell / nraf[0] / nraf[1] };

  // gauss point id where we compute the jacobian
  int p[3];
  {
    int ipg = get_local_id(0);
    p[0] = ipg % npg[0];
    p[1] = (ipg / npg[0]) % npg[1];
    p[2] = ipg / npg[0] / npg[1];
  }

  // ref coordinates
  real h[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  int offset[3] = {gauss_lob_offset[deg[0]] + p[0],
		   gauss_lob_offset[deg[1]] + p[1],
		   gauss_lob_offset[deg[2]] + p[2]};

  real x[3] = {h[0] * (ic[0] + gauss_lob_point[offset[0]]),
	       h[1] * (ic[1] + gauss_lob_point[offset[1]]),
	       h[2] * (ic[2] + gauss_lob_point[offset[2]]) };

  real wpg
    = h[0] * gauss_lob_weight[offset[0]]
    * h[1] * gauss_lob_weight[offset[1]]
    * h[2] * gauss_lob_weight[offset[2]];

  real codtau[3][3];
  {
    real dtau[3][3];
    get_dtau(x[0], x[1], x[2], physnode, dtau);
    compute_codtau(dtau, codtau);
  }

  real wL[_M];
  int ipgL = ipg(npg, p, icell);
  for(int iv = 0; iv < m; iv++) {
    int imemL = VARINDEX(param, ipgL, iv);
    wL[iv] = wn[imemL];
  }

  real flux[_M];
  for(int dim0 = 0; dim0 < 3; dim0++) {
    int q[3] = {p[0], p[1], p[2]};

    // Loop on the "cross" points
    for(int iq = 0; iq < npg[dim0]; iq++) {
      q[dim0] = (p[dim0] + iq) % npg[dim0];
      real dphiref[3] = {0, 0, 0};
      dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) * nraf[dim0];
      real dphi[3];
      compute_dphi(dphiref, codtau, dphi);
      NUMFLUX(wL, wL, dphi, flux);

      int ipgR = ipg(npg, q, icell);
      int imemR0 = VARINDEX(param, ipgR, 0);
      __global real *dtwn0 = dtwn + imemR0; 
      for(int iv = 0; iv < m; iv++) {
#ifdef FP_FAST_FMA
	dtwn0[iv] = fma(flux[iv], wpg, dtwn0[iv]);
#else
	dtwn0[iv] += flux[iv] * wpg;
#endif
      }
    }

  } // dim0 loop

}
  
// Compute the volume  terms inside  one macrocell
__kernel
void DGVolume(__constant int *param,     // interp param
	      __constant real *physnode, // macrocell nodes
              __global real *wn,         // field values
	      __global real *dtwn,       // time derivative
	      __local real *wnloc,       // cache for wn
	      __local real *dtwnloc      // cache for dtwn
	      )
{
  // Use __local memory in DGVolume kernel?
#define DGVolume_LOCAL 1

#if DGVolume_LOCAL

  prefetch_macrocell(wn, wnloc, param);
  zero_macrocell_buffer(dtwnloc, param);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  compute_volume(param, physnode, wnloc, dtwnloc);

  barrier(CLK_LOCAL_MEM_FENCE);

  postfetch_macrocell(dtwnloc, dtwn, param);

#else

  compute_volume_global(param, physnode, wn, dtwn);
  
#endif
}

inline void mass_division(__constant int *param,
			  const  __global real *mass,
			  __local real *dtwnloc)
{
  const int m = param[0];

  // The mass is __global
  int ipg = get_global_id(0);
  real overmassipg = 1.0 / mass[ipg];

  // dtwnloc is __local
  int ipgloc = get_local_id(0);
  __local real *dtwnloc0 =  dtwnloc + ipgloc * m;
  for(int iv = 0; iv < m; iv++) {
    dtwnloc0[iv] *= overmassipg;
  }
}

// Apply division by the mass matrix on one macrocell
__kernel
void DGMass(__constant int *param,      // interp param
            __constant real *physnode,  // macrocell nodes
	    const __global real *mass,  // macrocell masses
            __global real *dtwn,        // time derivative
	    __local real *dtwnloc       // cache for dtwn
	    )
{
  // TODO: if there's enough space, make mass __constant  

#if 0
  const int ipg = get_global_id(0);
  const int m = param[0];

  real overmass = 1.0 / mass[ipg];
  __global real *dtwn0 = dtwn + m * ipg;
  for(int iv = 0; iv < m; iv++) {
    dtwn0[iv] *= overmass;
  }

#else
  prefetch_macrocell(dtwn, dtwnloc, param);

  barrier(CLK_LOCAL_MEM_FENCE);
  
  mass_division(param, mass, dtwnloc);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  postfetch_macrocell_overwrite(dtwnloc, dtwn, param);
#endif
}

__kernel
void ExtractInterface(__constant int *param,   // interp param
		      const int ifa,
		      const __global real *wn, // volumic input
		      __global real *wface     // output
		      )
{
  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int raf[3] = {param[4], param[5], param[6]};
  
  const int ipgf = get_global_id(0);
  const int iv = get_global_id(1);

  int pic[3];
  int pix[3];
  
  ipgf_to_pxyz(deg, raf, ifa, ipgf, pic, pix);

  int ic[3];
  int ix[3];
  unpermute_indices(ic, pic, ifa);
  unpermute_indices(ix, pix, ifa);
  
  // FIXME: do something with this.

  int vpos = xyz_to_ipg(raf, deg, ic, ix) * m + iv;
  int fpos = ipgf * m + iv;
  
  wface[fpos] = wn[vpos];
}

__kernel
void InsertInterface(__constant int *param,   // interp param
		     const int ifa,
		     __global real *dtwn, // volumic input
		     const __global real *wface     // output
		     )
{
  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int raf[3] = {param[4], param[5], param[6]};
  
  const int ipgf = get_global_id(0);
  const int iv = get_global_id(1);

  int pic[3];
  int pix[3];
  
  ipgf_to_pxyz(deg, raf, ifa, ipgf, pic, pix);

  int ic[3];
  int ix[3];
  unpermute_indices(ic, pic, ifa);
  unpermute_indices(ix, pix, ifa);
  
  // FIXME: do something with this.

  int vpos = xyz_to_ipg(raf, deg, ic, ix) * m + iv;
  int fpos = ipgf * m + iv;
  
  dtwn[vpos] = wface[fpos];
}

__kernel
void ExtractedDGInterfaceFlux(__constant int *param,
			      int Rcorner,
			      __global real *faceL,
			      __global real *faceR)
{

  const int m = param[0];
  const int ipgf = get_global_id(0);
  
  for(int iv = 0; iv < m; ++iv) {
    faceL[ipgf * m + iv] = 0.0;
  }

  for(int iv = 0; iv < m; ++iv) {
    faceR[ipgf * m + iv] = 0.0;
  }

  // The kernel is launched with ND range given by {the number of
  // points in the first L facial direction, the number of points in
  // the second L facial direciton}

  // FIXME

  /* const int m = param[0]; */
  
  /* // Copy Intrafces into __local memory */
  
  /* for(int i = 0; i < m; ++i) { */
  /*   int pos = get_local_id(1) */
  /*     + get_local_id(0)* get_local_size(1) */
  /*     + i * get_local_size(1) * get_local_size(1); */
  /*   fL[pos] = faceL[pos]; */
  /*   fR[pos] = faceR[pos]; */
  /* } */
   
  /* barrier(CLK_LOCAL_MEM_FENCE); */

  /* // Compute flux */

  /* // The point with 2D face indices (d0, d1) is at memory location */
  /* // d1 * n0 * m  + d0 * m */
  /* const int d0 = get_global_id(0); // first dimension */
  /* const int d1 = get_global_id(1); // second dimension */
  /* const int n0 = get_global_size(0); */

  /* barrier(CLK_LOCAL_MEM_FENCE); */
  
  /* // Copy flux back into input buffers */
  /* for(int i = 0; i < m; ++i) { */
  /*   int pos = get_local_id(1) */
  /*     + get_local_id(0) * get_local_size(1) */
  /*     + i * get_local_size(1) * get_local_size(1); */
  /*   faceL[pos] = fL[pos]; */
  /*   faceR[pos] = fR[pos]; */
  /* } */

}

__kernel
void ExtractedDGBoundaryFlux(__constant int *param,
			     __constant real *physnode,
			     __global real *faceL,
			     int locfa,
			     real tnow)
{
  // FIXME
  const int raf[3] = {param[4], param[5], param[6]};
  const int deg[3] = {param[1], param[2], param[3]};
  const int m = param[0];
  
  int ipgf = get_global_id(0);

  for(int iv = 0; iv < m; ++iv) {
    faceL[ipgf * m + iv] = 0.0;
  }

  real xref[3];    // reference coordinates
  real wpg;        // weighting for the GL point
  int ipgL = ref_pg_face(deg, raf, locfa, ipgf, xref, &wpg);

  // Normal vector
  real vnds[3];
  // Physical coordinates
  real xphy[3];
  {
    real gradphi[20][4];
    compute_gradphi(xref, gradphi);
    
    compute_xphy(physnode, gradphi, xphy);

    real dtau[3][3];  
    compute_dtau(physnode, gradphi, dtau);

    real codtau[3][3];
    compute_codtau(dtau, codtau);

    ComputeNormal(codtau, locfa, vnds);
  }
  
  real wL[_M];
  
  __global real *wn0 = faceL + ipgf * m;;
  for(int iv = 0; iv < m; ++iv) {
    wL[iv] = wn0[iv];
  }

#ifndef BOUNDARYFLUX
#define BOUNDARYFLUX BoundaryFlux
#endif

  real flux[_M];
  BOUNDARYFLUX(xphy, tnow, wL, vnds, flux);

  for(int iv = 0; iv < m; ++iv) {
    wn0[iv] = flux[iv];
  }
  
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
__kernel
void DGMacroCellInterface(__constant int *param,        // interp param
                          int locfaL,                   // left face index
			  int locfaR,                   // right face index
                          __constant real *physnodeL,   // left physnode
			  __constant real *physnodeR,   // right physnode
                          __global real *wnL,           // field 
                          __global real *dtwnL,         // time derivative
                          __global real *wnR,           // field 
                          __global real *dtwnR          // time derivative
			  )
{
  // Index of the point on the face.
  int ipgfL = get_global_id(0);

  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int raf[3] = {param[4], param[5], param[6]};

  real xpgref[3]; // reference point for L
  real xpgref_in[3]; // reference point slightly in R
  real wpg;

  // Compute volumic index in the left MacroCell.
  int ipgL = ref_pg_face(deg, raf, locfaL, ipgfL, xpgref, &wpg);

  int icL[3];
  int ixL[3];
  ipg_to_xyz(raf, deg, icL, ixL, &ipgL);

  compute_xpgin(raf, deg, icL, ixL, locfaL, xpgref_in);
  
  // Normal vector at gauss point based on xpgref
  real vnds[3];
  {
    real gradphi[20][4];
    compute_gradphi(xpgref, gradphi);

    real dtau[3][3];
    compute_dtau(physnodeL, gradphi, dtau);

    real codtau[3][3];
    compute_codtau(dtau, codtau);

    ComputeNormal(codtau, locfaL, vnds);
  }

  // Find the volumic index for the opposite side:
  int ipgR;
  {
    real gradphi[20][4];
    compute_gradphi(xpgref_in, gradphi);

    real xphy[3];
    compute_xphy(physnodeL, gradphi, xphy);
    
#if defined(_PERIODX) || defined(_PERIODY) || defined(_PERIODZ)
#ifndef _PERIODX
#define _PERIODX -1
#endif
#ifndef _PERIODY
#define _PERIODY -1
#endif
#ifndef _PERIODZ
#define _PERIODZ -1
#endif
    real period[3] = {_PERIODX, _PERIODY, _PERIODZ};
    PeriodicCorrection(xphy, period);
#endif

    real xrefR[3];
    Phy2Ref(physnodeR, xphy, xrefR);
    ipgR = ref_ipg(raf, deg, xrefR);
  }
  
  real wL[_M];
  real wR[_M];

  int imemL0 = VARINDEX(param, ipgL, 0);
  int imemR0 = VARINDEX(param, ipgR, 0);

  __global real *wnL0 = wnL + imemL0;
  __global real *wnR0 = wnR + imemR0;
  for(int iv = 0; iv < m; iv++) {
    wL[iv] = wnL0[iv];
    wR[iv] = wnR0[iv];
  }

  real flux[_M];
  NUMFLUX(wL, wR, vnds, flux);

  __global real *dtwnL0 = dtwnL + imemL0;
  __global real *dtwnR0 = dtwnR + imemR0;
  for(int iv = 0; iv < m; ++iv) {
    real fluxwpg = flux[iv] * wpg;
    dtwnL0[iv] -= fluxwpg;
    dtwnR0[iv] += fluxwpg;
  }
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
__kernel
void DGBoundary(__constant int *param,      // interp param
		real tnow,                  // current time
		int locfa,                 // left face index
		__constant real *physnode, // geometry for all mcells
		__global real *wn,          // field 
		__global real *dtwn        // time derivative
		)
{
  // TODO: use __local real *cache.

  int ipgf = get_global_id(0);

  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int raf[3] = {param[4], param[5], param[6]};

  real xref[3];    // reference coordinates
  real wpg;        // weighting for the GL point
  int ipgL = ref_pg_face(deg, raf, locfa, ipgf, xref, &wpg);

  // Normal vector
  real vnds[3];
  // Physical coordinates
  real xphy[3];
  {
    real gradphi[20][4];
    compute_gradphi(xref, gradphi);
    
    compute_xphy(physnode, gradphi, xphy);

    real dtau[3][3];  
    compute_dtau(physnode, gradphi, dtau);

    real codtau[3][3];
    compute_codtau(dtau, codtau);

    ComputeNormal(codtau, locfa, vnds);
  }

  real wL[_M];
  
  int imemL0 = VARINDEX(param, ipgL, 0);
  __global real *wn0 = wn + imemL0;
  for(int iv = 0; iv < m; ++iv) {
    wL[iv] = wn0[iv];
  }

#ifndef BOUNDARYFLUX
#define BOUNDARYFLUX BoundaryFlux
#endif

  real flux[_M];
  BOUNDARYFLUX(xphy, tnow, wL, vnds, flux);

  __global real *dtwn0 = dtwn + imemL0; 
  for(int iv = 0; iv < m; ++iv) {
#ifdef FP_FAST_FMA
    dtwn0[iv] = fma(-flux[iv], wpg, dtwn0[iv]);
#else
    dtwn0[iv] -= flux[iv] * wpg;
#endif
  }
}

// Sample source function
void ZeroSource(const real *x, const real t, const real *w, real *source,
		int m)
{
  for(int i = 0; i < m; ++i) 
    source[i] = 0.0;
}

// Sample source function
void OneSource(const real *x, const real t, const real *w, real *source,
	       int m)
{
  for(int i = 0; i < m; ++i) 
    source[i] = 1.0;
}

void compute_source(__constant int *param,     // interp param
		    __constant real *physnode, // macrocell nodes
		    const  __global real *mass,   // collocation point weights
		    const real tnow,           // the current time
		    __local real *wnloc,       // cache for wn
		    __local real *dtwnloc      // cache for dtwn
		    )
{
  // TODO: if there's enough space, make mass __constant

  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int nraf[3] = {param[4], param[5], param[6]};

  int ipg = get_global_id(0);
  int ipgloc = get_local_id(0);
  
  // Compute xref
  real xref[3];
  real wpg;
  ref_pg_vol(deg, nraf, ipg, xref, &wpg);

  // Compute xphy
  real xphy[3];
  real gradphi[20][4];
  compute_gradphi(xref, gradphi);
  compute_xphy(physnode, gradphi, xphy);
  
  // Copy w
  real w[_M];
  __local real *wnloc0 = wnloc + ipgloc * m;
  for(int iv = 0; iv < m; iv++) {
    w[iv] = wnloc0[iv];
  }

#ifndef _SOURCE_FUNC
#define _SOURCE_FUNC ZeroSource
#endif
  
  // Compute source using w and xref
  real source[_M];
  _SOURCE_FUNC(xphy, tnow, w, source, _M);

  // The mass point is in __global memory, so get the global id.
  real massipg = mass[ipg];

  // Add the source buffer to dtw
  __local real *dtwnloc0 =  dtwnloc + ipgloc * m;;
  for(int iv = 0; iv < m; iv++) {
    dtwnloc0[iv] = source[iv] * massipg;
  }
}

// Compute the source terms inside  one macrocell
__kernel
void DGSource(__constant int *param,     // interp param
	      __constant real *physnode, // macrocell nodes
	      const __global real *mass, // collocation point weights
	      const real tnow,           // the current time
              __global real *wn,         // field values
	      __global real *dtwn,       // time derivative
	      __local real *wnloc,       // cache for wn
	      __local real *dtwnloc      // cache for dtwn
	      )
{
  prefetch_macrocell(wn, wnloc, param);
  prefetch_macrocell(dtwn, dtwnloc, param);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  compute_source(param, physnode, mass, tnow, wnloc, dtwnloc);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  postfetch_macrocell(dtwnloc, dtwn, param);
}

// Out-of-place RK stage
__kernel
void RK_out_CL(__global real *wnp1, 
	       __global const real *wn, 
	       __global const real *dtwn, 
	       const real dt)
{
  int ipg = get_global_id(0);

#ifdef FP_FAST_FMA
  wnp1[ipg] = fma(dt, dtwn[ipg], wn[ipg]);
#else
  wnp1[ipg] = wn[ipg] + dt * dtwn[ipg];
#endif
}

// In-place RK stage
__kernel
void RK_in_CL(__global real *wnp1, 
	      __global real *dtwn, 
	      const real dt)
{
  int ipg = get_global_id(0);
#ifdef FP_FAST_FMA
  wnp1[ipg] = fma(dt, dtwn[ipg], wnp1[ipg]);
#else
  wnp1[ipg] += dt * dtwn[ipg];
  #endif
}

// Out-of-place RK stage
__kernel
void RK4_first_stages(__global real *wnp1, 
		     __global const real *wn, 
		     __global const real *dtwn, 
		      const real dt)
{
  int ipg = get_global_id(0);

#ifdef FP_FAST_FMA
  wnp1[ipg] = fma(dt, dtwn[ipg], wn[ipg]);
#else
  wnp1[ipg] = wn[ipg] + dt * dtwn[ipg];
#endif
}

// RK4 final stage
__kernel
void RK4_final_stage(__global real *w,
		     __global real *l1,
		     __global real *l2,
		     __global real *l3,
		     __global real *dtw, 
		     const real dt)
{
  const real b = -1.0 / 3.0;
  const real a0 = 1.0 / 3.0;
  const real a1 = 2.0 / 3.0;
  const real a2 = 1.0 / 3.0;
  const real a3 = dt / 6.0;
  int i = get_global_id(0);

#ifdef FP_FAST_FMA
  w[i] = fma(b, w[i],
	     fma(a0, l1[i],
		 fma(a1, l2[i],
		     fma(a2, l3[i],
			 a3 * dtw[i])
		     )
		 )
	     );
#else
  w[i] 
    = b * w[i]
    + a0 * l1[i]
    + a1 * l2[i]
    + a2 * l3[i]
    + a3 * dtw[i];
#endif
}

// -*- mode: c; -*-

#pragma OPENCL EXTENSION cl_khr_fp64: enable

//#define _FPREC
// for single precision computations
//#define _FPREC f

#define NULL 0

// For each face, give the dimension index i
__constant int axis_permut[6][4] = {
  {0, 2, 1, 0},
  {1, 2, 0, 1},
  {2, 0, 1, 1},
  {2, 1, 0, 0},
  {0, 1, 2, 1},
  {1, 0, 2, 0}
};

__constant int h20_refnormal[6][3]={{0,-1,0},
				    {1,0,0},
				    {0,1,0},
				    {-1,0,0},
				    {0,0,1},
				    {0,0,-1}};

double dlag(int deg, int ib, int ipg);

// Return the 1d derivative of lagrange polynomial ib at glop ipg
double dlag(int deg, int ib, int ipg) {
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg] + ib * (deg + 1) + ipg];
}

int varindex(__constant int* param, int elem, int ipg, int iv);

// Memory location of w : component iv, macrocell elem and gauss point
// id in the macrocell ipg
int varindex(__constant int *param, int elem, int ipg, int iv) {
  int npg
    = (param[1] + 1) * (param[2] + 1) * (param[3] + 1)
    * param[4] * param[5] * param[6];
  return iv + param[0] * ( ipg + npg * elem);
}

int ref_ipg(__constant int* param, double* xref);

void Ref2Phy(__constant double* physnode,
             double xref[3],
             double dphiref[3],
             int ifa,
             double xphy[3],
             double dtau[3][3],
             double codtau[3][3],
             double dphi[3],
             double vnds[3]);

void Phy2Ref(__constant double* physnode,
             double xphy[3], double xref[3]);

int ref_pg_face(__constant int* param, int ifa, int ipg,
                double* xpg, double* wpg, double* xpgin);

int ref_pg_face(__constant int* param, int ifa, int ipg,
                double* xpg, double* wpg, double* xpgin) {
  //int ipgf=ipg;

  // approximation degree in each direction
  int deg[3] = {param[axis_permut[ifa][0]],
		param[axis_permut[ifa][1]],
		param[axis_permut[ifa][2]]};

  // number of subcells in each direction
  int nraf[3] = {param[3 + axis_permut[ifa][0]],
		 param[3 + axis_permut[ifa][1]],
		 param[3 + axis_permut[ifa][2]]};

  // Compute permuted indices
  int ix = ipg % (deg[0] + 1);
  ipg /= (deg[0] + 1);

  int iy = ipg % (deg[1] + 1);
  ipg /= (deg[1] + 1);

  // Equals 0 or d depending on the face
  int iz = axis_permut[ifa][3] * deg[2];

  // Compute permuted indices of the subface
  int ncx = ipg % nraf[0];

  double h[3];
  h[0] = 1.0 / (double) nraf[0];
  ipg /= nraf[0];

  int ncy = ipg;
  h[1] = 1.0 / (double) nraf[1];

  // Equals 0 or nraf-1 depending on the face
  int ncz = axis_permut[ifa][3] * (nraf[2] - 1);
  h[2] = 1.0 / (double) nraf[2];

  // Compute non permuted indices for points and subfaces
  int ipgxyz[3];
  ipgxyz[axis_permut[ifa][0]] = ix;
  ipgxyz[axis_permut[ifa][1]] = iy;
  ipgxyz[axis_permut[ifa][2]] = iz;

  int ncpgxyz[3];
  ncpgxyz[axis_permut[ifa][0]] = ncx;
  ncpgxyz[axis_permut[ifa][1]] = ncy;
  ncpgxyz[axis_permut[ifa][2]] = ncz;

  // Compute the global index of the Gauss-Lobatto point in the volume
  int ipgv
    = ipgxyz[0]
    + (param[0] + 1)
    * (ipgxyz[1] + (param[1] + 1)
       * (ipgxyz[2] + (param[2] + 1)
	  * (ncpgxyz[0] + param[3]
	     * (ncpgxyz[1] + param[4]
		* ncpgxyz[2])
	     )
	  )
       );

  // Compute the reference coordinates of the Gauss-Lobatto point in
  // the volume
  int offset[2] = {gauss_lob_offset[deg[0]] + ix,
		   gauss_lob_offset[deg[1]] + iy};
  //printf("offset=%d\n",offset);

  xpg[axis_permut[ifa][0]] = h[0] * (ncx + gauss_lob_point[offset[0]]);
  xpg[axis_permut[ifa][1]] = h[1] * (ncy + gauss_lob_point[offset[1]]);
  xpg[axis_permut[ifa][2]] = axis_permut[ifa][3];

  *wpg = h[0] * h[1] *
    gauss_lob_weight[offset[0]] * gauss_lob_weight[offset[1]];

  // If xpgin exists, compute a point slightly INSIDE the opposite
  // subcell along the face.
  if(xpgin != NULL) {
    double small = 1e-3;//0.001
    double vsmall = 1e-6;//0.000001;

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

  return ipgv;
};


void NumFlux(double wL[], double wR[], double* vnorm, double* flux);

void NumFlux(double wL[], double wR[], double* vnorm, double* flux) {
  double s2 = 0.707106781186547524400844362105;
  double vn = s2 * (vnorm[0] + vnorm[1]);

  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;

  flux[0] = vnp * wL[0] + vnm * wR[0];
};

void BoundaryFlux(double x[3], double t, double wL[], double* vnorm,
                  double* flux);

void BoundaryFlux(double x[3], double t, double wL[], double* vnorm,
                  double* flux) {
  double wR[1];
  double s2 = 0.707106781186547524400844362105;
  double vx = s2 * (x[0] + x[1]);
  wR[0] = cos(vx - t);

  NumFlux(wL, wR, vnorm, flux);
};

//! \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
double wglop(int deg, int i);

double wglop(int deg, int i) {
  return gauss_lob_weight[gauss_lob_offset[deg] + i];
}

void get_dtau(double x, double y, double z,
	      __constant double* physnode, double dtau[][3]);

// Get the logical index of the gaussian point given the coordinate
// p[] of the point in the subcell and the index of the subcell icell.
int ipg(const int npg[], const int p[], const int icell) {
  return npg[0] * npg[1] * npg[2] * icell
    + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
}

// Compute the volume and subcell-interface terms on one macrocell
__kernel
void DGVolume(__constant int* param,        // interp param
	      int ie,            // macrocel index
	      __constant double* physnode,  // macrocell nodes
              __global double* wn,       // field values
	      __global double* dtwn) {       // time derivative

  const int m = param[0];
  const int deg[3] = {param[1],param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};

  // subcell id
  int icL[3];
  int icell = get_group_id(0);
  icL[0] = icell % nraf[0];
  icL[1] = (icell / nraf[0]) % nraf[1];
  icL[2]= icell / nraf[0] / nraf[1];

  // gauss point id where we compute the jacobian
  int p[3];
  {
    int ipg = get_local_id(0);
    p[0] = ipg % npg[0];
    p[1] = (ipg / npg[0]) % npg[1];
    p[2] = ipg / npg[0] / npg[1];
  }

  // ref coordinates
  double hx = 1. / (double) nraf[0];
  double hy = 1. / (double) nraf[1];
  double hz = 1. / (double) nraf[2];

  int offset[3] = {gauss_lob_offset[deg[0]] + p[0],
		   gauss_lob_offset[deg[1]] + p[1],
		   gauss_lob_offset[deg[2]] + p[2]};

  double x = hx * (icL[0] + gauss_lob_point[offset[0]]);
  double y = hy * (icL[1] + gauss_lob_point[offset[1]]);
  double z = hz * (icL[2] + gauss_lob_point[offset[2]]);

  double wpg = hx * hy * hz
    * gauss_lob_weight[offset[0]]
    * gauss_lob_weight[offset[1]]
    * gauss_lob_weight[offset[2]];

  double dtau[3][3];
  get_dtau(x, y, z, physnode, dtau);

  double codtau[3][3];
  codtau[0][0] =  dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
  codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
  codtau[0][2] =  dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
  codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
  codtau[1][1] =  dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
  codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
  codtau[2][0] =  dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
  codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
  codtau[2][2] =  dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];

  /* double det=dtau[0][0]*dtau[1][1]*dtau[2][2]-dtau[0][0]*dtau[1][2]*dtau[2][1]-dtau[1][0]*dtau[0][1]*dtau[2][2]+ */
  /*   dtau[1][0]*dtau[0][2]*dtau[2][1]+dtau[2][0]*dtau[0][1]*dtau[1][2]-dtau[2][0]*dtau[0][2]*dtau[1][1]; */

#define _M 1  /// TODO let schnaps specify m !!!!!!!!

  double wL[_M];
  for(int iv = 0; iv < m; iv++) {
    // gauss point id in the macrocell
    int ipgL = ipg(npg, p, icell);
    int imemL = varindex(param, ie, ipgL, iv);
    //int imemL= iv + m * ( get_global_id(0) + nnpg * *ie);
    wL[iv] = wn[imemL];
  }

  for(int dim0 = 0; dim0 < 3; dim0++) {
    int q[3] = {p[0], p[1], p[2]};

    // Loop on the "cross" points
    for(int iq = 0; iq < npg[dim0]; iq++) {
      q[dim0] = (p[dim0] + iq) % npg[dim0];
      double dphiref[3] ={0, 0, 0};
      dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) * nraf[dim0];
      double dphi[3];
      for(int ii = 0; ii < 3; ii++) {
	dphi[ii] = 0;
	for(int jj = 0; jj < 3; jj++) {
	  dphi[ii] += codtau[ii][jj] * dphiref[jj];
	}
      }

      double flux[_M];
      NumFlux(wL, wL, dphi, flux); // TODO: let schnaps gives fluxnum

      int ipgR = ipg(npg, q, icell);
      for(int iv=0; iv < m; iv++) {
	int imemR = varindex(param, ie, ipgR, iv); // TODO !
    	dtwn[imemR] += flux[iv] * wpg;
      }
    }

    // Compute the inter-subcell fluxes if needed
    if(p[dim0] == 0 || p[dim0] == npg[dim0] - 1) {
      int sgn = (p[dim0] > 0) ? 1 : -1;
      int dim1 = (dim0 + 1) % 3;
      int dim2 = (dim1 + 1) % 3;

      // Logical coordinates of the right subcell
      int icR[3] = {icL[0], icL[1], icL[2]};
      icR[dim0] += sgn;

      // if we are not at the boundary of the macrocell
      if(icR[dim0] >= 0 && icR[dim0] < nraf[dim0]) {
	double vnds[3];
	double h1h2 = 1. / nraf[dim1] / nraf[dim2];
	vnds[0] = sgn * codtau[0][dim0] * h1h2;
	vnds[1] = sgn * codtau[1][dim0] * h1h2;
	vnds[2] = sgn * codtau[2][dim0] * h1h2;

        int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
        int q[3] = {p[0], p[1], p[2]};
        q[dim0] = (sgn == -1) ? npg[dim0] - 1 : 0;
	int ipgR = ipg(npg, q, ncR);

	double wR[_M];
        for(int iv = 0; iv < m; iv++) {
          int imemR = varindex(param, ie, ipgR, iv);
          wR[iv] = wn[imemR];
        }

	double wpgs = wglop(deg[dim1], p[dim1]) * wglop(deg[dim2], p[dim2]);
        double flux[_M];
        NumFlux(wL, wR, vnds, flux); // TODO: let schnaps gives fluxnum
        for(int iv = 0; iv < m; iv++) {
          int ipgL = ipg(npg, p, icell);
          int imemL = varindex(param, ie, ipgL, iv);
          dtwn[imemL] -= flux[iv] * wpgs;
        }
      }
    }

  } // dim0 loop

}

// Apply division by the mass matrix on one macrocell
__kernel
void DGMass(__constant int *param,        // interp param
            int ie,            // macrocel index
            __constant double *physnode,  // macrocell nodes
            __global double *dtwn)        // time derivative
{

  int ipg = get_global_id(0);
  const int m = param[0];
  const int npg[3] = {param[1] + 1, param[2] + 1, param[3] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};

  const int npgie = npg[0] * npg[1] * npg[2] * nraf[0] * nraf[1] * nraf[2];

  //ref_pg_vol(param+1, ipg,xpgref,&wpg,NULL);
  int ix = ipg % npg[0];
  ipg /= npg[0];
  int iy = ipg % npg[1];
  ipg /= npg[1];
  int iz = ipg % npg[2];
  ipg /= npg[2];

  int ncx = ipg % nraf[0];
  ipg /= nraf[0];
  int ncy = ipg % nraf[1];
  ipg /= nraf[1];
  int ncz = ipg;

  double hx = 1.0 / (double) nraf[0];
  double hy = 1.0 / (double) nraf[1];
  double hz = 1.0 / (double) nraf[2];

  int offset[3] = {gauss_lob_offset[param[1]] + ix,
		   gauss_lob_offset[param[2]] + iy,
		   gauss_lob_offset[param[3]] + iz};

  double x = hx * (ncx + gauss_lob_point[offset[0]]);
  double y = hy * (ncy + gauss_lob_point[offset[1]]);
  double z = hz * (ncz + gauss_lob_point[offset[2]]);

  double wpg = hx * hy * hz
    * gauss_lob_weight[offset[0]]
    * gauss_lob_weight[offset[1]]
    * gauss_lob_weight[offset[2]];

  // end of ref_pg_vol
  //////////////////////////////////////////////

  //Ref2Phy(physnode, // phys. nodes
  //        xpgref,  // xref
  //        NULL,-1, // dpsiref, ifa
  //        NULL, dtau,  // xphy, dtau
  //        codtau,NULL,NULL); // codtau, dpsi,vnds

  double dtau[3][3];
  get_dtau(x, y, z, physnode, dtau);

  //codtau[0][0] = dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
  //codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
  //codtau[0][2] = dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
  //codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
  //codtau[1][1] = dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
  //codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
  //codtau[2][0] = dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
  //codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
  //codtau[2][2] = dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
  //double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
  // dtau[0][2]*codtau[0][2];

  // end of Ref2Phy
  //////////////////////////////////////////////////////////

  double det
    = dtau[0][0] * dtau[1][1] * dtau[2][2]
    - dtau[0][0] * dtau[1][2] * dtau[2][1]
    - dtau[1][0] * dtau[0][1] * dtau[2][2]
    + dtau[1][0] * dtau[0][2] * dtau[2][1]
    + dtau[2][0] * dtau[0][1] * dtau[1][2]
    - dtau[2][0] * dtau[0][2] * dtau[1][1];

  for(int iv = 0; iv < m; iv++) {
    // varindex
    int imem = iv + m * (get_global_id(0) + npgie * ie);
    // end of varindex
    /////////////////////////////////////
    //printf("imem=%d dtw=%f\n", imem, dtwn[imem]);
    //printf("det=%f wpg=%f imem=%d h=%f %f %f\n", det,wpg, imem,hx,hy,hz);

    dtwn[imem] /= (wpg * det);
  }
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
__kernel
void DGMacroCellInterface(__constant int *param,        // interp param
                          double tnow,  // current time
                          int ieL, int ieR,  // left and right elem ids
                          int locfaL, int locfaR, // current face local indices
                          __constant double *physnodeL, // left macrocell nodes
                          __constant double *physnodeR, // right macrocell nodes
                          __global double *wn,
                          __global double *dtwn) // time derivative
{
  int ipgfL=get_global_id(0);

  double xpgref[3], xpgref_in[3], wpg;
  // Get the coordinates of the Gauss point and coordinates of a
  // point slightly inside the opposite element in xref_in
  int ipgL = ref_pg_face(param + 1, locfaL, ipgfL,
			 xpgref, &wpg, xpgref_in);

  // Normal vector at gauss point ipg
  double vnds[3], xpg[3];
  {
    double dtau[3][3], codtau[3][3];
    Ref2Phy(physnodeL,
            xpgref,
            NULL, locfaL, // dpsiref, ifa
            xpg, dtau,
            codtau, NULL, vnds); // codtau, dpsi,vnds
  }

  double wL[_M];
  
  if (ieR >= 0) {  // The right element exists
    double xrefL[3];
    {
      double xpg_in[3];
      Ref2Phy(physnodeL,
              xpgref_in,
              NULL, -1, // dpsiref, ifa
              xpg_in, NULL,
              NULL, NULL, NULL); // codtau, dpsi,vnds
      Phy2Ref(physnodeR, xpg_in, xrefL);
    }

    int ipgR = ref_ipg(param + 1, xrefL);

    // Test code
    /* { */
    /*   double xpgR[3], xrefR[3], wpgR; */
    /*   ref_pg_vol(param + 1, ipgR, xrefR, &wpgR, NULL); */
    /*   Ref2Phy(physnodeR, */
    /* 	  xrefR, */
    /* 	  NULL, -1, // dphiref, ifa */
    /* 	  xpgR, NULL,   */
    /* 	  NULL, NULL, NULL); // codtau, dphi,vnds */
    /*   assert(Dist(xpgR, xpg) < 1e-10); */
    /* }	 */

    double wR[_M];
    for(int iv = 0; iv < _M; iv++) {
      int imemL = varindex(param, ieL, ipgL, iv);
      wL[iv] = wn[imemL];
      int imemR = varindex(param, ieR, ipgR, iv);
      wR[iv] = wn[imemR];
    }

    // int_dL F(wL, wR, grad phi_ib)
    double flux[_M];
    NumFlux(wL, wR, vnds, flux);

    // Add flux to both sides
    for(int iv = 0; iv < _M; iv++) {
      // The basis functions is also the gauss point index
      int imemL = varindex(param, ieL, ipgL, iv);
      int imemR = varindex(param, ieR, ipgR, iv);
      dtwn[imemL] -= flux[iv] * wpg;
      dtwn[imemR] += flux[iv] * wpg;
    }

  } else { // The point is on the boundary.
    for(int iv = 0; iv < _M; iv++) {
      int imemL = varindex(param, ieL, ipgL, iv);
      wL[iv] = wn[imemL];
    }

    double flux[_M];
    BoundaryFlux(xpg, tnow, wL, vnds, flux);

    //printf("ipgfL=%d ipgL=%d tnow=%f wL=%f flux=%f\n", ipgfL, ipgL,tnow,wL[0],flux[0]);

    for(int iv = 0; iv < _M; iv++) {
      // The basis functions is also the gauss point index
      int imemL = varindex(param, ieL, ipgL, iv);
      dtwn[imemL] -= flux[iv] * wpg;
    }
  }
}

void get_dtau(double x, double y, double z,
	      __constant double *p, double dtau[][3]) {

  // gradient of the shape functions and value (4th component)
  // of the shape functions
  /* double gradphi[20][3]; */
  /* //double x,y,z; */
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

  dtau[0][0]=2*(-1+z)*(-1+y)*(-1+x)*p[0]+2*x*(-1+z)*(-1+y)*p[3]-2*x*y*(-1+z)*p[6]-2*y*(-1+z)*(-1+x)*p[9]-2*z*(-1+y)*(-1+x)*p[12]-2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]+2*y*z*(-1+x)*p[21]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[0]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[3]-y*(-1+z)*(2*y-2*z-3+2*x)*p[6]-y*(-1+z)*(2*x+2*z+1-2*y)*p[9]-z*(-1+y)*(2*x+2*y-2*z+1)*p[12]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[15]+y*z*(2*y+2*z-5+2*x)*p[18]+y*z*(2*x-2*z+3-2*y)*p[21]-4*(-1+z)*(-1+y)*(-1+x)*p[24]-4*x*(-1+z)*(-1+y)*p[24]-4*y*(-1+z)*(-1+y)*p[27]-4*z*(-1+z)*(-1+y)*p[30]+4*y*(-1+z)*(-1+y)*p[33]+4*z*(-1+z)*(-1+y)*p[36]+4*y*(-1+z)*(-1+x)*p[39]+4*x*y*(-1+z)*p[39]-4*y*z*(-1+z)*p[42]+4*y*z*(-1+z)*p[45]+4*z*(-1+y)*(-1+x)*p[48]+4*x*z*(-1+y)*p[48]+4*y*z*(-1+y)*p[51]-4*y*z*(-1+y)*p[54]-4*y*z*(-1+x)*p[57]-4*x*y*z*p[57];

  dtau[0][1]=2*(-1+z)*(-1+y)*(-1+x)*p[0]-2*x*(-1+z)*(-1+y)*p[3]-2*x*y*(-1+z)*p[6]+2*y*(-1+z)*(-1+x)*p[9]-2*z*(-1+y)*(-1+x)*p[12]+2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]-2*y*z*(-1+x)*p[21]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[0]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[3]-x*(-1+z)*(2*y-2*z-3+2*x)*p[6]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[9]-z*(-1+x)*(2*x+2*y-2*z+1)*p[12]-x*z*(-2*y+2*z-3+2*x)*p[15]+x*z*(2*y+2*z-5+2*x)*p[18]+z*(-1+x)*(2*x-2*z+3-2*y)*p[21]-4*x*(-1+z)*(-1+x)*p[24]-4*(-1+z)*(-1+y)*(-1+x)*p[27]-4*y*(-1+z)*(-1+x)*p[27]-4*z*(-1+z)*(-1+x)*p[30]+4*x*(-1+z)*(-1+y)*p[33]+4*x*y*(-1+z)*p[33]+4*x*z*(-1+z)*p[36]+4*x*(-1+z)*(-1+x)*p[39]-4*x*z*(-1+z)*p[42]+4*z*(-1+z)*(-1+x)*p[45]+4*x*z*(-1+x)*p[48]+4*z*(-1+y)*(-1+x)*p[51]+4*y*z*(-1+x)*p[51]-4*x*z*(-1+y)*p[54]-4*x*y*z*p[54]-4*x*z*(-1+x)*p[57];

  dtau[0][2]=2*(-1+z)*(-1+y)*(-1+x)*p[0]-2*x*(-1+z)*(-1+y)*p[3]+2*x*y*(-1+z)*p[6]-2*y*(-1+z)*(-1+x)*p[9]+2*z*(-1+y)*(-1+x)*p[12]-2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]-2*y*z*(-1+x)*p[21]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[0]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[3]-x*y*(2*y-2*z-3+2*x)*p[6]-y*(-1+x)*(2*x+2*z+1-2*y)*p[9]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[12]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[15]+x*y*(2*y+2*z-5+2*x)*p[18]+y*(-1+x)*(2*x-2*z+3-2*y)*p[21]-4*x*(-1+y)*(-1+x)*p[24]-4*y*(-1+y)*(-1+x)*p[27]-4*(-1+z)*(-1+y)*(-1+x)*p[30]-4*z*(-1+y)*(-1+x)*p[30]+4*x*y*(-1+y)*p[33]+4*x*(-1+z)*(-1+y)*p[36]+4*x*z*(-1+y)*p[36]+4*x*y*(-1+x)*p[39]-4*x*y*(-1+z)*p[42]-4*x*y*z*p[42]+4*y*(-1+z)*(-1+x)*p[45]+4*y*z*(-1+x)*p[45]+4*x*(-1+y)*(-1+x)*p[48]+4*y*(-1+y)*(-1+x)*p[51]-4*x*y*(-1+y)*p[54]-4*x*y*(-1+x)*p[57];

  dtau[1][0]=2*(-1+z)*(-1+y)*(-1+x)*p[1]+2*x*(-1+z)*(-1+y)*p[4]-2*x*y*(-1+z)*p[7]-2*y*(-1+z)*(-1+x)*p[10]-2*z*(-1+y)*(-1+x)*p[13]-2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]+2*y*z*(-1+x)*p[22]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[1]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[4]-y*(-1+z)*(2*y-2*z-3+2*x)*p[7]-y*(-1+z)*(2*x+2*z+1-2*y)*p[10]-z*(-1+y)*(2*x+2*y-2*z+1)*p[13]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[16]+y*z*(2*y+2*z-5+2*x)*p[19]+y*z*(2*x-2*z+3-2*y)*p[22]-4*(-1+z)*(-1+y)*(-1+x)*p[25]-4*x*(-1+z)*(-1+y)*p[25]-4*y*(-1+z)*(-1+y)*p[28]-4*z*(-1+z)*(-1+y)*p[31]+4*y*(-1+z)*(-1+y)*p[34]+4*z*(-1+z)*(-1+y)*p[37]+4*y*(-1+z)*(-1+x)*p[40]+4*x*y*(-1+z)*p[40]-4*y*z*(-1+z)*p[43]+4*y*z*(-1+z)*p[46]+4*z*(-1+y)*(-1+x)*p[49]+4*x*z*(-1+y)*p[49]+4*y*z*(-1+y)*p[52]-4*y*z*(-1+y)*p[55]-4*y*z*(-1+x)*p[58]-4*x*y*z*p[58];

  dtau[1][1]=2*(-1+z)*(-1+y)*(-1+x)*p[1]-2*x*(-1+z)*(-1+y)*p[4]-2*x*y*(-1+z)*p[7]+2*y*(-1+z)*(-1+x)*p[10]-2*z*(-1+y)*(-1+x)*p[13]+2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]-2*y*z*(-1+x)*p[22]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[1]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[4]-x*(-1+z)*(2*y-2*z-3+2*x)*p[7]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[10]-z*(-1+x)*(2*x+2*y-2*z+1)*p[13]-x*z*(-2*y+2*z-3+2*x)*p[16]+x*z*(2*y+2*z-5+2*x)*p[19]+z*(-1+x)*(2*x-2*z+3-2*y)*p[22]-4*x*(-1+z)*(-1+x)*p[25]-4*(-1+z)*(-1+y)*(-1+x)*p[28]-4*y*(-1+z)*(-1+x)*p[28]-4*z*(-1+z)*(-1+x)*p[31]+4*x*(-1+z)*(-1+y)*p[34]+4*x*y*(-1+z)*p[34]+4*x*z*(-1+z)*p[37]+4*x*(-1+z)*(-1+x)*p[40]-4*x*z*(-1+z)*p[43]+4*z*(-1+z)*(-1+x)*p[46]+4*x*z*(-1+x)*p[49]+4*z*(-1+y)*(-1+x)*p[52]+4*y*z*(-1+x)*p[52]-4*x*z*(-1+y)*p[55]-4*x*y*z*p[55]-4*x*z*(-1+x)*p[58];

  dtau[1][2]=2*(-1+z)*(-1+y)*(-1+x)*p[1]-2*x*(-1+z)*(-1+y)*p[4]+2*x*y*(-1+z)*p[7]-2*y*(-1+z)*(-1+x)*p[10]+2*z*(-1+y)*(-1+x)*p[13]-2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]-2*y*z*(-1+x)*p[22]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[1]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[4]-x*y*(2*y-2*z-3+2*x)*p[7]-y*(-1+x)*(2*x+2*z+1-2*y)*p[10]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[13]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[16]+x*y*(2*y+2*z-5+2*x)*p[19]+y*(-1+x)*(2*x-2*z+3-2*y)*p[22]-4*x*(-1+y)*(-1+x)*p[25]-4*y*(-1+y)*(-1+x)*p[28]-4*(-1+z)*(-1+y)*(-1+x)*p[31]-4*z*(-1+y)*(-1+x)*p[31]+4*x*y*(-1+y)*p[34]+4*x*(-1+z)*(-1+y)*p[37]+4*x*z*(-1+y)*p[37]+4*x*y*(-1+x)*p[40]-4*x*y*(-1+z)*p[43]-4*x*y*z*p[43]+4*y*(-1+z)*(-1+x)*p[46]+4*y*z*(-1+x)*p[46]+4*x*(-1+y)*(-1+x)*p[49]+4*y*(-1+y)*(-1+x)*p[52]-4*x*y*(-1+y)*p[55]-4*x*y*(-1+x)*p[58];

  dtau[2][0]=2*(-1+z)*(-1+y)*(-1+x)*p[2]+2*x*(-1+z)*(-1+y)*p[5]-2*x*y*(-1+z)*p[8]-2*y*(-1+z)*(-1+x)*p[11]-4*y*(-1+z)*(-1+y)*p[29]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[2]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[5]-y*(-1+z)*(2*y-2*z-3+2*x)*p[8]-y*(-1+z)*(2*x+2*z+1-2*y)*p[11]-z*(-1+y)*(2*x+2*y-2*z+1)*p[14]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[17]+y*z*(2*y+2*z-5+2*x)*p[20]+y*z*(2*x-2*z+3-2*y)*p[23]-4*(-1+z)*(-1+y)*(-1+x)*p[26]-4*x*(-1+z)*(-1+y)*p[26]-2*z*(-1+y)*(-1+x)*p[14]-2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]+2*y*z*(-1+x)*p[23]-4*z*(-1+z)*(-1+y)*p[32]+4*y*(-1+z)*(-1+y)*p[35]+4*z*(-1+z)*(-1+y)*p[38]+4*y*(-1+z)*(-1+x)*p[41]+4*x*y*(-1+z)*p[41]-4*y*z*(-1+z)*p[44]+4*y*z*(-1+z)*p[47]+4*z*(-1+y)*(-1+x)*p[50]+4*x*z*(-1+y)*p[50]+4*y*z*(-1+y)*p[53]-4*y*z*(-1+y)*p[56]-4*y*z*(-1+x)*p[59]-4*x*y*z*p[59];

  dtau[2][1]=2*(-1+z)*(-1+y)*(-1+x)*p[2]-2*x*(-1+z)*(-1+y)*p[5]-2*x*y*(-1+z)*p[8]+2*y*(-1+z)*(-1+x)*p[11]-2*z*(-1+y)*(-1+x)*p[14]+2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]-2*y*z*(-1+x)*p[23]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[2]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[5]-x*(-1+z)*(2*y-2*z-3+2*x)*p[8]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[11]-z*(-1+x)*(2*x+2*y-2*z+1)*p[14]-x*z*(-2*y+2*z-3+2*x)*p[17]+x*z*(2*y+2*z-5+2*x)*p[20]+z*(-1+x)*(2*x-2*z+3-2*y)*p[23]-4*x*(-1+z)*(-1+x)*p[26]-4*(-1+z)*(-1+y)*(-1+x)*p[29]-4*y*(-1+z)*(-1+x)*p[29]-4*z*(-1+z)*(-1+x)*p[32]+4*x*(-1+z)*(-1+y)*p[35]+4*x*y*(-1+z)*p[35]+4*x*z*(-1+z)*p[38]+4*x*(-1+z)*(-1+x)*p[41]-4*x*z*(-1+z)*p[44]+4*z*(-1+z)*(-1+x)*p[47]+4*x*z*(-1+x)*p[50]+4*z*(-1+y)*(-1+x)*p[53]+4*y*z*(-1+x)*p[53]-4*x*z*(-1+y)*p[56]-4*x*y*z*p[56]-4*x*z*(-1+x)*p[59];

  dtau[2][2]=2*(-1+z)*(-1+y)*(-1+x)*p[2]-2*x*(-1+z)*(-1+y)*p[5]+2*x*y*(-1+z)*p[8]-2*y*(-1+z)*(-1+x)*p[11]+2*z*(-1+y)*(-1+x)*p[14]-2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]-2*y*z*(-1+x)*p[23]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[2]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[5]-x*y*(2*y-2*z-3+2*x)*p[8]-y*(-1+x)*(2*x+2*z+1-2*y)*p[11]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[14]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[17]+x*y*(2*y+2*z-5+2*x)*p[20]+y*(-1+x)*(2*x-2*z+3-2*y)*p[23]-4*x*(-1+y)*(-1+x)*p[26]-4*y*(-1+y)*(-1+x)*p[29]-4*(-1+z)*(-1+y)*(-1+x)*p[32]-4*z*(-1+y)*(-1+x)*p[32]+4*x*y*(-1+y)*p[35]+4*x*(-1+z)*(-1+y)*p[38]+4*x*z*(-1+y)*p[38]+4*x*y*(-1+x)*p[41]-4*x*y*(-1+z)*p[44]-4*x*y*z*p[44]+4*y*(-1+z)*(-1+x)*p[47]+4*y*z*(-1+x)*p[47]+4*x*(-1+y)*(-1+x)*p[50]+4*y*(-1+y)*(-1+x)*p[53]-4*x*y*(-1+y)*p[56]-4*x*y*(-1+x)*p[59];

}

void Ref2Phy(__constant double* physnode,
             double xref[3],
             double dphiref[3],
             int ifa,
             double xphy[3],
             double dtau[3][3],
             double codtau[3][3],
             double dphi[3],
             double vnds[3]) {

  // compute the mapping and its jacobian
  double x = xref[0];
  double y = xref[1];
  double z = xref[2];

  // gradient of the shape functions and value (4th component)
  // of the shape functions
  double gradphi[20][4];
  {
    double t1 = -1 + z;
    double t2 = -1 + y;
    double t3 = t1 * t2;
    double t4 = 2 * y;
    double t5 = 2 * z;
    double t6 = 4 * x;
    double t9 = -1 + x;
    double t10 = t1 * t9;
    double t11 = 2 * x;
    double t12 = 4 * y;
    double t15 = t2 * t9;
    double t16 = 4 * z;
    double t24 = x * t1;
    double t27 = x * t2;
    double t33 = y * t1;
    double t38 = x * y;
    double t48 = y * t9;
    double t54 = z * t2;
    double t57 = z * t9;
    double t67 = x * z;
    double t75 = y * z;
    double t94 = t11 - 1;
    double t98 = 4 * t24 * t9;
    double t100 = 4 * t27 * t9;
    double t104 = 4 * t33 * t2;
    double t105 = t4 - 1;
    double t111 = 4 * y * t2 * t9;
    double t114 = z * t1;
    double t116 = 4 * t114 * t2;
    double t118 = 4 * t114 * t9;
    double t119 = t5 - 1;
    double t128 = 4 * t38 * t2;
    double t132 = 4 * t67 * t1;
    double t141 = 4 * t38 * t9;
    double t145 = 4 * t75 * t1;
    double t158 = 4 * t67 * t9;
    double t162 = 4 * t75 * t2;

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

  if (xphy != NULL) {
    for(int ii = 0; ii < 3; ++ii) {
      xphy[ii] = 0;
      for(int i = 0; i < 20; ++i) {
	xphy[ii] += physnode[ 3 * i + ii] * gradphi[i][3];
      }
    }
  }

  if (dtau != NULL) {
    for(int ii = 0; ii < 3; ii++) {
      for(int jj = 0; jj < 3; jj++) {
	dtau[ii][jj] = 0;
      }
      for(int i = 0; i < 20; i++) {
	for(int jj = 0; jj < 3; jj++) {
	  dtau[ii][jj] += physnode[3 * i + ii] * gradphi[i][jj];;
	}
      }
    }
  }

  if (codtau != NULL) {
    //assert(dtau != NULL);
    /* codtau[0] =  dtau[4] * dtau[8] - dtau[5] * dtau[7]; */
    /* codtau[1] = -dtau[3] * dtau[8] + dtau[5] * dtau[6]; */
    /* codtau[2] =  dtau[3] * dtau[7] - dtau[4] * dtau[6]; */
    /* codtau[3] = -dtau[1] * dtau[8] + dtau[2] * dtau[7]; */
    /* codtau[4] =  dtau[0] * dtau[8] - dtau[2] * dtau[6]; */
    /* codtau[5] = -dtau[0] * dtau[7] + dtau[1] * dtau[6]; */
    /* codtau[6] =  dtau[1] * dtau[5] - dtau[2] * dtau[4]; */
    /* codtau[7] = -dtau[0] * dtau[5] + dtau[2] * dtau[3]; */
    /* codtau[8] =  dtau[0] * dtau[4] - dtau[1] * dtau[3]; */
    codtau[0][0] =  dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
    codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
    codtau[0][2] =  dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
    codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
    codtau[1][1] =  dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
    codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
    codtau[2][0] =  dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
    codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
    codtau[2][2] =  dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
  }

  if (dphi != NULL) {
    //assert(codtau != NULL);
    for(int ii = 0; ii < 3; ii++) {
      dphi[ii]=0;
      for(int jj = 0; jj < 3; jj++) {
        dphi[ii] += codtau[ii][jj] * dphiref[jj];
      }
    }
  }

  if (vnds !=NULL) {
    //assert(codtau != NULL);
    //assert(ifa >=0);
    for(int ii = 0; ii < 3; ii++) {
      vnds[ii]=0;
      for(int jj = 0; jj < 3; jj++) {
        vnds[ii] += codtau[ii][jj] * h20_refnormal[ifa][jj];
      }
    }
  }

}

void Phy2Ref(__constant double *physnode, double xphy[3], double xref[3]) 
{
#define ITERNEWTON 10
  double dxref[3], dxphy[3];
  xref[0] = 0.5;
  xref[1] = 0.5;
  xref[2] = 0.5;
  for(int iter = 0; iter < ITERNEWTON; ++iter ) {
    double dtau[3][3], codtau[3][3];
    int ifa =- 1;
    Ref2Phy(physnode, xref, 0, ifa, dxphy, dtau, codtau, 0,0);
    dxphy[0] -= (xphy)[0];
    dxphy[1] -= (xphy)[1];
    dxphy[2] -= (xphy)[2];
    double det
      = dtau[0][0] * codtau[0][0]
      + dtau[0][1] * codtau[0][1]
      + dtau[0][2] * codtau[0][2];
    //assert(det>0);
    for(int ii = 0; ii < 3; ++ii) {
      dxref[ii] = 0;
      for(int jj = 0; jj < 3; ++jj) {
        dxref[ii] += codtau[jj][ii] * dxphy[jj];
      }
      xref[ii] -= dxref[ii] / det;
    }
  }
  //double eps=1e-2;  // may be to constraining...
  //assert(xref[0]<1+eps && xref[0]>-eps);
  //assert(xref[1]<1+eps && xref[1]>-eps);
  //assert(xref[2]<1+eps && xref[2]>-eps);
}

// From a reference point find the nearest gauss point
// Warning: works only  degree 1, 2, or 3 (FIXME: why?)
int ref_ipg(__constant int *param, double *xref) 
{
  // approximation degree in each direction
  int deg[3] = {param[0], param[1], param[2]};

  // number of subcells in each direction
  int nraf[3] = {param[3], param[4], param[5]};

  double hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  // get the subcell id
  int ncx = floor(xref[0] * nraf[0]);
  int ncy = floor(xref[1] * nraf[1]);
  int ncz = floor(xref[2] * nraf[2]);

  //printf("x=%f ncx=%d nrafx=%d\n",xref[0], ncx,nraf[0]);
  //printf("y=%f ncy=%d nrafy=%d\n",xref[1], ncy,nraf[1]);
  //printf("z=%f ncz=%d nrafz=%d\n",xref[2], ncz,nraf[2]);
  //assert(ncx >=0 && ncx<nraf[0]);
  //assert(ncy >=0 && ncy<nraf[1]);
  //assert(ncz >=0 && ncz<nraf[2]);

  // subcell index in the macrocell
  int nc = ncx + nraf[0] * (ncy + nraf[1] * ncz);
  int offset = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1) * nc;

  // round to the nearest integer
  int ix = floor((xref[0] - ncx * hh[0]) / hh[0] * deg[0] + 0.5);
  int iy = floor((xref[1] - ncy * hh[1]) / hh[1] * deg[1] + 0.5);
  int iz = floor((xref[2] - ncz * hh[2]) / hh[2] * deg[2] + 0.5);
  //int iz=floor(xref[2]*deg[2]+0.5);

  //printf("xref %f %f %f ix=%d iy=%d iz=%d\n",
  //	 xref[0],xref[1],xref[2], ix, iy, iz);

  return ix + (deg[0] + 1) * (iy + (deg[1] + 1) * iz) + offset;
};

// Out-of-place RK stage
__kernel
void RK_out_CL(__global double *wnp1, 
	       __global double *wn, 
	       __global double *dtwn, 
	       const double dt)
{
  int ipg = get_global_id(0);
  wnp1[ipg] = wn[ipg] + dt * dtwn[ipg];
}

// In-place RK stage
__kernel
void RK_in_CL(__global double *wnp1, 
	      __global double *dtwn, 
	      const double dt)
{
  int ipg = get_global_id(0);
  wnp1[ipg] += dt * dtwn[ipg];
}

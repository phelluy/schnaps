#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

#include "geometry.h"

// utilitie functions for computing interpolation
// on a macrocell
//! \brief a struct for managing geometric mapping
typedef struct Interpolation{
  //! \brief interpolation parameters
  //! \details generally the convention is <BR>
  //! param[0] = M number of variables <BR>
  //! param[1] = deg x <BR>
  //! param[2] = deg y <BR>
  //! param[3] = deg z <BR>
  //! param[4] = raf x <BR>
  //! param[5] = raf y <BR>
  //! param[6] = raf z <BR>
  //! param[7..] = others param or return from interp
  int interp_param[8];

  //! \brief underlying geometry mapping
  Geom geo;
  //! \brief number of Gauss points in the volume
  int npgv;
  //! \brief current Gauss point index
  int ipgv;
  //! \brief current face id
  int ifa;
  //! \brief number of Gauss points on the current face
  int npgf;
  //! \brief current face Gauss point index
  int ipgf;
  //! \brief basis function index
  int ib;
  //! \brief current volume Gauss weight
  real wpgv;
  //! \brief  current face Gauss weight
  real wpgf;
  //! \brief Gauss point ref location
  real xpgref[3];
  //! \brief Gauss point physical location
  real xpg[3];
  //! \brief basis function values
  real phi;
  //! \brief basis function reference gradient
  real dphiref[3];
  //! \brief basis function physical gradient
  real dphi[3];

} Interpolation;

// in all these functions
// param is an input integer array
// param[0..2] : approximation degree in the 3 directions on
// the reference cube
// param[3..5] : number of subcells in the 3 directions on
// the reference cube
// param[6..nparam]: can be used for returning
// parameters (the user has to ensure that enough memory is reserved

//! \brief number of Gauss-LObatto Points (GLOPs) on the macro cell
//! \param[in] param the param list
int NPG(int param[]);

//! \brief number of GLOPs on the face ifa of the macrocell
//! \param[in] param the param list
//! \param[in] ifa face index
int NPGF(int param[],int ifa);

//! \brief compute 3d glop and subcell indices from the index
//! of the glop in the macrocell
//! \param[in] deg the degrees list (size [3])
//! \param[in] raf the refinments list (size [3])
//! \param[in] ipg the glop index in the macrocell
//! \param[out] ic the 3 subcell indices in x,y,z directions
//! \param[out] ix the 3 glop indices in the subcell
#pragma start_opencl
void ipg_to_xyz(const int *raf, const int *deg, int *ic, int *ix, const 
		int *ipg);
#pragma end_opencl

//! \brief compute the index of the glop in the macrocell
//!  from 3d glop and subcell indices
//! \param[in] deg the degrees list (size [3])
//! \param[in] raf the refinments list (size [3])
//! \param[out] ipg the glop index in the macrocell
//! \param[in] ic the 3 subcell indices in x,y,z directions
//! \param[in] ix the 3 glop indices in the subcell
#pragma start_opencl
void xyz_to_ipg(int* raf,int* deg,int* ic,int* ix,int *ipg);
#pragma end_opencl

//! \brief return the reference coordinates xpg[3] and weight wpg of the GLOP ipg
//! \param[in] param interp. params list
//! \param[in] ipg Gauss point index
//! \param[out] xpg  reference Gauss point coordinates
//! \param[out] wpg reference Gauss weight
//! \param[in] xpg_in same as xpg but slightly moved such
//! that the resulting point is in the interior of the ref. element
void ref_pg_vol(int* param,int ipg,
		real* xpg,real* wpg,real* xpg_in);

//! \brief from a reference point find the nearest
//! gauss point
//! \param[in] param interp. params list
//! \param[in] xref  reference Gauss point coordinates
//! \return Gauss point index
#pragma start_opencl
int ref_ipg(__constant int* param,real* xref);
#pragma end_opencl

//! \brief compute the position xpg of glop ipg in the local
//! numbering on face ifa. If xpgin is not NULL also compute
//! the position of point slightly inside the opposite subcell.
//! \param[inout] param interp. params list. param[6] also contains the volume Gauss point index
//! \param[in] ifa local face index (0..5)
//! \param[in] ipgf local 2d Gauss point index on the face ifa
//! \param[out] xpg Gauss point coordinates
//! \param[out] wpg Gauss point weight.
//! \param[out] xpgin same as xpg but slightly moved such
//! that the resulting point is in the interior of the ref. element
void ref_pg_face(int* param,int ifa,int ipgf,real* xpg,real* wpg,
		 real* xpgin);
//! \brief return the value and the gradient of the basis
//! functions.
//! Warning: the value of the gradient is
//! not reliable if xref is on the boundary
//! of a subcell (because the gradient is discontinuous)
//! \param[in] param interpolation parammeters (degrees and refinements)
//! \param[in] ib basis function index
//! \param[in] xref position of a point in the reference element
//! \param[out] psi value of the basis function
//! \param[out] dpsiref gradient of the basis function in the reference element
void psi_ref(int* param, int ib, real* xref, real* psi, real* dpsiref);

//! \brief  gradient of a basis function at a given glop 
//! \param[in]  param interpolation parammeters (degrees and refinements)
//! \param[in] ib basis function index
//! \param[in] ipg glop index
//! \param[out] dpsiref gradient of the basis function in the reference element
void grad_psi_pg(int* param,int ib,int ipg,real* dpsiref);

//! \brief  gradient of a basis function at a given glop in a given subcell
//! \param[in]  param interpolation parammeters (degrees and refinements)
//! \param[in] is a vector of three integer indices locating the subcell
//! \param[in] ib basis function index
//! \param[in] xref position of a point in the reference element
//! \param[out] psi value of the basis function
//! \param[out] dpsiref gradient of the basis function in the reference element
void psi_ref_subcell(int* param, int* is,int ib, real* xref, real* psi, real* dpsiref);

//! \brief return the 1d ith GLOP position for degree deg
//! \param[in] deg degree
//! \param[in] i GLOP 1D index
//! \returns the position in [0,1]
real glop(int deg,int i);

//! \brief return the 1d ith GLOP weight for degree deg
//! \param[in] deg degree
//! \param[in] i GLOP 1D index
//! \returns the weight
#pragma start_opencl
real wglop(int deg,int i);
#pragma end_opencl

//! \brief return the 1d derivative of lagrange polynomial ib at glop ipg
//! \param[in] deg degree
//! \param[in] ib basis function index
//! \param[in] ipg index of the Gauss point where the derivative is computed
//! \returns the value of the derivative
#pragma start_opencl
real dlag(int deg,int ib,int ipg);
#pragma end_opencl


//! \brief return the value of a 1D lagrange polynomial
//! \param[in] p value of the Lagrange polynomial
//! \param[in] subdiv the list of interpolation points (of size deg+1)
//! \param[in] deg polynomial degree
//! \param[in] ii index of the Lagrange polynomial
//! \param[in] x position where to compute p
void lagrange_polynomial(real* p,const real* subdiv,
			 int deg,int ii,real x);


//! \brief return the derivative of a 1D lagrange polynomial
//! \param[in] dp value of the derivative of the Lagrange polynomial
//! \param[in] subdiv the list of interpolation points (of size deg+1)
//! \param[in] deg polynomial degree
//! \param[in] ii index of the Lagrange polynomial
//! \param[in] x position where to compute dp
void dlagrange_polynomial(real* dp,const real* subdiv,
			  int deg,int ii,real x);



#endif

#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "interpolation.h"
#include "model.h"

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif

//! \brief A simple struct for packing a field
//! and a faces range.  To be passed to a thread
//! as a void* pointer.
typedef struct MacroFace {
  int first; //!< first cell/face index
  int last_p1;  //!< last cell/face index + 1
  //field *field; //! pointer to a  field
} MacroFace;

//! \brief A simple struct for packing a field
//! and a cells range.  To be passed to a thread
//! as a void* pointer.
typedef struct MacroCell {
  int first; //!< first cell/face index
  int last_p1;  //!< last cell/face index + 1
  //  field *field; //! pointer to a  field
} MacroCell;

//! \brief Data structure for managing a  discrete vector field
//! solution of a DG approximation
typedef struct field {
  //! Underlying mesh
  MacroMesh macromesh;
  //! Physical and numerical model
  Model model;
  //! Interpolation used for each component of the field
  Interpolation interp;
  //! A copy of the interpolation parameters
  int interp_param[8];
  //! Current time
  double tnow;
  //! CFL parameter min_i (vol_i / surf_i)
  double hmin;
  //! Time step
  //! dt has to be smaller than hmin / vmax
  double dt;
  
  int itermax;

  //! Activate or not 2D computations
  bool is2d;

  //! Size of the field buffers
  int wsize;
  //! fields at time steps n
  double *wn;
  //! Time derivative of the field
  double *dtwn;
  //! vmax
  double vmax;

  //! \brief Memory arrangement of field components
  //! \param[in] param interpolation parameters
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  int (*varindex)(int* param, int elem, int ipg, int iv);

  // Array of pointers to MacroFaces used for DGMacroCellInterface
  MacroFace *mface;
  // Array of pointers to MacroCell used for various DG routines
  MacroCell *mcell;

#ifdef _WITH_OPENCL
  //! \brief opencl data
  CLInfo cli;
  //! \brief copy of the dtwn array
  cl_mem wn_cl;
  cl_mem dtwn_cl;
  //! \brief copy of the params
  cl_mem param_cl;
  //! \brief copy physnode
  cl_mem physnode_cl;
  cl_double *physnode;

  cl_mem physnodeR_cl;
  cl_double *physnodeR;

  //! opencl kernels
  cl_kernel dgmass;
  cl_kernel dgflux;
  cl_kernel dgvolume;
  cl_kernel dginterface;
  cl_kernel RK_out_CL;
  cl_kernel RK_in_CL;
  cl_kernel RK4_final_stage;
  cl_kernel zero_buf;

  // OpenCL events

  // used in update_physnode_cl
  cl_event clv_mapdone; 
  
  // set_buf_to_zero event
  cl_event clv_zbuf; 
  
  cl_event clv_physnodeupdate;

  // Subcell mass events
  cl_event clv_mass; 

  // Subcell flux events
  cl_event *clv_flux;

  // Subcell volume events
  cl_event clv_volume; 

  // Macrocell interface events
  cl_event clv_mci;
  cl_event clv_interkernel; 
  cl_event clv_interupdate;
  cl_event clv_interupdateR;

  // OpenCL timing
  cl_ulong zbuf_time, mass_time, vol_time, flux_time, minter_time, rk_time;
#endif

} field;

#pragma start_opencl
//! \brief memory arrangement of field components.
//! Generic implementation.
//! \param[in] param interpolation parameters
//! param[0] = M
//! param[1] = deg x
//! param[2] = deg y
//! param[3] = deg z
//! param[4] = raf x
//! param[5] = raf y
//! param[6] = raf z
//! \param[in] elem macro element index
//! \param[in] ipg glop index
//! \param[in] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn.
int GenericVarindex(__constant int *param, int elem, int ipg, int iv);
#pragma end_opencl

#pragma start_opencl
//! \brief memory arrangement of field components.
//! with 3D components
//! \param[in] param interpolation parameters
//! param[0] = M
//! param[1] = deg x
//! param[2] = deg y
//! param[3] = deg z
//! param[4] = raf x
//! param[5] = raf y
//! param[6] = raf z
//! \param[in] elem macro element index
//! \param[in] ix components of the glop in its subcell
//! \param[in] ic components of the subcell in the macrocell
//! \param[in] iv component of the conservative variable
/* //! \returns the memory position in the arrays wn wnp1 or dtwn. */
/* int GenericVarindex3d(int m, */
/* 		      int nx0, int nx1, int nx2,  */
/* 		      int nc0, int nc1, int nc2,  */
/* 		      int elem, int *ix,int *ic, int iv); */
#pragma end_opencl

//! field initialization. Computation of the initial at each glop.
//! \param[inout] f a field
void Initfield(field *f);

//! free the buffers created in Initfield.
//! \param[inout] f a field
void Freefield(field *f);

//! \brief compute the Discontinuous Galerkin volume terms
//! slow version (no optimization using the tensors products)
//! \param[inout] f a field
void DGVolumeSlow(field *f);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. One subcell implementation.
//! \param[inout] f a field
void dtfieldSlow(field *f);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. Works with several subcells.
//! Fast version: multithreaded and with tensor products optimizations
//! \param[inout] f a field
void dtfield(field *f, double *w, double *dtw);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void DGMacroCellInterfaceSlow(void *mcell, field *f, double *w, double *dtw);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms second implementation with a loop on the faces
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void DGMacroCellInterface(void *mface, field *f, double *w, double *dtw);

//! \brief compute the Discontinuous Galerkin volume terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void DGVolume(void *mcell, field *f, double *w, double *dtw);

//! \brief compute the Discontinuous Galerkin inter-subcells terms
//! \param[inout] mcell a MacroCell
void DGSubCellInterface(void *mcell, field *f, double *w, double *dtw);

//! \brief  apply the DG mass term
//! \param[inout] mcell a MacroCell
void DGMass(void *mcell, field *f, double *w, double *dtw);

//! \brief An out-of-place RK stage
//! \param[out] fwnp1 field at time n+1
//! \param[in] fwn field at time n
//! \param[in] fdtwn time derivative of the field
//! \param[in] time step
//! \param[in] size of the field buffer
void RK_out(double *fwnp1, double *fwn, double *fdtwn, const double dt, 
	    const int sizew);

//! \brief An in-place RK stage
//! \param[inout] fwnp1 field at time n+1
//! \param[in] fdtwn time derivative of the field
//! \param[in] time step
//! \param[in] size of the field buffer
void RK_in(double *fwnp1, double *fdtwn, const double dt, const int sizew);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2(field *f, double tmax);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK4(field *f, double tmax);

#ifdef _WITH_OPENCL
//! \brief OpenCL version of RK2
//! time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2_CL(field *f, double tmax, 
	    cl_uint nwait, cl_event *wait, cl_event *done);
void RK4_CL(field *f, double tmax, 
	    cl_uint nwait, cl_event *wait, cl_event *done);
#endif

//! \brief save the results in the gmsh format
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! \param[in] f a field
//! with the analytical solution
//! \param[in] filename the path to the gmsh visualization file.
void Plotfield(int typplot, int compare, field *f, char *fieldname, 
	       char *filename);

//! \brief  display the field on screen
//! \param[in] f the field.
void Displayfield(field *f);

//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] f the field.
//! \returns the error.
double L2error(field *f);

#endif

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
  real tnow;
  //! CFL parameter min_i (vol_i / surf_i)
  real hmin;


  // TODO: once the output of the diagnostics is done by appending,
  // remove dt, ieter_time, itermax, nb_diags, and Diagnostics.
  int iter_time;
  //! final time iter
  int itermax;
  //! nb of diagnostics
  int nb_diags;
  //! table for diagnostics
  real *Diagnostics;

  //! Size of the field buffers
  int wsize;
  //! fields at time steps n
  real *wn;
  //! Time derivative of the field
  real *dtwn;
  //! vmax
  real vmax;

  //! \brief Pointer to a generic function called before computing dtfield. 
  //! \param[inout] f a field (to be converted from void*)
  void (*pre_dtfield)(void *f, real *w);

  //! \brief Pointer to a generic function called after computing dtfield. 
  //! \param[inout] f a field (to be converted from void*)
  void (*post_dtfield)(void *f, real *w);

  //! \brief generic update function called 
  //! \brief called at each runge-kutta sustep
  //! \param[inout] f a field (to be converted from void*)
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  void (*update_after_rk)(void *f, real *w);

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
  real *physnode;

  cl_mem physnodeR_cl;
  real *physnodeR;

  bool use_source_cl;
  char *sourcename_cl;

  //! opencl kernels
  cl_kernel dgmass;
  cl_kernel dgflux;
  cl_kernel dgvolume;
  cl_kernel dgsource;
  cl_kernel dginterface;
  cl_kernel dgboundary;
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

  // Subcell volume events
  cl_event clv_source; 

  // Macrocell interface events
  cl_event clv_mci;
  cl_event clv_interkernel; 
  cl_event clv_interupdate;
  cl_event clv_interupdateR;

  // OpenCL timing
  cl_ulong zbuf_time;
  cl_ulong mass_time;
  cl_ulong vol_time;
  cl_ulong flux_time;
  cl_ulong minter_time;
  cl_ulong boundary_time;
  cl_ulong source_time;
  cl_ulong rk_time;

  // OpenCL roofline measurements
  unsigned long int flops_vol, flops_flux, flops_mass; 
  unsigned long int reads_vol, reads_flux, reads_mass; 
#endif
} field;

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
#pragma start_opencl
int GenericVarindex(__constant int *param, int elem, int ipg, int iv);
#pragma end_opencl

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
#pragma start_opencl
/* int GenericVarindex3d(int m, */
/* 		      int nx0, int nx1, int nx2, */
/* 		      int nc0, int nc1, int nc2, */
/* 		      int elem, int *ix,int *ic, int iv); */
#pragma end_opencl

//! field initialization. Computation of the initial at each glop.
//! \param[inout] f a field
void Initfield(field *f);

void init_empty_field(field *f);

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
void dtfield(field *f, real *w, real *dtw);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void DGMacroCellInterfaceSlow(void *mcell, field *f, real *w, real *dtw);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms second implementation with a loop on the faces
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mface a MacroFace
void DGMacroCellInterface(void *mface, field *f, real *w, real *dtw);

//! \brief compute the Discontinuous Galerkin volume terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void DGVolume(void *mcell, field *f, real *w, real *dtw);

//! \brief compute the Discontinuous Galerkin inter-subcells terms
//! \param[inout] mcell a MacroCell
void DGSubCellInterface(void *mcell, field *f, real *w, real *dtw);

//! \brief  apply the DG mass term
//! \param[inout] mcell a MacroCell
void DGMass(void *mcell, field *f, real *dtw);

//! \brief Add the source term
//! \param[inout] mcell a MacroCell
//! \param[in] w: the field
//! \param[out] dtw: the derivative
void DGSource(void *mcell, field *f, real *w, real *dtw);

//! \brief An out-of-place RK stage
//! \param[out] fwnp1 field at time n+1
//! \param[in] fwn field at time n
//! \param[in] fdtwn time derivative of the field
//! \param[in] time step
//! \param[in] size of the field buffer
void RK_out(real *fwnp1, real *fwn, real *fdtwn, const real dt, 
	    const int sizew);

//! \brief An in-place RK stage
//! \param[inout] fwnp1 field at time n+1
//! \param[in] fdtwn time derivative of the field
//! \param[in] time step
//! \param[in] size of the field buffer
void RK_in(real *fwnp1, real *fdtwn, const real dt, const int sizew);

real set_dt(field *f);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2(field *f, real tmax, real dt);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK4(field *f, real tmax, real dt);

#ifdef _WITH_OPENCL
//! \brief OpenCL version of RK2
//! time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done);
void RK4_CL(field *f, real tmax, real dt,
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

//! \brief interpolate field at a reference point a macrocell
//! \param[in] f a field
//! \param[in] ie the macrocell index
//! \param[in] xref reference coordinates
//! \param[out] w the m field values
void InterpField(field *f,int ie,real* xref,real* w);

//! \brief  display the field on screen
//! \param[in] f the field.
void Displayfield(field *f);

//! \brief Save 1D results in a text file
//! \param[in] f the field.
//! \param[in] dir fixed direction to plot
//! \param[in] fixval fixed value to plot
//! \param[in] filename the path to the gmsh visualization file.
void Gnuplot(field* f,int dir, real fixval,char* filename);

//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] f the field.
//! \returns the error.
real L2error(field *f);

#endif

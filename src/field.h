#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "interpolation.h"
#include "model.h"

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif

//! \brief Data structure for managing a  discrete vector field
//! solution of a DG approximation
typedef struct Field{
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

  //! Activate or not 2D computations
  bool is2d;

  //! Size of the field buffers
  int wsize;
  //! Fields at time steps n
  double* wn;
  //! Fields at time steps n+1
  double* wnp1;
  //! Time derivative of the field
  double* dtwn;

  //! \brief Memory arrangement of field components
  //! \param[in] param interpolation parameters
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  int (*varindex)(int* param, int elem, int ipg, int iv);

#ifdef _WITH_OPENCL
  //! \brief opencl data
  CLInfo cli;
  //! \brief copy of the dtwn array
  cl_mem wn_cl;
  cl_mem wnp1_cl;
  cl_mem dtwn_cl;
  //! \brief copy of the params
  cl_mem param_cl;

  //! opencl kernels for mass inversion
  cl_kernel dgmass;
  cl_kernel dgvolume;
  cl_kernel dginterface;
  cl_kernel RK_out_CL;
  cl_kernel RK_in_CL;

#endif

} Field;

//! \brief A simple struct for packing a field
//! and a cells range.  To be passed to a thread
//! as a void* pointer.
typedef struct MacroCell {
  int first; //!< first cell/face index
  int last_p1;  //!< last cell/face index + 1
  Field *field; //! pointer to a  field
} MacroCell;

//! \brief A simple struct for packing a field
//! and a faces range.  To be passed to a thread
//! as a void* pointer.
typedef struct MacroFace {
  int first; //!< first cell/face index
  int last_p1;  //!< last cell/face index + 1
  Field *field; //! pointer to a  field
} MacroFace;

#pragma start_opencl
//! \brief memory arrangement of field components.
//! Generic implementation.
//! \param[in] param interpolation parameters
//! \param[in] elem macro element index
//! \param[in] ipg glop index
//! \param[in] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn.
int GenericVarindex(int *param, int elem, int ipg, int iv);
#pragma end_opencl

//! field initialization. Computation of the initial
//! at each glop.
//! \param[inout] f a field
void InitField(Field *f);

//! copy back the field to host memory
//! \param[inout] f a field
void CopyFieldtoCPU(Field *f);

//! \brief compute the Discontinuous Galerkin volume terms
//! slow version (no optimization using the tensors products)
//! \param[inout] f a field
void DGVolumeSlow(Field *f);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. One subcell implementation.
//! \param[inout] f a field
void dtFieldSlow(Field *f);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. Works with several subcells.
//! Fast version: multithreaded and with tensor products optimizations
//! \param[inout] f a field
void dtField(Field *f);

//! \brief OpenCL version of dtField : 
//! apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the field. Works with several subcells.
//! Fast version: multithreaded and with tensor products optimizations
//! \param[inout] f a field
void dtField_CL(Field *f);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void* DGMacroCellInterfaceSlow(void *mcell);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms second implementation with a loop on the faces
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void* DGMacroCellInterface(void *mface);

void* DGMacroCellInterface_CL(void *mface);

//! \brief compute the Discontinuous Galerkin volume terms
//! The argument has to be void* (for compatibility with pthread)
//! but it is logically a MacroCell*
//! \param[inout] mcell a MacroCell
void* DGVolume(void *mcell);

void* DGVolume_CL(void *mcell);

//! \brief compute the Discontinuous Galerkin inter-subcells terms
//! \param[inout] mcell a MacroCell
void* DGSubCellInterface(void *mcell);

//! \brief  apply the DG mass term
//! \param[inout] mcell a MacroCell
void* DGMass(void *mcell);

//! \brief  apply the DG mass term OpenCL version
//! \param[inout] mcell a MacroCell
void* DGMass_CL(void *mcell);

//! \brief exchange two pointers
//! \param[inout] a first pointer
//! \param[inout] b second pointer
void swap_pdoubles(double **a, double **b);

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
void RK2(Field *f, double tmax);

void RK2_step1_CL(Field *f);

//! \brief OpenCL version of RK2
//! time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2_CL(Field *f,double tmax);

//! \brief time integration by a second order Runge-Kutta algorithm.
//! slow version
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2Copy(Field *f,double tmax);

//! \brief save the results in the gmsh format
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! \param[in] f a field
//! with the analytical solution
//! \param[in] filename the path to the gmsh visualization file.
void PlotField(int typplot, int compare, Field *f, char *fieldname, 
	       char *filename);

//! \brief  display the field on screen
//! \param[in] f the field.
void DisplayField(Field *f);

//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] f the field.
//! \returns the error.
double L2error(Field *f);

#endif

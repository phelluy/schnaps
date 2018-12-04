#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "interpolation.h"
#include "model.h"
#include "linear_solver.h"
#include "skyline_spu.h"

#include <starpu.h>

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif


typedef int (*varindexptr)(int*, int*, int, int, int);


//! \brief Data structure for managing a  discrete vector field
//! solution of a DG approximation
typedef struct field {
  //! Physical nodes of the macrocell
  schnaps_real physnode[20][3];
  
  //! Physical and numerical model
  Model model;
  //! Interpolation used for each component of the field
  Interpolation interp;

  //! Refinement of the macrocell in each direction
  int raf[3];

  //! Degrees of interpolation in each direction
  int deg[3];

  //! number og Gauss points in each direction
  int npg[3];

  //! id of the field in the macromesh
  int id;
  
  //! Current time and time steps
  schnaps_real tnow;
  schnaps_real dt;

  //! Crank-Nicholson parameter
  schnaps_real theta;

  //! ref length of the mesh subcells
  schnaps_real hmin;

  //! two arrays for storing jacobians and det at glops
  //! NULL if store_det is false
  schnaps_real *jacobian;
  schnaps_real *det;
  bool store_det;

  //! period in each direction
  //! if negative: non-periodic computation (default)
  schnaps_real period[3];


  //! PIC struct pointer (=NULL if not used)
  //void *pic;

  //! a solver for the locally implicit scheme
  LinearSolver* solver;
  Skyline_SPU* solver_spu;
  
  //! a vector handle for the solver rhs
  starpu_data_handle_t rhs_handle;
  bool local_source_cl_init;
  
  //! another solver for computing the rhs
  //! the matrix is not factorized
  LinearSolver* rmat;
  Skyline_SPU* rmat_spu;
  
  //! starpu registering flag
  bool starpu_registered;

  
  //! Size of the field buffers
  int wsize;
  //! fields at current time step
  schnaps_real *wn;
  starpu_data_handle_t wn_handle;
  
  //! Time derivative of the field
  //! is equal to res
  //! only used in explicit schemes
  schnaps_real *dtwn;
  starpu_data_handle_t dtwn_handle;

  //! residual
  //! is equal to dtwn
  //! only used in implicit schemes
  schnaps_real *res;
  starpu_data_handle_t res_handle;

  //! \brief Pointer to a generic function called before computing dtfield. 
  //! \param[inout] f a field (to be converted from void*)
  void (*pre_dtfield)(void *f, schnaps_real *w);

  //! \brief Pointer to a generic function called after computing dtfield. 
  //! \param[inout] f a field (to be converted from void*)
  void (*post_dtfield)(void *f, schnaps_real *w);

  //! \brief generic update function called 
  //! \brief called at each runge-kutta sustep
  //! \param[inout] f a field (to be converted from void*)
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  void (*update_after_rk)(void *f, schnaps_real *w);

  //! \brief Memory arrangement of field components
  //! \param[in] param interpolation parameters  WRONG 
  //! \param[in] elem macro element index WRONG
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  int (*varindex)(int* deg, int *raf, int m, int ipg, int iv);

  //! \brief Memory arrangement of field components/Reverse lookup 
  //! \param[in] deg interpolation  degree vector
  //! \param[in] deg interpolation subcell refinement level vector
  //! \param[in] m model number of field components
  //! \param[in] imem memory position
  //! \param[out] ipg glop index in the arrays wn wnp1 or dtwn.
  //! \param[out] iv field component index
  void (*varindex_to_ipgiv)(int *def,int *raf,int m, int imem, int *ipg,int *iv);
} field;

//! \brief memory arrangement of field components.
//! Generic implementation.
//! \param[in] deg degrees parameters
//! \param[in] raf refinement parameters
//! \param[in] m number of conservative variables
//! \param[in] ipg glop index
//! \param[in] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn.
#pragma start_opencl
int GenericVarindex(__constant int *deg, __constant int *raf, int m,
		    int ipg, int iv);
#pragma end_opencl

//! \brief memory arrangement of field components/reverse lookup function 
//! CAUTION : must be consistent with GenericVarindex
//! Generic implementation. 
//! \param[in] deg degrees parameters
//! \param[in] raf refinement parameters
//! \param[in] m number of conservative variables
//! \param[in] imem memory position in the arrays wn wnp1 or dtwn.
//! \param[out] ipg glop index
//! \param[out] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn.

void GenericVarindex_toipgiv(int *deg, int *raf, int m, int imem,int *ipg,int *iv);

//! \brief memory arrangement of field components.
//! Generic implementation continuous case
//! \param[in] deg degrees parameters
//! \param[in] raf refinement parameters
//! \param[in] m number of conservative variables
//! \param[in] ipg glop index
//! \param[in] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn.
#pragma start_opencl
int GenericVarindex_CG(__constant int *deg, __constant int *raf, int m,
		    int ipg, int iv);
#pragma end_opencl

//! \brief field initialization. Computation of the initial at each glop.
//! \param[inout] f a field
//! \param[in] m a model
//! \param[in] physnode list of geometrical nodes of the macroelement
//! \param[in] deg degrees parameters 
//! \param[in] raf refinements parameters 
//! \param[in] w a pointer to field value (if NULL memory will be allocated)
//! \param[inout] dtw a pointer to derivatives (if NULL memory will be allocated)
void Initfield(field *f, Model m, schnaps_real physnode[][3], int *deg, int *raf, schnaps_real *w, schnaps_real* dtw);

//! \brief clean empty field initialization
//! \param[inout] f a field
void init_empty_field(field *f);

//! \brief init of starpu data inside a field (called at the end of Initfield)
//! \param[inout] f a field
void Initfield_SPU(field *f);

//! free the buffers created in Initfield.
//! \param[inout] f a field
void Freefield(field *f);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms second implementation with a loop on the faces
//! \param[in] locfaL local index of the face in the left element
//! \param[inout] fL a left field
//! \param[in] offsetL left offset for accessing data in w and dtw
//! \param[inout] fR a right field
//! \param[in] offsetR right offset for accessing data in w and dtw
//! \param[in] w field data 
//! \param[out] dtw time derivative of the field data w
void DGMacroCellInterface(int locfaL,
			  field *fL, int offsetL, field *fR, int offsetR,
			  schnaps_real *w, schnaps_real *dtw);

//! \brief compute the Discontinuous Galerkin volume terms
//! \param[in] f a field
void DGVolume(field *f, schnaps_real *w, schnaps_real *dtw);


//! \brief compute the average value on a cell for each cell of the mesh
//! \param[in] f a field
void DGAverage(field *f, schnaps_real *w);

//! \brief compute the Discontinuous Galerkin inter-subcells terms
//! \param[in] f a field
void DGSubCellInterface(field *f, schnaps_real *w, schnaps_real *dtw);

//! \brief  apply the DG mass term
//! \param[in] f a field
void DGMass(field *f, schnaps_real *w, schnaps_real *dtw);

//! \brief Add the source term
//! \param[in] f a field
void DGSource(field *f, schnaps_real *w, schnaps_real *dtw);

/* //! \brief save the results in the gmsh format */
/* //! \param[in] typplot index of the field variable to plot. */
/* //! \param[in] compare if true, the numerical solution is compared */
/* //! with the analytical solution */
/* //! \param[in] f a field */
/* //! \param[in] fieldname name of the plotted data */
/* //! \param[in] filename the path to the gmsh visualization file. */
/* void Plotfield(int typplot, int compare, field *f, char *fieldname,  */
/* 	       char *filename); */

//! \brief interpolate field at a reference point a macrocell
//! \param[in] f a field
//! \param[in] ie the macrocell index
//! \param[in] xref reference coordinates
//! \param[out] w the m field values
void InterpField(field *f,schnaps_real* xref,schnaps_real* w);

//! \brief  display the field on screen
//! \param[in] f the field.
void Displayfield(field *f);

//! \brief  register field data for starpuy
//! \param[inout] f the field.
void Registerfield_SPU(field *f);

//! \brief  unregister field data for starpuy
//! \param[inout] f the field.
void Unregisterfield_SPU(field *f);


#endif

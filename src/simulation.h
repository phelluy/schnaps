#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "macromesh.h"
#include "field.h"
#include "interface.h"

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif


//! \brief Data structure for managing a schnaps numerical simulation
typedef struct Simulation {
  //! Underlying mesh
  MacroMesh macromesh;
  //! List of fields for each macrocell
  field *fd;

  //! interfaces list
  Interface *interface;

  //! memory spaces for w and dtw
  schnaps_real *w;
  schnaps_real *dtw;
  schnaps_real *res;


  //! array of starpu data handles
  starpu_data_handle_t* w_handle;
  starpu_data_handle_t* dtw_handle;
  starpu_data_handle_t* res_handle;


  //! sum of sizes of field data
  int wsize;

  //!
  int interp_param[7];

  //! Current time
  schnaps_real tnow;
  //! CFL parameter min_i (vol_i / surf_i)
  schnaps_real hmin;
  //! Maximal CFL velocity
  schnaps_real vcfl;

  //! PIC struct pointer (=NULL if not used)
  void *pic;


  //! final time of simulation
  schnaps_real tmax;

  //! time step, theta (Crank-Nicholson) and cfl
  schnaps_real dt,theta,cfl;

  //! current iteration of the RK algorithm
  int iter_time_rk;
  //! maximal number of iterations
  int itermax_rk;
  //! nb of diagnostics
  int nb_diags;
  //! table for diagnostics
  schnaps_real *Diagnostics;


  //! \brief Pointer to a generic function called before computing dtfield.
  //! \param[inout] si a simulation (to be converted from void*)
  void (*pre_dtfields)(void *si);

  //! \brief Pointer to a generic function called after computing dtfield.
  //! \param[inout] si a simulation (to be converted from void*)
  void (*post_dtfields)(void *si);

  //! \brief generic update function called
  //! \brief called at each runge-kutta sustep
  //! \param[inout] si a simulation (to be converted from void*)
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  void (*update_after_rk)(void *si, schnaps_real *w);


  //! vmax
  schnaps_real vmax;

#ifdef _WITH_OPENCL
  //! \brief opencl data
  CLInfo cli;
  //! \brief copy of the dtwn array
  cl_mem w_cl;
  cl_mem dtw_cl;
  //! \brief copy of the params
  cl_mem param_cl;
  //! \brief copy physnode
  cl_mem physnode_cl;
  cl_mem physnodes_cl; // The physnodes for all the macrocells

  cl_mem physnodeR_cl;
  schnaps_real *physnodeR;

  bool use_source_cl;
  char *sourcename_cl;

  //! opencl kernels
  cl_kernel dgmass;
  cl_kernel dgflux;
  cl_kernel dgvolume;
  cl_kernel dgsource;
  cl_kernel dgcharge;
  cl_kernel dginterface;
  cl_kernel dgboundary;
  cl_kernel RK_out_CL;
  cl_kernel RK_in_CL;
  cl_kernel RK4_final_stage;
  cl_kernel zero_buf;

  // OpenCL events

  // set_buf_to_zero event
  cl_event clv_zbuf;

  // Subcell mass events
  cl_event *clv_mass;

  // Subcell flux events
  cl_event *clv_flux0, *clv_flux1, *clv_flux2;

  // Subcell volume events
  cl_event *clv_volume;

  // Subcell volume events
  cl_event *clv_source;

  // Subcell charge events
  cl_event *clv_charge;

  // Macrocell poisson event
  cl_event *clv_poisson;
  
  // Macrocell interface events
  cl_event *clv_mci;
  // Boundary term events
  cl_event *clv_boundary;

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


} Simulation;


//! \brief create an empty simulation in a coherent state
//! \param[inout] simu a simulation
void EmptySimulation(Simulation *simu);

//! \brief simulation initialization.
//! Computation of the initial data at each glop.
//! \param[inout] simu a simulation
//! \param[in] mesh a macromesh
//! \param[in] deg degrees parameters
//! \param[in] raf refinements parameters
//! \param[in] model a model
void InitSimulation(Simulation *simu, MacroMesh *mesh,
		    int *deg, int *raf, Model *model);

//! \brief init interface structs
//! \param[inout] simu a simulation
void InitInterfaces(Simulation *simu);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! \param[inout] simu A simulation
void DtFields(Simulation *simu);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! \param[inout] simu A simulation
void DtFields_old(Simulation *simu, schnaps_real *w, schnaps_real *dtw);

//! \brief compute the time step of the RK scheme
//! respecting a cfl condition
//! \param[inout] simu A simulation
schnaps_real Get_Dt_RK(Simulation *simu);

//! \brief An out-of-place RK stage
//! \param[out] fwnp1 field at time n+1
//! \param[in] fwn field at time n
//! \param[in] fdtwn time derivative of the field
//! \param[in] dt time step
//! \param[in] sizew size of the field buffer
void RK_out_old(schnaps_real *fwnp1, schnaps_real *fwn, schnaps_real *fdtwn, const schnaps_real dt,
	    const int sizew);

//! \brief An in-place RK stage
//! \param[inout] fwnp1 field at time n+1
//! \param[in] fdtwn time derivative of the field
//! \param[in] dt time step
//! \param[in] sizew size of the field buffer
void RK_in_old(schnaps_real *fwnp1, schnaps_real *fdtwn, const schnaps_real dt, const int sizew);

//! \brief An out-of-place RK stage
void RK_out(Simulation *simu, schnaps_real * wn);

//! \brief An-in-place RK stage
void RK_in(Simulation *simu);

void RK_Copy(Simulation * simu,schnaps_real * w, schnaps_real * w_temp);

void RK4_final_inplace(Simulation *simu, schnaps_real *l1, schnaps_real *l2, schnaps_real *l3, schnaps_real *l4);


//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] simu a simulation
//! \param[in] tmax physical duration of the simulation
void RK4_new(Simulation *simu, schnaps_real tmax);

//! \brief Final in-place RK stage
//! \param[out] w field at time n+1
//! \param[in] l1 first rk4 vector
//! \param[in] l2 second rk4 vector
//! \param[in] l3 third rk4 vector
//! \param[in] dtw last rk4 vector
//! \param[in] dt time step
//! \param[in] sizew size of the field buffer
void RK4_final_inplace_old(schnaps_real *w, schnaps_real *l1, schnaps_real *l2, schnaps_real *l3,
		       schnaps_real *dtw, const schnaps_real dt, const int sizew);


//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] simu a simulation
//! \param[in] tmax physical duration of the simulation
void RK1(Simulation *simu, schnaps_real tmax);


//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] simu a simulation
//! \param[in] tmax physical duration of the simulation
void RK2(Simulation *simu, schnaps_real tmax);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] simu a simulation
//! \param[in] tmax physical duration of the simulation
void RK4(Simulation *simu, schnaps_real tmax);


// TODO: see how to manage opencl...
/* #ifdef _WITH_OPENCL */
/* //! \brief OpenCL version of RK2 */
/* //! time integration by a second order Runge-Kutta algorithm */
/* //! \param[inout] f a field */
/* //! \param[in] tmax physical duration of the simulation */
/* void RK2_CL(field *f, real tmax, real dt, */
/* 	    cl_uint nwait, cl_event *wait, cl_event *done); */
/* void RK4_CL(field *f, real tmax, real dt, */
/* 	    cl_uint nwait, cl_event *wait, cl_event *done); */
/* #endif */

//! \brief save the results in the gmsh format
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! with the analytical solution
//! \param[in] simu a simulation
//! \param[in] fieldname name of the plotted data
//! \param[in] filename the path to the gmsh visualization file.
void PlotFields(int typplot, int compare, Simulation *simu, char *fieldname,
		char *filename);

//! \brief  list the valeus of the simulation
//! \param[in] simu a simulation.
void DisplaySimulation(Simulation *simu);

//! \brief Save 1D results in a text file
//! \param[in] simu a simulation.
//! \param[in] dir fixed direction to plot
//! \param[in] fixval fixed value to plot
//! \param[in] filename the path to the gmsh visualization file.
void Gnuplot(Simulation *simu,int dir, schnaps_real fixval,char* filename);

//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] simu a simulation.
//! \returns the error.
schnaps_real L2error(Simulation *simu);


//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] simu a simulation.
//! \param[in] nbvar index of one variable.
//! \returns the error.
schnaps_real L2error_onefield(Simulation *simu, int nbvar);

//! \brief frees any Simulation object
//! \param[inout] simu: a Simulation object
void freeSimulation(Simulation* simu);

//! \brief register simulation data into starpu
//! \param[inout] simu: a Simulation object
void RegisterSimulation_SPU(Simulation* simu);

//! \brief unregister simulation data into starpu
//! \param[inout] simu: a Simulation object
void UnregisterSimulation_SPU(Simulation* simu);

//! \brief Display an array
//! \param[in] array a real array
//! \param[in] size array size
//! \param[in] name an array name (appears on every line: make it short)
void DisplayArray(schnaps_real* array,
                  size_t size,
                  const char* name);
//! \brief Compute the gradient of a scalar field [NEED VALIDATION]
//! \param[in]  simu  a simulation object
//! \param[out] wd pointer to a real array of size (3 * total_nb_gauss_points)
//! \param[in] nbfield index of the field  
void Compute_derivative(Simulation *simu, schnaps_real * wd, int nbfield);
//! \brief save the results in the gmsh binary format / Sparser projection [NEED VALIDATION]
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! with the analytical solution
//! \param[in] simu a simulation
//! \param[in] fieldname name of the plotted data
//! \param[in] filename the path to the gmsh visualization file.
void PlotFieldsBinSparse(int typplot, int compare, Simulation *simu, char *fieldname,
		char *filename);
//! \brief save the results in the gmsh ascii format / Sparser projection [NEED VALIDATION]
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! with the analytical solution
//! \param[in] simu a simulation
//! \param[in] fieldname name of the plotted data
//! \param[in] filename the path to the gmsh visualization file.
void PlotFieldsAsciiSparse(int typplot, int compare, Simulation *simu, char *fieldname,
		char *filename);
#endif

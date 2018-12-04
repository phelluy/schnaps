#ifndef _LBMGENERIC_H
#define _LBMGENERIC_H
#include <string.h>
#include <stdbool.h>
#include <global.h>
#include "lbm_models_data.h"
#include "model.h"
#include "simulation.h"
#include "linear_solver.h"
#include "klu_csr.h"
//
static const int LBM_MAX_POLYSTRSIZE = 1024;
static const int LBM_MAX_MACROVARNAMESIZE = 32;
static const char LBM_POLYVARNAMES[3] = { 'X', 'Y', 'Z' };

//!\brief compact structure for encoding multivariate polynomials describing velocity moments. 
typedef struct MPolyDescriptor {
  //! space dimension i.e number of variables
  int ndim;
  //! number of multivariate monomials
  int nc;
  //! [nc] array monomial coefficients
  schnaps_real *c;
  //! [nc][nd] array exponents multiplets
  int **e;
  //! human readable name
  char *name;
  char *poly_string; 
} MPolyDescriptor;
//
static const MPolyDescriptor MPolyDescriptor_NULL;	// zeroed out structure
//!\brief Constructor 
//!\param[inout] mom a MPolyDescriptor
//!\param[in] ndim dimension of velocity space
//!\param[in] nc number of monomials
//!\param[in] name a human readable name
void NewMPolyDescriptor(MPolyDescriptor * mpd, int ndim, int nc,
			char *name);
//!\brief MPolyDescriptor Destructor
//!\param[inout] mom a MPolyDescriptor
void DestroyMPolyDescriptor(MPolyDescriptor * mpd);
//!\brief Sets the actual values of coefficients and exponents. Creates the string polynomial representation.
//!\param[inout] mom a MPolyDescriptor
//!\param[in] c [nc] array of monomial coefficients
//!\param[in] e [nc][ndim] array of exponents 
void SetMPolyDescriptor(MPolyDescriptor * mpd, schnaps_real c_in[mpd->nc],
			int e_in[mpd->nc][mpd->ndim]);
//!\brief Output name and Polynomial form to standard output
//!\param[in] mom a MPolyDescriptor
void DisplayMPolyDescriptor(MPolyDescriptor * mpd);
//!\brief Computes the value of multivariate polynomial
//!\param[in] mom a MPolyDescriptor
//!\param[in] V [ndim] array of values
//!\returns value of the polynomial at V 
schnaps_real EvaluatePoly(MPolyDescriptor * mpd, schnaps_real * V);
//
//!\brief create MPolyDescriptor from commonly used hard-coded models
//!\param[inout] mpd MPolyDescriptor
//!\param[in] pointer to predefined model
//!\param[in] name a human readable name, if set to NULL the default name encoded in the model is used.
void NewMPolyDescriptorFromModel(MPolyDescriptor * mpd,
				 const MPolyData * model, char *name);
//!\brief A structure describing a LBM model. 
typedef struct LBModelDescriptor {
  int d;			//!velocity space dimension
  int q;			//!number of velocity nodes
  int nb_macro;			//!number of conserved macro quantities NB : for multi-lattice models nc>q may be true.
  char **macro_names;		//! human readable names of macroscopic quantities
//  int nr;//!number of macro quantities to relax;
  schnaps_real cref;		//! velocity scale to apply to reference nodes
  schnaps_real **vi;		//! [q][d} array of velocity nodes;
  int *iopposite;		//! [q] array index of node with opposite velocity
  schnaps_real vmax;		//! maximum velocity 
  MPolyDescriptor *Moments;	//![q] array of moments descriptors for multi-relaxation models;
  int *inode_min;		//! [q] array ; for each moment starting index of node range for moment computation
  int *inode_max;		//! [q] array ; for each moment end index of node range for moment computation
  schnaps_real **M;		//! [q][q] moment matrix allowing to compute moments from f 
  LinearSolver *Msolv; //! linear solver object associated with moment matrix
  //schnaps_real **Minv; //! [q][q] inverse of moment matrix
  schnaps_real *nu;		//! [q] array of relaxation frequencies, will apply to either f or moments
  schnaps_real *s;		//! [q] array of relaxation parameters, will apply to either f or moments
  //!\brief compute macro quantities vector from f
  //!\param[in] vector f[size q] values of the distribution functions
  //!\param[out] wmac vector of nb_macro conserved macro quantities
  void (*f_to_macro) (const schnaps_real * f, schnaps_real * wmac);	// computation of macro-quantities from f
  //!\brief computation of equilibrium function for one velocity node from macro quantities vector
  //!\param[in] inode velocity node index in [0,q-1]
  //!\param[in] nb_macro size of macro quantities vector
  //!\param[in] wmac vector of macro quantities values
  //!\returns f_i
  schnaps_real(*feq_i) (int inode, int nb_macro, schnaps_real * wmac);
  //!\brief computation of a component of the jacobian equilibrium function for an (i,j) pair of velocity
  //! nodes from macro quantities vector.
  //!\param[in] inode velocity node index in [0,q-1]
  //!\param[in] jnode velocity node index in [0,q-1] 
  //!\param[in] nb_macro size of macro quantities vector
  //!\param[in] wmac vector of macro quantities values
  //!\returns df_i/df_j 
  schnaps_real(*Jfeq_ij) (int inode, int jnode, int nb_macro, schnaps_real * wmac);
  //!\brief computation of all equilibrium functions for all velocity nodes from macro quantities vector 
  //!\param[in] q number of velocity nodes, size of vector f
  //!\param[in] nb_macro number of macro quantities size of vector wmac
  //!\param[in] wmac vector of macro quantities values
  //!\param[out] f values of the equilibrium functions
  void (*feq) (int q, int nb_macro, schnaps_real *wmac, schnaps_real *f);
  //!\brief computation of the full jacobian matrix from macro quantities vector 
  //!\param[in] q number of velocity nodes, size of vector f
  //!\param[in] nb_macro number of macro quantities size of vector wmac
  //!\param[in] wmac vector of macro quantities values
  //!\param[out] Jf [q * q] array of values of the jacobian matrix dfi/dfj = J[j*q+i]
  void (*Jfeq) (int q, int nb_macro, schnaps_real *wmac, schnaps_real *J);
  //
  void *(model_spec_params);	// pointer to model specific parameters each model needing it should have its specific struc to hold
  // such parameters and init them in the corresponding set routine; 
} LBModelDescriptor;
//
static const LBModelDescriptor LBModelDescriptor_NULL;	//! zeroed out stucture
//
//!\brief LBModelDescriptor Constructor / Initializer
//!\param[inout] lb a LBModelDescriptor
//!\param[in] d dimension of velocity space
//!\param[in] q number of velocity nodes
void NewLBModelDescriptor(LBModelDescriptor * lb, int d, int nb_macro,
			  int q);
//!\brief LBModelDescriptor Destructor
//!\param[inout] lb a LBModelDescriptor
//!\paral[in] d dimension of velocity space
//!\paral[in] q number of velocity nodes
void DestroyLBModelDescriptor(LBModelDescriptor * lb);
//!\brief Moment matric computation (assuming Moments,cref and nodes are set)+ associated linear solver
void ComputeLBModelDescriptorMomentMatrix(LBModelDescriptor * lb);
//!brief print moment basis associated polynomials
void DisplayLBModelDescriptorMomentPoly(LBModelDescriptor * lb);
//!\brief print moment matrix M to standard output
void DisplayLBModelDescriptorMomentMatrix(LBModelDescriptor * lb);
//!\brief compute max velocity of node grid
void ComputeLBModelDescriptorVmax(LBModelDescriptor * lb);
//!\brief check that macro quantities are conserved;
void CheckLBModelDescriptorMacroConservation(LBModelDescriptor * lb,
					     bool verbose);
void CheckLBMMomentMatrixInversion(LBModelDescriptor *lb,bool verbose);
//
//******************************************************************************//
void LBM_Dummy_InitMacroData(schnaps_real x[3], schnaps_real w[]);
void LBM_Dummy_InitData(schnaps_real x[3], schnaps_real w[]);
void LBM_Dummy_InitData_OneNode(schnaps_real x[3], schnaps_real w[]);
//******************************************************************************//
//******************************************************************************//
typedef struct LBMSimulation {
  int d;			//!velocity space dimension
  int nb_macro_fields;		//! number of macro quantities
  int q;			//! total number of velocity nodes 
  schnaps_real vmax;		//!
  schnaps_real cfl;		//!
  schnaps_real dt;		//!
  schnaps_real tmax;		//!
  schnaps_real itermax;		//!
  schnaps_real tnow;
  schnaps_real iter_time;
  Model macro_model;		// Model object to intiate macroscopic simulation
  Simulation macro_simu;	//! a simulation object containing only macroscopic fields
  Model micro_model;		//! array of nb_lattices Models object to initiate corresponding simulations.
  Simulation micro_simu;	//! array of nb_lattices simulation object, one for each lattice.
  Model model_advec;		//! model for advection of 1 velocity node, used in the spliited per node implicit scheme ;
  //
  schnaps_real *wmic_buffer;	//! small [q] buffer to store the current micro state at the current ipg during assembly.
  schnaps_real *wmac_buffer;	//! small [q] buffer to store the current macro state at the current ipg during assembly.
  //
  void (*pre_advec) (void *lbsimu);	//!
  void (*post_advec_one_node) (void*lbsimu);  //! executed after each node individual advection step
  void (*post_advec) (void *lbsimu);	//! 
  void (*post_tstep) (void *lbsimu, schnaps_real * w);	//! w to stay interface compatible with update_after_rk of simu
  //
  schnaps_real diag_2d_period;	//! period for dumping of 2D diagnostics.
  void (*collect_diags) (void *lbsimu, schnaps_real * wmac, schnaps_real * wmic);	//! routine for 1D time traces collection
  //
  //
} LBMSimulation;
//*******************************************************************************//
void InitLBMSimulation(LBMSimulation * lbsimu,
		       LatticeBoltzmannSimData * lsd, MacroMesh * mesh,
		       int deg[3], int raf[3]);
void FreeBMSimulation(LBMSimulation * lbsimu);	//
void LB_Relaxation_bgk_f(void *lbsimu);	//
void LB_Relaxation_bgk_f_full(void *lbsimu);	//
/*void LB_Relaxation_Moments( LBMSimulation *lbsimu);*/
void LB_ComputeMacroFromMicro(void *lbsimu);	//
// wrappers for compatibility with RK schemes operating on 1 simulation //
// in that case only the micsimu object is passed to rk. The global lbsim object, which is 
// necessary for pre/post advection steps is passed through a global pointer (in schnaps_lbm_simdata) to the current simu.
// the following wrappers get a pointer to the lbsimu object and call the corresponding routines (pre_advec,post_advec,post_tstep)
void LBM_pre_dtfields_wrapper(void *simu);	// redirected to pre_advec;
void LBM_post_dtfields_wrapper(void *simu);	// redirected to post_advec;
void LBM_update_after_rk_wrapper(void *simu, schnaps_real * w);	// redirected to post_tstep;
//*******************************************************************************//
//******************************************************************************//
//******************************************************************************//
// flux functions for LB models
//******************************************************************************//
//******************************************************************************//
//! \brief flux for a lattice model;
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void LBM_OneLatticeNumFlux(schnaps_real * wL, schnaps_real * wR,
			   schnaps_real * vnorm, schnaps_real * flux);
//! \brief flux for a lattice model restricted to one velocity node whose index is set in the global lattice simdata structure
//! by the (current_advec_node_index) field.
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void LBM_OneLatticeNumFlux_SelectVel(schnaps_real * wL, schnaps_real * wR,
			   schnaps_real * vnorm, schnaps_real * flux);
//! \brief flux for a lattice model and a specific velocity nodes model and node index specified by the global shared
//! LatticeBoltzmannSimData structure
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void LBM_OneNodeNumFlux(schnaps_real * wL, schnaps_real * wR,
			schnaps_real * vnorm, schnaps_real * flux);
// basic source functions for relaxation
//! \brief standard source function for BGK relaxation of a lbm model
//! \param[in] x : spatial position 
//! \param[in] t : time
//! \param[in] w : local data value
//! \param[out] s: source term value 
void LBM_BGK_Relaxation_Source(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *s);
void LBM_BGK_Relaxation_Source_Jacobian(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *J);
void LBM_BGK_Relaxation_Source_TestDecay(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *s);
void LBM_BGK_Relaxation_Source_Jacobian_TestDecay(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *J);
//
// per model specific routines
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//!\brief zero equilibrium function for tests - component wise version
schnaps_real LBM_dummy_zero_feq_i(int inode, int nb_macro,
				 schnaps_real * w);
//!\brief zero equilibrium function for tests - vector version
void LBM_dummy_zero_feq(int q,int nb_macro, schnaps_real *w, schnaps_real *f);
//!\brief zero equilibrium function jacobian for tests - component wise version
schnaps_real LBM_dummy_zero_Jfeq_ij(int inode, int jnode, int nb_macro,
				 schnaps_real * w);
//!\brief zero equilibrium function jacobian for tests - matrix version
void LBM_dummy_zero_Jfeq(int q, int nb_macro, schnaps_real *w, schnaps_real *J);

// *********** artficial toy model 2DQ3_TEST for tests *************************//
// this is a 1DQ3 extended to 2D for light transport tests with a small number of velocities
// the velocity grid is a 1DQ3 (-1,0,1) model on an axis
// a parameter alpha sets the angle (in degrees 0 to 180) between the axis and the x reference axis
// the equilibrium conserve density and parallel momentum/velocity (ie aligned with the axis)
// caution: obviously this model is anisotropic
typedef struct params_D2Q3_TEST{
  schnaps_real alpha_deg;
  schnaps_real cosa;
  schnaps_real sina;
  schnaps_real theta;
  schnaps_real invtheta;
} params_D2Q3_TEST;
void LBM_Set_D2Q3_TEST_model(LBModelDescriptor *lb,schnaps_real alpha,schnaps_real cref);
void LBM_f_to_macro_D2Q3_TEST(const schnaps_real *f,schnaps_real *w);
schnaps_real LBM_feq_i_D2Q3_TEST(int inode,int nb_macro, schnaps_real *w);
void LBM_feq_D2Q3_TEST(int q,int nb_macro,schnaps_real *w,schnaps_real *f);
//******* Hydrodynamic models (Euler, Navier-Stokes) ***************************//
typedef struct params_D2Q9_ISOTH{
  schnaps_real theta;
  schnaps_real invtheta;
} params_D2Q9_ISOTH; 
// D2Q9 isothermal
void LBM_Set_D2Q9_ISOTH_model(LBModelDescriptor * lb, schnaps_real cref);
void LBM_f_to_macro_D2Q9_ISOTH(const schnaps_real * f, schnaps_real * w);
schnaps_real LBM_feq_i_D2Q9_ISOTH(int inode, int nb_macro, schnaps_real * w);
void LBM_feq_D2Q9_ISOTH(int q, int nb_macro, schnaps_real * w, schnaps_real *f);
// D2Q9 isothermal incompressible modif (with reference density rho0=1)
// see He and Luo Journal of Statistical Physics vol 88 (1997)
typedef struct params_D2Q9_ISOTH_INC{
  schnaps_real theta;
  schnaps_real invtheta;
  schnaps_real rho0;
} params_D2Q9_ISOTH_INC; 
void LBM_Set_D2Q9_ISOTH_INC_model(LBModelDescriptor * lb, schnaps_real cref);
void LBM_f_to_macro_D2Q9_ISOTH_INC(const schnaps_real * f, schnaps_real * w);
schnaps_real LBM_feq_i_D2Q9_ISOTH_INC(int inode, int nb_macro, schnaps_real * w);
void LBM_feq_D2Q9_ISOTH_INC(int q, int nb_macro, schnaps_real * w, schnaps_real *f);
// D2Q9 isothermal linearized (2D wave equation)
typedef struct params_D2Q9_ISOTH_LINEARIZED{
  schnaps_real theta;
  schnaps_real invtheta;
} params_D2Q9_ISOTH_LINEARIZED; 
void LBM_Set_D2Q9_ISOTH_LINEARIZED_model(LBModelDescriptor * lb,
					 schnaps_real cref);
void LBM_f_to_macro_D2Q9_ISOTH_LINEARIZED(const schnaps_real * f,
					  schnaps_real * w);
schnaps_real LBM_feq_i_D2Q9_ISOTH_LINEARIZED(int inode, int nb_macro,
					   schnaps_real * w);
void LBM_feq_D2Q9_ISOTH_LINEARIZED(int q, int nb_macro, schnaps_real * w, schnaps_real *f);
//
//*******************************************************************************//
//***************************** MDH models **************************************//
//*******************************************************************************//
// 2D D2Q9 + 3x 2DQ5 basic mhd model
typedef struct params_MHD_D2Q9_2D2Q5{
  schnaps_real theta;
  schnaps_real invtheta;
  schnaps_real theta_mag;
  schnaps_real invtheta_mag;
} params_MHD_D2Q9_2D2Q5; 
void LBM_Set_MHD_D2Q9_2D2Q5_model(LBModelDescriptor * lb, schnaps_real cref);
void LBM_f_to_macro_MHD_D2Q9_2D2Q5(const schnaps_real *f, schnaps_real *w);
schnaps_real LBM_feq_i_MHD_D2Q9_2D2Q5(int inode, int nb_macro, schnaps_real *w);
void LBM_feq_MHD_D2Q9_2D2Q5(int q, int nb_macro, schnaps_real * w, schnaps_real *f);
#endif

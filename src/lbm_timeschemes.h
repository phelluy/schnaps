#ifndef LBM_TIMESCHEMES_H
#define LBM_TIMESCHEMES_H
#include "lbm_generic.h"

//!\brief small struct for multi-index revovery
typedef struct LBM_glop_multindex{
  //! element index (macrocell)
  int ie; 
  //! field varindex
  int ivarindex;
  //! gauss-lobatto point  index in the element  
  int ipg;
  //! field component index
  int iv;
  //! subcell cell index
  int ic[3];
  //! subcell glop local index
  int ix[3];
} LBM_glop_multindex;
//!\brief recover multi index from global memory location
//!\param[in] simu
//!\param[in] imem memory location in wn
//!\param[out] glob a LBM_glop_multindex 
void LBM_imem_to_glopmulti(Simulation *simu,int imem,LBM_glop_multindex *glop);
//!\brief recover global memory location from multiindex
//!\param[in] simu
//!\param[in] glob a LBM_glop_multindex 
//!\returns imem memory location in wn
int LBM_glopmuti_to_imem(Simulation * simu,LBM_glop_multindex *glop);
//!\brief multiindex/memory mapping test
int CheckLBM_imem_multivarindex_constistency(Simulation *simu);
// LBM specific time schemes
//!\brief Theta (Impl/Expl) Time scheme for LBM with splitted advection/relaxation
//! advections are operated per velocity node using a one node temporary Simulation  object
//!\param[inout] lbmsimu a LBMSimulation object 
//!\param[in] theta , explicit/implicit weight (1/2 yields Crank-Nicholson)
//!\param[in] tmax end time of simulation
//!\param[in] dt time step duration
void LBMThetaTimeScheme(LBMSimulation * lbsimu, schnaps_real theta,
			schnaps_real tmax, schnaps_real dt);
// coupling routines for per velocity node splitted theta scheme
void LBM_MassCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky);
void LBM_InternalCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky);
void LBM_FluxCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky);
void LBM_InterfaceCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky);
void LBM_InitLinearSolver_OneNode(Simulation *simu, LinearSolver *solver, MatrixStorage ms);
//
void LBM_FluxAssembly_OneNode(Simulation * simu, LinearSolver * solver,
			   schnaps_real theta, schnaps_real dt);
void LBM_InterfaceAssembly_OneNode(Simulation * simu, LinearSolver * solver,
			   schnaps_real theta, schnaps_real dt);

//!\brief Multistep time scheme based on splitted per-node advection and CN relaxation 
void LBMThetaTimeScheme_OneNode_Multistep(LBMSimulation * lbsimu, schnaps_real theta,
			schnaps_real tmax, schnaps_real dt);

// modified assembly routines for per velocity node splitted theta scheme
// the modification was made necessary to allow for a selection of the velocity node index
// when performing per node advections in the LBMThetaTimeScheme routine 
//!\brief builds the (linear) operator for the weighted implicit/timescheme
//!\param[in] simu Simulation object
//!\param[inout] solver a linear solver
//!\param[in] theta weight factor of the scheme
//!\param[in] dt timestep duration
void LBM_AssemblyImplicitLinearSolver(Simulation * simu,
				      LinearSolver * solver,
				      schnaps_real theta, schnaps_real dt);
void LBM_InterfaceAssembly(Simulation * simu, LinearSolver * solver,
			   schnaps_real theta, schnaps_real dt);
void LBM_SourceAssembly(Simulation * simu, LinearSolver * solver,
			schnaps_real theta, schnaps_real dt);
//************************************************************************************//
void LBM_Test_Coupling_BlockCell(Simulation *simu, LinearSolver *solver,int isky);
void LBM_Mass_SourceJac_Coupling(Simulation *simu, LinearSolver *solver, int isky);
//! \brief Assembly of the DG operator into a sparse matrix with a constant velocity
//! the velocity node index is recovered from the global schnaps_lbm_simdata->current_advec_node_index
//! the velocity is recovered from schnaps_lbm_simdata->lb_model->vi
//! the routine takes into account the velocity norm and velocity index to remove uncessary couplings,
//! independently of actual cell geometry.
//! prepare the matrix structure of the fluxes inside the fields / volumic terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void LBM_InternalCoupling_SelectVel(Simulation *simu,  LinearSolver *solver);

//! \brief Assembly of the DG operator into a sparse matrix with a constant velocity
//! the velocity node index is recovered from the global schnaps_lbm_simdata->current_advec_node_index
//! the velocity is recovered from schnaps_lbm_simdata->lb_model->vi
//! the routine takes into account the velocity norm/direction and velocity index to remove uncessary couplings
//! prepare the matrix structure of the fluxes inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void LBM_FluxCoupling_SelectVel(Simulation *simu,  LinearSolver *solver);

void LBM_InterfaceCoupling_SelectVel(Simulation *simu,  LinearSolver *solver);

//!\brief Nonlinear weighted implicit/explicit time scheme
//! the scheme uses a symmetrized sequence of steps 
//! for each step 
//! -advection operates on one node while relaxation to equilibrium on all nodes.
//! -the block/triangular structure of the operators, detection of the block structure is done using klu 
//! - nonlinear resolution is done blockwise by a newton method.(hence the term local as a block maps to a cell)
//!\param[inout] lbsimu a LBMSimulation object
//!\param[in] theta weight facor of the scheme
//!\param[in] tmax end time of simulation
//!\param[in] dt time step duration

void LBMThetaTimeScheme_NL_LocalNewton(LBMSimulation *lbsimu,schnaps_real theta, schnaps_real tmax,schnaps_real dt);

void LBMThetaTimeScheme_NL_LocalNewton_solve(Simulation *simu,LinearSolver *impsolv,schnaps_real *exp_rhs,schnaps_real theta,schnaps_real dt);
int CheckKLUBLockCellConstitency(Simulation *simu,LinearSolver *solv,bool verbose);

void LBM_NL_LocalNewton_Assembly_Blockwise(Simulation *simu,LinearSolver *impsolv, const schnaps_real theta,  
  schnaps_real *exp_rhs,bool exp_is_assembled, int iblock);
#endif

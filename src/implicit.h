#ifndef _IMPLICIT_LS_H
#define _IMPLICIT_LS_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "simulation.h"
#include "linear_solver.h"
#include "advanced_linear_solver.h"


//! \brief Construct the profile of the linear solver
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] ms matrix storage type 
void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver, MatrixStorage ms);

//! \brief Assembly of the DG operator into a sparse matrix
//! computations of all the terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! computations of all the terms
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void AssemblyFieldImplicitSolver(field *fd,schnaps_real theta, schnaps_real dt);

//! \brief Assembly and solving of the implicit DG operator into
//! one macrocell. No global matrix allocation
//! thanks to downwind numbering
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FieldDownwindSolve(field *fd,schnaps_real theta, schnaps_real dt);

//! \brief residual of the DG operator into
//! one macrocell.
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FieldResidual(field *fd,schnaps_real theta, schnaps_real dt);


//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the differential terms inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] itest should be 0 (1 is used for debugging purposes)
void InternalCoupling(Simulation *simu,  LinearSolver *solver, int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the differential terms inside the fields
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] itest should be 0 (1 is used for debugging purposes)
void InternalLocalCoupling(field *fd, int itest);




//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the fluxes inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] itest should be 0 (1 for debugging purposes)
void FluxCoupling(Simulation *simu,  LinearSolver *solver,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the fluxes inside the fields
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] itest should be 0 (1 for debugging purposes)
void FluxLocalCoupling(field *fd,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the interface fluxes between fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] itest should be 0 (1 for debugging purposes)
void InterfaceCoupling(Simulation *simu,  LinearSolver *solver,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the differential terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InternalAssembly(Simulation *simu,  LinearSolver *solver,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of DG differential terms
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InternalLocalAssembly(field *fd, schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the internal flxes of the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FluxAssembly(Simulation *simu,  LinearSolver *solver,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the internal flxes of the fields
//! case of the local implicit solver
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FluxLocalAssembly(field* fd,schnaps_real theta, schnaps_real dt);

//! \brief time-stepping by the crank nicholson scheme
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void ThetaTimeScheme(Simulation *simu, schnaps_real tmax, schnaps_real dt);

//! \brief time-stepping by the crank nicholson scheme
//! macrocell local version
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void LocalThetaTimeScheme(Simulation *simu, schnaps_real tmax, schnaps_real dt);

//! \brief time-stepping by the crank nicholson scheme StarPU version
//! macrocell local version 
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void LocalThetaTimeScheme_SPU(Simulation *simu, schnaps_real tmax, schnaps_real dt);

//! \brief Crank-Nicholson time-stepping for transport solver 
//! graph_based and StarPU version
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void GraphThetaTimeScheme_SPU(Simulation *simu, schnaps_real tmax, schnaps_real dt);

//! \brief Crank-Nicholson time-stepping for transport solver 
//! graph_based and StarPU version without matrix assembly
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void GraphThetaTimeSchemeSubCell_SPU(Simulation *simu, schnaps_real tmax, schnaps_real dt);

///! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the interface fluxes between the neighboring fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InterfaceAssembly(Simulation *simu,  LinearSolver *solver,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the right hand side of the linear system:
//! volume terms and boundary terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void SourceAssembly(Simulation *simu,  LinearSolver *solver,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the right hand side of the linear system:
//! volume terms
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void SourceLocalAssembly(field *fd,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the right hand side of the linear system:
//! volume terms
//! case of the locally implicit scheme
//! StarPU version
//! \param[inout] field a field
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void SourceLocalAssembly_SPU(field *fd,schnaps_real theta, schnaps_real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the mass terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void MassAssembly(Simulation *simu,  LinearSolver *solver);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the mass terms
//! case of the locally implicit scheme
//! \param[inout] field a field
void MassLocalAssembly(field *fd);


//! \brief initialize the struct JFLinear solver for our DG scheme
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a Free jacobian linear solver
void InitImplicitJFLinearSolver(Simulation *simu, JFLinearSolver *solver);

//! \brief Assembly of the DG vectors for the JF solver
//! computations of all the terms
//! \param[inout] simu a simulation
//! \param[inout] solver a jacobian free linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void AssemblyImplicitJFLinearSolver(Simulation *simu, JFLinearSolver *solver, schnaps_real dt);


//! ADD DESCRIPTION
void ThetaTimeScheme_WithJF(Simulation *simu, schnaps_real tmax, schnaps_real dt);

//! \brief function which compute the implicit nonlinear vector for the Jacobian free for the sol solvector
//! \param[i] simu a simulation
//! \param[in] solver a linear solver
//! \param[in] solvector a vector of unknown used to construct the nonlinear vector
//! \param[inout] nlvector a nonlinear vector obtained at the end
void ImplicitNonlinearVector_computation(Simulation * simu,void* lsol,schnaps_real * solvector,schnaps_real *nlvector);

//! \brief function which compute the nonlinear vector for the theta scheme
//! \param[i] simu a simulation
//! \param[in] solvector a vector of unknown used to construct the nonlinear vector
//! \param[inout] nlvector a nonlinear vector obtained at the end
//! \param[in] theta a coefficeny
//! \param[in] dt a time step
void NonlinearThetaVector_assembly(Simulation * simu,schnaps_real * solvector,schnaps_real *nlvector,schnaps_real theta, schnaps_real dt);



#endif

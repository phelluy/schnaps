#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H

#include "global.h"
#include <stdbool.h>
#include <starpu.h>

typedef enum MatrixStorage{SKYLINE,CSR,SKYLINE_SPU,KLU_CSR} MatrixStorage;
typedef enum Solver{LU,GMRES} Solver;
typedef enum PC{NONE,JACOBI} PC;

//! class for managing linear solvers
typedef struct LinearSolver{

  //! \brief number of equations
  int neq;
  //! \brief storage struct for the matrix
  //! the actual type depends on the chosen format
  void* matrix;
  //! name of the storage method;
  MatrixStorage storage_type; 
  //! \brief true if the matrix is symmetric
  bool is_sym;
  //! \brief true if the struct is initialized
  bool is_init;
  //! \brief true if the arrays are allocated
  bool is_alloc;
   //! \brief true if the matrix is assembly
  bool mat_is_assembly;
   //! \brief true if the matrix is assembly
  bool rhs_is_assembly;
  //! solver type;
  Solver solver_type;
  //! solver type;
  bool is_CG;
  //! name of the storage method;
  PC pc_type; 
  //! \brief solution of the linear system
  schnaps_real* sol;
  //! \brief rhs of the linear system
  schnaps_real* rhs;
  //! tolerance iterative solver
  schnaps_real tol;
  //! restart for gmres
  int restart_gmres;
  //! number max of iteration
  int iter_max;

  //! \brief compute a matrix vector product
  //! \param[in] lsol the LinearSolver object containing matrix A
  //! \param[in] x a vector
  //! \param[out] prod Ax
  void (*MatVecProduct)(void* lsol,schnaps_real x[],schnaps_real prod[]);

} LinearSolver;


//! \brief init the LinearSolver structure with an empty matrix
//! \param[inout] lsol the LinearSolver object
//! \param[in] n number of equations
//! \param[in] matstor storage type (optional)
//! \param[in] solvtyp solver type (optional)
void InitLinearSolver(LinearSolver* lsol,int n,
		      MatrixStorage* matstor,
		      Solver* solvtyp);

//! \brief free the allocated arrays
//! \param[inout] lsol the LinearSolver object
//! \param[in] freeAll: an integer whose purpose is to choose whether we free everything. (testing purposes)
void FreeLinearSolver(LinearSolver* lsol);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
void IsNonZero(LinearSolver* lsol,int i,int j); 

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] lsol the LinearSolver object
void AllocateLinearSolver(LinearSolver* lsol);

//! \brief add  to elem (i,j)  value val
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddLinearSolver(LinearSolver* lsol,int i,int j,schnaps_real val);

//! \brief set  to elem (i,j)  value val
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetLinearSolver(LinearSolver* lsol,int i,int j,schnaps_real val); 

//! \brief get elem (i,j)
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
schnaps_real GetLinearSolver(LinearSolver* lsol,int i,int j); 


//! \brief display the matrix
//! \param[inout] lsol the LinearSolver object
void DisplayLinearSolver(LinearSolver* lsol); 

//! \brief compute a matrix vector product
//! \param[in] system the LinearSolver object containing matrix A
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVect(void * system,schnaps_real x[],schnaps_real prod[]);

//! \brief compute a matrix vector product
//! use the rhs and sol starpu data handles
//! \param[in] system the LinearSolver object containing matrix A
//! \param[in] sol_handle starpu handle to the multiplied vector
//! \param[out] rhs_handle starpu handle to the result
void MatVect_SPU(void * system, starpu_data_handle_t sol_handle, starpu_data_handle_t rhs_handle);

//! \brief compute a matrix vector product
//! use the rhs and sol stored inside the structure
//! \param[inout] system the LinearSolver object containing matrix A
void MatVectIn(void * system);

//! \brief compute a matrix vector product
//! \param[in] system the LinearSolver object containing matrix A
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVect_slow(void * system,schnaps_real x[],schnaps_real prod[]);

//! \brief compute the inplace LU decomposition
//! \param[inout] lsol the LinearSolver object
void LUDecompLinearSolver(LinearSolver* lsol);

//! \brief solve the linear system
//! \param[inout] lsol the LinearSolver object
void SolveLinearSolver(LinearSolver* lsol);



//! \brief copy vector
//! \param[in] x vector
//! \param[inout] prod is a copy of x
//! \param[in] N size
void Vector_copy(schnaps_real x[],schnaps_real prod[],int N);

//! \brief return the dot product
//! \param[in] x vector
//! \param[in] y vector
//! \param[in] N size
schnaps_real Vector_prodot(schnaps_real x[],schnaps_real y[],int N);

//! \brief return the l2 norm
//! \param[in] x vector
//! \param[in] N size
schnaps_real Vector_norm2(schnaps_real x[],int  N);


//! \brief solve the linear system with the GMREs of the cerfacs
//! \param[inout] lsol contains the matrices rhs and sol
void GMRESSolver(LinearSolver* lsol);

//! \brief Jacobi preconditioner
//! \param[in] lsol contains the matrices rhs and sol
void Jacobi_PC(LinearSolver* lsol, schnaps_real* sol, schnaps_real* rhs);

#endif

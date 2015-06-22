#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H

#include "global.h"
#include <stdbool.h>


typedef enum MatrixStorage{SKYLINE,CSR} MatrixStorage;
typedef enum Solver{LU,GMRES,PAR_GMRES,PAR_FGMRES,PAR_CG,PAR_BICGSTAB,PAR_AMG,PAR_LU,PAR_QR} Solver;
typedef enum PC{NONE,PAR_JACOBI,PAR_ILU,PAR_MULTICOLOREDSGS,PAR_MULTICOLOREDGS,PAR_MULTICOLOREDILU,PAR_AMG_PC} PC;

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

  //! solver type;
  Solver solver_type;

  //! name of the storage method;
  PC pc_type; 

  //! \brief solution of the linear system
  real* sol;
  //! \brief rhs of the linear system
  real* rhs;


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
void AddLinearSolver(LinearSolver* lsol,int i,int j,real val); 

//! \brief get elem (i,j)
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
real GetLinearSolver(LinearSolver* lsol,int i,int j); 


//! \brief display the matrix
//! \param[inout] lsol the LinearSolver object
void DisplayLinearSolver(LinearSolver* lsol); 

//! \brief compute a matrix vector product
//! \param[in] lsol the LinearSolver object containing matrix A
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVecLinearSolver(LinearSolver* lsol,real x[],real prod[]);

//! \brief compute the inplace LU decomposition
//! \param[inout] lsol the LinearSolver object
void LUDecompLinearSolver(LinearSolver* lsol);

//! \brief solve the linear system
//! \param[in] lsol the LinearSolver object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
void SolveLinearSolver(LinearSolver* lsol);

//! \brief copy vector
//! \param[in] x vector
//! \param[in] N size
//! \param[in] copy x in prod
void Vector_copy(real x[],real prod[],int N);

  //! \brief dot product
//! \param[in] x vector
//! \param[in] y vector
//! \param[in] N size
//! \param[in] prod dot product between x and y
void Vector_prodot(real x[],real y[],real prod[],int N);

//! \brief solve the linear system with paralution
//! \param[in] lsol contains the matrices rhs and sol
void Solver_Paralution(LinearSolver* lsol);




void GMRESSolver(LinearSolver* lsol);


#endif

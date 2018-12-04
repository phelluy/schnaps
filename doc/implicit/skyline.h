#ifndef _SKYLINE_H
#define _SKYLINE_H

#include <stdbool.h>
#include <complex.h>
#define dcmplx  double _Complex



//! \brief a struct for managing skyline linear system
typedef struct Skyline{

  //! \brief number of equations
  int neq;

  //! \brief size in memory
  int nmem;

  //! \brief array for the upper part of  the matrix
  dcmplx* vkgs;

  //! \brief array for the diagonal part of  the matrix (size neq)
  dcmplx* vkgd;

  //! \brief array for the lower part of  the matrix
  dcmplx* vkgi;

  //! \brief array for the upper part of  the matrix (copy used if factorized)
  dcmplx* copy_vkgs;

  //! \brief array for the diagonal part of  the matrix (size neq) (copy used if factorized)
  dcmplx* copy_vkgd;

  //! \brief array for the lower part of  the matrix (copy used if factorized)
  dcmplx* copy_vkgi;

  //! \brief profile of the matrix (size neq)
  int* prof;

  //! \brief array of indices of column start (size neq+1)
  int* kld;

  //! \brief true if the matrix is symmetric
  bool is_sym;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

  //! \brief true if the arrays of the copied matrix are allocated
  bool copy_is_alloc;

  //! \brief true if the matrix is factorized
  bool is_lu;

} Skyline;

//! \brief init the skyline structure with an empty matrix
//! \param[inout] sky the skyline object
//! \param[in] n number of equations
void InitSkyline(Skyline* sky,int n);

//! \brief free the allocated arrays
//! \param[inout] sky the skyline object
void FreeSkyline(Skyline* sky);

//! \brief fill allocated arrays with zeros
//! \param[inout] sky the skyline object
void ZeroSkyline(Skyline* sky);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
void SwitchOn(Skyline* sky,int i,int j); 

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the skyline object
void AllocateSkyline(Skyline* sky);

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the skyline object
void AllocateCopySkyline(Skyline* sky);

//! \brief set elem (i,j) to value val
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetSkyline(Skyline* sky,int i,int j,dcmplx val);

//! \brief add elem (i,j) to value val
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddSkyline(Skyline* sky,int i,int j,dcmplx val); 

//! \brief get elem (i,j)
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
dcmplx GetSkyline(Skyline* sky,int i,int j); 


//! \brief display the matrix
//! \param[inout] sky the skyline object
void DisplaySkyline(Skyline* sky);

//! \brief compute a matrix vector product
//! \param[in] sky a skyline matrix
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVectSkyline(Skyline * sky, dcmplx * x, dcmplx * prod);


//! \brief compute the inplace LU decomposition
//! \param[inout] sky the skyline object
void FactoLU(Skyline* sky);

//! \brief solve the linear system
//! \param[in] sky the skyline object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
void SolveSkyline(Skyline* sky,dcmplx* rhs,dcmplx* sol);

//! \brief solve the linear system
//! \param[in] sky the skyline object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
void FastSolveSkyline(Skyline* sky,dcmplx* rhs,dcmplx* sol);

//! \brief compute the inverse of a square matrix
//! \param[in] n size of the matrix
//! \param[in] a the square matrix ( aij = a[i*n+j] )
//! \param[out] b inverse of a
void InvertSquare(int n, dcmplx *a, dcmplx *b);

//! \brief compute the PLU decomposition of a square matrix
//! with partial row pivoting (P is a permutation matrix)
//! \param[in] n size of the matrix
//! \param[in] a the square matrix ( aij = a[i*n+j] )
//! \param[out] b contains LU decomposition of a
//! \param[out] sigma contains permutation of rows of a
void PLU_Square(int n, dcmplx *a, int *sigma);

//! \brief use the PLU decomposition for solving a linear system
//! \param[in] n size of the matrix
//! \param[in] a contains LU decomposition of a
//! \param[in] sigma contains permutation of rows of a
//! \param[in] b right hand side
//! \param[out] x solution
void PLU_Solve(int n, dcmplx *a, int *sigma, dcmplx *b,
		 dcmplx *x);

//! \brief small test of PLU functions
void TestPLU(void);


#endif


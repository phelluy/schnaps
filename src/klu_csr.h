#ifndef _KLU_CSR_H
#define _KLU_CSR_H
#include "global.h"
#include "cs.h"
#include <stdbool.h>
#include "klu.h"


//! \brief a struct for managing KLU linear system
typedef struct KLU {

  //! \brief number of equations
  int neq;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

  //! \brief true if the arrays of the copied matrix are allocated
  bool copy_is_alloc;

  //! \brief true if the matrix is factorized
  bool is_lu;

  //! \brief symbolic and numeric objects
  klu_symbolic *symbolic;
  klu_numeric *numeric;
  klu_common common;

  //! csparse struct for triplet storage
  cs_di *T;

  //! csparse struct for csr storage
  cs_di *A;

  //! csparse struct for csr copy storage
  cs_di *Acopy;
  
  //! reverse line lookup table for symbolic permutation 
  //int *symb_Qinv;
  
  //! reverse column lookup table for symbolic permutation 
  //int *symb_Pinv;
  //!
  //! display switches
  //! the matrix in its original order
  bool display_orig;
  //! true : print actual values in A false: put symbols in nonzero entries 
  bool display_orig_values;
  // show klu block structure obtained from symbolic analysis (only with * symbols and blocks bounds)
  bool display_klu_blocks;
  //!
} KLU;

//! \brief init the KLU structure with an empty matrix
//! \param[inout] klumat the KLU object
//! \param[in] n number of equations
void InitKLU(KLU * klumat, int n);

//! \brief free the allocated arrays
//! \param[inout] klumat the KLU object
void FreeKLU(KLU * klumat);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] klumat the KLU object
//! \param[in] i row index
//! \param[in] j column index
void SwitchOnKLU(KLU * klumat, int i, int j);

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] klumat the KLU object
void AllocateKLU(KLU * klumat);

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] klumat the KLU object
void AllocateCopyKLU(KLU * klumat);

//! \brief set elem (i,j) to value val
//! \param[inout] klumat the KLU object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetKLU(KLU * klumat, int i, int j, schnaps_real val);

//! \brief add elem (i,j) to value val
//! \param[inout] klumat the KLU object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddKLU(KLU * klumat, int i, int j, schnaps_real val);

//! \brief get elem (i,j)
//! \param[inout] klumat the KLU object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
schnaps_real GetKLU(KLU * klumat, int i, int j);


//! \brief display the matrix
//! \param[inout] klumat the KLU object
void DisplayKLU(KLU * klumat);

//! \brief compute a matrix vector product
//! \param[in] klumat a KLU matrix
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVectKLU(KLU * klumat, schnaps_real * x, schnaps_real * prod);


//! \brief compute the inplace LU decomposition
//! \param[inout] klumat the KLU object
void FactoKLU(KLU * klumat);

//! \brief compute the inplace LU decomposition - 1 block restriction
//! \param[inout] klumat the KLU object
void FactoKLU_Blockwise(KLU * klumat, const int iblock);

//! \brief recompute the inplace LU decomposition
//! \param[inout] klumat the KLU object
void ReFactoKLU(KLU * klumat);

//! \brief recompute the inplace LU decomposition - 1 block restriction
//! \param[inout] klumat the KLU object
//! \param[in] the block
void ReFactoKLU_Blockwise(KLU * klumat,const int iblock);

//! \brief solve the linear system
//! \param[in] klumat the KLU object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
void SolveKLU(KLU * klumat, schnaps_real * rhs, schnaps_real * sol);

//! \brief solve the linear system - restriction to one block
//! \param[in] klumat the KLU object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
//! \param[in] iblock block index
void SolveKLU_Blockwise(KLU * klumat, schnaps_real * rhs, schnaps_real * sol,const int iblock);

#endif

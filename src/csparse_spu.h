#ifndef _CSPARSESPU_H
#define _CSPARSESPU_H

#include "global.h"
#include <starpu.h>
#include <csparse.h>


//! \brief a struct for managing csparse linear system
//! wrapper and StarPU mem management
typedef struct CSparse_SPU{

  //! \brief number of equations
  int neq;

  //! \brief size in memory
  int nmem;

  //! \brief rhs of the linear solver
  schnaps_real* rhs;
  starpu_data_handle_t rhs_handle;
  
  //! \brief solution of the linear solver
  schnaps_real* sol;
  starpu_data_handle_t sol_handle;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

  //! \brief true if the matrix is factorized
  bool is_lu;

  //! \brief true if the arrays are registerd in StarPU
  bool is_registered;

} CSparse_SPU;

//! \brief init the csparse structure with an empty matrix
//! \param[inout] spm the csparse object
//! \param[in] n number of equations
void Init_CSparse_SPU(CSparse_SPU* spm,int n);

//! \brief register the csparse arrays in starPU
//! \param[inout] spm the csparse object
void Register_CSparse_SPU(CSparse_SPU* spm);

//! \brief unregister the csparse arrays in starPU
//! \param[inout] spm the csparse object
void UnRegister_CSparse_SPU(CSparse_SPU* spm);

//! \brief free the allocated arrays
//! \param[inout] spm the csparse object
void Free_CSparse_SPU(CSparse_SPU* spm);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] spm the csparse object
//! \param[in] i row index
//! \param[in] j column index
void SwitchOn_SPU(CSparse_SPU* spm,int i,int j); 

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] spm the csparse object
void Allocate_CSparse_SPU(CSparse_SPU* spm);

//! \brief set elem (i,j) to value val
//! \param[inout] spm the csparse object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void Set_CSparse_SPU(CSparse_SPU* spm,int i,int j,schnaps_real val);

//! \brief add elem (i,j) to value val
//! \param[inout] spm the csparse object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void Add_CSparse_SPU(CSparse_SPU* spm,int i,int j,schnaps_real val); 

//! \brief get elem (i,j)
//! \param[inout] spm the csparse object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
schnaps_real Get_CSparse_SPU(CSparse_SPU* spm,int i,int j); 


//! \brief display the matrix
//! \param[inout] spm the csparse object
void Display_CSparse_SPU(CSparse_SPU* spm);


//! \brief compute the inplace LU decomposition
//! \param[inout] spm the csparse object
void LU_CSparse_SPU(CSparse_SPU* spm);

//! \brief solve the linear system
//! \param[in] spm the csparse object
void Solve_CSparse_SPU(CSparse_SPU* spm);

//! \brief compute a matrix vector product
//! \param[in] spm a csparse matrix
//! \param[in] sol_handle starpu handle to the multiplied vector
//! \param[out] rhs_handle starpu handle to the result
void MatVect_CSparse_SPU(CSparse_SPU * spm,
			starpu_data_handle_t sol_handle,
			starpu_data_handle_t rhs_handle);


#endif

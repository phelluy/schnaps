#ifndef _SKYLINESPU_H
#define _SKYLINESPU_H

#include "global.h"
#include <starpu.h>


//! \brief a struct for managing skyline linear system
//! StarPU version
typedef struct Skyline_SPU{

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

  //! \brief array for the upper part of  the matrix
  schnaps_real* vkgs;
  starpu_data_handle_t vkgs_handle;

  //! \brief array for the diagonal part of  the matrix (size neq)
  schnaps_real* vkgd;
  starpu_data_handle_t vkgd_handle;

  //! \brief array for the lower part of  the matrix
  schnaps_real* vkgi;
  starpu_data_handle_t vkgi_handle;

  //! \brief profile of the matrix (size neq)
  int* prof;

  //! \brief array of indices of column start (size neq+1)
  int* kld;
  starpu_data_handle_t kld_handle;

  //! \brief true if the matrix is symmetric
  bool is_sym;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

  //! \brief true if the matrix is factorized
  bool is_lu;

  //! \brief true if the arrays are registerd in StarPU
  bool is_registered;

} Skyline_SPU;

//! \brief init the skyline structure with an empty matrix
//! \param[inout] sky the skyline object
//! \param[in] n number of equations
void InitSkyline_SPU(Skyline_SPU* sky,int n);

//! \brief register the skyline arrays in starPU
//! \param[inout] sky the skyline object
void RegisterSkyline_SPU(Skyline_SPU* sky);

//! \brief unregister the skyline arrays in starPU
//! \param[inout] sky the skyline object
void UnRegisterSkyline_SPU(Skyline_SPU* sky);

//! \brief free the allocated arrays
//! \param[inout] sky the skyline object
void FreeSkyline_SPU(Skyline_SPU* sky);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
void SwitchOn_SPU(Skyline_SPU* sky,int i,int j); 

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the skyline object
void AllocateSkyline_SPU(Skyline_SPU* sky);

//! \brief set elem (i,j) to value val
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetSkyline_SPU(Skyline_SPU* sky,int i,int j,schnaps_real val);

//! \brief add elem (i,j) to value val
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddSkyline_SPU(Skyline_SPU* sky,int i,int j,schnaps_real val); 

//! \brief get elem (i,j)
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
schnaps_real GetSkyline_SPU(Skyline_SPU* sky,int i,int j); 


//! \brief display the matrix
//! \param[inout] sky the skyline object
void DisplaySkyline_SPU(Skyline_SPU* sky);


//! \brief compute the inplace LU decomposition
//! \param[inout] sky the skyline object
void FactoLU_SPU(Skyline_SPU* sky);

//! \brief solve the linear system
//! \param[in] sky the skyline object
void SolveSkyline_SPU(Skyline_SPU* sky);

//! \brief compute a matrix vector product
//! \param[in] sky a skyline matrix
//! \param[in] sol_handle starpu handle to the multiplied vector
//! \param[out] rhs_handle starpu handle to the result
void MatVectSkyline_SPU(Skyline_SPU * sky,
			starpu_data_handle_t sol_handle,
			starpu_data_handle_t rhs_handle);


#endif

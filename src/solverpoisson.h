#ifndef _SOLVERPOISSON_H
#define _SOLVERPOISSON_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "skyline.h"

#define _Dirichlet_Poisson_BC (1)
#define _Periodic_Poisson_BC (2)


//! \brief a struct for sorting and pasting the
//! nodes of the DG mesh for obtaining a FE mesh
typedef struct FatNode{

  //! \brief index in the dg mesh
  int dg_index;
  //! \brief index in the fe mesh
  int fe_index;

  //! \brief physical coordinates of the node
  real x[3];
  //! \brief int converted coordinates for sorting and searching
  int x_int[3];


} FatNode;

//! \brief a struct for the resolution of the poisson equation:
//! conversion between a DG and Finite Element FE mesh
//! FE assembly and resolution, etc. 
typedef struct PoissonSolver{

  //! \brief a field (gives the mesh and the charge)
  field* fd;

  //! \brief charge index in the conservative variables vector
  int charge_index;

  //! \brief vector containing the charge (right hand side)
  real* rhs;

  //! \brief vector containing the potential (solution)
  real* sol;

  //! \brief number of FE nodes
  int nb_fe_nodes;

  //! \brief number of DG nodes
  int nb_dg_nodes;

  //! \brief node list with coordinates converted to integer
  //! and DG/FE indices
  FatNode* fn_list;

  //! \brief connectivity DG node -> FE node
  int* dg_to_fe_index;

  //! \brief list that marks boundary nodes
  int* is_boundary_node;

} PoissonSolver;



//! \brief compare two nodes (used by quicksort).
//! Lexicographic order on the coordinates converted to integers.
//! \param[in] a first node
//! \param[in] b second node
//! \returns a value v, v<0 if a<b, v=0 if a==b, v>0 if a>b
int CompareFatNode(const void* a,const void* b);

//! \brief build the fat nodes list from a field
//! \param[in] a initialized field
//! \param[out] an allocated, prepared and sorted list of fat nodes
//! \returns the size of the list
int BuildFatNodeList(field* f,FatNode* fn_list);

//! \brief init a poisson solver
//! \param[inout] ps a PoissonSolver struct
//! \param[in] fd a Field
//! \param[in] charge_index charge index in the field variables
void InitPoissonSolver(PoissonSolver* ps, field* fd,int charge_index);

//! \brief solve a 1D poisson problem
//! \param[in] f a field (contains the mesh)
//! \param[in] w the field values (for computing the charge
//! , returning the potential and the electric field)
//! \param[in] type_bc the boundary condition type
//!  (1->dirichlet ; 2-> periodic)
//! \param[in] bc_l left boundary value (dirichlet case)
//! \param[in] bc_r right boundary value (dirichlet case)
void SolvePoisson1D(field *f,real * w,
		    int type_bc, real bc_l, real bc_r);



//! \brief solve a 2D poisson problem
//! \param[in] f a field (contains the mesh)
//! \param[in] w the field values (for computing the charge
//! , returning the potential and the electric field)
//! \param[in] type_bc the boundary condition type
//!  (1->dirichlet ; 2-> periodic)
void SolvePoisson2D(PoissonSolver* ps,int type_bc);

#endif

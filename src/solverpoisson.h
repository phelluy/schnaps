#ifndef _SOLVERPOISSON_H
#define _SOLVERPOISSON_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "skyline.h"

//! \brief a struct for the resolution of the poisson equation:
//! conversion between a DG and Finite Element FE mesh
//! FE assembly and resolution, etc. 
typedef struct PoissonSolver{

  //! \brief a field (gives the mesh and the charge)
  Field* fd;

  //! \brief charge index in the conservative variables vector
  int charge_index;

  //! \brief number of FE nodes
  int nb_fe_nodes;

  //! \brief number of DG nodes
  int nb_dg_nodes;

  //! \brief connectivity DG node -> FE node
  int* dg_to_fe_index;


} PoissonSolver;




void SolvePoisson(field *f,real * w,int type_bc, real bc_l, real bc_r);

#endif

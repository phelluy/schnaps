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
void SolvePoisson2D(field *f,real * w,int type_bc);

#endif

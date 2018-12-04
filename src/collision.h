#ifndef _COLLISION_H
#define _COLLISION_H

//! \brief number of conservative variables
//! values of  the distribution function at the velocity glops
//! and value of the potential 

#include "model.h"
#include "field.h"
// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void VlasovP_Lagrangian_NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);

//! \brief  compute the source term of the collision
//! model: electric force + true collisions
//! \param[in] x space position
//! \param[in] t time
//! \param[in] w the distribution function
//! \param[out] source the source
void VlasovP_Lagrangian_Source(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);

void BGK_Source(const schnaps_real* w, schnaps_real* source);
		    
#endif

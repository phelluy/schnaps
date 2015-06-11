#ifndef _COLLISION_H
#define _COLLISION_H

#define _NB_ELEM_V 2
#define _DEG_V 1

//! \brief number of conservative variables
//! values of  the distribution function at the velocity glops
//! and value of the potential 
#define _MV (_NB_ELEM_V *  _DEG_V + 1 )

#define _INDEX_MAX_KIN (_MV-1)
#define _INDEX_PHI (_MV)
#define _INDEX_EX (_MV+1)
#define _INDEX_RHO (_MV+2)
#define _INDEX_VELOCITY (_MV+3)
#define _INDEX_PRESSURE (_MV+4)
#define _INDEX_TEMP (_MV+5)
#define _INDEX_MAX (_MV+6)

#define _VMAX 6.
#define _DV (2*_VMAX / _NB_ELEM_V)

#include "model.h"
#include "field.h"
// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void VlasovP_Lagrangian_NumFlux(real *wL, real *wR, real *vn, real *flux);

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
//! \param[in] w the distribution function
//! \param[in] f the force 
//! \param[out] source the source
void VlasovP_Lagrangian_Source(const real *x, const real t, const real *w, 
			       real *source);


#endif

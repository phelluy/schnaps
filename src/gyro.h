#ifndef _GYRO_H
#define _GYRO_H

#define _NB_ELEM_V 4
#define _DEG_V 1

#define _MV (_NB_ELEM_V *  _DEG_V + 1)
#define _VMAX 6
#define _DV (2*_VMAX / _NB_ELEM_V)


#include "model.h"
#include "field.h"
// gyrokinetic models

//! \brief particular flux for the gyro model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Gyro_Lagrangian_NumFlux(real wL[],real wR[],real vn[3],real* flux);

//! \brief particular boundary flux for the gyro model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Gyro_Lagrangian_BoundaryFlux(real* x,real t,real* wL,real* vn,
			   real* flux);

//! \brief particular init data for the gyro model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void GyroInitData(real* x,real* w);

//! \brief particular imposed data for the  gyro model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void GyroImposedData(real* x,real t,real* w);

//! \brief particular imposed data for the  gyro model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
real Gyro_ImposedKinetic_Data(real* x,real t,real v);

//! \brief compute gyro L2 error in x and v
//! \param[in] f : a field
real GyroL2_Kinetic_error(field* f);

//! \brief compute square of velocity L2 error
//! \param[in] x,t : space and time position
//! \param[in] w : values of f at glops
real GyroL2VelError(real* x,real t,real *w);

//! \brief compute compute the source term of the gyro
//! model: electric force + true gyros
//! \param[in] w : 
//! \param[in] f : to be removed TODO !!!
//! \param[out] source : source terms
void GyroSource(real* force, real* w, real* source);

#endif

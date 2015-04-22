#ifndef _MAXWELL2D_H
#define _MAXWELL2D_H

#include "model.h"

//! \brief The particular flux for the maxwell2d model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Maxwell2DNumFlux(double wL[], double wR[], double vn[3], double* flux);


//! \brief The particular imposed data for the maxwell2d model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void Maxwell2DImposedData(double* x, double t, double* w);



//! \brief The particular boundary flux for the maxwell2d model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Maxwell2DBoundaryFlux(double* x, double t, double* wL, double* vn,
		       double* flux);


//! \brief The particular init data for the maxwell2d model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void Maxwell2DInitData(double* x, double* w);




#endif

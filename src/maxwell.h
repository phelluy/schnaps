#ifndef _MAXWELL_H
#define _MAXWELL_H

#include "model.h"

//! \brief The particular flux for the maxwell2d model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void Maxwell2DNumFlux_uncentered(real *wL, real *wR, real *vn, real *flux);
#pragma end_opencl

//! \brief The particular flux for the maxwell2d model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void Maxwell2DNumFlux_centered(real *wL, real *wR, real *vn, real *flux);
#pragma end_opencl

void Maxwell2DNumFlux_unoptimised(real *wL, real *wR, real *vnorm, real *flux);
//! \brief The particular imposed data for the maxwell2d model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void Maxwell2DImposedData(const real * x, const real t, real *w);
#pragma end_opencl

//! \brief The particular boundary flux for the maxwell2d model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void Maxwell2DBoundaryFlux_uncentered(real *x, real t, real *wL, 
				    real *vn, real *flux);
#pragma end_opencl

//! \brief The particular init data for the maxwell2d model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void Maxwell2DInitData(real *x, real *w);

#pragma start_opencl
void Maxwell2DSource(const real *x, const real t, const real *w, real *source);
#pragma end_opencl



#pragma start_opencl
void Maxwell3DNumFlux_uncentered(real *wL, real *wR, real *vnorm, real *flux);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DNumFluxClean_uncentered(real *wL, real *wR, real *vnorm, 
					real *flux);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DImposedData(const real *x, const real t, real *w);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DInitData(real *x, real *w);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DBoundaryFlux_uncentered(real *x, real t, 
					real *wL, real *vnorm, real *flux);
#pragma end_opencl

#endif

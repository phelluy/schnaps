#ifndef _MAXWELL_H
#define _MAXWELL_H

#include "model.h"

//! \brief The particular flux for the maxwell2d model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void Maxwell2DNumFlux_upwind(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

//! \brief The particular flux for the maxwell2d model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void Maxwell2DNumFlux_centered(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

void Maxwell2DNumFlux_unoptimised(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux);
//! \brief The particular imposed data for the maxwell2d model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void Maxwell2DImposedData(const schnaps_real * x, const schnaps_real t, schnaps_real *w);
#pragma end_opencl

//! \brief The particular boundary flux for the maxwell2d model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void Maxwell2DBoundaryFlux_upwind(schnaps_real *x, schnaps_real t, schnaps_real *wL,
				    schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

//! \brief The particular init data for the maxwell2d model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void Maxwell2DInitData(schnaps_real *x, schnaps_real *w);

#pragma start_opencl
void Maxwell2DSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);
#pragma end_opencl



#pragma start_opencl
void Maxwell3DNumFlux_upwind(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DNumFluxClean_upwind(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm,
					schnaps_real *flux);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DInitData(schnaps_real *x, schnaps_real *w);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DBoundaryFlux_upwind(schnaps_real *x, schnaps_real t,
					schnaps_real *wL, schnaps_real *vnorm, schnaps_real *flux);
#pragma end_opencl

#pragma start_opencl
void Maxwell3DBoundaryFluxClean_upwind(schnaps_real *x, schnaps_real t,
					schnaps_real *wL, schnaps_real *vnorm, schnaps_real *flux);
#pragma end_opencl

#endif

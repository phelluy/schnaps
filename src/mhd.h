#ifndef _MHD_H
#define _MHD_H
#include "model.h"

#pragma start_opencl
//! \brief computes the conservatives states from the primitives
//! \param[in] y : primitives states
//! \param[out] w : conservatives states
void conservatives(real *y, real *w);

//! \brief Numerical flux for the MHD model
//! \param[in] w : states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void fluxnum(real *w, real *vn, real *flux);

//! \brief particular flux for the MHD model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void MHDNumFlux(real *wL, real *wR, real *vn, real *flux);

//! \brief particular boundary flux for the MHD model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void MHDBoundaryFlux(real *x, real t, real *wL, real *vn, real *flux);

//! \brief particular init data for the MHD model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void MHDInitData(real *x, real *w);

//! \brief particular imposed data for the MHD model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void MHDImposedData(const real *x, const real t, real *w);

void primitives(real *W, real *Y);

// FIXME: [][] breaks OpenCL.  Please use * instead.
void jacobmhd(real* W,real* vn, real M[9][9]);
void matrix_vector(real A[9][9], real B[9], real* C);
void matrix_matrix(real A[9][9],real B[9][9],real C[9][9]);
void write_matrix(real A[9][9],real *second, real B[9][9+1]);
void gauss(real A[9][9], real b[9], real *x);
void MHDNumFlux_2(real *wL, real *wR, real *vn, real *flux);

void MHDNumFlux1D(real wL[],real wR[],real* vn, real* flux);
#pragma end_opencl

#endif

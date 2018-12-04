#ifndef _MHD_H
#define _MHD_H
#include "model.h"

//! \brief computes the conservatives states from the primitives
//! \param[in] y : primitives states
//! \param[out] w : conservatives states
#pragma start_opencl
void conservatives(schnaps_real *y, schnaps_real *w);
void primitives(schnaps_real *W, schnaps_real *Y);
#pragma end_opencl


//! \brief Numerical flux for the MHD model
//! \param[in] w : states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void fluxnum(schnaps_real *w, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

//! \brief particular flux for the MHD model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void MHDNumFluxRusanov(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
void MHDNumFluxP2(schnaps_real *wL,schnaps_real *wR,schnaps_real *vn, schnaps_real *flux);
void MHDNumFlux1D(schnaps_real *wL,schnaps_real *wR,schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

//! \brief particular boundary flux for the MHD model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void MHDBoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);
void MHDBoundaryFluxChocFort(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);
void MHDBoundaryFluxOrszagTang(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);
void MHDBoundaryFluxReconnexion(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);
void MHDBoundaryFluxKelvinHelmotz(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);
void MHDBoundaryFluxDoubleTearing(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

//! \brief particular init data for the MHD model
//! \param[in] x : space position
//! \param[out] w : init state at point x
#pragma start_opencl
void MHDInitData(schnaps_real *x, schnaps_real *w);
void MHDInitDataChocFort(schnaps_real *x, schnaps_real *w);
void MHDInitDataOrszagTang(schnaps_real *x, schnaps_real *w);
void MHDInitDataReconnexion(schnaps_real *x, schnaps_real *w);
void MHDInitDataKelvinHelmotz(schnaps_real *x, schnaps_real *w);
void MHDInitDataDoubleTearing(schnaps_real *x, schnaps_real *w);
#pragma end_opencl

//! \brief particular imposed data for the MHD model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void MHDImposedData(const schnaps_real *x,const schnaps_real t, schnaps_real *w);
void MHDImposedDataChocFort(const schnaps_real *x,const schnaps_real t, schnaps_real *w);
void MHDImposedDataOrszagTang(const schnaps_real *x,const  schnaps_real t, schnaps_real *w);
void MHDImposedDataReconnexion(const schnaps_real *x,const  schnaps_real t, schnaps_real *w);
void MHDImposedDataKelvinHelmotz(const schnaps_real *x,const  schnaps_real t, schnaps_real *w);
void MHDImposedDataDoubleTearing(const schnaps_real *x,const  schnaps_real t, schnaps_real *w);
#pragma end_opencl

// FIXME: using "real var[]" instead of "real *var" breaks OpenCL on
// certain platforms.
#pragma start_opencl
void jacobmhd(schnaps_real* W,schnaps_real* vn, schnaps_real *M);
void matrix_vector(schnaps_real *A, schnaps_real *B, schnaps_real* C);
#pragma end_opencl



#endif

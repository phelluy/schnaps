#ifndef _GYRO_H
#define _GYRO_H

// FIXME:

// The use of defines to set varialbes is going to cause massive
// problems with maintainability.  Please stop doing this.

/* #define _NB_ELEM_V 1 */
/* #define _DEG_V 1 */

/* #define _MV (_NB_ELEM_V *  _DEG_V + 1) */
/* #define _INDEX_MAX_KIN (_MV-1) */
/* #define _INDEX_PHI (_MV) */
/* #define _INDEX_EX (_MV+1) */
/* #define _INDEX_EY (_MV+2) */
/* #define _INDEX_EZ (_MV+3) */
/* #define _INDEX_MAX (_MV+4) */
/* #define _VMAX 6. */
/* #define _DV (2*_VMAX / _NB_ELEM_V) */


#define OCL 1
#include "model.h"
#include "field.h"
#include "diagnostics_vp.h"
// gyrokinetic models

//! \brief compute the maximal advection velocity
//! \param[inout] simu : the simulation
void GyroCFLVelocity(Simulation* simu);

//! \brief particular flux for the gyro model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void GyroUpwindNumFlux(schnaps_real wL[],schnaps_real wR[],
			       schnaps_real* vnorm,schnaps_real* flux);
#pragma end_opencl

#pragma start_opencl
void GC_OCLUpwindNumFlux(schnaps_real wL[],schnaps_real wR[],
			       schnaps_real* vnorm,schnaps_real* flux);
#pragma end_opencl

//! \brief particular flux for the gyro model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void GyroCenteredNumFlux(schnaps_real wL[],schnaps_real wR[],
			      schnaps_real* vnorm,schnaps_real* flux);
#pragma end_opencl

//! \brief particular flux for the gyro model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void GyroZeroNumFlux(schnaps_real wL[],schnaps_real wR[],
			     schnaps_real* vnorm,schnaps_real* flux);
#pragma end_opencl

//! \brief particular boundary flux for the gyro model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void GyroBoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux);
#pragma end_opencl

#pragma start_opencl
void ChargeOCLBoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux);
#pragma end_opencl

#pragma start_opencl
void GC_OCLBoundaryFlux(schnaps_real x[3],schnaps_real t,
		      schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux);
#pragma end_opencl

//! \brief particular init data for the gyro model
//! \param[in] x : space position
//! \param[out] w : init state at point x
#pragma start_opencl
void GyroInitData(schnaps_real x[3],schnaps_real w[]);
#pragma end_opencl


#pragma start_opencl
void ChargeOCLInitData(schnaps_real x[3],schnaps_real w[]);
#pragma end_opencl

//! GC_OCL : Guiding center OpenCL test case
#pragma start_opencl
void GC_OCLInitData(schnaps_real x[3],schnaps_real w[]);
#pragma end_opencl

//! \brief particular imposed data for the  gyro model
//! \param[in] x  space position
//! \param[in] t time
//! \param[out] w  init state at point x
#pragma start_opencl
void GyroImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[]);
#pragma end_opencl

#pragma start_opencl
void ChargeOCLImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[]);
#pragma end_opencl

#pragma start_opencl
void GC_OCLImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[]);
#pragma end_opencl

//! \brief particular imposed data for the  gyro model
//! \param[in] x  space
//! \param[in] t time
//! \param[in] v  velocity
//! \returns value of the distribution function
schnaps_real GyroImposedKineticData(const schnaps_real x[3], const schnaps_real t, schnaps_real v);


//! \brief compute gyro L2 error in x and v
//! \param[in] f : a field
schnaps_real GyroL2_Kinetic_error(field* f);

//! \brief compute square of velocity L2 error
//! \param[in] x,t : space and time position
//! \param[in] w : values of f at glops
schnaps_real GyroL2VelError(schnaps_real* x,schnaps_real t,schnaps_real *w);

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void GyroSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);

//! \brief gnuplot file for the distribution function
//! \param[in] w : values of f at glops
void Velocity_distribution_plot(schnaps_real *w);


//! \brief charge computation, quasineutrality and electric field
//! \param[inout] si : a simulation (as a void pointer)
//! \param[inout] w : a field
void UpdateGyroPoisson(void *si);

void PlotVP(void* field, schnaps_real * w);

#endif

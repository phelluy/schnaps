#ifndef _MODEL_H
#define _MODEL_H

#include "global.h"

// #include "user_model.h"

// typedefs for function pointers
// Numerical flux function pointer
typedef void (*fluxptr)(schnaps_real*, schnaps_real*, schnaps_real*, schnaps_real*);
// Boundary flux function pointer
typedef void (*bfluxptr)(schnaps_real*, schnaps_real, schnaps_real*, schnaps_real*, schnaps_real*);
// Init data function pointer
typedef void (*initdataptr)(schnaps_real*, schnaps_real*);
// Imposed data function pointer
typedef void (*imposeddataptr)(const schnaps_real*, const schnaps_real, schnaps_real* w); 
// Source function pointer
typedef void (*sourceptr)(const schnaps_real*, const schnaps_real, const schnaps_real*, schnaps_real *); 
// Source jacobian function pointer
typedef void (*sourcejacptr)(const schnaps_real*, const schnaps_real, const schnaps_real*, schnaps_real *); 
// Return a pointer to a numflux function based on the string given.
fluxptr numflux(const char *numfluxname);

// Return a pointer to a boundary blux function based on the string given.
bfluxptr bflux(const char *bfluxname);

// Return a pointer to an init data function based on the string given.
initdataptr initdata(const char *name);

// Return a pointer to an imposed data function based on the string given.
imposeddataptr imposeddata(const char *name);

//! \brief a unified framework for all physical models
typedef struct Model {
  //! Model name
  char name[50];

  //! Number of conservative variables
  int m;

  //! CFL coefficient for time-stepping
  schnaps_real cfl; // 0.05

  //! Number of conservative variables in each dimension (NB: their
  //! product must equal m).
  int vlasov_mx, vlasov_my, vlasov_mz;

  //! The conservative variables have velocity in [-vlasov_vmax, vlasov_vmax]
  schnaps_real vlasov_vmax;

  //! \brief A pointer to the numflux function
  //! \param[in] wL, wR : left and right states
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*NumFlux)(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);

  //! \brief A pointer to the boundary flux function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] wL : left state
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*BoundaryFlux)(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux);

   //! \brief a pointer to the source function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
  void (*Source)(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);

  //! \brief a pointer to the source jacobian function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] jacsource : the jacobian as an m*m 1d array dS_i/dw_j=jacsource[j*q+i]  
  void (*SourceJac)(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *sourcejac);
  
  //! \brief A pointer to the init data function
  // !\param[in] x : space position
  //! \param[out] w : init state at point x
  void (*InitData)(schnaps_real *x, schnaps_real *w);

  //! \brief A pointer to the imposed data function
  //!\param[in] x, t : space and time position
  //! \param[out] w : imposed state at point x and time t
  void (*ImposedData)(const schnaps_real *x, const schnaps_real t, schnaps_real *w);

} Model;

//! \brief The particular flux for the transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TransNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl
//! \brief The particular flux for the 2d transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TransNumFlux2d(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

#pragma start_opencl
void VecTransNumFlux2d(__private schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
#pragma end_opencl

//! \brief The particular boundary flux for the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TransBoundaryFlux(schnaps_real* x, schnaps_real t, schnaps_real* wL, schnaps_real* vn, schnaps_real* flux);
#pragma end_opencl

//! \brief The particular boundary flux for the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TransBoundaryFlux2d(schnaps_real* x, schnaps_real t, schnaps_real* wL, schnaps_real* vn, schnaps_real* flux);
#pragma end_opencl

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
#pragma start_opencl
void TransInitData(schnaps_real* x, schnaps_real* w);
#pragma end_opencl
//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData2d(schnaps_real *x, schnaps_real *w);

void VecTransInitData2d(schnaps_real *x, schnaps_real *w);

//void vTransImposedData2d(real *x, real t, real *w);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void TransImposedData(const schnaps_real* x, const schnaps_real t, schnaps_real* w);
#pragma end_opencl

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void TransImposedData2d(const schnaps_real *x, const schnaps_real t, schnaps_real *w);
#pragma end_opencl

#pragma start_opencl
void VecTransImposedData2d(const schnaps_real* x, const schnaps_real t, schnaps_real *w);
#pragma end_opencl

//! \brief The particular flux for testing the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TestTransBoundaryFlux(schnaps_real* x, schnaps_real t, schnaps_real* wL, schnaps_real* vn,
			   schnaps_real* flux);
#pragma end_opencl

//! \brief The particular flux for testing the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux2d(schnaps_real* x, schnaps_real t, schnaps_real* wL, schnaps_real* vn,
			     schnaps_real* flux);

#pragma start_opencl
void VecTransBoundaryFlux2d(schnaps_real* x, schnaps_real t, schnaps_real* wL, schnaps_real* vn,
			    schnaps_real* flux);
#pragma end_opencl

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData(schnaps_real* x, schnaps_real* w);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData2d(schnaps_real* x, schnaps_real* w);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void TestTransImposedData(const schnaps_real* x, const schnaps_real t, schnaps_real* w);
#pragma end_opencl

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData2d(const schnaps_real* x, const schnaps_real t, schnaps_real* w);

//! Parameters used for computing the velocity for the Vlasov equation.
static int m;
static int vlasov_mx;
static int vlasov_my;
static int vlasov_mz;
static schnaps_real vlasov_vmax;

void set_global_m(int m0);
void set_vlasov_params(Model *mod);
schnaps_real vlasov_vel(const int id, const int md, schnaps_real vlasov_vmax);
void vlaTransInitData2d(schnaps_real *x, schnaps_real *w);
void vlaTransNumFlux2d(schnaps_real *wL, schnaps_real *wR, schnaps_real* vnorm, schnaps_real* flux);
void vlaTransBoundaryFlux2d(schnaps_real *x, schnaps_real t, 
			    schnaps_real *wL, schnaps_real* vnorm,
			    schnaps_real* flux);
void vlaTransImposedData2d(const schnaps_real *x, const schnaps_real t, schnaps_real* w);
schnaps_real compact_bump(schnaps_real r);
schnaps_real icgaussian(schnaps_real r, schnaps_real c);

void cemracs2014_TransBoundaryFlux(schnaps_real *x, schnaps_real t, 
				   schnaps_real *wL, schnaps_real *vnorm,
				   schnaps_real *flux);
void cemracs2014_TransInitData(schnaps_real *x, schnaps_real *w);
void cemcracs2014_imposed_data(const schnaps_real *x, const schnaps_real t, schnaps_real *w); 
void cemracs2014a_TransInitData(schnaps_real *x, schnaps_real *w);
void cemcracs2014a_imposed_data(const schnaps_real *x, const schnaps_real t, schnaps_real *w); 

void OneSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);


void schnaps_model_load(char* model_source, Model *m);

#endif

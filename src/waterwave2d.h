#ifndef _WATERWAVE2D_H
#define _WATERWAVE2D_H

#include "model.h"


#define _SPEED_WAVE (1)

#define _LENGTH_DOMAIN (1.0)

#define _GRAVITY (1.0)


//! \brief boundardy flux based on the upwind scheme for wave
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Upwind_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,schnaps_real *flux);


//! \brief upwind flux for wave
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Upwind_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);

//! \brief rusanov flux for wave
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Rusanov_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);

//! \brief centered flux for wave
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Centered_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);

//! \brief compute exact solution for x and t
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestPeriodic_Wave_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);


//! \brief init periodic solution for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestPeriodic_Wave_InitData(schnaps_real *x, schnaps_real *w);


//! \brief roe flux for shallows water
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void ShallowWater_Roe_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);


//! \brief HLL flux for shallow water
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void ShallowWater_HLL_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);

//! \brief Rusanov flux for Shallow water
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void ShallowWater_Rusanov_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);


 //! \brief a pointer to the source function for the equilibrium case
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_classical_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);


 //! \brief a pointer to the source function for the periodic case
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_periodic_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);

 //! \brief a pointer to the source function for the Steady State case with u imposed
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_SteadyState_U_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);

 //! \brief a pointer to the source function for the Steady State case with p imposed
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_SteadyState_P_SourceTerm(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);


//! \brief init solution for x
//! \param[in] x : current position
//! \param[in] t : current time
//! \param[out] w : solution exact
void TestSH_equilibrium_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);


//! \brief init solution for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestSH_equilibrium_InitData(schnaps_real *x, schnaps_real *w);

//! \brief init solution for x
//! \param[in] x : current position
//! \param[in] t : current time
//! \param[out] w : solution exact
void TestSH_periodic_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);

//! \brief imposed SteadyState with u imposed solution for x
//! \param[in] x : current position
//! \param[in] t : current time
//! \param[out] w : solution exact
void TestSH_SteadyState_U_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);

//! \brief imposed SteadyState with p imposed solution for x
//! \param[in] x : current position
//! \param[in] t : current time
//! \param[out] w : solution exact
void TestSH_SteadyState_P_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);

//! \brief init solution for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestSH_periodic_InitData(schnaps_real *x, schnaps_real *w);

//! \brief init a steady solution with u imposed for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestSH_SteadyState_U_InitData(schnaps_real *x, schnaps_real *w);

//! \brief init a steady solution with p imposed for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestSH_SteadyState_P_InitData(schnaps_real *x, schnaps_real *w);


#endif

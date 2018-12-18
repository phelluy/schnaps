#ifndef _SIMULATION_SPU_H
#define _SIMULATION_SPU_H

#include "simulation.h"
#include <starpu.h>

enum simulation_spu_task_prio{
    ZeroBuffer_SPU_PRIO,
    ExtractInterface_SPU_PRIO,
    DGMacroCellInterface_SPU_PRIO,
    DGMacroCellBoundaryFlux_SPU_PRIO,
    AddBuffer_SPU_PRIO,
    DGSubCellInterface_SPU_PRIO,
    DGSource_SPU_PRIO,
    DGVolume_SPU_PRIO,
    DGMass_SPU_PRIO,
    SSPU_NB_PRIO // Should be always last
};

schnaps_real seconds();

//! \brief Display a starpu data handle on runtime
//! The function calls starpu_task_wait_for_all before display.
//! StarPU version
//! \param[in] handle a starpu handle
//! \param[in] name a starpu handle name (appears on every line: make it short)
void DisplayHandle_SPU(starpu_data_handle_t handle,
                       const char* name);


//! \brief Init and get ZeroBuffer codelet.
struct starpu_codelet* ZeroBuffer_codelet();

//! \brief Create a ZeroBuffer task and submit it.
//! StarPU version
//! \param[in,out] handle handle of the buffer to nullifiy
void ZeroBuffer_SPU2(starpu_data_handle_t handle, int ie);

void ZeroBuffer_SPU(starpu_data_handle_t handle);


//! \brief Init and get AddBuffer codelet.
struct starpu_codelet* AddBuffer_codelet();

//! \brief Add a scaled buffer to another buffer.
//! handle_out = handle_out + alpha * handle_in
//! StarPU version
//! \param[in] alpha the scaling factor
//! \param[in] handle_in handle to the scaled buffer
//! \param[in,out] handle_out handle to the result buffer
void AddBuffer_SPU(schnaps_real alpha,
                   starpu_data_handle_t handle_in,
                   starpu_data_handle_t handle_out);

void AddBuffer_SPU2(schnaps_real alpha,
                   starpu_data_handle_t handle_in,
                    starpu_data_handle_t handle_out,
                    int ie);



//! \brief Init and get DGSubCellInterface codelet.
struct starpu_codelet* DGSubCellInterface_codelet();

//! \brief Apply the flux terms inside a macrocell
//! StarPU version
//! \param[in,out] f field
void DGSubCellInterface_SPU(field* f);

void DGSubCellInterface_SPU2(field* f, int ie);


//! \brief Init and get DGVolume codelet.
struct starpu_codelet* DGVolume_codelet();

//! \brief Apply the "cross" derivative terms inside a macrocell
//! StarPU version
//! \param[in,out] f field
void DGVolume_SPU2(field* f, int ie);

void DGVolume_SPU(field* f);


//! \brief Init and get DGSource codelet.
struct starpu_codelet* DGSource_codelet();

//! \brief Apply the source terms inside a macrocell
//! StarPU version
//! \param[in,out] f field
void DGSource_SPU(field* f);

void DGSource_SPU2(field* f, int ie);


//! \brief Init and get DGMass codelet.
struct starpu_codelet* DGMass_codelet();

//! \brief Apply the inverse of the mass matrix in a macrocell
//! StarPU version
//! \param[in,out] f field
void DGMass_SPU(field* f);

void DGMass_SPU2(field* f, int ie);


//! \brief Init and get DGMacroCellInterface codelet.
struct starpu_codelet* DGMacroCellInterface_codelet();

//! \brief Apply the interface fluxes to a neighbouring field
//! StarPU version
//! \param[in] inter interface
//! \param[in,out] side side: left if == 0 right if == 1
void DGMacroCellInterface_SPU(Interface* inter, int side);

void DGMacroCellInterface_SPU2(Interface* inter, int side);


//! \brief Init and get DGMacroCellBoundaryFlux codelet.
struct starpu_codelet* DGMacroCellBoundaryFlux_codelet();

//! \brief Apply the boundary flux to a field
//! StarPU version
//! \param[in,out] inter interface
void DGMacroCellBoundaryFlux_SPU(Interface* inter);

void DGMacroCellBoundaryFlux_SPU2(Interface* inter);



//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! starpu version
//! \param[inout] simu a simulation
//! \param[inout] w a starpu handle to the field value
//! \param[out] dtw a starpu handle to the time derivatives
void DtFields_SPU(Simulation *simu,
		  starpu_data_handle_t* w_handle,
		  starpu_data_handle_t* dtw_handle);

//! \brief RK2 integration of the DG approximation
//! starpu version
//! \param[inout] simu a simulation
//! \param[in] tmax tmax
void RK2_SPU(Simulation *simu, schnaps_real tmax);


//! \brief  apply the interface fluxes to a neighbouring field
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void InterfaceExplicitFlux_bis(Interface* inter, int side);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! test version before starpu implementation
//! \param[in] simu a simulation
//! \param[inout] w an array to the field value
//! \param[out] dtw an array to the time derivatives
void DtFields_bis(Simulation *simu,
		  schnaps_real* w,
		  schnaps_real* dtw);

void SmartPrefetch_SPU(Simulation *simu);


#endif

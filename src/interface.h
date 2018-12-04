#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "field.h"


typedef struct Interface{

  //! \brief Left field
  field *fL;

  //! \brief Right field
  field *fR;

  //! \brief left macrocell index
  int ieL;

  //! \brief right macrocell index
  int ieR;


  //! \brief local interface index in Left field
  int locfaL;

  //! \brief local interface index in Left field
  int locfaR;

  //! \brief number of left glops
  int npgL;

  //! \brief number of right glops
  int npgR;

  //! \brief Left glops volume indices from the interface glop indices
  int* vol_indexL;
  starpu_data_handle_t vol_indexL_handle;

  //! \brief Right glops volume indices from the interface glop indices
  int * vol_indexR;
  starpu_data_handle_t vol_indexR_handle;

  //! \brief Size in memory of the stored conservative data at Left interface glops
  int wsizeL;

  //! \brief Size in memory of the stored conservative data at Left interface glops
  int wsizeR;

  //! \brief Left conservative glops data
  schnaps_real *wL;
  starpu_data_handle_t wL_handle;

  //! \brief Right conservative glops data
  schnaps_real *wR;
  starpu_data_handle_t wR_handle;

  //! \brief weighted normal vectors at interface glops
  schnaps_real *vnds;
  starpu_data_handle_t vnds_handle;

  //! \brief glops coordinates
  schnaps_real *xpg;
  starpu_data_handle_t xpg_handle;

  //! \brief glops weights
  schnaps_real *wpg;
  starpu_data_handle_t wpg_handle;

  //! \brief starpu registering flag
  bool starpu_registered;


} Interface;

//! \brief  registration of starpu data for an interface
//! \param[inout] inter an Interface
void RegisterInterface_SPU(Interface* inter);

//! \brief Unregister interface data from starpu
//! \param[in,out] inter interface
void UnregisterInterface_SPU(Interface* inter);


//! \brief  extract the values of the neighbouring fields to the interface
//! \param[inout] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void ExtractInterface(Interface* inter, int side);

//! \brief  extract the values of the neighbouring fields to the interface
//! \param[inout] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
//! \param[in] w a vector with the volume data to extract
void ExtractInterface_bis(Interface* inter, int side, schnaps_real* w);


//! \brief Init and get ExtractInterface codelet.
struct starpu_codelet* ExtractInterface_codelet();

//! \brief Extract the values of the neighbouring fields to the interface
//! StarPU version
//! \param[in,out] inter interface
//! \param[in] side side: left if == 0 right if == 1
void ExtractInterface_SPU(Interface* inter, int side);

void ExtractInterface_SPU2(Interface* inter, int side);


//! \brief  apply the interface fluxes to a neighbouring field
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void InterfaceExplicitFlux(Interface* inter, int side);


//! \brief  apply the interface fluxes to a neighbouring field
//! \brief StarPU version
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void InterfaceExplicitFlux_SPU(Interface* inter, int side);

//! \brief  apply the boundary fluxes to the Left field
//! \brief StarPU version
//! \param[in] inter an Interface
void InterfaceBoundaryFlux_SPU(Interface* inter);

//! \brief  assembly of the inter-fields matrix
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
//! \param[in] theta Cranck-Nicolson parameter
//! \param[in] dt time step
void InterfaceLocalAssembly(Interface *inter,  schnaps_real theta, schnaps_real dt);

//! \brief  varindex function of the interface
//! \param[in] npg number of interface Gauss points
//! \param[in] m number of conservative variables
//! \param[in] ipgf Gauss point index
//! \param[in] iv conservative variable index
//! \returns the memory position of the variable
#pragma start_opencl
int VarindexFace(int npg, int m, int ipgf, int iv);
#pragma end_opencl

#endif

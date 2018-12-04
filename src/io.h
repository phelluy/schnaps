#ifndef _SCHNAPS_IO_H
#define _SCHNAPS_IO_H

#include "macromesh.h"
#include "field.h"
#include "interface.h"
#include "simulation.h"
#include "model.h"


//!\brief xmdf fields dump - pure xml version containing mesh and data - each macrocell is duped as a structured grid
//! with all information at gauss-lobatto points dumped - duplicate points at subcell interface are removed
//!\param[in] simu simulation object
//!\param[in] filename name of the xdmf file (be sure to change it if you write multiple files)
//!\param[in] data_select : selector of  macro (_KIRSCH_SIMU_DIAG_DATA_MACRO) or micro fields (_KIRSCH_SIMU_DIAG_DATA_MICRO) 
//! to output
void schnaps_simulation_xdmf_plot_fields_xml_structured_glops(Simulation *simu, char *filename);

#ifdef _WITH_HDF5
#include <hdf5.h>
//!\brief write macro and micro fields data to hdf5 file ( for restart purpose)
//! BEWARE : this is NOT a full simulation state dump;  
//!\param[in] simu a simulation object
//!\param[in] filename_pref prefix for the hdf5 data file . 
//!\ The routine appends both  mpi process (%03d) and iteration number (%06d)  to build the actual filename.
void schnaps_simulation_dump_field_state(Simulation *simu,char* filename_pref);

//!\brief load macro and micro fields data from hdf5 file ( for restart purpose)
//! BEWARE : this is NOT a full simulation state dump; 
//!\param[in] simu a simulation object
//!\param[in] filename_pref prefix for the hdf5 data file . 
//!\ The routine appends both  mpi process (%03d) and iteration number (%06d)  to build the actual filename.
void schnaps_simulation_load_field_state(Simulation *simu,char* filename_pref);
#endif
#endif
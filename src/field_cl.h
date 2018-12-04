#ifndef _FIELD_CL_H
#define _FIELD_CL_H

#ifdef _WITH_OPENCL

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif

#include "simulation.h"
#include "clinfo.h"

//! copy back the field to host memory
//! \param[inout] f a field
void CopyfieldtoCPU(Simulation *simu);

//! copy back the field to GPU
void CopyCPU2GPU(Simulation *simu);

void update_physnode_cl(Simulation *simu, int ie, cl_mem physnode_cl, schnaps_real *physnode,
			cl_ulong *time,
			cl_uint nwait, cl_event *wait, cl_event *done);
void set_source_CL(Simulation *simu, char *sourcename_cl);
void set_buf_to_zero_cl(cl_mem *buf, int size, Simulation *simu,
			cl_uint nwait, cl_event *wait,  cl_event *done);
void dtfield_CL(Simulation *simu, cl_mem *dtwn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done);
void DGFlux_CL(Simulation *simu, int d, int ie, cl_mem *wn_cl,
	       cl_uint nwait, cl_event *wait, cl_event *done);
void DGVolume_CL(int ie, Simulation *simu, cl_mem *dtwn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done);
void DGMacroCellInterface_CL(int ifa, Simulation *simu, cl_mem *wn_cl,
			     cl_uint nwait, cl_event *wait, cl_event *done);
void DGBoundary_CL(int ifa, Simulation *simu, cl_mem *wn_cl,
		   cl_uint nwait, cl_event *wait, cl_event *done);
void DGMass_CL(int ie, Simulation *simu,
	       cl_uint nwait, cl_event *wait, cl_event *done);

void show_cl_timing(Simulation *simu);

void init_field_cl(Simulation *simu);
void set_physnodes_cl(Simulation *simu);


//! \brief OpenCL version of RK2
//! time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2_CL(Simulation *simu, schnaps_real tmax, schnaps_real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done);
void RK4_CL(Simulation *simu, schnaps_real tmax, schnaps_real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done);

#endif // _WITH_OPENCL

#endif // _FIELD_CL_H

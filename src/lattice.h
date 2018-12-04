#ifndef _LATTICE_H
#define _LATTICE_H

#include "model.h"
#include "field.h"
#include "simulation.h"
#include "global.h"
#include "implicit.h"
//!
//! \brief flux for the lattice model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Lattice_NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux);
void Lattice_OneNodeNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux);

//! \brief computation of equilibrium distribution (global wrapper)
//! \para[in] simu a simulation object
//! \para[inout] w_eq real array  
void Compute_distribution_eq(Simulation * simu,schnaps_real * w_eq);

void Compute_relaxation(Simulation * simu,schnaps_real * w_eq);

void Compute_moments(Simulation * simu);

// init routine for one node model, as simu init does not cope with NULL InitData
void Lattice_Dummy_InitData(schnaps_real x[3],schnaps_real w[]);


// equilibrium functions


//! \brief equiilibrium distribution for the euler/navier stokes isothermal model
//! \paral[in] i_node index of velocity node 
//! \para[in]  lattice object 
//! \para[in] rho, ux, uy, iz, temp, p macroscopic variabbles
//! \returns
schnaps_real feq_isothermal_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);
schnaps_real feq_isothermal_linearwave_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);


// time schemes


//! \brief implicit scheme for lattice model ; creates m=1 simulation defined by model_advec for each velocity nodes
//! \param[inout] simu Simulation with the full lattice model (with all f and macro quantities)
//! \param[in] model_advec a model with m=1 describing the advection of a single velocity node 
//! \param[in] tmax end time of simulation
//! \param[in] dt time step
void LatticeThetaTimeScheme(Simulation *simu, Model *model_advec,schnaps_real tmax, schnaps_real dt);



// Plotting and diagnostics routines (WIP) 



void PlotExtFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void PlotExtFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void PlotFieldsBinSparseMultitime(int typplot, int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep);
void PlotExtScalarFieldBinMultitime(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename,int create_file, schnaps_real t,int istep);
void PlotVecFieldsBinSparseMultitime(int typplot[3], int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep);
//
void Store_Lattice_diags(Simulation *simu);
void Dump_Lattice_Diagnostics(Simulation *simu,char simtag[3]);

// \brief output the vorticy z component along with the u vector field in a gmsh file (meant for 2D simulations)
void Compute_and_dump_Vorticity_2d(Simulation *simu,char *filename, int create_file, schnaps_real t, int istep);
//
#endif

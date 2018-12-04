#ifndef _LBM_DIAGNOSTICS_H
#define _LBM_DIAGNOSTICS_H
#include "simulation.h"
#include "lbm_generic.h"
void LBM_PlotExtFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void LBM_PlotExtFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void LBM_PlotFieldsBinSparseMultitime(int typplot, int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep);
void LBM_PlotExtScalarFieldBinMultitime(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename,int create_file, schnaps_real t,int istep);
void LBM_PlotVecFieldsBinSparseMultitime(int typplot[3], int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep);
//
void LBM_Store_Lattice_diags(LBMSimulation *lbsimu);
void LBM_Dump_Lattice_Diagnostics(LBMSimulation *lbsimu,char simtag[3]);

// \brief output the vorticy z component along with the u vector field in a gmsh file (meant for 2D simulations)
//void LBM_Compute_and_dump_Vorticity_2d(LBMSimulation *lbsimu,char *filename, int create_file, schnaps_real t, int istep);
#endif

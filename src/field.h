#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "interpolation.h"
#include "model.h"

typedef struct Field{
  //underlying mesh and model
  MacroMesh macromesh;
  Model model;
  Interpolation interp;
  int interp_param[8];
  // current time
  double tnow;
  // cfl parameter min_i (vol_i / surf_i)
  double hmin;
  // time step
  // dt has to be smaller than hmin / vmax
  double dt;

  // activate or not 2D computations
  bool is2d;

  // fields at time steps n and n+1 
  double* wn;
  // time derivative of the field
  double* wnp1;
  double* dtwn;

  // memory location of field component iv at Gauss point ipg and
  // element elem
  int (*varindex)(int* param, int elem, int ipg, int iv);

} Field;

// memory location of field component iv at Gauss point ipg and
// element elem (generic access)
int GenericVarindex(int* param, int elem, int ipg, int iv);

void InitField(Field* f);

// apply the Discontinuous Galerkin approximation for computing
// the time derivative of the field (one subcell version)
void dtFieldSlow(Field* f);
// same function but works with subcells
void dtField(Field* f);
// compute the Discontinuous Galerkin volume terms
void DGVolume(Field* f);
// compute the Discontinuous Galerkin volume terms
// slow version
void DGVolumeSlow(Field* f);
// compute the Discontinuous Galerkin inter-subcells terms
void DGSubCellInterface(Field* f);
// compute the Discontinuous Galerkin inter-macrocells boundary terms
void DGMacroCellInterface(Field* f);
// apply the DG mass term
void DGMass(Field* f);

// time integration by a second order Runge-Kutta algorithm 
void RK2(Field* f,double tmax);
// time integration by a second order Runge-Kutta algorithm 
// slow version
void RK2Copy(Field* f,double tmax);

// save the results in the gmsh format
void PlotField(int typplot,int compare,Field* f,char* filename);

// display the field on screen
// typplot: index of the plotted variable
// int compare == true -> compare with the exact value
void DisplayField(Field* f);

// compute the normalized L2 distance with the imposed data
double L2error(Field* f);


#endif

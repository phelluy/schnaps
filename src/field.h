#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "model.h"

typedef struct Field{
  //underlying mesh and model
  MacroMesh macromesh;
  Model model;
  
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
// the time derivative of the field
void dtField(Field* f);

// time integration by a second order Runge-Kutta algorithm 
void RK2(Field* f,double tmax);

// save the results in the gmsh format
void DisplayField(Field* f,char* filename);


#endif

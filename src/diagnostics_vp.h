#ifndef _DIAGSNOTICS_VP_H
#define _DIAGSNOTICS_VP_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"


//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
//! \param[in] x : point of the mesh
//! \param[in] t : time
//! \param[in] t : type of L2norm. if type_norm=0 this is the numerical solution if type_norm=1 this is the error 
double L2VelError(field * f,double* x,double *w);

double L2_Kinetic_error(field* f);

double local_kinetic_energy(field * f,double* x,double *w);

void Energies(field* f,double k_energy, double e_energy, double t_energy);

void Plot_Energies(field* f);


#endif

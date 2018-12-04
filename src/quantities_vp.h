#ifndef _QUANTITIES_VP_H
#define _QUANTITIES_VP_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "simulation.h"

void Computation_charge_density(Simulation *simu);

//void Compute_electric_field(field * f, real * w);
void ComputeElectricField(field* f);

schnaps_real Computation_charge_average(Simulation *simu);

schnaps_real Computation_Maxwellian(schnaps_real rho, schnaps_real U, schnaps_real T, schnaps_real v);

void Computation_Fluid_Quantities(Simulation *simu);

void Computation_Fluid_Quantities_loc(Simulation *simu, schnaps_real *w);

void Collision_Source(Simulation *simu);

#endif

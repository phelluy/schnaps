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
#include "simulation.h"


//! \brief compute square of velocity L2 error
schnaps_real L2VelError(field * f,schnaps_real* x,schnaps_real *w);

schnaps_real L2_Kinetic_error(Simulation * simu);

schnaps_real local_kinetic_energy(field * f,schnaps_real* x,schnaps_real *w);

void Energies(Simulation* simu,schnaps_real * w, schnaps_real k_energy,
	      schnaps_real e_energy, schnaps_real t_energy,int first_diag);

void Charge_total(Simulation * simu, schnaps_real * w, schnaps_real t_charge,
		  int first_diag);

void Taux_instability(Simulation * simu, schnaps_real * w, schnaps_real mode,
		      schnaps_real taux_ins,int first_diag);

void Plot_Energies(Simulation *simu, schnaps_real dt);


#endif

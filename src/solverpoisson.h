#ifndef _SOLVERPOISSON_H
#define _SOLVERPOISSON_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "skyline.h"


void SolvePoisson(field *f,int type_bc, double bc_l, double bc_r);

#endif

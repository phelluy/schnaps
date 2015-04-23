#ifndef _SOLVERPOISSON_H
#define _SOLVERPOISSON_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "skyline.h"


void SolvePoisson(field *f,real * w,int type_bc, real bc_l, real bc_r);

#endif

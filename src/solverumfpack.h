#ifndef _SOLVERUMFPACK_H
#define _SOLVERUMFPACK_H

#ifdef __APPLE__
#include "umfpack.h"
#else
#include "suitesparse/umfpack.h"
#endif
#include "linear_solver.h"

int smalltestumfpack(void);

#endif

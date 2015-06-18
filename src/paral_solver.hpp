#ifndef _PARALSOLVER_HPP
#define _PARALSOLVER_HPP

#include "global.h"
#include "linearsolver.h"
#include <stdbool.h>

#include <paralution.hpp>

void Solver_Paral_CG_Jacobi(LinearSolver* lsol);

void MatrixVector_Skyline_to_Paralution(LinearSolver* lsol,LocalMatrix<double> mat,LocalVector<double> RHS);

void Vector_Paralution_to_Skyline(LinearSolver* lsol,LocalVector<double> Sol);

#endif

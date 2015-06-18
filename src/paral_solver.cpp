#include "global.h"
#include "paral_solver.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <paralution.hpp>

using namespace paralution;

void Solver_Paral_CG_Jacobi(LinearSolver* lsol){

  
  init_paralution();
  info_paralution();
  
LocalVector<double> RHS;
LocalVector<double> Sol;
LocalMatrix<double> mat;

 Sol.Allocate(”Sol”,lsol->neq);
 RHS.Allocate(”RHS”,lsol->neq);
 
 switch(lsol->storage_type) {
  case SKYLINE :
    MatrixVector_Skyline_to_Paralution(lsol,mat,RHS);
    break;

  case CSR :
    assert(1==2);
    break;

  default : 
    assert(1==2);
  }



 switch(lsol->storage_type) {
  case SKYLINE :
    Vector_Paralution_to_Skyline(lsol,Sol);
    break;
    
  case CSR :
    assert(1==2);
    break;

  default : 
    assert(1==2);
  }

}

#include "linear_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "skyline.h"


void InitLinearSolver(LinearSolver* lsol,int n,
		      MatrixStorage* matstor,
		      Solver* solvtyp){

  lsol->neq=n;
  lsol->storage_type = SKYLINE;
  lsol->solver_type = LU;
  lsol->is_sym = false;
  lsol->is_init = false;
  lsol->is_alloc = false;
  lsol->rhs=NULL;
  lsol->sol=NULL;

  if (matstor != NULL) lsol->storage_type = *matstor;
  if (solvtyp != NULL) lsol->solver_type = *solvtyp;

  switch(lsol->storage_type) {

  Skyline* sky;

  case SKYLINE :
    sky = malloc(sizeof(Skyline));
    assert(sky);
    lsol->matrix = (void*) sky;
    InitSkyline(sky,n);
    lsol->is_init = true;
    break;

  case CSR :
    assert(1==2);
    break;

  default : 
    assert(1==2);

  }


}

void FreeLinearSolver(LinearSolver* lsol){

  assert(lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    FreeSkyline((Skyline*)lsol->matrix);
    free(lsol->matrix);
    break;

  default : 
    assert(1==2);
   
  }

  lsol->is_alloc= false;
    

}


void IsNonZero(LinearSolver* lsol,int i,int j){

  assert(lsol->is_init);
  assert(!lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    SwitchOn((Skyline*)lsol->matrix,i,j);
    break;

  default : 
    assert(1==2);
   
  }
  

} 

void AllocateLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(!lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    AllocateSkyline((Skyline*)lsol->matrix);
    break;

  default : 
    assert(1==2);
   
  }
  lsol->is_alloc=true;
}

void AddLinearSolver(LinearSolver* lsol,int i,int j,real val){

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    SetSkyline((Skyline*)lsol->matrix,i,j,val);
    break;

  default : 
    assert(1==2);
   
  }

} 

real GetLinearSolver(LinearSolver* lsol,int i,int j){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  real val;

  switch(lsol->storage_type) {

  case SKYLINE :
    val=GetSkyline((Skyline*)lsol->matrix,i,j);
    break;

  default : 
    assert(1==2);
   
  }
 
  return val;

} 


void DisplayLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  switch(lsol->storage_type) {

  case SKYLINE :
    DisplaySkyline((Skyline*)lsol->matrix);
    break;

  default : 
    assert(1==2);
  }

} 

void MatVecLinearSolver(LinearSolver* lsol,real* x,real* prod){

}

void LUDecompLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  switch(lsol->storage_type) {

  case SKYLINE :
    FactoLU((Skyline*)lsol->matrix);
    break;

  default : 
    assert(1==2);
  }

}

void SolveLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);
  
  switch(lsol->storage_type) {

    Skyline* sky;

  case SKYLINE :
    sky=(Skyline*)lsol->matrix;
    if (!sky->is_lu) FactoLU(sky);
    SolveSkyline(sky,lsol->rhs,lsol->sol);
    break;

  default : 
    assert(1==2);
  }

}

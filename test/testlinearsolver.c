#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  
  // unit tests
    
  int resu=TestLinearSolver();
	 

  if (resu) printf("Linear solver test OK !\n");
  else printf("Linear solver test failed !\n");

  return !resu;
} 




int TestLinearSolver(void){

  int test=1;

  LinearSolver sky;


#define _NN 5
  
  // preliminary work on the skyline struct
  // _NN is the size of the linear system to be solved
  InitLinearSolver(&sky,_NN,NULL,NULL);

  real A[_NN][_NN];
  real vf[_NN],sol[_NN];

  A[0][0] = 0.2e1;
  A[0][1] = -0.1e1;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = 0;
  A[1][0] = -0.1e1;
  A[1][1] = 0.2e1;
  A[1][2] = -0.1e1;
  A[1][3] = 0;
  A[1][4] = 0;
  A[2][0] = 0;
  A[2][1] = -0.1e1;
  A[2][2] = 0.2e1;
  A[2][3] = -0.1e1;
  A[2][4] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -0.1e1;
  A[3][3] = 0.2e1;
  A[3][4] = -0.1e1;
  A[4][0] = 0;
  A[4][1] = 0;
  A[4][2] = 0;
  A[4][3] = -0.1e1;
  A[4][4] = 0.2e1;
  vf[0] = 0;
  vf[1] = 0;
  vf[2] = 0;
  vf[3] = 0;
  vf[4] = 0.6e1;



  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	AddLinearSolver(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetLinearSolver(&sky,i,j,2); */
      /* } */
    }
  }


  // printf for checking...
  DisplayLinearSolver(&sky);


  // LU decomposition
  LUDecompLinearSolver(&sky);


  // solve from a decomposed matrix
  // vf: rhs
  // sol: solution

  sky.rhs=vf;
  sky.sol=sol;
  SolveLinearSolver(&sky);


  // checking
  real verr=0;
  printf("sol=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);
  

  test= (verr<1e-10);

  InitLinearSolver(&sky,_NN,NULL,NULL);
  // now test a symmetric matrix
  A[0][0] = 0.3e1;
  A[0][1] = -0.1e1;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = -0.1e1;
  A[1][0] = -0.1e1;
  A[1][1] = 0.2e1;
  A[1][2] = -0.1e1;
  A[1][3] = 0;
  A[1][4] = 0;
  A[2][0] = 0;
  A[2][1] = -0.1e1;
  A[2][2] = 0.2e1;
  A[2][3] = -0.1e1;
  A[2][4] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -0.1e1;
  A[3][3] = 0.2e1;
  A[3][4] = -0.1e1;
  A[4][0] = -0.1e1;
  A[4][1] = 0;
  A[4][2] = 0;
  A[4][3] = -0.1e1;
  A[4][4] = 0.3e1;
  vf[0] = -0.4e1;
  vf[1] = 0;
  vf[2] = 0;
  vf[3] = 0;
  vf[4] = 0.10e2;

  sky.is_sym=false;

  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
      //if ((A[i][j] != 0) && (j >= i)) IsNonZero(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	AddLinearSolver(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetLinearSolver(&sky,i,j,2); */
      /* } */
    }
  }

  // printf for checking...
  DisplayLinearSolver(&sky);

  // LU decomposition
  LUDecompLinearSolver(&sky);

  // printf for checking...
  DisplayLinearSolver(&sky);

  // solve from a decomposed matrix
  // vf: rhs
  // sol: solution
  sky.rhs=vf;
  sky.sol=sol;
  SolveLinearSolver(&sky);


  // checking
  verr=0;
  printf("sol=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);

  test= test && (verr<1e-10);


  return test;

}

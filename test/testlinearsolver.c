#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

#define _NN 5

int TestLinearSolver(void);

int main(void) {
#ifndef _DOUBLE_PRECISION
  printf("Test not available in single precision\n");
  return 0;
#endif
  
  int resu=TestLinearSolver();	 

  if (resu) printf("Linear solver test OK !\n");
  else printf("Linear solver test failed !\n");

  return !resu;
} 



int TestLinearSolver(void){

  int test=0,test1=1,test2=1,test3=1;
  Simulation simu;

  LinearSolver sky;

  //MatrixStorage ms = SKYLINE_SPU;
  MatrixStorage ms = SKYLINE;
  
  //InitLinearSolver(&sky,_NN,NULL,NULL);
  InitLinearSolver(&sky,_NN,&ms,NULL);

  sky.solver_type = LU;
  sky.pc_type=NONE;

  schnaps_real A[_NN][_NN];
  schnaps_real vf[_NN],sol[_NN];

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
  vf[0] = 0.0;
  vf[1] = 0.0;
  vf[2] = 0.0;
  vf[3] = 0.0;
  vf[4] = 0.6e1;

  sol[0] = 0.0;
  sol[1] = 0.0;
  sol[2] = 0.0;
  sol[3] = 0.0;
  sol[4] = 0.0;


  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
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
    }
  }
 
  for(int i=0;i<_NN;i++){
    sky.rhs[i]=vf[i];
    sky.sol[i]=sol[i];
  }
  
  // printf for checking...
  DisplayLinearSolver(&sky);
  
  Advanced_SolveLinearSolver(&sky,&simu);

  // checking
  schnaps_real verr=0;
  for(int i=0;i<_NN;i++){
    printf("%.8e ",sky.sol[i]);
    verr+=fabs(sky.sol[i]-i-1);
  }
  printf("\n");
  printf("\n");
  // deallocate memory
  FreeLinearSolver(&sky);

  test1 = test1 && (verr<1e-10);
  printf("Error =%.12e\n",verr);

 
  InitLinearSolver(&sky,_NN,&ms,NULL);

  sky.solver_type = GMRES;
  sky.pc_type=NONE;
  
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

  sol[0] = 0;
  sol[1] = 0;
  sol[2] = 0;
  sol[3] = 0;
  sol[4] = 0;

  sky.is_sym=false;
 
  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
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
    }
  }

  for(int i=0;i<_NN;i++){
    sky.rhs[i]=vf[i];
    sky.sol[i]=sol[i];
  }

  Advanced_SolveLinearSolver(&sky,&simu);
  
  // checking
  verr=0;
  printf("sol of gmres=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sky.sol[i]);
    verr+=fabs(sky.sol[i]-i-1.0);
  }
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);

  test2 = test2 && (verr<1e-6);
  printf("Error =%.12e\n",verr);

  int NPoisson=60;
  schnaps_real h=1.0/NPoisson;
  

  InitLinearSolver(&sky,NPoisson,&ms,NULL);

  sky.solver_type = GMRES;
  sky.pc_type=JACOBI;
  sky.iter_max=20000;
  sky.tol=1.e-7;
  
  for(int i=0;i<NPoisson;i++){
    if (i==0){ 
      IsNonZero(&sky,0,0);
      IsNonZero(&sky,0,1);
    } 
    else if (i==NPoisson-1){ 
      IsNonZero(&sky,NPoisson-1,NPoisson-1);
      IsNonZero(&sky,NPoisson-1,NPoisson-2);
    } 
    else { 
    IsNonZero(&sky,i,i);
    IsNonZero(&sky,i,i+1);
    IsNonZero(&sky,i,i-1);
    } 
  }
  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  for(int i=0;i<NPoisson;i++){
    if (i==0){ 
      AddLinearSolver(&sky,0,0,2.0/(h*h));
      AddLinearSolver(&sky,0,1,-1.0/(h*h));
    } 
    else if (i==NPoisson-1){ 
      AddLinearSolver(&sky,NPoisson-1,NPoisson-1,2.0/(h*h));
      AddLinearSolver(&sky,NPoisson-1,NPoisson-2,-1.0/(h*h));
    } 
    else { 
      AddLinearSolver(&sky,i,i,2.0/(h*h));
      AddLinearSolver(&sky,i,i+1,-1.0/(h*h));
      AddLinearSolver(&sky,i,i-1,-1.0/(h*h));
    } 
    sky.sol[i]=0.0;
    sky.rhs[i]=2.0;
  }

  schnaps_real bigval=1.e+15;
  AddLinearSolver(&sky,NPoisson-1,NPoisson-1,bigval);
  AddLinearSolver(&sky,0,0,bigval);
 
  Advanced_SolveLinearSolver(&sky,&simu);


  // checking
  verr=0;
  printf("sol of laplacien with gmres=");
  for(int i=0;i<NPoisson;i++){
    printf("%.5e ",sky.sol[i]);
    verr+=h*fabs(sky.sol[i]-(i*h)*(1-i*h))*fabs(sky.sol[i]-(i*h)*(1-i*h));
  }
  verr=sqrt(verr);
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);

  test3 = test3 && (verr<5.e-2);
  printf("Error =%.12e\n",verr);

  if(test1==1 && test2==1 && test3==1) test=1;

  

  return test;

}

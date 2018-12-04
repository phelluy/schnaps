#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int TestSkyline_SPU(void);

int main(void) {
  
  // unit tests
    
  int resu=TestSkyline_SPU();
	 

  if (resu) printf("Skyline_SPU test OK !\n");
  else printf("Skyline_SPU test failed !\n");

  return !resu;
} 




int TestSkyline_SPU(void){

  int test=1;

  int nbtests = 5;

  starpu_use = true;

  Skyline_SPU sky[nbtests];


#define _NN 5
  
  schnaps_real A[_NN][_NN];
  schnaps_real vf[_NN];
  schnaps_real sol[_NN]={1,2,3,4,5};
  schnaps_real vf2[_NN]={0,0,0,0,6};

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

  for (int is = 0; is < nbtests; is++){
    // preliminary work on the skyline struct
    // _NN is the size of the linear system to be solved
    InitSkyline_SPU(sky + is,_NN);
    for(int i = 0 ; i < _NN; i++) sky[is].rhs[i] = 0; // xxx
    for(int i = 0 ; i < _NN; i++) sky[is].sol[i] = 0; // xxx
  }

  for (int is = 0; is < nbtests; is++){
    // first mark the nonzero values in A
    for(int i=0;i<_NN;i++){
      for(int j=0;j<_NN;j++){
	if (A[i][j] != 0) SwitchOn_SPU(sky + is,i,j);
	//if (i==j) SwitchOn(&sky,i,j);
      }
    }
    
    // once the nonzero positions are known allocate memory
    AllocateSkyline_SPU(sky + is);

    // now set the nonzero terms
    for(int i=0;i<_NN;i++){
      for(int j=0;j<_NN;j++){
	if (A[i][j] != 0){
	  SetSkyline_SPU(sky + is,i,j,A[i][j]);
	}
	/* if (i==j){ */
	/* 	SetSkyline_SPU(&sky,i,j,2); */
	/* } */
      }
    }
  }
  
  schnaps_real vf3[_NN];
    
  // test the product non symmetric case
  for(int i=0; i < _NN; i++){
    vf3[i]=0;
    for(int j=0; j< _NN; j++){
      vf3[i] += A[i][j] * sol[j];
    }
  }


  for(int is = 0; is < nbtests; is++){
  
    for(int i = 0 ; i < _NN; i++) sky[is].sol[i] = sol[i];

    MatVectSkyline_SPU(sky + is, NULL, NULL);
    UnRegisterSkyline_SPU(sky + is);

    printf("test num. %d\n",is);
    for(int i=0; i < _NN; i++) {
      printf("%d vf=%f vf3=%f\n",i,sky[is].rhs[i],vf3[i]);
      test = test && fabs(sky[is].rhs[i]-vf3[i]) < _SMALL;
    }
  }

  


  
  
  // LU decomposition
  for(int is = 0; is < nbtests; is++){
    FactoLU_SPU(sky + is);

    // printf for checking...
    DisplaySkyline_SPU(sky + is);
    
    // solve from a decomposed matrix
    UnRegisterSkyline_SPU(sky + is);

    for(int i = 0 ; i < _NN; i++) sky[is].rhs[i] = vf[i];
    for(int i = 0 ; i < _NN; i++) sky[is].sol[i] = 0;
    
    SolveSkyline_SPU(sky + is);
    UnRegisterSkyline_SPU(sky + is);


    // checking
    schnaps_real verr=0;
    printf("sol test num. %d\n", is);
    for(int i=0;i<_NN;i++){
      printf("%f ",sky[is].sol[i]);
      verr+=fabs(sky[is].sol[i]-i-1);
    }
    test = test && (verr < _SMALL);
    printf("\n");
    
    // deallocate memory
    FreeSkyline_SPU(sky + is);

  }




  return test;

}

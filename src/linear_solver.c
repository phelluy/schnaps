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
  lsol->pc_type = None;
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

void MatVecLinearSolver(LinearSolver* lsol,real x[],real prod[]){
  int i,j;
  real aij;
  
  switch(lsol->storage_type) {

  case SKYLINE :
  
    for(i=0;i<lsol->neq;i++)
      {
	prod[i]=0;
	for(j=0;j<lsol->neq;j++) {
	  aij=GetLinearSolver(lsol,i,j);
	  prod[i] += aij*x[j];
	}
      }
    
    break;

  default : 
    assert(1==2);
  }

  
}

void Vector_copy(real x[],real prod[],int N){
  int i;
 
    for(i=0;i<N;i++)
      {
	prod[i] = x[i];
      }     
  
}


void Vector_prodot(real x[],real y[],real prod[],int N){
  int i;
 
    for(i=0;i<N;i++)
      {
	prod[i] = x[i]*y[i];
      }     
  
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
  
  switch(lsol->solver_type) {
    
    Skyline* sky;
    
  case LU :
    
    switch(lsol->storage_type) {
    case SKYLINE :
         sky=(Skyline*)lsol->matrix;
      if (!sky->is_lu) FactoLU(sky);
      SolveSkyline(sky,lsol->rhs,lsol->sol);
      break;
      
    default : 
      assert(1==2);      
    }
    break;

  case GMRES_CERFACS:
    GMRESSolver(lsol);
    break;  

  case GMRES :
    Solver_Paralution(lsol);
    break;   
    
  default : 
    assert(1==2);      
  }

}


void Solver_Paralution(LinearSolver* lsol){
  int * rows=NULL;
  int * cols=NULL;
  real * coefs=NULL;
  
  double * mat_coefs=NULL;
  double * RHS=NULL;
  double * Sol=NULL;
  char * solver;
  char * pc;
  char * storage;
  double residu=0; 
  int nnz=0,n=0;
  
  solver="GMRES";
  storage="CSR";
  pc="None";
  
  int basis_size_gmres=30, ILU_p=0,ILU_q=0;
  int iter_final=0,ierr=0,maxit=10000;
  double a_tol=1.e-13,r_tol=1.e-8,div_tol=1.e+2;
 

  n=lsol->neq;
  RHS = calloc(n, sizeof(double));
  Sol = calloc(n, sizeof(double));

  for(int i=1;i<n;i++){
    RHS[i] = (double) lsol->rhs[i];
    Sol[i] = (double )lsol->sol[i];
  }

   
 switch(lsol->storage_type) {
  case SKYLINE :
    nnz=Matrix_Skyline_to_COO(lsol->matrix,rows,cols,coefs);

    mat_coefs = calloc(n, sizeof(double));
    for(int i=1;i<nnz;i++){
      mat_coefs[i] = (double) coefs[i];
    }
    
#ifdef PARALUTION
    paralution_fortran_solve_coo(n,n,nnz,solver,storage,pc,storage,
    				 rows,cols,mat_coefs,RHS,a_tol,r_tol,div_tol,maxit,
    				 basis_size_gmres,ILU_p,ILU_q,Sol,iter_final,residu,ierr);
#endif /* PARALUTION */
    break;

  case CSR :
    assert(1==2);
    break;

  default : 
    assert(1==2);
  }
  
  for(int i=1;i<n;i++){
    lsol->sol[i] = (real) Sol[i];
  }
  
}



void GMRESSolver(LinearSolver* lsol){
  int revcom, colx, coly, colz, nbscal, lwork, li_maxiter;
  int irc[5], icntl[8], info[3];
  int ierr,N,m,c,i;
  real cntl[5];
  real rinfo[2];
  real sum,err,sum_rhs,lr_tol;
  real* work;
  real*loc_x;
  real*loc_z;
  real*loc_y;
  int matvec=1, precondLeft=2, precondRight=3, dotProd=4;
 
    
  assert(lsol->is_init);
  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);
  
  //call init_dgmres(icntl,cntl) DECALER INDICE TABLEAU ICNTL
     
  icntl[2]  = 6 ;           // output unit
  icntl[6]  = 2000; // Maximum number of iterations
  icntl[3]  = 0; //!1            // preconditioner (1) = left preconditioner
  icntl[4]  = 1; ////3            // orthogonalization scheme
  icntl[5]  = 1; //1            // initial guess  (1) = user supplied guess
  icntl[7]  = 1; //1            
   
     
  cntl[0]  = 0.00000001; //       ! stopping tolerance
  cntl[1]  = 1.0;
  cntl[2]  = 1.0;
  cntl[2]  = 1.0;         
  cntl[4]  = 1.0;

  m = 30; 
  N = lsol->neq;
  lwork = m*m + m*(N+5) + 5*N + m + 1;

  work = calloc(lwork, sizeof(real));
  loc_x = calloc(N, sizeof(real));
  loc_z = calloc(N, sizeof(real));
  loc_y = calloc(N, sizeof(real));
  
  for(int ivec = 0; ivec < N; ivec++) {
    work[ivec]     = lsol->sol[ivec];                    
    work[N+ivec]    = lsol->rhs[ivec];
  }

  
  //*****************************************
  //** Reverse communication implementation
  //*****************************************
  //L10:    drive_dgmres(N,N,m,lwork,work,irc,icntl,cntl,info,rinfo)
   
  //revcom = irc(1);
  //colx   = irc(2);
  //coly   = irc(3);
  //colz   = irc(4);
  //nbscal = irc(5);
  
  colx=N;
  colz=0;
  coly=N;
  revcom=4;
  
  for(int ivec = 0; ivec < N; ivec++) {
    loc_z[ivec]=work[colz+ivec];                    
    loc_x[ivec]=work[colx+ivec];
    loc_y[ivec]=work[coly+ivec];
  }

  if (revcom == matvec) {                 // perform the matrix vector product
    // work(colz) <-- A * work(colx)  
    MatVecLinearSolver(lsol,loc_x,loc_z);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]=loc_z[ivec];                    
      work[colx+ivec]=loc_x[ivec]; 
    }
    // goto L10;
  }
  else if(revcom == precondLeft)  {                 // perform the matrix vector product
    // work(colz) <-- A * work(colx)  
    Vector_copy(loc_x,loc_z,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]=loc_z[ivec];                    
      work[colx+ivec]=loc_x[ivec]; 
    }
    // goto L10;
  }

  else if(revcom == precondRight)  {                 // perform the matrix vector product
    // work(colz) <-- A * work(colx)  
    Vector_copy(loc_x,loc_z,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]=loc_z[ivec];                    
      work[colx+ivec]=loc_x[ivec]; 
    }
    //  goto L10;
  }

  else if(revcom == dotProd)  {                 // perform the matrix vector product
    // work(colz) <-- A * work(colx)  
    Vector_prodot(loc_x,loc_y,loc_z,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]=loc_z[ivec];                    
      work[colx+ivec]=loc_x[ivec];
      work[coly+ivec]=loc_y[ivec]; 
    }
    //	 goto L10;
  }

  //******************************** end of GMRES reverse communication
  

  for(int ivec = 0; ivec < N; ivec++) {
    lsol->sol[ivec] = work[ivec];                    
  }
  

}

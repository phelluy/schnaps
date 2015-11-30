#include "linear_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "skyline.h"
#include "paralution_c.h"
#include "dpackfgmres.h"

void InitLinearSolver(LinearSolver* lsol,int n,
		      MatrixStorage* matstor,
		      Solver* solvtyp)
{
  lsol->neq=n;
  lsol->storage_type = SKYLINE;
  lsol->solver_type = LU;
  lsol->pc_type = NONE;
  lsol->is_sym = false;
  lsol->is_init = false;
  lsol->is_alloc = false;
  lsol->is_assembly = false;
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

void FreeLinearSolver(LinearSolver* lsol)
{
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

void IsNonZero(LinearSolver* lsol,int i,int j)
{
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

void AllocateLinearSolver(LinearSolver* lsol)
{
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

void AddLinearSolver(LinearSolver* lsol,int i,int j,real val)
{
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

real GetLinearSolver(LinearSolver* lsol,int i,int j)
{
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

void DisplayLinearSolver(LinearSolver* lsol)
{
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

void MatVecLinearSolver(LinearSolver* lsol,real x[],real prod[])
{
  int i,j;
  real aij;
  
  switch(lsol->storage_type) {
  case SKYLINE :
    for(i = 0; i < lsol->neq; i++) {
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

void Vector_copy(real x[], real prod[],int N)
{
  for(int i = 0; i < N; i++) {
    prod[i] = x[i];
  }      
}

real Vector_norm2(real x[],int  N)
{
  real norm=0;
  for(int i = 0; i < N; i++) {
    norm += x[i] * x[i];
  }
  return sqrt(norm);
}

real Vector_prodot(real x[],real y[],int N)
{
  real prod = 0;
  for(int i = 0; i <N; i++) {
    prod += x[i]*y[i];
  }     
  return prod;
}

void LUDecompLinearSolver(LinearSolver* lsol)
{
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

void SolveLinearSolver(LinearSolver* lsol)
{
  assert(lsol->is_init);
  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);

  if(lsol->solver_type == LU) {
    Skyline* sky;
    switch(lsol->storage_type) {
    case SKYLINE :
      sky=(Skyline*)lsol->matrix;
      if (!sky->is_lu) FactoLU(sky);
      SolveSkyline(sky,lsol->rhs,lsol->sol);
      break;
    default : 
      assert(1==2);      
    }
  }
  else if(lsol->solver_type == GMRES) {
    GMRESSolver(lsol);
  }
  else {
#ifdef PARALUTION
    Solver_Paralution(lsol);
#else
    printf("paralution is not build this solver is not accessible.");
    assert(1==2);
#endif
  }  
}

void Solver_Paralution(LinearSolver* lsol)
{
  int * rows=NULL;
  int * cols=NULL;
  real * coefs=NULL;
  
  real * mat_coefs=NULL;
  real * RHS=NULL;
  real * Sol=NULL;
  char * solver;
  char * pc;
  char * storage;
  real * residu=0; 
  int nnz=0,n=0,c=0;
  Skyline * mat;
  
  int basis_size_gmres=30, ILU_p=2,ILU_q=2;
  int* iter_final=0;
  int* ierr=0;
  int maxit=100000;
  real norm_rhs=0;
  real a_tol=0,r_tol=0,div_tol=1.e+8;

  storage="CSR";
  norm_rhs=Vector_norm2(lsol->rhs,lsol->neq);
  a_tol=1.e-8*(1.0+1.e-20*norm_rhs);

  switch(lsol->solver_type){
  case PAR_CG :
    solver="CG";
    break;
  case PAR_GMRES :
    solver="GMRES";
    break;
  case PAR_FGMRES :
    solver="FGMRES";
    break;
  case PAR_BICGSTAB :
    solver="BiCGStab";
    break;
  case PAR_AMG :
    solver="AMG";
    break;
  case PAR_LU :
    solver="LU";
    storage="DENSE";
    break;
  case PAR_QR :
    solver="QR";
    storage="DENSE";
    break;
  default : 
    assert(1==2);   
  }

  switch(lsol->pc_type){
  case NONE :
    pc="None";
    break;   
  case PAR_JACOBI :
    pc="Jacobi";
    break;
  case PAR_ILU :
    pc="ILU";
    break;
  case PAR_MULTICOLOREDGS :
    pc="MultiColoredGS";
    break;
  case PAR_MULTICOLOREDSGS :
    pc="MultiColoredSGS";
    break;
  case PAR_MULTICOLOREDILU :
    pc="MultiColoredILU";
    break;
  case PAR_AMG_PC :
    pc="AMG";
    break;
  case PAR_ELIMI :
    pc="ELIMI";
    break;  
  default : 
    assert(1==2);   
  }

  n=lsol->neq;
  RHS = calloc(n,sizeof(real));
  Sol = calloc(n,sizeof(real));

  for(int i=0;i<n;i++){
    RHS[i] = (real) lsol->rhs[i];
    Sol[i] = (real )lsol->sol[i];   
  }
  
  switch(lsol->storage_type) {
  case SKYLINE :

    mat = lsol->matrix;
    
    nnz= mat->neq+2*mat->nmem;
  
    rows = (int*) malloc(nnz*sizeof(int)); 
    cols = (int*) malloc(nnz*sizeof(int));
    coefs = (real*) malloc(nnz*sizeof(real));
    assert(rows);
  
    for (int i=0;i< mat->neq; i++) {
      for (int j=0;j< mat->neq; j++) {
	if (j-i <= mat->prof[j] && i-j <= mat->prof[i]){
	  if (i==j){
	    coefs[c]=mat->vkgd[i];
	    rows[c]=i;
	    cols[c]=j;
	    c++;
	  }
	  else if ( j>i){
	    int k=mat->kld[j+1]-j+i;
	    coefs[c]=mat->vkgs[k];
	    rows[c]=i;
	    cols[c]=j;
	    c++; 
	  }
	  else {
	    int k=mat->kld[i+1]-i+j;
	    coefs[c]=mat->vkgi[k];
	    rows[c]=i;
	    cols[c]=j;
	    c++; 
	  }
	}
      }
    }    
    
    mat_coefs = malloc(nnz*sizeof(real));
    for(int i=0;i<nnz;i++){
      mat_coefs[i] = (real) coefs[i];
    }
    
#ifdef PARALUTION
    paralution_fortran_solve_coo(n,n,nnz,solver,storage,pc,storage,
    				 rows,cols,mat_coefs,RHS,a_tol,r_tol,div_tol,
				 maxit,
    				 basis_size_gmres,ILU_p,ILU_q,Sol); 
#endif /* PARALUTION */
    break;
  case CSR :
    assert(1==2);
    break;
  default : 
    assert(1==2);
  }
  
  for(int i = 0; i < n; i++) {
    lsol->sol[i] = (real) Sol[i];
  }
}

void GMRESSolver(LinearSolver* lsol)
{
  // FIXME: why are we using GOTOs here?  This should be cleaned up.
  
  int revcom, colx, coly, colz, nbscal;
  int li_maxiter;
  int m,lwork,N;
  int * pt_m;
  int * pt_lwork;
  int * pt_Size;

  int irc[5+1];
  int icntl[8+1];
  int info[3+1];
  real cntl[5+1];
  real rinfo[2+1];
  real sum,err,sum_rhs,lr_tol;
  real * work;
  real *loc_x;
  real *loc_y;
  real *loc_z;
  real prodot=0.0;
  int res=0;
  int matvec=1, precondLeft=2, precondRight=3, dotProd=4;

  res=init_dgmres(icntl,cntl);
  
  icntl[3]  = 6 ;           // output unit
  icntl[7]  = 2000; // Maximum number of iterations
  icntl[4]  = 0; //!1            // preconditioner (1) = left preconditioner
  icntl[5]  = 1; ////3            // orthogonalization scheme
  icntl[6]  = 1; //1            // initial guess  (1) = user supplied guess
  icntl[8]  = 1; //1            
     
  cntl[1]  = 0.000000001; //       ! stopping tolerance
  cntl[2]  = 1.0;
  cntl[3]  = 1.0;
  cntl[4]  = 1.0;         
  cntl[5]  = 1.0;

  N = lsol->neq;
  
  if(N<61) {
    m = (int) (N/2)+1;
  }
  else {
    m = 30;
  }
  lwork = m*m + m*(N+5) + 5*N + m + 1 +1; //(+ one because  ) ??

  pt_m = &m;
  pt_Size = &N;
  pt_lwork = &lwork;

  work = calloc(lwork, sizeof(real));
  loc_x = calloc(N, sizeof(real));
  loc_y = calloc(N, sizeof(real));
  loc_z = calloc(N, sizeof(real));
  
  for(int ivec = 0; ivec < N; ivec++) {
    work[ivec+1]     = lsol->sol[ivec];                    
    work[N+ivec+1]    = lsol->rhs[ivec];
  }
  
  //*****************************************
  //** Reverse communication implementation
  //*****************************************
 L10:    res=drive_dgmres(pt_Size,pt_Size,pt_m,pt_lwork,&work[1],&irc[1],&icntl[1],&cntl[1],&info[1],&rinfo[1]);

  revcom = irc[1];
  colx   = irc[2];
  coly   = irc[3];
  colz   = irc[4];
  nbscal = irc[5];
  
  for(int ivec = 0; ivec < N; ivec++) {
    loc_z[ivec]= work[colz+ivec];                    
    loc_x[ivec]= work[colx+ivec];
    loc_y[ivec]= work[coly+ivec];    
  }

  if (revcom == matvec) {                 // perform the matrix vector product
    // work(colz) <-- A * work(colx)
    MatVecLinearSolver(lsol,loc_x,loc_z);//MatVecLinearSolver(lsol,loc_x,loc_z);
    for(int ivec = 0; ivec < N; ivec++) {       
      work[colz+ivec]= loc_z[ivec];                    
      work[colx+ivec]= loc_x[ivec]; 
    }

    // https://xkcd.com/292/
    goto L10;
  }
  else if(revcom == precondLeft)  {                 // perform the matrix vector product
    // work(colz) <-- M-1 * work(colx)  
    Vector_copy(loc_x,loc_z,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]= loc_z[ivec];                    
      work[colx+ivec]= loc_x[ivec]; 
    }

    // https://xkcd.com/292/
    goto L10;
  }

  else if(revcom == precondRight)  {                 // perform the matrix vector product
    // work(colz) <-- M-1 * work(colx)  
    Vector_copy(loc_x,loc_z,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]= loc_z[ivec];                    
      work[colx+ivec]= loc_x[ivec]; 
    }

    // https://xkcd.com/292/
    goto L10;
  }

  else if(revcom == dotProd)  {// perform the matrix vector product
    // work(colz) <-- work(colx) work(coly)
    prodot=Vector_prodot(loc_x,loc_y,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colx+ivec]= loc_x[ivec];
      work[coly+ivec]= loc_y[ivec]; 
    }
    work[colz]= prodot;  

    // https://xkcd.com/292/
    goto L10;
  }

  //******************************** end of GMRES reverse communication
  
  for(int ivec = 0; ivec < N; ivec++) {
    lsol->sol[ivec] = work[ivec+1];                    
  }
}

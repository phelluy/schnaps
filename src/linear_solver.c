#include "linear_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "skyline.h"
#include "skyline_spu.h"
#include "dpackfgmres.h"
#include "klu_csr.h"

void InitLinearSolver(LinearSolver * lsol, int n,
		      MatrixStorage * matstor, Solver * solvtyp)
{

  lsol->neq = n;
  lsol->storage_type = SKYLINE;
  lsol->solver_type = LU;
  lsol->pc_type = NONE;
  lsol->is_sym = false;
  lsol->is_init = false;
  lsol->is_alloc = false;
  lsol->rhs_is_assembly = false;
  lsol->mat_is_assembly = false;
  lsol->rhs = NULL;
  lsol->sol = NULL;
  lsol->MatVecProduct = NULL;
  lsol->tol = 1.e-9;
  lsol->restart_gmres = 1;
  lsol->iter_max = 100;
  lsol->is_CG = false;

  if (matstor != NULL)
    lsol->storage_type = *matstor;
  if (solvtyp != NULL)
    lsol->solver_type = *solvtyp;

  if (matstor != NULL) {
    if (*matstor != SKYLINE_SPU) {
      lsol->rhs = calloc(n, sizeof(schnaps_real));
      lsol->sol = calloc(n, sizeof(schnaps_real));
    }
  }

  switch (lsol->storage_type) {

    Skyline *sky;
    Skyline_SPU *sky_spu;
    KLU *klumat;
  case SKYLINE:
    sky = malloc(sizeof(Skyline));
    assert(sky);
    lsol->matrix = (void *) sky;
    InitSkyline(sky, n);
    lsol->is_init = true;
    break;

  case SKYLINE_SPU:
    sky_spu = malloc(sizeof(Skyline_SPU));
    assert(sky_spu);
    lsol->matrix = (void *) sky_spu;
    InitSkyline_SPU(sky_spu, n);
    lsol->rhs = sky_spu->rhs;
    lsol->sol = sky_spu->sol;
    lsol->is_init = true;
    break;

  case KLU_CSR:
    klumat = (KLU *) malloc(sizeof(KLU));
    assert(klumat);
    lsol->matrix = (void *) klumat;
    InitKLU(klumat, n);
    lsol->is_init = true;
    break;

  case CSR:
    assert(1 == 2);
    break;

  default:
    assert(1 == 2);

  }

}

void FreeLinearSolver(LinearSolver * lsol)
{

  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);
  // Free rhs
  if (lsol->rhs != NULL) {
    free(lsol->rhs);
  }
  // Free sol
  if (lsol->sol != NULL) {
    free(lsol->sol);
  }

  switch (lsol->storage_type) {

  case SKYLINE:
    FreeSkyline((Skyline *) lsol->matrix);
    // Free matrix
    if (lsol->matrix != NULL) {
      free(lsol->matrix);
    }
    break;

  case SKYLINE_SPU:
    FreeSkyline_SPU((Skyline_SPU *) lsol->matrix);
    free(lsol->matrix);
    break;

  case KLU_CSR:
    FreeKLU((KLU *) lsol->matrix);
    free(lsol->matrix);
    lsol->matrix = NULL;
    break;

  default:
    assert(1 == 2);

  }

  lsol->is_alloc = false;


}


void IsNonZero(LinearSolver * lsol, int i, int j)
{

  assert(lsol->is_init);
  assert(!lsol->is_alloc);

  switch (lsol->storage_type) {

  case SKYLINE:
    SwitchOn((Skyline *) lsol->matrix, i, j);
    break;

  case SKYLINE_SPU:
    SwitchOn_SPU((Skyline_SPU *) lsol->matrix, i, j);
    break;

  case KLU_CSR:
    SwitchOnKLU((KLU *) lsol->matrix, i, j);
    break;

  default:
    assert(1 == 2);

  }


}

void AllocateLinearSolver(LinearSolver * lsol)
{

  assert(lsol->is_init);
  assert(!lsol->is_alloc);

  switch (lsol->storage_type) {

  case SKYLINE:
    AllocateSkyline((Skyline *) lsol->matrix);
    break;

  case SKYLINE_SPU:
    AllocateSkyline_SPU((Skyline_SPU *) lsol->matrix);
    break;

  case KLU_CSR:
    AllocateKLU((KLU *) lsol->matrix);
    break;

  default:
    assert(1 == 2);

  }
  lsol->is_alloc = true;
}

void AddLinearSolver(LinearSolver * lsol, int i, int j, schnaps_real val)
{

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch (lsol->storage_type) {

  case SKYLINE:
    AddSkyline((Skyline *) lsol->matrix, i, j, val);
    break;

  case SKYLINE_SPU:
    AddSkyline_SPU((Skyline_SPU *) lsol->matrix, i, j, val);
    break;

  case KLU_CSR:
    AddKLU((KLU *) lsol->matrix, i, j, val);
    break;

  default:
    assert(1 == 2);

  }

}

void SetLinearSolver(LinearSolver * lsol, int i, int j, schnaps_real val)
{

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch (lsol->storage_type) {

  case SKYLINE:
    SetSkyline((Skyline *) lsol->matrix, i, j, val);
    break;

  case SKYLINE_SPU:
    SetSkyline_SPU((Skyline_SPU *) lsol->matrix, i, j, val);
    break;

  case KLU_CSR:
    SetKLU((KLU *) lsol->matrix, i, j, val);
    break;

  default:
    assert(1 == 2);

  }

}

schnaps_real GetLinearSolver(LinearSolver * lsol, int i, int j)
{

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  schnaps_real val;

  switch (lsol->storage_type) {

  case SKYLINE:
    val = GetSkyline((Skyline *) lsol->matrix, i, j);
    break;

  case KLU_CSR:
    val = GetKLU((KLU *) lsol->matrix, i, j);
    break;
  default:
    assert(1 == 2);

  }

  return val;

}


void DisplayLinearSolver(LinearSolver * lsol)
{

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch (lsol->storage_type) {

  case SKYLINE:
    DisplaySkyline((Skyline *) lsol->matrix);
    break;

  case SKYLINE_SPU:
    DisplaySkyline_SPU((Skyline_SPU *) lsol->matrix);
    break;

  case KLU_CSR:
    DisplayKLU((KLU *) lsol->matrix);
    break;

  default:
    assert(1 == 2);
  }

  printf("rhs=");
  for (int i = 0; i < lsol->neq; i++) {
    printf("%f ", lsol->rhs[i]);
  }
  printf("\n");

}

void MatVect(void *system, schnaps_real x[], schnaps_real prod[])
{
  int i, j;
  schnaps_real aij;
  LinearSolver *lsol = system;

  switch (lsol->storage_type) {

  case SKYLINE:

    //if (!((Skyline*)lsol->matrix)->is_lu){
    MatVectSkyline((Skyline *) lsol->matrix, x, prod);
    //}
    /*else
       {
       for(i=0;i<lsol->neq;i++)
       { 
       prod[i]=0; 
       for(j=0;j<lsol->neq;j++) { 
       aij=GetLinearSolver(lsol,i,j); 
       prod[i] += aij*x[j]; 
       } 
       } 
       } */
    break;

  case SKYLINE_SPU:
    assert(1 == 2);
    break;

  case KLU_CSR:
    MatVectKLU((KLU *) lsol->matrix, x, prod);
    break;

  default:
    assert(1 == 2);
  }


}

void MatVectIn(void *system)
{
  int i, j;
  schnaps_real aij;
  LinearSolver *lsol = system;

  switch (lsol->storage_type) {

  case SKYLINE:
    assert(1 == 2);		// not yet implemented and verified
    break;
    /* for(i=0;i<lsol->neq;i++) */
    /*   { */
    /*  prod[i]=0; */
    /*  for(j=0;j<lsol->neq;j++) { */
    /*    aij=GetLinearSolver(lsol,i,j); */
    /*    prod[i] += aij*x[j]; */
    /*  } */
    /*   } */

  case SKYLINE_SPU:
    MatVectSkyline_SPU((Skyline_SPU *) lsol->matrix, NULL, NULL);
    break;

  case KLU_CSR:
    assert(1 == 2);
    break;

  default:
    assert(1 == 2);
  }


}

void MatVect_SPU(void *system, starpu_data_handle_t sol_handle,
		 starpu_data_handle_t rhs_handle)
{
  int i, j;
  schnaps_real aij;
  LinearSolver *lsol = system;

  switch (lsol->storage_type) {

  case SKYLINE:
    assert(1 == 2);		// not yet implemented and verified
    break;
    /* for(i=0;i<lsol->neq;i++) */
    /*   { */
    /*  prod[i]=0; */
    /*  for(j=0;j<lsol->neq;j++) { */
    /*    aij=GetLinearSolver(lsol,i,j); */
    /*    prod[i] += aij*x[j]; */
    /*  } */
    /*   } */

  case SKYLINE_SPU:
    MatVectSkyline_SPU((Skyline_SPU *) lsol->matrix, sol_handle,
		       rhs_handle);
    break;

  case KLU_CSR:
    assert(1 == 2);
    break;

  default:
    assert(1 == 2);
  }


}


void Vector_copy(schnaps_real x[], schnaps_real prod[], int N)
{
  int i;

  for (i = 0; i < N; i++) {
    prod[i] = x[i];
  }
}


schnaps_real Vector_norm2(schnaps_real x[], int N)
{
  int i;
  schnaps_real norm = 0;

  for (i = 0; i < N; i++) {
    norm += x[i] * x[i];
  }
  norm = sqrt(norm);
  return norm;
}




schnaps_real Vector_prodot(schnaps_real x[], schnaps_real y[], int N)
{
  int i;
  schnaps_real prod;

  prod = 0;
  for (i = 0; i < N; i++) {
    prod += x[i] * y[i];
  }
  return prod;
}



void LUDecompLinearSolver(LinearSolver * lsol)
{

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch (lsol->storage_type) {

  case SKYLINE:
    FactoLU((Skyline *) lsol->matrix);
    break;

  case SKYLINE_SPU:
    FactoLU_SPU((Skyline_SPU *) lsol->matrix);
    break;

  case KLU_CSR:
    FactoKLU((KLU *) lsol->matrix);
    break;

  default:
    assert(1 == 2);
  }

}

void SolveLinearSolver(LinearSolver * lsol)
{

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);

  if (lsol->solver_type == LU) {
    Skyline *sky;
    Skyline_SPU *sky_spu;
    KLU *klumat;
    switch (lsol->storage_type) {

    case SKYLINE:
      sky = (Skyline *) lsol->matrix;
      if (!sky->is_lu) {
	FactoLU(sky);
      }
      SolveSkyline(sky, lsol->rhs, lsol->sol);
      break;

    case SKYLINE_SPU:
      sky_spu = (Skyline_SPU *) lsol->matrix;
      if (!sky_spu->is_lu) {
	//UnRegisterSkyline_SPU(sky_spu);
	//FactoLU_SPU(sky_spu);
	printf("to do: insert factolu in the starpu system !!!!!\n");
	//RegisterSkyline_SPU(sky_spu);
      }
      SolveSkyline_SPU(sky_spu);
      break;

    case KLU_CSR:
      klumat = (KLU *) lsol->matrix;
      if (!klumat->is_lu) {
	FactoKLU(klumat);
      }
      SolveKLU(klumat, lsol->rhs, lsol->sol);
      break;

    default:
      assert(1 == 2);
    }
  } else if (lsol->solver_type == GMRES) {
    GMRESSolver(lsol);
  }

}



void GMRESSolver(LinearSolver * lsol)
{
  int revcom, colx, coly, colz, nbscal;
  int li_maxiter;
  int m, lwork, N;
  int *pt_m;
  int *pt_lwork;
  int *pt_Size;

  int irc[5 + 1];
  int icntl[8 + 1];
  int info[3 + 1];
  schnaps_real cntl[5 + 1];
  schnaps_real rinfo[2 + 1];
  schnaps_real sum, err, sum_rhs, lr_tol;
  schnaps_real *work;
  schnaps_real *loc_x;
  schnaps_real *loc_y;
  schnaps_real *loc_z;
  schnaps_real prodot = 0.0;
  int res = 0;
  int matvec = 1, precondLeft = 2, precondRight = 3, dotProd = 4;


  lsol->MatVecProduct = MatVect;

  res = init_dgmres(icntl, cntl);

  icntl[3] = 6;			// output unit
  icntl[7] = lsol->iter_max;	// Maximum number of iterations
  icntl[4] = 0;			//!1            // preconditioner (1) = left preconditioner
  icntl[5] = 1;			////3            // orthogonalization scheme
  icntl[6] = 1;			//1            // initial guess  (1) = user supplied guess
  icntl[8] = 1;			//1


  if (lsol->pc_type == JACOBI) {
    icntl[4] = 2;
  }

  cntl[1] = lsol->tol;		//       ! stopping tolerance
  cntl[2] = 1.0;
  cntl[3] = 1.0;
  cntl[4] = 1.0;
  cntl[5] = 1.0;

  N = lsol->neq;

  if (lsol->restart_gmres == 1) {
    if (N < 61) {
      m = (int) (N / 2) + 1;
    } else {
      m = 30;
    }
  } else {
    m = lsol->restart_gmres;
  }
  lwork = m * m + m * (N + 5) + 5 * N + m + 1;	//(+ one because  ) ??

  pt_m = &m;
  pt_Size = &N;
  pt_lwork = &lwork;

  work = calloc(lwork, sizeof(schnaps_real));
  loc_x = calloc(N, sizeof(schnaps_real));
  loc_y = calloc(N, sizeof(schnaps_real));
  loc_z = calloc(N, sizeof(schnaps_real));

  for (int ivec = 0; ivec < N; ivec++) {
    work[ivec + 1] = lsol->sol[ivec];
    work[N + ivec + 1] = lsol->rhs[ivec];
  }


  schnaps_real *Ax2 = calloc(lsol->neq, sizeof(schnaps_real));
  lsol->MatVecProduct(lsol, lsol->sol, Ax2);
  schnaps_real errorb = 0;
  schnaps_real errorb2 = 0;
  for (int i = 0; i < N; i++) {
    errorb =
	errorb + fabs((Ax2[i] - lsol->rhs[i]) * (Ax2[i] - lsol->rhs[i]));
    errorb2 = errorb2 + fabs(lsol->sol[i] * lsol->sol[i]);
  }
  printf(" error gmres begin %.5e \n",
	 sqrt(errorb) / (sqrt(errorb2) + 1.0));

  //*****************************************
  //** Reverse communication implementation
  //*****************************************
L10:res =
      drive_dgmres(pt_Size, pt_Size, pt_m, pt_lwork, &work[1], &irc[1],
		   &icntl[1], &cntl[1], &info[1], &rinfo[1]);

  revcom = irc[1];
  colx = irc[2];
  coly = irc[3];
  colz = irc[4];
  nbscal = irc[5];

  for (int ivec = 0; ivec < N; ivec++) {

    loc_z[ivec] = work[colz + ivec];
    loc_x[ivec] = work[colx + ivec];
    loc_y[ivec] = work[coly + ivec];
  }

  if (revcom == matvec) {	// perform the matrix vector product
    // work(colz) <-- A * work(colx)   
    lsol->MatVecProduct(lsol, loc_x, loc_z);
    for (int ivec = 0; ivec < N; ivec++) {
      work[colz + ivec] = loc_z[ivec];
      work[colx + ivec] = loc_x[ivec];
    }
    goto L10;
  } else if (revcom == precondLeft) {	// perform the matrix vector product
    // work(colz) <-- M-1 * work(colx)
    Vector_copy(loc_x, loc_z, N);
    for (int ivec = 0; ivec < N; ivec++) {
      work[colz + ivec] = loc_z[ivec];
      work[colx + ivec] = loc_x[ivec];
    }
    goto L10;
  }

  else if (revcom == precondRight) {
    if (lsol->pc_type == JACOBI) {
      Jacobi_PC(lsol, loc_z, loc_x);
    } else {
      Vector_copy(loc_x, loc_z, N);
    }


    for (int ivec = 0; ivec < N; ivec++) {
      work[colz + ivec] = loc_z[ivec];
      work[colx + ivec] = loc_x[ivec];
    }
    goto L10;
  }

  else if (revcom == dotProd) {	// perform the matrix vector product
    // work(colz) <-- work(colx) work(coly)
    prodot = Vector_prodot(loc_x, loc_y, N);
    for (int ivec = 0; ivec < N; ivec++) {
      work[colx + ivec] = loc_x[ivec];
      work[coly + ivec] = loc_y[ivec];
    }
    work[colz] = prodot;
    goto L10;
  }
  //******************************** end of GMRES reverse communication



  for (int ivec = 0; ivec < N; ivec++) {
    lsol->sol[ivec] = work[ivec + 1];
  }

  schnaps_real *Ax = calloc(lsol->neq, sizeof(schnaps_real));
  lsol->MatVecProduct(lsol, lsol->sol, Ax);
  schnaps_real error = 0;
  schnaps_real error2 = 0;
  for (int i = 0; i < N; i++) {
    error = error + fabs((Ax[i] - lsol->rhs[i]) * (Ax[i] - lsol->rhs[i]));
    error2 = error2 + fabs(lsol->sol[i] * lsol->sol[i]);
  }
  printf(" error gmres end %.5e \n", sqrt(error) / (sqrt(error2) + 1.0));


  free(work);
  free(loc_x);
  free(loc_y);
  free(loc_z);

}



void Jacobi_PC(LinearSolver * lsol, schnaps_real * sol, schnaps_real * rhs)
{

  for (int i = 0; i < lsol->neq; i++) {
    assert(GetLinearSolver(lsol, i, i) != 0);
    sol[i] = rhs[i] / GetLinearSolver(lsol, i, i);
  }
}

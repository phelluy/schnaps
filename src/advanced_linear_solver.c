#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "skyline.h"
#include "skyline_spu.h"
#include "dpackfgmres.h"
#include "advanced_linear_solver.h"
#include "linear_solver.h"
#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "klu_csr.h"

void Advanced_SolveLinearSolver(LinearSolver * lsol, Simulation * simu)
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

    case KLU_CSR :
      klumat = (KLU*)lsol->matrix;
      if (!klumat->is_lu) {
        FactoKLU(klumat);
      }
      SolveKLU(klumat,lsol->rhs,lsol->sol);
      break;

    default:
      assert(1 == 2);
    }
  } else if (lsol->solver_type == GMRES) {
    GMRESSolver(lsol);
  }
}



void InitJFLinearSolver(JFLinearSolver * lsol, int n, Solver * solvtyp)
{

  lsol->neq = n;
  lsol->solver_type = GMRES;
  lsol->pc_type = NONE;
  lsol->rhs = NULL;
  lsol->sol = NULL;
  lsol->soln = NULL;
  lsol->MatVecProduct = NULL;
  lsol->NonlinearVector_computation = NULL;
  lsol->tol = 1.e-6;
  lsol->restart_gmres = 1;
  lsol->iter_max = 10000;
  lsol->eps = 0.000001;

  if (solvtyp != NULL)
    lsol->solver_type = *solvtyp;

  lsol->rhs = calloc(n, sizeof(schnaps_real));
  lsol->sol = calloc(n, sizeof(schnaps_real));
  lsol->soln = calloc(n, sizeof(schnaps_real));

}

void FreeJFLinearSolver(JFLinearSolver * lsol)
{

  free(lsol->rhs);
  free(lsol->sol);
  free(lsol->soln);

}

void MatVecJacobianFree(Simulation * simu, void *system, schnaps_real x[],
			schnaps_real prod[])
{
  int i, j;
  schnaps_real aij;
  JFLinearSolver *lsol = system;
  schnaps_real *U;
  schnaps_real *Up;
  schnaps_real *solnp;

  solnp = calloc(lsol->neq, sizeof(schnaps_real));
  U = calloc(lsol->neq, sizeof(schnaps_real));
  Up = calloc(lsol->neq, sizeof(schnaps_real));


  for (i = 0; i < lsol->neq; i++) {
    solnp[i] = lsol->soln[i] + lsol->eps * x[i];
  }

  lsol->NonlinearVector_computation(simu, system, lsol->soln, U);
  lsol->NonlinearVector_computation(simu, system, solnp, Up);

  for (i = 0; i < lsol->neq; i++) {
    prod[i] = (Up[i] - U[i]) / lsol->eps;
  }

  free(solnp);
  free(U);
  free(Up);


}


void SolveJFLinearSolver(JFLinearSolver * lsol, Simulation * simu)
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


  lsol->MatVecProduct = MatVecJacobianFree;

  res = init_dgmres(icntl, cntl);

  icntl[3] = 6;			// output unit
  icntl[7] = lsol->iter_max;	// Maximum number of iterations
  icntl[4] = 0;			//!1            // preconditioner (1) = left preconditioner
  icntl[5] = 1;			////3            // orthogonalization scheme
  icntl[6] = 1;			//1            // initial guess  (1) = user supplied guess
  icntl[8] = 1;			//1            


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
  lwork = m * m + m * (N + 5) + 5 * N + m + 1 + 1;	//(+ one because  ) ??

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
    lsol->MatVecProduct(simu, lsol, loc_x, loc_z);
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

  else if (revcom == precondRight) {	// perform the matrix vector product
    // work(colz) <-- M-1 * work(colx)  
    Vector_copy(loc_x, loc_z, N);
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

  free(work);
  free(loc_x);
  free(loc_y);
  free(loc_z);
}

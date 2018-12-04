/* Sample program for the UMFPACK sparse matrix solver
*/

#include <stdio.h>
#include <stdlib.h>	/* for EXIT_FAILURE */
#include "solverumfpack.h"

void error_report(int status, const char *file, const char *func, int line)
{
	fprintf(stderr, "in %s: file %s, line %d: ", func, file, line);

	switch (status) {
		case UMFPACK_ERROR_out_of_memory:
			fprintf(stderr, "out of memory!\n");
			break;
		case UMFPACK_WARNING_singular_matrix:
			fprintf(stderr, "matrix is singular!\n");
			break;
		default:
			fprintf(stderr, "UMFPACK error code %d\n", status);
	}
}



int smalltestumfpack(void)
/* { */
/* 	int n = 5; */
/* 	double x[5]; */
/* 	void *Symbolic, *Numeric; */
/* 	int i; */

/* 	/\* cumulative count of entries, as matrix is scanned columnwise *\/ */
/* 	int Ap[] = { 0, 2, 5, 9, 10, 12 }; */

/* 	/\* row indices of entries, as matrix is scanned columnwise *\/ */
/* 	int Ai[] =    { 0, 1, 0,  2, 4, 1,  2, 3, 4, 2, 1, 4 }; */

/* 	/\* matrix entries *\/ */
/* 	double Ax[] = { 2, 3, 3, -1, 4, 4, -3, 1, 2, 2, 6, 1 }; */

/* 	/\* the right hand side *\/ */
/* 	double b[] = { 8, 45, -3, 3, 19 }; */

/* 	/\* symbolic analysis *\/ */
/* 	umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL); */

/* 	/\* LU factorization *\/ */
/* 	umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL); */
/* 	umfpack_di_free_symbolic(&Symbolic); */

/* 	/\* solve system *\/ */
/* 	umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL); */
/* 	umfpack_di_free_numeric(&Numeric); */

/* 	for (i = 0; i < n; i++) */
/* 		printf("x[%d] = %g\n", i, x[i]); */

/* 	return 0; */
/* } */

/* umfpack-ex2.c
   
   Sample program for the UMFPACK sparse matrix solver

   This is just like umfpack-ex1.c, but we create the matrix in the
   "triplet form" then call umfpack_di_triplet_to_col() to convert
   to the compress-column for that was used in umfpack-ex1.c.

   The triplet form is suitable for assembling the mass matrix in
   finite elements.

   I have added error-checking.

   Solves the system Ax=b, where:

   A = 
       2   3   0   0   0
       3   0   4   0   6
       0  -1  -3   2   0
       0   0   1   0   0
       0   4   2   0   1

   b = (8, 45, -3, 3, 19)

   The solution x is:

   x = (1, 2, 3, 4, 5)

   RR, November 2003
*/


//int main(void)
{

#define n   4		/* matrix is nxn */
#define nz 6		/* number of (non-zero) entries */

	void *Symbolic, *Numeric;
	double x[n];
	int Ap[n+1];
	int Ai[nz];
	double Ax[nz];
	int status;
	int i;
	int permut[n];
	double Control[UMFPACK_CONTROL];

	/* In the example below, vectors Ti, Tj, Tx are of length nz each
	   because each entry of the matrix is referred to exaclty once.
	   It is possible to have multiple references to an entry.  In that
	   case their values will be summed.  This is useful when assembling
	   a finite element mass matrix.
	*/

	/* Ti[k] is row index of entry k, as matrix is scanned columnwise */
	int Ti[] =    {0,1,1,2,3,3};

	/* Tj[k] is column index of entry k, as matrix is scanned columnwise */
	int Tj[] =    {0,0,1,2,2,3};

	/* value of entry k, as matrix is scanned columnwise */
	double Tx[] = {1,1,1,1,1,1};

	/* the right hand side */
	double b[] = {1,3,3,7};

	/* convert matrix from triplet form to compressed-column form */
	status = umfpack_di_triplet_to_col(n, n, nz, Ti, Tj, Tx,
			Ap, Ai, Ax, NULL);

	if (status != UMFPACK_OK) {
		error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}

	/* symbolic analysis */
	status = umfpack_di_symbolic(n, n, Ap, Ai, Ax, &Symbolic, NULL, NULL);

	if (status != UMFPACK_OK) {
		error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}

	/* LU factorization */
	umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &Numeric, NULL, NULL);

	if (status != UMFPACK_OK) {
		error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}

	umfpack_di_free_symbolic(&Symbolic);

	/* solve system */
	umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, NULL, NULL);

	if (status != UMFPACK_OK) {
		error_report(status, __FILE__, __func__, __LINE__);
		return EXIT_FAILURE;
	}


	Control[UMFPACK_PRL]=4;
	umfpack_di_report_symbolic(Symbolic, Control);


	umfpack_di_free_numeric(&Numeric);
	/* for (i = 0; i < n; i++) */
	/* 	printf("permut[%d] = %d\n", i, permut[i]); */

	for (i = 0; i < n; i++)
		printf("x[%d] = %g\n", i, x[i]);

	return EXIT_SUCCESS;
}



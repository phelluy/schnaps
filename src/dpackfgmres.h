#ifndef _DPACKGMRES_H
#define _DPACKGMRES_H

int drive_dgmres(int *n, int *nloc, int *m, 
	int *lwork, real *work, int *irc, int *icntl, 
		   real *cntl, int *info, real *rinfo);

int dgmres_(int *n, int *m, real *b, 
	real *x, real *h__, real *w, real *r0, 
	real *v, real *dot, real *ycurrent, real *
	xcurrent, real *rotsin, real *rotcos, int *irc, 
	    int *icntl, real *cntl, int *info, real *rinfo);

int init_dgmres(int *icntl, real *cntl);

int gcopy(int *n, real *dx, int *incx, 
	  real *dy, int *incy);

int gconstantmut(int *n, real *da, real *dx, 
	   int *incx, real *dy, int *incy);

int MatVectorProduct(char *trans, int *m, int *n, real *
	alpha, real *a, int *lda, real *x, int *incx, 
	   real *beta, real *y, int *incy);

real gnorm2(int *n, real *x, int *incx);

int UpperLowerSolver(char *uplo, char *trans, char *diag, int *n, 
		     real *a, int *lda, real *x, int *incx);

int PlaneRotation(int *n, real *dx, int *incx, 
		  real *dy, int *incy, real *c, real *s);

int GivensRot(real *da, real *db, real *c, 
	  real *s);

real Sign(real *a, real *b);

int Compare_char(char *ca, char *cb);

int ErrorBlas(char *srname, int *info);

#endif

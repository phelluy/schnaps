#ifndef _DPACKGMRES_H
#define _DPACKGMRES_H

#include "global.h"

int drive_dgmres(int *n, int *nloc, int *m, 
	int *lwork, schnaps_real *work, int *irc, int *icntl, 
		   schnaps_real *cntl, int *info, schnaps_real *rinfo);

int dgmres_(int *n, int *m, schnaps_real *b, 
	schnaps_real *x, schnaps_real *h__, schnaps_real *w, schnaps_real *r0, 
	schnaps_real *v, schnaps_real *dot, schnaps_real *ycurrent, schnaps_real *
	xcurrent, schnaps_real *rotsin, schnaps_real *rotcos, int *irc, 
	    int *icntl, schnaps_real *cntl, int *info, schnaps_real *rinfo);

int init_dgmres(int *icntl, schnaps_real *cntl);

int gcopy(int *n, schnaps_real *dx, int *incx, 
	  schnaps_real *dy, int *incy);

int gconstantmut(int *n, schnaps_real *da, schnaps_real *dx, 
	   int *incx, schnaps_real *dy, int *incy);

int MatVectorProduct(char *trans, int *m, int *n, schnaps_real *
	alpha, schnaps_real *a, int *lda, schnaps_real *x, int *incx, 
	   schnaps_real *beta, schnaps_real *y, int *incy);

schnaps_real gnorm2(int *n, schnaps_real *x, int *incx);

int UpperLowerSolver(char *uplo, char *trans, char *diag, int *n, 
		     schnaps_real *a, int *lda, schnaps_real *x, int *incx);

int PlaneRotation(int *n, schnaps_real *dx, int *incx, 
		  schnaps_real *dy, int *incy, schnaps_real *c, schnaps_real *s);

int GivensRot(schnaps_real *da, schnaps_real *db, schnaps_real *c, 
	  schnaps_real *s);

schnaps_real Sign(schnaps_real *a, schnaps_real *b);

int Compare_char(char *ca, char *cb);

int ErrorBlas(char *srname, int *info);

#endif

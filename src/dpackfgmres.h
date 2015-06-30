#ifndef _DPACKGMRES_H
#define _DPACKGMRES_H

int drive_dgmres(int *n, int *nloc, int *m, 
	int *lwork, double *work, int *irc, int *icntl, 
		   double *cntl, int *info, double *rinfo);

int dgmres_(int *n, int *m, double *b, 
	double *x, double *h__, double *w, double *r0, 
	double *v, double *dot, double *ycurrent, double *
	xcurrent, double *rotsin, double *rotcos, int *irc, 
	    int *icntl, double *cntl, int *info, double *rinfo);

int init_dgmres(int *icntl, double *cntl);

int gcopy(int *n, double *dx, int *incx, 
	  double *dy, int *incy);

int gconstantmut(int *n, double *da, double *dx, 
	   int *incx, double *dy, int *incy);

int MatVectorProduct(char *trans, int *m, int *n, double *
	alpha, double *a, int *lda, double *x, int *incx, 
	   double *beta, double *y, int *incy);

double gnorm2(int *n, double *x, int *incx);

int UpperLowerSolver(char *uplo, char *trans, char *diag, int *n, 
		     double *a, int *lda, double *x, int *incx);

int PlaneRotation(int *n, double *dx, int *incx, 
		  double *dy, int *incy, double *c, double *s);

int GivensRot(double *da, double *db, double *c, 
	  double *s);

double Sign(double *a, double *b);

int Compare_char(char *ca, char *cb);

int ErrorBlas(char *srname, int *info);

#endif

#ifndef _FE_IMPLICIT_H
#define _FE_IMPLICIT_H
#include "skyline.h"


///////////////// physical data ////////////////////

// nb of kinetic equations
#define M 3

// sound speed
#define CSON (0.6)

// relaxation parameter
//#define RELAX (1e6)
#define RELAX (1e10)
//#define RELAX (4000)

// Riemann data
#define RL 2
#define UL -0.

#define RR 1
#define UR 0.

// estimated maximal wave speed
#define VMAX (1.)

// theta of the theta time scheme
#define THETA (0.5)

// physical flux
void flux(dcmplx* w, dcmplx* flux);

// jacobian of the physical flux
void dflux(dcmplx* w, dcmplx* dflux);

// exact or reference solution 
void solexacte(double x, dcmplx t, dcmplx* w);

///////////// interp. data ////////////////////

// number of finite elements
#define NB_ELEMS (200)
// polynomial order
#define DEG (2)


// nb of finite element nodes
#define NB_NODES ((DEG+1) * NB_ELEMS)

// nb of unknowns on the mesh
#define WSIZE (M * NB_NODES)

/////////// main container ////////////

// a struct for managing all this...
typedef struct galerkin{

  dcmplx wn[WSIZE];
  dcmplx wnm1[WSIZE];

  // nodes positions
  double xnode[NB_NODES];
  
  // mesh bounds
  double xmin;
  double xmax;
  double dx;

  // final time, time step and cfl
  double tmax, cfl, t;
  dcmplx dt;
  
  // time step multiplicator (1/2 or (1+I)/2)
  dcmplx smul;

} galerkin;

// constructor: parameters and init. condition.
void gal_construct(galerkin *gal,
		   double xmin,
		   double xmax,
		   double cfl,
		   double tmax,
		   dcmplx smult);

// data access functions
// memory location of var iv at gauss point ipg
int vindex(int ipg, int iv);

// connectivity array
// return global node index of local node iloc in elem ie
int connec(int ie, int iloc);

// compute and invert the local implicit matrix
void loc_assembly(galerkin *gal, dcmplx *mat, double vit);

// advance one time step
void gal_step(galerkin *gal, dcmplx tnow);

// advance one time step non linear scheme
void gal_step_nonlin(galerkin *gal, dcmplx tnow);

// outputs
void gal_plot(galerkin *gal);

// error measure
double gal_L2_error(galerkin *gal, int numvar);

// interpolate another function from a different mesh
// xi: interp points
// w_sav: stored data
// deg_sav, nbnodes_sav: quantities of saved data
// x: point where is computed the interp values
// ws: computed interpolation
void interpolate(double *xi, dcmplx *w_sav,
		 int deg_sav, int nbnodes_sav,
		 double x, dcmplx *ws);
void lagrange(double *p, double *subdiv,
	      int deg, int ii, double x);
// BGK relaxation
void bgk_relax(dcmplx *w, dcmplx dt);

// BGK projection: infinite relaxation
void bgk_project(dcmplx *w);

// local dg update from boundary data
void local_dg(galerkin *gal, dcmplx *wloc, dcmplx wvnm1, dcmplx wvn, int ivel);


// BGK projection: i=0 or 2
// projection conserving w0 and w2 (momentum)
// but not the total mass
// returns the loss of mass  
dcmplx bgk_project_i(dcmplx *w, int i);

#endif

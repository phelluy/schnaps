#ifndef _ADERDG_H
#define _ADERDG_H

// number of conservative variables
#define _M 2
const double velocity[_M]={1,-1};
// polynomial degree
#define _D 4



typedef struct ADERDG{

  int nbelems;
  int nbfaces; // = nbelems +1

  double xmin,xmax;

  double* face;

  // two arrays for maintaining previous and next
  // values of the _M conservatives variables at the _D+1
  // Gauss points of each cell -> size = _M * (_D+1) * nbelems
  double* wnow;
  double* wnext;
  // one array for the current predicted values
  double* wpred;

  // current time evolution on the bigger cells
  double tnow;
  // macro time step
  double dt;
  // maximal cell size
  double dx;

  // number of CFL levels
  int ncfl;

  int* cell_level;  
  int* face_level;  
  

} ADERDG;

void InitADERDG(ADERDG* adg,int nbelems,double xmin,double xmax);

void Predictor(ADERDG* adg,int ie,double s);

void VolumeTerms(ADERDG* adg,int ie);


// perform a macro time step of the ADER method
void BigStep(ADERDG* adg);

double stretching(double xh);


void NumFlux(double* wL,double* wR,double* flux);


void ExactSol(double x,double t,double* w);


// plotting in gnuplot
//plot 'adgplot.dat' using 1:2 w l , 'adgplot.dat' using 1:4  
//plot 'adgplot.dat' using 1:3 w l , 'adgplot.dat' using 1:5 
void Plot(ADERDG* adg);




//! Gauss LObatto Points (GLOP) up to order 4
static const double gauss_lob_point[] = {
  0.5,
  0,
  1,
  0,
  0.5,
  1,
  0,
  0.276393202250021030359082633127,
  0.723606797749978969640917366873,
  1,
  0,
  0.172673164646011428100853771877,
  0.5,
  0.827326835353988571899146228123,
  1
};

//! GLOP weights up to order 4
static const double gauss_lob_weight[] = {
  1,
  0.5,
  0.5,
  0.166666666666666666666666666667,
  0.666666666666666666666666666668,
  0.166666666666666666666666666667,
  0.0833333333333333333333333333333,
  0.416666666666666666666666666666,
  0.416666666666666666666666666666,
  0.0833333333333333333333333333333,
  0.05,
  0.272222222222222222222222222223,
  0.355555555555555555555555555556,
  0.272222222222222222222222222219,
  0.05
};

//! indirection for finding the GLOP
//! data for a given degree in the previous arrays
static const int gauss_lob_offset[] = {0, 1, 3, 6, 10};

static const double gauss_lob_dpsi[] = {
  0.0,
  -1.,
  -1.,
  1.,
  1.,
  -3.,
  -1.,
  1.,
  4.,
  0.,
  -4.,
  -1.,
  1.,
  3.,
  -6,
  -1.61803398874989484820458683436,
  .618033988749894848204586834362,
  -1,
  8.09016994374947424102293417177,
  0,
  -2.23606797749978969640917366872,
  3.09016994374947424102293417184,
  -3.09016994374947424102293417182,
  2.23606797749978969640917366872,
  0,
  -8.09016994374947424102293417177,
  1,
  -.618033988749894848204586834362,
  1.61803398874989484820458683436,
  6,
  -10,
  -2.48198050606196571569743868436,
  .75,
  -.518019493938034284302561315632,
  1,
  13.5130049774484800076860550594,
  0,
  -2.67316915539090667050969419631,
  1.52752523165194666886268239794,
  -2.82032835588485332564727827404,
  -5.33333333333333333333333333336,
  3.49148624377587810025755976667,
  0,
  -3.49148624377587810025755976662,
  5.33333333333333333333333333336,
  2.82032835588485332564727827399,
  -1.52752523165194666886268239791,
  2.67316915539090667050969419635,
  0,
  -13.5130049774484800076860550594,
  -1,
  .518019493938034284302561315631,
  -.75,
  2.48198050606196571569743868437,
  10
};

//! indirection for finding the GLOP
//! data for a given degree in the previous arrays
static const int gauss_lob_dpsi_offset[] = {0, 1, 5, 14, 30};


//return glop weight i
double wglop(int deg,int i);

// returns glop i
double glop(int deg,int i);
// return the 1d derivative of lagrange polynomial ib at glop ipg
double dlag(int deg,int ib,int ipg);

#endif

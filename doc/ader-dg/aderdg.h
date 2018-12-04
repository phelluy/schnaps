#ifndef _ADERDG_H
#define _ADERDG_H
#include <stdbool.h>
// Number of conservative variables
#define _M 2

const double velocity[_M] = {1, -1};

// Polynomial degree
#define _D 1

// Number of interior elements (without fictional boundary elements)
#define _NBELEMS_IN 1000

// CFL
#if (_D == 0)
#define _CFL 0.56
#elif (_D == 1)
#define _CFL 0.5
#elif (_D == 2)
#define _CFL 0.17
#elif (_D == 3)
#define _CFL 0.08
#elif (_D == 4)
#define _CFL 0.05
#endif

// Number of gauss lobatto points in one element
#define _NGLOPS (_D + 1)
// Total number of elements (counting the two boundary elements)
#define _NBELEMS (_NBELEMS_IN + 2)
// Number of interfaces (counting the two boundary interfaces)
#define _NBFACES (_NBELEMS_IN + 1)



//! \brief 1-D ADER-DG data
typedef struct ADERDG {

  //! \brief 1-D mesh boundaries (without fictional boundary elements)
  double xmin, xmax;

  //! \brief Interface coordinates
  double face[_NBFACES];

  //! \brief Maximal time step (of the biggest cells)
  double dt;

  //! \brief Minimal time step (of the smallest cells)
  double dt_small;

  //! \brief Number of levels
  int ncfl;

  //! \brief Cell levels
  int cell_level[_NBELEMS];

  //! \brief Interface levels
  int face_level[_NBFACES];

  //! \brief Interface ids ordered according to their level
  int face_order[_NBFACES];

  //! \brief Current time (last global update)
  double tnow;

  //! \brief Cell current time (last local update)
  double cell_tnow[_NBELEMS];


  // Three arrays for maintaining previous and next values of the
  // _M conservatives variables at the _D + 1 gauss points and their derivatives
  // ---------------------------------------------------------------------------
  //! \brief Field at current local time
  double wnow[_NBELEMS][_D + 1][_M];

  //! \brief Field at future local time
  double wnext[_NBELEMS][_D + 1][_M];

  //! \brief Derivatives
  double dtw[_NBELEMS][_D + 1][_M];


  //! \brief Predicted values
  double wpred[_NBELEMS][_D + 1][_M];

  //! \brief Indicate if predicted values have been computed for current step
  bool pred_done[_NBELEMS];


} ADERDG;


//! \brief ADER-DG initialization function
//! \param[in] xmin mesh left boundary
//! \param[in] xmax mesh right boundary
//! \todo put the velocity
void ADERDG_Init(ADERDG*, double xmin, double xmax);

//! \brief Plot ADERDG data in gnuplot
// $ plot 'adgplot.dat' using 1:2 w l , 'adgplot.dat' using 1:4
// $ plot 'adgplot.dat' using 1:3 w l , 'adgplot.dat' using 1:5
void ADERDG_Plot(ADERDG*);

//! \brief Perform a resolution by the ADER-DG method
//! \param[in] tmax final time
void ADERDG_Solve(ADERDG*, double tmax);

//! \brief Predictor
//! \param[in] ie element id
//! \param[in] tpred prediction time
void ADERDG_Predictor(ADERDG*, int ie, double tpred);

//! \brief Perform a standard time step with ADER-DG method
void ADERDG_TimeStep(ADERDG*);

//! \brief Perform a macro time step with ADER-DG method
void ADERDG_BigStep(ADERDG*);



// Model
// -----------------------------------------------------------------------------

//! \brief Exact solution
//! \param[in] x coordinate
//! \param[in] t time
//! \param[out] w exact solution
void ExactSol(double x, double t, double w[_M]);

//! \brief Numerical flux
//! \param[in] wL left field
//! \param[in] wR right field
//! \param[out] flux flux
void NumFlux(double* wL, double* wR, double* flux);

//! \brief Prediction
//! \param[in] wnow current field
//! \param[in] pw prediction weight
//! \param[out] wpred predicted field
void Prediction(double* wnow, double pw, double* wpred);



// Interpolation
// -----------------------------------------------------------------------------

//! \brief GLOP ipg
//! \param[in] deg interpolation degree
//! \param[in] ipg gauss point
//! \return local GLOP coordinate
double glop(int deg, int ipg);

//! \brief Weight of GLOP ipg
//! \param[in] deg interpolation degree
//! \param[in] ipg gauss point
//! \return GLOP weight
double wglop(int deg, int ipg);

//! \brief 1-D derivative of lagrange polynomial ib at GLOP ipg
//! \param[in] deg interpolation degree
//! \param[in] ib lagrange polynomial
//! \param[in] ipg gauss point
//! \return derivative
double dlag(int deg, int ib, int ipg);

#endif

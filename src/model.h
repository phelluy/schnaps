#ifndef _MODEL_H
#define _MODEL_H

//! \brief a unified framework for all physical models
typedef struct Model {
  //! Number of conservative variables
  int m;

  //! Number of conservative variables in each dimension (NB: their
  //! prouduct must equal m).
  int mx, my, mz;

  //! The conservative variables have velocity in [-vmax, vmax]
  double vmax;

  //! \brief A pointer to the numflux function
  //! \param[in] wL, wR : left and right states
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*NumFlux)(double wL[], double wR[], double vn[3], double flux[],
		  int mx, int my, int mz, double vmax);

  //! \briefA pointer to the boundary flux function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] wL : left state
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*BoundaryFlux)(double x[3], double t, double wL[], double vn[3],
		       double flux[], 
		       int mx, int my, int mz, double vmax);

  //! \brief A pointer to the init data function
  // !\param[in] x : space position
  //! \param[out] w : init state at point x
  void (*InitData)(double x[3], double w[], 
		   int mx, int my, int mz, double vmax);

  //! \brief A pointer to the imposed data function
  //!\param[in] x, t : space and time position
  //! \param[out] w : imposed state at point x and time t
  void (*ImposedData)(double x[3], double t, double w[], 
		      int mx, int my, int mz, double vmax);

} Model;

//! \brief The particular flux for the transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransNumFlux(double wL[], double wR[], double vn[3], double* flux, 
		      int mx, int my, int mz, double vmax);

//! \brief The particular flux for the 2d transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransNumFlux2d(double wL[], double wR[], double vn[3], double* flux, 
			int mx, int my, int mz, double vmax);

void VecTransNumFlux2d(double wL[], double wR[], double vn[3], double* flux, 
		       int mx, int my, int mz, double vmax);

//! \brief The particular boundary flux for the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransBoundaryFlux(double* x, double t, double* wL, double* vn,
			   double* flux,
			   int mx, int my, int mz, double vmax);

//! \brief The particular boundary flux for the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransBoundaryFlux2d(double* x, double t, double* wL, double* vn,
			     double* flux,
			     int mx, int my, int mz, double vmax);

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData(double* x, double* w,
		       int mx, int my, int mz, double vmax);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData2d(double *x, double *w,
			 int mx, int my, int mz, double vmax);

void VecTransInitData2d(double *x, double *w,
			int mx, int my, int mz, double vmax);

void vTransImposedData2d(double *x, double t, double *w, 
			 int mx, int my, int mz, double vmax);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TransImposedData(double* x, double t, double* w,
			  int mx, int my, int mz, double vmax);

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TransImposedData2d(double* x, double t, double* w, 
			    int mx, int my, int mz, double vmax);

void VecTransImposedData2d(double* x, double t, double* w, 
			   int mx, int my, int mz, double vmax);

//! \brief The particular flux for testing the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux(double* x, double t, double* wL, double* vn,
			       double* flux,
			       int mx, int my, int mz, double vmax);

//! \brief The particular flux for testing the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux2d(double* x, double t, double* wL, double* vn,
				 double* flux,
				 int mx, int my, int mz, double vmax);

void VecTransBoundaryFlux2d(double* x, double t, double* wL, double* vn,
				 double* flux,
			    int mx, int my, int mz, double vmax);

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData(double* x, double* w,
			   int mx, int my, int mz, double vmax);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData2d(double* x, double* w,
			     int mx, int my, int mz, double vmax);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData(double* x, double t, double* w, 
			      int mx, int my, int mz, double vmax);

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData2d(double* x, double t, double* w, 
				int mx, int my, int mz, double vmax);

#endif

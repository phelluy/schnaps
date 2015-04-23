// {{{ MHD
#pragma start_opencl
//! \brief computes the conservatives states from the primitives
//! \param[in] y : primitives states
//! \param[out] w : conservatives states
void conservatives(double* y, double* w);
//! \brief Numerical flux for the MHD model
//! \param[in] w : states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void fluxnum(double* w,double* vn, double* flux);
//! \brief particular flux for the MHD model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void MHDNumFlux(double wL[],double wR[],double vn[3],double* flux);
//! \brief particular boundary flux for the MHD model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void MHDBoundaryFlux(double* x,double t,double* wL,double* vn,
		     double* flux);
//! \brief particular init data for the MHD model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void MHDInitData(double* x,double* w);
//! \brief particular imposed data for the MHD model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void MHDImposedData(double* x,double t,double* w);

void primitives(double* W, double* Y);
void jacobmhd(double* W,double* vn, double M[9][9]);
void matrix_vector(double A[9][9], double B[9], double* C);
void matrix_matrix(double A[9][9],double B[9][9],double C[9][9]);
void write_matrix(double A[9][9],double *second, double B[9][9+1]);
void gauss(double A[9][9], double b[9], double *x);
void MHDNumFlux_2(double wL[],double wR[],double* vn, double* flux);

void MHDNumFlux1D(double wL[],double wR[],double* vn, double* flux);
#pragma end_opencl
// }}}


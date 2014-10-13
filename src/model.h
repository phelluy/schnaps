#ifndef _MODEL_H
#define _MODEL_H

//! \brief a unified framework for all physical models
typedef struct Model{
  //! number of conservative variables
  int m;
  //! \brief a pointer to the numflux function
  void (*NumFlux)(double* wL,double* wR,double* vn,double* flux);
  //! \brief a pointer to the boundary flux function
  //! \param[in] wL[m]: left state
  // !\param[in] x[3],t : space and time positions
  //! \param[in] vn[3]: normal vector
  //! \param[out] flux[m]: the flux
  void (*BoundaryFlux)(double* x,double t,double* wL,double* vn,double* flux);
  //! \brief a pointer to the init data function
  // !\param[in] x[3] : space position
  //! \param[out] w[m]: init state at point x
  void (*InitData)(double* x,double* w);
  //! \brief a pointer to the imposed data function
  // !\param[in] x[3],t : space and time position
  //! \param[out] w[m]: imposed state at point x and time t
  void (*ImposedData)(double* x,double t,double* w);

} Model;


//! particular flux for the transport model
//! \param[in] wL[m],wR[m]: left and right states
//! \param[in] vn[3]: normal vector
//! \param[out] flux[m]: the flux
void TransportNumFlux(double* wL,double* wR,double* vn,double* flux);
void TransportBoundaryFlux(double* x,double t,double* wL,double* vn,
			     double* flux);
  void TransportInitData(double* x,double* w);
  void TransportImposedData(double* x,double t,double* w);
  void TransportImposedDataLinear(double* x,double t,double* w);

// particular functions for unit testing
  void TestTransportBoundaryFlux(double* x,double t,double* wL,double* vn,
			     double* flux);
  void TestTransportInitData(double* x,double* w);
  void TestTransportImposedData(double* x,double t,double* w);
  void TestTransportImposedDataLinear(double* x,double t,double* w);


#endif

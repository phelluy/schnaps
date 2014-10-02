#ifndef _MODEL_H
#define _MODEL_H

// model management
// impose a unified framework for all modeling
typedef struct Model{
  //number of conservative variables
  int m;
  void (*NumFlux)(double* wL,double* wR,double* vn,double* flux);
  void (*BoundaryFlux)(double* x,double t,double* wL,double* vn,double* flux);
  void (*InitData)(double* x,double* w);
  void (*ImposedData)(double* x,double t,double* w);

} Model;


// particular functions for the transport model
  void TransportNumFlux(double* wL,double* wR,double* vn,double* flux);
  void TransportBoundaryFlux(double* x,double t,double* wL,double* vn,
			     double* flux);
  void TransportInitData(double* x,double* w);
  void TransportImposedData(double* x,double t,double* w);

#endif

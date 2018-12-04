#ifndef _LBMMODELSDATA_H
#define _LBMMODELSDATA_H
#include "global.h"
//!\file database of lattice boltzmann models static data (nodes and polynomials describing moments)

// Moments Database (Uncentered and non-normalized moments) macro quantities non-linearly derived from those moments should
// be computed at higher level if necessary 
typedef struct MPolyData {
  int ndim;			//! number of dimension
  int nc;			//! number of monomials
  void *c;			//! opaque pointer to schnaps_real array of nc coefficients
  void *e;			//! opaque pointer to int[nc][ndim] array of exponent multiplets
  void *s;			//! opaque pointer to default name 
} MPolyData;
// density 
static const MPolyData model_rho1D = {
  .ndim = 1,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[1][1]) {{0}}),.s = &("rho")
};

static const MPolyData model_rho2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{0}, {0}}),.s = &("rho")
};

static const MPolyData model_rho3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{0}, {0}, {0}}),.s = &("rho")
};

// momentum flux
static const MPolyData model_jx1D = {
  .ndim = 1,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[1][1]) {{1}}),.s = &("jx")
};

//
static const MPolyData model_jx2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),

  .e = &((int[2][1]) {{1}, {0}}),.s = &("jx")
};

static const MPolyData model_jy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{0}, {1}}),.s = &("jy")
};

//
static const MPolyData model_jx3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{1}, {0}, {0}}),.s = &("jx")
};

static const MPolyData model_jy3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{0}, {1}, {0}}),.s = &("jy")
};

static const MPolyData model_jz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{0}, {0}, {1}}),.s = &("jz")
};

// stress tensor  
static const MPolyData model_Mxx1D = {
  .ndim = 1,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[1][1]) {{2}}),.s = &("Mxx")
};

//
static const MPolyData model_Mxx2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{2}, {0}}),.s = &("Mxx")
};

static const MPolyData model_Myy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{0}, {2}}),.s = &("Myy")
};

static const MPolyData model_Mxy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{1}, {1}}),.s = &("Mxy")
};

//
static const MPolyData model_Mxx3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{2}, {0}, {0}}),.s = &("Mxx")
};

static const MPolyData model_Myy3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{0}, {2}, {0}}),.s = &("Myy")
};

static const MPolyData model_Mzz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{0}, {0}, {2}}),.s = &("Mzz")
};

static const MPolyData model_Mxy3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{1}, {1}, {0}}),.s = &("Mxy")
};

static const MPolyData model_Mxz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{1}, {0}, {1}}),.s = &("Mxz")
};

static const MPolyData model_Myz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[3][1]) {{0}, {1}, {1}}),.s = &("Myz")
};

//
// trace and diff of stress tensor in 2D
static const MPolyData model_trMab2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){1.0, 1.0}),
  .e = &((int[2][2]) {{2, 0}, {0, 2}}),.s = &("Mxx+Myy")
};

static const MPolyData model_DifMab2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){1.0, -1.0}),
  .e = &((int[2][2]) {{2, 0}, {0, 2}}),.s = &("Mxx-Myy")
};

// kinetic energy 
static const MPolyData model_ekin1D = {
  .ndim = 1,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[1][1]) {{2}}),.s = &("ekin")
};

static const MPolyData model_ekin2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){.5, .5}),
  .e = &((int[2][2]) {{2, 0}, {0, 2}}),.s = &("ekin")
};

static const MPolyData model_ekin3D = {
  .ndim = 3,.nc = 3,
  .c = &((schnaps_real[]){.5, .5, .5}),
  .e = &((int[3][3]) {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}}),.s = &("ekin")
};

// Third order Moments (without 0.5 factors)
// Energy flux Individual 3rd (Qabc) order Components and contraction to kinetic energy flux vector components qa= (Qaaa+Qabb+Qacc) 
static const MPolyData model_Mxxx2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{3}, {0}}),.s = &("Mxxx")
};

static const MPolyData model_Mxyy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{1}, {2}}),.s = &("Mxyy")
};

static const MPolyData model_Myyy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{0}, {3}}),.s = &("Myyy")
};

static const MPolyData model_Myxx2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){1.}),
  .e = &((int[2][1]) {{2}, {1}}),.s = &("Myxx")
};

static const MPolyData model_Mabcx2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){1., 1.}),
  .e = &((int[2][2]) {{3, 0}, {1, 2}}),.s = &("Mxxx+Mxyy")
};

static const MPolyData model_Mabcy2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){1., 1.}),
  .e = &((int[2][2]) {{0, 3}, {2, 1}}),.s = &("Myxx+Myyy")
};

//
static const MPolyData model_Qxxx1D = {
  .ndim = 1,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[1][1]) {{3}}),.s = &("Qxxx")
};

static const MPolyData model_qx1D = {
  .ndim = 1,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[1][1]) {{3}}),.s = &("qx")
};

//
static const MPolyData model_Qxxx2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[2][1]) {{3}, {0}}),.s = &("Qxxx")
};

static const MPolyData model_Qxyy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[2][1]) {{1}, {2}}),.s = &("Qxyy")
};

static const MPolyData model_Qyyy2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[2][1]) {{0}, {3}}),.s = &("Qyyy")
};

static const MPolyData model_Qyxx2D = {
  .ndim = 2,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[2][1]) {{2}, {1}}),.s = &("Qyxx")
};

static const MPolyData model_qx2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){.5, .5}),
  .e = &((int[2][2]) {{3, 0}, {1, 2}}),.s = &("qx")
};

static const MPolyData model_qy2D = {
  .ndim = 2,.nc = 2,
  .c = &((schnaps_real[]){.5, .5}),
  .e = &((int[2][2]) {{0, 3}, {2, 1}}),.s = &("qy")
};

//
static const MPolyData model_Qxxx3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{3}, {0}, {0}}),.s = &("Qxxx")
};

static const MPolyData model_Qxyy3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{1}, {2}, {0}}),.s = &("Qxyy")
};

static const MPolyData model_Qxzz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{1}, {0}, {2}}),.s = &("Qxzz")
};

static const MPolyData model_Qyxx3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{2}, {1}, {0}}),.s = &("Qyxx")
};

static const MPolyData model_Qyzz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{0}, {1}, {2}}),.s = &("Qyzz")
};

static const MPolyData model_Qzxx3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{2}, {0}, {1}}),.s = &("Qzxx")
};

static const MPolyData model_Qzyy3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{0}, {2}, {1}}),.s = &("Qzyy")
};

static const MPolyData model_Qzzz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{0}, {0}, {3}}),.s = &("Qzzz")
};

static const MPolyData model_Qxyz3D = {
  .ndim = 3,.nc = 1,
  .c = &((schnaps_real[]){.5}),
  .e = &((int[3][1]) {{1}, {1}, {1}}),.s = &("Qxyz")
};

//
static const MPolyData model_qx3D = {
  .ndim = 3,.nc = 3,
  .c = &((schnaps_real[]){.5, .5, .5}),
  .e = &((int[3][3]) {{3, 0, 0}, {1, 2, 0}, {1, 0, 2}}),.s = &("qx")
};

static const MPolyData model_qy3D = {
  .ndim = 3,.nc = 3,
  .c = &((schnaps_real[]){.5, .5, .5}),
  .e = &((int[3][3]) {{0, 3, 0}, {2, 1, 0}, {0, 1, 2}}),.s = &("qy")
};

static const MPolyData model_qz3D = {
  .ndim = 3,.nc = 3,
  .c = &((schnaps_real[]){.5, .5, .5}),
  .e = &((int[3][3]) {{0, 0, 3}, {2, 0, 1}, {0, 2, 1}}),.s = &("qz")
};

//
// fourth order
static const MPolyData model_V42D = {
  .ndim = 2,.nc = 3,
  .c = &((schnaps_real[]){1.0, 1.0, 2.0}),
  .e = &((int[2][3]) {{4, 0, 2}, {0, 4, 2}}),.s = &("Mxxxx+Myyyy+2MxxMyy")
};

//
// 1D nodes
static const schnaps_real LBM_D1Q3_nodes[3][1] = { {-1.}, {0.}, {1.} };
static const int LBM_D1Q3_iopposite[3] = { 2, 1, 0 };

// 2D nodes
static const schnaps_real LBM_D2Q5_nodes[5][2] = {
  {0., 0.},
  {1., 0.}, {0., 1.}, {-1., 0.}, {0., -1.},
};
static const int LBM_D2Q5_iopposite[5] = { 0, 3, 4, 1, 2 };

static const schnaps_real LBM_D2Q9_nodes[9][2] = {
  {0., 0.},
  {1., 0.}, {0., 1.}, {-1., 0.}, {0., -1.},
  {1., 1.}, {-1., 1.}, {-1., -1.}, {1., -1.}
};
static const int LBM_D2Q9_iopposite[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

// 3D nodes
// TODO fill up 

// weights collections for equilibrium functions
static const schnaps_real LBM_WEIGHTS_D1Q3_ISOTH[3]={4./6.,1./6.,1./6.};

static const schnaps_real LBM_WEIGHTS_D2Q9_ISOTH[9] =
    { 4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36.,
  1. / 36., 1. / 36.
};
static const schnaps_real LBM_WEIGHTS_D2Q5_ISOTH[5] =
    { 1. / 3., 1. / 6., 1. / 6., 1. / 6., 1. / 6. };
#endif

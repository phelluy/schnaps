#include "geometry.h"
#include "global.h"
#include<stdio.h>
#include <math.h>
#include <assert.h>
const int h20_refnormal[6][3]={{0,-1,0},
			       {1,0,0},
			       {0,1,0},
			       {-1,0,0},
			       {0,0,1},
			       {0,0,-1} };

// 20 nodes of the reference element
const double h20_ref_node[20][3]={
  0  ,0  ,0  ,
  1  ,0  ,0  ,
  1  ,1  ,0  ,
  0  ,1  ,0  ,
  0  ,0  ,1  ,
  1  ,0  ,1  ,
  1  ,1  ,1  ,
  0  ,1  ,1  ,
  0.5,0  ,0  ,
  0  ,0.5,0  ,
  0  ,0  ,0.5,
  1  ,0.5,0  ,
  1  ,0  ,0.5,
  0.5,1  ,0  ,
  1  ,1  ,0.5,
  0  ,1  ,0.5,
  0.5,0  ,1  ,
  0  ,0.5,1  ,
  1  ,0.5,1  ,
  0.5,1  ,1
};

// Return the dot-product of the doubles a[3] and b[3]
double dot_product(double a[3], double b[3])
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Return the dot-product of the doubles a[3] and b[3]
double norm(double a[3])
{
  return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

double Dist(double a[3], double b[3]) 
{
  double d[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  return norm(d);
}

void PrintPoint(double x[3]) 
{
  printf("%f %f %f\n", x[0], x[1], x[2]);
}

void GeomRef2Phy(Geom* g) 
{
  Ref2Phy(g->physnode, g->xref, g->dphiref,
          g->ifa, g->xphy, g->dtau, g->codtau,
          g->dphi, g->vnds);
  g->det 
    = g->codtau[0][0] * g->dtau[0][0] 
    + g->codtau[0][1] * g->dtau[0][1]
    + g->codtau[0][2] * g->dtau[0][2];
};

void Ref2Phy(double physnode[20][3],
             double xref[3],
             double dphiref[3],
             int ifa,
             double xphy[3],
             double dtau[3][3],
             double codtau[3][3],
             double dphi[3],
             double vnds[3]) 
{
  // compute the mapping and its jacobian
  double x = xref[0];
  double y = xref[1];
  double z = xref[2];

  // Gradient of the shape functions and value (4th component) of the
  // shape functions
  double gradphi[20][4];
#include "h20phi.h"  // this file fills the values of phi and gradphi

  if (xphy != NULL) {
    for(int ii = 0; ii < 3; ii++) {
      xphy[ii]=0;
      for(int i=0;i<20;i++) {
	xphy[ii] += physnode[i][ii] * gradphi[i][3];
      }
    }
  }

  if (dtau != NULL) {
    for(int ii = 0; ii < 3; ii++) {
      for(int jj = 0; jj < 3; jj++) {
	dtau[ii][jj]=0;
      }
      for(int i = 0; i < 20; i++) {
	for(int jj = 0; jj < 3; jj++) {
	  dtau[ii][jj] += physnode[i][ii] * gradphi[i][jj];;
	}
      }
    }
  }

  if (codtau != NULL) {
    assert(dtau != NULL);
    /* codtau[0] =  dtau[4] * dtau[8] - dtau[5] * dtau[7]; */
    /* codtau[1] = -dtau[3] * dtau[8] + dtau[5] * dtau[6]; */
    /* codtau[2] =  dtau[3] * dtau[7] - dtau[4] * dtau[6]; */
    /* codtau[3] = -dtau[1] * dtau[8] + dtau[2] * dtau[7]; */
    /* codtau[4] =  dtau[0] * dtau[8] - dtau[2] * dtau[6]; */
    /* codtau[5] = -dtau[0] * dtau[7] + dtau[1] * dtau[6]; */
    /* codtau[6] =  dtau[1] * dtau[5] - dtau[2] * dtau[4]; */
    /* codtau[7] = -dtau[0] * dtau[5] + dtau[2] * dtau[3]; */
    /* codtau[8] =  dtau[0] * dtau[4] - dtau[1] * dtau[3]; */
    codtau[0][0] =  dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
    codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
    codtau[0][2] =  dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
    codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
    codtau[1][1] =  dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
    codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
    codtau[2][0] =  dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
    codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
    codtau[2][2] =  dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
  }

  if (dphi != NULL) {
    assert(codtau != NULL);
    for(int ii = 0; ii < 3; ii++) {
      dphi[ii]=0;
      for(int jj = 0; jj < 3; jj++) {
        dphi[ii]+=codtau[ii][jj]*dphiref[jj];
      }
    }
  }

  if (vnds !=NULL) {
    assert(codtau != NULL);
    assert(ifa >= 0);
    for(int ii = 0; ii < 3; ii++) {
      vnds[ii] = 0.0;
      for(int jj = 0; jj < 3; jj++) {
        vnds[ii] += codtau[ii][jj] * h20_refnormal[ifa][jj];
      }
    }

  }

}

void GeomPhy2Ref(Geom* g) 
{
  Phy2Ref(g->physnode,g->xphy,g->xref);
}

void Phy2Ref(double physnode[20][3], double xphy[3], double xref[3]) 
{
#define ITERNEWTON 10

  double dtau[3][3], codtau[3][3];
  double dxref[3], dxphy[3];
  double det;
  int ifa =- 1;
  xref[0] = 0.5;
  xref[1] = 0.5;
  xref[2] = 0.5;

  for(int iter = 0; iter < ITERNEWTON; ++iter) {
    Ref2Phy(physnode, xref, 0,ifa, dxphy, dtau, codtau, 0,0);
    dxphy[0] -= (xphy)[0];
    dxphy[1] -= (xphy)[1];
    dxphy[2] -= (xphy)[2];
    det = dot_product(dtau[0], codtau[0]);
    assert(det > 0);

    for(int ii = 0; ii < 3; ii ++ ) {
      dxref[ii] = 0;
      for(int jj = 0; jj < 3; jj ++ ) {
        dxref[ii] += codtau[jj][ii] * dxphy[jj];
      }
      xref[ii] -= dxref[ii] / det;
    }
  }

  double eps = 1e-2;  // may be to constraining...
  assert(xref[0] < 1 + eps && xref[0] > -eps);
  assert(xref[1] < 1 + eps && xref[1] > -eps);
  assert(xref[2] < 1 + eps && xref[2] > -eps);
}

void RobustPhy2Ref(double physnode[20][3], double xphy[3], double xref[3]) 
{
#define _ITERNEWTON 5
#define _NTHETA 3

  double dtau[3][3], codtau[3][3];
  double dxref[3], dxphy[3];
  int ifa =- 1;
  xref[0] = 0.5;
  xref[1] = 0.5;
  xref[2] = 0.5;

  // compute a linear transformation
  // based on nodes 0,1,3,4
  double M[3][3]={
    physnode[1][0]-physnode[0][0],
    physnode[3][0]-physnode[0][0],
    physnode[4][0]-physnode[0][0],

    physnode[1][1]-physnode[0][1],
    physnode[3][1]-physnode[0][1],
    physnode[4][1]-physnode[0][1],

    physnode[1][2]-physnode[0][2],
    physnode[3][2]-physnode[0][2],
    physnode[4][2]-physnode[0][2],

  };
  

  // initial affine hexaedron 
  double physnode0[20][3];
  for(int ino=0;ino<20;ino++){
    for(int ii=0;ii<3;ii++){
      physnode0[ino][ii]=physnode[0][ii];
      for (int jj=0;jj<3;jj++){
	physnode0[ino][ii]+=M[ii][jj]*h20_ref_node[ino][jj];
      }
    }
  }

  // debug: test the initial points.
  //int itest;
  //itest=0;
  //itest=1;
  //itest=3;
  //itest=4;
  //for(int ii=0;ii<3;ii++){
  //  assert(fabs(physnode0[itest][ii]-physnode[itest][ii])<1e-11);
  //}


  // homotopy path
  // theta=0 -> affine mapping
  // theta=1 -> full nonlinear mapping
  double physnode1[20][3];
  double dtheta=1./_NTHETA;

  for(int itheta=0;itheta<=_NTHETA;itheta++){
    double theta=itheta*dtheta;
    // intermediate curved hexaedron
    for(int ino=0;ino<20;ino++){
      for(int ii=0;ii<3;ii++){
	physnode1[ino][ii]=theta*physnode[ino][ii]+(1-theta)*physnode0[ino][ii];
      }
    }
  

    for(int iter = 0; iter < _ITERNEWTON; ++iter) {
      Ref2Phy(physnode1, xref, 0,ifa, dxphy, dtau, codtau, 0,0);
      dxphy[0] -= (xphy)[0];
      dxphy[1] -= (xphy)[1];
      dxphy[2] -= (xphy)[2];
      double det = dot_product(dtau[0], codtau[0]);
      assert(det > 0);

      for(int ii = 0; ii < 3; ii ++ ) {
	dxref[ii] = 0;
	for(int jj = 0; jj < 3; jj ++ ) {
	  dxref[ii] += codtau[jj][ii] * dxphy[jj];
	}
	xref[ii] -= dxref[ii] / det;
      }
      //printf("iter= %d dxref=%f %f %f\n",iter,dxref[0],dxref[1],dxref[2]);
    }
  }

}

#include "mhd.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#pragma start_opencl
// FIXME: Do not set variables with #define.
#define _CH (5)
#define _GAM (1.666666666666)
#pragma end_opencl

// FIXME: documentation???
// FIXME: Do not set variables with #define.

#pragma start_opencl
void conservatives(double* y, double* w) {
  double gam = _GAM; // FIXME: this should not be a #define

  w[0] = y[0];
  w[1] = y[0]*y[1];
  w[2] = y[2]/(gam-1) + y[0]*(y[1]*y[1]+y[3]*y[3]+y[4]*y[4])/2
    + (y[7]*y[7]+y[5]*y[5]+y[6]*y[6])/2;
  w[3] = y[0]*y[3];
  w[4] = y[0]*y[4];
  w[5] = y[5];
  w[6] = y[6];
  w[7] = y[7];        // Bx
  w[8] = y[8];        // psi
}
#pragma end_opencl

#pragma start_opencl
// FIXME: documentation????
void primitives(double *W, double *Y) {
  double gam = _GAM;

  Y[0] = W[0];
  Y[1] = W[1]/W[0];
  Y[2] = (gam - 1) * (W[2] - W[0] * (W[1]/W[0]*W[1]/W[0]
				     + W[3]/W[0]*W[3]/W[0]
				     + W[4]/W[0]*W[4]/W[0])/2
		      - (W[7]*W[7]+W[5]*W[5]+W[6]*W[6])/2);
  Y[3] = W[3] / W[0];
  Y[4] = W[4] / W[0];
  Y[5] = W[5];
  Y[6] = W[6];
  Y[7] = W[7]; // Bx
  Y[8] = W[8]; // psi
}
#pragma end_opencl

// FIXME: documentation
void jacobmhd(double *W, double *vn, double M[9][9]) {
  int i, j;
  double gam = _GAM;
  double Y[9];

  double rho, ux, uy, uz, by, bz, p, bx;

  for(i = 0; i < 9; i++){
    for(j = 0; j < 9; j++) {
      M[i][j] = 0;
    }
  }

  primitives(W,Y);

  rho = Y[0];
  ux = Y[1];
  p = Y[2];
  uy = Y[3];
  uz = Y[4];
  by = Y[5];
  bz = Y[6];
  bx = Y[7];

  M[0][0] = ux * vn[0] + uy * vn[1] + uz * vn[2];
  M[0][1] = rho * vn[0];
  M[0][3] = rho * vn[1];
  M[0][4] = rho * vn[2];

  M[1][1] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[1][2] = 1/rho*vn[0];
  M[1][5] = -(vn[1]*bx-by*vn[0])/rho;
  M[1][6] = -(vn[2]*bx-bz*vn[0])/rho;
  M[1][7] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;

  M[2][1] = gam*p*vn[0];
  M[2][2] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[2][3] = gam*p*vn[1];
  M[2][4] = gam*p*vn[2];
  M[2][5] = -ux*vn[1]*bx-uz*bz*vn[1]-uy*by*vn[1]+uy*by*vn[1]*gam+\
    ux*gam*vn[1]*bx+uz*gam*bz*vn[1];
  M[2][6] = -uz*vn[2]*bz+ux*gam*bx*vn[2]+uy*gam*by*vn[2]+uz*bz*vn[2]*gam-\
    ux*vn[2]*bx-uy*by*vn[2];
  M[2][7] = -ux*vn[0]*bx-uy*by*vn[0]-uz*bz*vn[0]+ux*bx*vn[0]*gam+\
    uy*gam*vn[0]*by+uz*gam*bz*vn[0];
  M[2][8] = -bx*vn[0]*gam+bx*vn[0]-by*vn[1]*gam+by*vn[1]\
    -bz*vn[2]*gam+bz*vn[2];

  M[3][2] = 1/rho*vn[1];
  M[3][3] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[3][5] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;
  M[3][6] = -(vn[2]*by-bz*vn[1])/rho;
  M[3][7] = (vn[1]*bx-by*vn[0])/rho;

  M[4][2] = 1/rho*vn[2];
  M[4][4] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[4][5] = (vn[2]*by-bz*vn[1])/rho;
  M[4][6] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;
  M[4][7] = (vn[2]*bx-bz*vn[0])/rho;

  M[5][1] = by*vn[0];
  M[5][3] = -bx*vn[0]-bz*vn[2];
  M[5][4] = vn[2]*by;
  M[5][5] = ux*vn[0]+uz*vn[2];
  M[5][6] = -vn[2]*uy;
  M[5][7] = -vn[0]*uy;
  M[5][8] = vn[1];

  M[6][1] = bz*vn[0];
  M[6][3] = bz*vn[1];
  M[6][4] = -vn[0]*bx-by*vn[1];
  M[6][5] = -vn[1]*uz;
  M[6][6] = ux*vn[0]+uy*vn[1];
  M[6][7] = -vn[0]*uz;
  M[6][8] = vn[2];

  M[7][1] = -vn[1]*by-bz*vn[2];
  M[7][3] = bx*vn[1];
  M[7][4] = bx*vn[2];
  M[7][5] = -vn[1]*ux;
  M[7][6] = -vn[2]*ux;
  M[7][7] = uy*vn[1]+uz*vn[2];
  M[7][8] = vn[0];

  M[8][5] = _CH*_CH*vn[1];
  M[8][6] = _CH*_CH*vn[2];
  M[8][7] = _CH*_CH*vn[0];

  for(i = 0; i < 9; i++) {
    for(j = 0; j < 9; j++) {
      M[i][j] /= _CH;
    }
  }
}

// Matrix-Vector multiplication
// FIXME: [] is not suitable for OpenCL
void matrix_vector(double A[9][9], double B[9], double *C) {
  int i, j;
  for(i = 0; i < 9; i++) {
    C[i] = 0;
    for(j = 0; j < 9; j++) {
      C[i] += A[i][j]*B[j];
    }
  }
}

#pragma start_opencl
void fluxnum(double *W,double *vn, double *flux) {

  double gam = _GAM;

  double un = W[1]/W[0]*vn[0]+W[3]/W[0]*vn[1]+W[4]/W[0]*vn[2];
  double bn = W[7]*vn[0]+W[5]*vn[1]+W[6]*vn[2];

  double p = (gam-1)*(W[2] - W[0]*(W[1]/W[0]*W[1]/W[0] + W[3]/W[0]*W[3]/W[0]\
				   + W[4]/W[0]*W[4]/W[0])/2 - (W[7]*W[7]+W[5]*W[5]+W[6]*W[6])/2);

  flux[0] = W[0]*un;
  flux[1] = W[0]*un*W[1]/W[0] + (p + (W[7]*W[7] + W[5]*W[5] + W[6]*W[6])/2)\
    *vn[0] - bn*W[7];
  flux[2] = (W[2] + p + (W[7]*W[7] + W[5]*W[5] + W[6]*W[6])/2)*un\
    - (W[7]*W[1]/W[0] + W[5]*W[3]/W[0] + W[6]*W[4]/W[0])*bn;
  flux[3] = W[0]*un*W[3]/W[0] + (p + (W[7]*W[7] + W[5]*W[5]\
				      + W[6]*W[6])/2)*vn[1] - bn*W[5];
  flux[4] = W[0]*un*W[4]/W[0] + (p + (W[7]*W[7] + W[5]*W[5]\
				      + W[6]*W[6])/2)*vn[2] - bn*W[6];

  flux[5] = -bn*W[3]/W[0] + un*W[5] + W[8]*vn[1];
  flux[6] = -bn*W[4]/W[0] + un*W[6] + W[8]*vn[2];
  flux[7] = -bn*W[1]/W[0] + un*W[7] + W[8]*vn[0];

  flux[8] = _CH*_CH*bn;
}
#pragma end_opencl

#pragma start_opencl
void MHDNumFluxRusanov(double *wL, double *wR,double *vnorm, double *flux) {
  int i;
  double fluxL[9];
  double fluxR[9];
  fluxnum(wL, vnorm, fluxL);
  fluxnum(wR, vnorm, fluxR);

  for(i = 0; i < 9; i++){
    flux[i] = (fluxL[i] + fluxR[i]) / 2 - _CH * (wR[i] - wL[i]) / 2;
  }
}
#pragma end_opencl

void MHDNumFluxP2(double *wL, double *wR, double *vn, double *flux) {

  int i, j;
  double wmil[9];
  double wRmwL[9];
  double Z[9];

  double M[9][9];
  double M2[9][9];

  // Initialize matrix to 0
  for(i = 0; i < 9 ; i++) {
    for(j = 0; j < 9 ; j++) {
      M[i][j] = 0;
      M2[i][j] = 0;
    }
  }

  // Compute the middle state
  for (i = 0; i < 9; i++) {
    wmil[i] = (wL[i] + wR[i])/2;
  }

  // calcul de la matrice M
  jacobmhd(wmil, vn, M);

  double coef[3] = {1./2, 0., 1./2};

  // calcul de (wR-wL)
  for (i = 0; i < 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  double dabs[9];

  for (i = 0; i < 9 ; i++){
    dabs[i] = coef[2] * wRmwL[i];
  }

  double Mw[9];
  for(i = 1; i >= 0; i--){
    matrix_vector(M, dabs, Mw);
    for(j = 0; j < 9; j++){
      dabs[j] = coef[i] * wRmwL[j] + Mw[j];
    }
  }

  double fluxL[9];
  double fluxR[9];

  fluxnum(wL, vn, fluxL);
  fluxnum(wR, vn, fluxR);

  for(i = 0; i < 9; i++){
    flux[i] = (fluxL[i] + fluxR[i]) / 2 - _CH * dabs[i] / 2;
  }
}

void MHDNumFlux1D(double *WL, double *WR, double *vn, double *flux) {
  
  double gam = _GAM;

  double piL, piR, piyL, piyR, pizL, pizR;
  double a, aR, aL, al0, ar0;
  double cf, cfL, cfR, cL, cR;
  double Xl, Xr, pxl, pxr;
  double alpha;
  double cA, cB, b2;
  double us, pis, uys, uzs, piys, pizs;
  double pi, piy, piz;
  double sigma1, sigma2, sigma3;

  double YL[9];
  double YR[9];
  double ymil[9];
  double ystar[9];

  int i;

  primitives(WL, YL);
  primitives(WR, YR);

  double b = YL[7]; // En 1D BX est constant dans YL = YR

  // calcul des parametres issus de la relaxation
  piL = YL[2] + 0.5 * fabs(YL[5]*YL[5] + YL[6]*YL[6]) - 0.5*b*b;
  piyL = -b*YL[5];
  pizL = -b*YL[6];

  piR = YR[2] + 0.5 * fabs(YR[5]*YR[5] + YR[6]*YR[6]) - 0.5*b*b;
  piyR = -b*YR[5];
  pizR = -b*YR[6];


  for(i=0; i<9; i++){
    ymil[i] = 0.5*(YL[i]+YR[i]);
  }

  a = sqrt(gam*ymil[2]/ymil[0]);
  cf = sqrt(\
            0.5*((b*b + ymil[5]*ymil[5] + ymil[6]*ymil[6])/(ymil[0])+a*a)\
            + sqrt(\
                   0.25*(pow((b*b + ymil[5]*ymil[5] + ymil[6]*ymil[6])/(ymil[0])+a*a,2))-(a*a*b*b)/(ymil[0])\
		   )\
	    );


  aL = sqrt(gam*YL[2]/YL[0]);
  cfL = sqrt(\
	     0.5*((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(YL[0])+aL*aL)\
	     + sqrt(\
		    0.25*(pow((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(YL[0])+aL*aL,2))-(aL*aL*b*b)/(YL[0])\
		    )\
	     );


  aR = sqrt(gam*YR[2]/YR[0]);
  cfR = sqrt(\
	     0.5*((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(YR[0])+aR*aR)\
	     + sqrt(\
		    0.25*(pow((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(YR[0])+aR*aR,2))-(aR*aR*b*b)/(YR[0])\
		    )\
	     );


  // calcul des vitesses relaxees a gauche et a droite
  alpha = (gam-1)/2.;

  Xl = (fmax(YL[1]-YR[1], 0.0) + (fmax(piR-piL, 0.0))/(YL[0]*cfL+YR[0]*cfR))/cfL;
  Xr = (fmax(YL[1]-YR[1], 0.0) + (fmax(piL-piR, 0.0))/(YL[0]*cfL+YR[0]*cfR))/cfR;

  pxl = 1 - Xl/(1+alpha*Xl);
  pxr = 1 - Xr/(1+alpha*Xr);

  al0 = sqrt(\
	     0.5*((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(YL[0]*pxl)+aL*aL)\
	     + sqrt(\
		    0.25*(pow((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(YL[0]*pxl)+aL*aL,2))-(aL*aL*b*b)/(YL[0]*pxl)\
		    )\
	     );


  ar0 = sqrt(\
	     0.5*((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(YR[0]*pxr)+aR*aR)\
	     + sqrt(\
		    0.25*(pow((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(YR[0]*pxr)+aR*aR,2))-(aR*aR*b*b)/(YR[0]*pxr)\
		    )\
	     );

  cL = al0*YL[0] + alpha*YL[0]*(fmax(YL[1]-YR[1],0.0) + (fmax(piR-piL,0.0))/(YL[0]*cfL+YR[0]*cfR));
  cR = ar0*YR[0] + alpha*YR[0]*(fmax(YL[1]-YR[1],0.0) + (fmax(piL-piR,0.0))/(YL[0]*cfL+YR[0]*cfR));

  // pour le 3-ondes ondes on prend des vitesses simples
  cA = cf;
  cB = cf;
  b2 = 0.0;

  // calcul des etats intermediaires
  us = (cL*YL[1] + cR*YR[1] + piL-piR)/(cL+cR);
  pis = (cR*piL + cL*piR - cL*cR*(YR[1]-YL[1]))/(cL+cR);

  uys = (cL*YL[3] + cR*YR[3] + piyL-piyR)/(cL+cR);
  uzs = (cL*YL[4] + cR*YR[4] + pizL-pizR)/(cL+cR);

  piys = (cR*piyL + cL*piyR - cL*cR*(YR[3]-YL[3]))/(cL+cR);
  pizs = (cR*pizL + cL*pizR - cL*cR*(YR[4]-YL[4]))/(cL+cR);

  // calcul des vitesses caracteristiques
  sigma1 = YL[1] - cL/YL[0];
  sigma2 = us;
  sigma3 = YR[2] + cR/YR[0];

  // decentrement
  if(sigma1 > 0.0){
    ystar[0] = YL[0];
    ystar[1] = YL[1];
    ystar[2] = YL[2];
    ystar[3] = YL[3];
    ystar[4] = YL[4];
    ystar[5] = YL[5];
    ystar[6] = YL[6];
    pi = piL;
    piy = piyL;
    piz = pizL;
  }
  else
    if(sigma2 > 0.0){
      ystar[0] = 1.0/(1.0/YL[0] + (piL-piR+cR*(YR[1]-YL[1]))/(cL*(cL+cR)));
      ystar[5] = ystar[0]*(YL[5]/YL[0] + b/(cL*cL)*piyL - b/(cL*cL)*piys);
      ystar[6] = ystar[0]*(YL[6]/YL[0] + b/(cL*cL)*pizL - b/(cL*cL)*pizs);
      ystar[1] = us;
      ystar[3] = uys;
      ystar[4] = uzs;
      pi = pis;
      piy = piys;
      piz = pizs;
      ystar[2] = ystar[0]*(gam-1)*(YL[2]/((gam-1)*YL[0])\
				   + ((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(2.0*YL[0]))\
				   - piL*piL/(2.0*cL*cL) - (piyL*piyL+pizL*pizL)/(2.0*cL*cL)\
				   - ((b*b+ystar[5]*ystar[5]+ystar[6]*ystar[6])/(2.0*ystar[0]))\
				   + pis*pis/(2.0*cL*cL) + (piys*piys+pizs*pizs)/(2.0*cL*cL));
    }
    else
      if(sigma3 > 0.0){
	ystar[0] = 1.0/(1.0/YR[0] + (piR-piL+cL*(YR[1]-YL[1]))/(cR*(cL+cR)));
	ystar[5] = ystar[0]*(YR[5]/YR[0] + b/(cR*cR)*piyR - b/(cR*cR)*piys);
	ystar[6] = ystar[0]*(YR[6]/YR[0] + b/(cR*cR)*pizR - b/(cR*cR)*pizs);
	ystar[1] = us;
	ystar[3] = uys;
	ystar[4] = uzs;
	pi = pis;
	piy = piys;
	piz = pizs;
	ystar[2] = ystar[0]*(gam-1)*(YR[2]/((gam-1)*YR[0])\
				     + ((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(2.0*YR[0]))\
				     - piR*piR/(2.0*cR*cR) - (piyR*piyR+pizR*pizR)/(2.0*cR*cR)\
				     - ((b*b+ystar[5]*ystar[5]+ystar[6]*ystar[6])/(2.0*ystar[0]))\
				     + pis*pis/(2.0*cR*cR) + (piys*piys+pizs*pizs)/(2.0*cR*cR));
      }
      else{
	ystar[0] = YR[0];
	ystar[1] = YR[1];
	ystar[2] = YR[2];
	ystar[3] = YR[3];
	ystar[4] = YR[4];
	ystar[5] = YR[5];
	ystar[6] = YR[6];
	pi = piR;
	piy = piyR;
	piz = pizR;
      }

  flux[0] = ystar[0]*ystar[1];
  flux[1] = ystar[0]*(ystar[1]*ystar[1]) + pi;
  flux[2] = ystar[1]*(0.5*ystar[0]*(ystar[1]*ystar[1]+ystar[3]*ystar[3]+ystar[4]*ystar[4])\
                      + ystar[2]/(gam-1) + 0.5*(b*b+ystar[5]*ystar[5]+ystar[6]*ystar[6]) + pi)\
    + piy*ystar[3] + piz*ystar[4];
  flux[3]  = ystar[0]*ystar[1]*ystar[3] + piy;
  flux[4]  = ystar[0]*ystar[1]*ystar[4] + piz;
  flux[5]  = ystar[1]*ystar[5] - b*ystar[3];
  flux[6]  = ystar[1]*ystar[6] - b*ystar[4];
  flux[7]  = 0;
  flux[8]  = 0;
}

#pragma start_opencl
void MHDBoundaryFlux(double *x, double t, double *wL, double *vnorm,
		     double *flux) {
  int i;
  double wR[9];

  if(vnorm[1] > 0.0001 || vnorm[1] < -0.0001){
    MHDImposedData(x,t,wR);
  }
  else if(vnorm[0] > 0.0001 || vnorm[0] < -0.0001){
    for(i=0; i<9; i++){
      wR[i] = wL[i];
    }
  }
  else{
    //printf("Error in MHDBoundaryFlux !\n");
    //printf("vnorm = %f %f %f\n", vnorm[0], vnorm[1], vnorm[2]);
    //assert(1==2);
  }
  MHDNumFluxRusanov(wL,wR,vnorm,flux);
}
#pragma end_opencl

#pragma start_opencl
void MHDInitData(double *x, double *w) {

  double t=0;
  MHDImposedData(x,t,w);

}
#pragma end_opencl

#pragma start_opencl
void MHDImposedData(double *x, double t, double *w) {
  int i;
  double gam = _GAM;

  double yL[9], yR[9];
  double wL[9], wR[9];

  
//  yL[0] = 3.;
//  yL[1] = 1.3;
//  yL[3] = 0.;
//  yL[4] = 0.;
//  yL[2] = 3.;
//  yL[5] = 1.;
//  yL[7] = 1.;
//  yL[7] = 1.5;
//  yL[8] = 0.;
//
//  yR[0] = 1.;
//  yR[1] = 1.3;
//  yR[3] = 0.;
//  yR[4] = 0.;
//  yR[2] = 1.;
//  yR[5] = 0.0707372016677029;
//  yR[6] = 0.9974949866040544;
//  yR[7] = 1.5;
//  yR[8] = 0.;
//
//  conservatives(yL, wL);
//  conservatives(yR, wR);
//
//  if(x[0] < 5)
//    for(i=0; i<9; i++){
//      w[i] = wL[i];
//    }
//  else
//    for(i=0; i<9; i++){
//      w[i] = wR[i];
//    }

  yL[0] = gam*gam;
  yL[1] = -sin(x[1]);
  yL[2] = gam;
  yL[3] = sin(x[0]);
  yL[4] = 0.0;
  yL[5] = sin(2*(x[0]));
  yL[6] = 0.0;
  yL[7] = -sin(x[1]);
  yL[8] = 0.0;

  conservatives(yL, w);
  
}
#pragma end_opencl

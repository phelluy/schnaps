#include "mhd.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

//#pragma start_opencl
//// FIXME: Do not set variables with #define.
//#define _CH (6)
//#define _GAM (1.666666666666)
//#pragma end_opencl

#pragma start_opencl
void conservatives(schnaps_real *y, schnaps_real *w) {
  schnaps_real gam = 1.6666666666;

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
void primitives(schnaps_real *W, schnaps_real *Y) {
  schnaps_real gam = 1.6666666666;

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

#pragma start_opencl
void jacobmhd(schnaps_real* W,schnaps_real* vn, schnaps_real *M){
  schnaps_real gam = 1.6666666666;
  schnaps_real Y[9];

  schnaps_real rho, ux, uy, uz, by, bz, p, bx;

  for(int i = 0; i < 9; i++){
    for(int j = 0; j < 9; j++) {
      M[i*9+j] = 0;
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

  M[0*9+0] = ux * vn[0] + uy * vn[1] + uz * vn[2];
  M[0*9+1] = rho * vn[0];
  M[0*9+3] = rho * vn[1];
  M[0*9+4] = rho * vn[2];

  M[1*9+1] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[1*9+2] = 1/rho*vn[0];
  M[1*9+5] = -(vn[1]*bx-by*vn[0])/rho;
  M[1*9+6] = -(vn[2]*bx-bz*vn[0])/rho;
  M[1*9+7] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;

  M[2*9+1] = gam*p*vn[0];
  M[2*9+2] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[2*9+3] = gam*p*vn[1];
  M[2*9+4] = gam*p*vn[2];
  M[2*9+5] = -ux*vn[1]*bx-uz*bz*vn[1]-uy*by*vn[1]+uy*by*vn[1]*gam+\
    ux*gam*vn[1]*bx+uz*gam*bz*vn[1];
  M[2*9+6] = -uz*vn[2]*bz+ux*gam*bx*vn[2]+uy*gam*by*vn[2]+uz*bz*vn[2]*gam-\
    ux*vn[2]*bx-uy*by*vn[2];
  M[2*9+7] = -ux*vn[0]*bx-uy*by*vn[0]-uz*bz*vn[0]+ux*bx*vn[0]*gam+\
    uy*gam*vn[0]*by+uz*gam*bz*vn[0];
  M[2*9+8] = -bx*vn[0]*gam+bx*vn[0]-by*vn[1]*gam+by*vn[1]\
    -bz*vn[2]*gam+bz*vn[2];

  M[3*9+2] = 1/rho*vn[1];
  M[3*9+3] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[3*9+5] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;
  M[3*9+6] = -(vn[2]*by-bz*vn[1])/rho;
  M[3*9+7] = (vn[1]*bx-by*vn[0])/rho;

  M[4*9+2] = 1/rho*vn[2];
  M[4*9+4] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[4*9+5] = (vn[2]*by-bz*vn[1])/rho;
  M[4*9+6] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;
  M[4*9+7] = (vn[2]*bx-bz*vn[0])/rho;

  M[5*9+1] = by*vn[0];
  M[5*9+3] = -bx*vn[0]-bz*vn[2];
  M[5*9+4] = vn[2]*by;
  M[5*9+5] = ux*vn[0]+uz*vn[2];
  M[5*9+6] = -vn[2]*uy;
  M[5*9+7] = -vn[0]*uy;
  M[5*9+8] = vn[1];

  M[6*9+1] = bz*vn[0];
  M[6*9+3] = bz*vn[1];
  M[6*9+4] = -vn[0]*bx-by*vn[1];
  M[6*9+5] = -vn[1]*uz;
  M[6*9+6] = ux*vn[0]+uy*vn[1];
  M[6*9+7] = -vn[0]*uz;
  M[6*9+8] = vn[2];

  M[7*9+1] = -vn[1]*by-bz*vn[2];
  M[7*9+3] = bx*vn[1];
  M[7*9+4] = bx*vn[2];
  M[7*9+5] = -vn[1]*ux;
  M[7*9+6] = -vn[2]*ux;
  M[7*9+7] = uy*vn[1]+uz*vn[2];
  M[7*9+8] = vn[0];

  M[8*9+5] = 6*6*vn[1];
  M[8*9+6] = 6*6*vn[2];
  M[8*9+7] = 6*6*vn[0];

  for(int i = 0; i < 9; i++) {
    for(int j = 0; j < 9; j++) {
      M[i*9+j] /= 6;
    }
  }
  
}
#pragma end_opencl

// Matrix-Vector multiplication
// FIXME: [] is not suitable for OpenCL
#pragma start_opencl
void matrix_vector(schnaps_real *A, schnaps_real *B, schnaps_real *C) {
  for(int i = 0; i < 9; i++) {
    C[i] = 0;
    for(int j = 0; j < 9; j++) {
      C[i] += A[i*9+j]*B[j];
    }
  }
}
#pragma end_opencl

#pragma start_opencl
void fluxnum(schnaps_real *W,schnaps_real *vn, schnaps_real *flux) {

  schnaps_real gam = 1.6666666666;

  schnaps_real un = W[1]/W[0]*vn[0]+W[3]/W[0]*vn[1]+W[4]/W[0]*vn[2];
  schnaps_real bn = W[7]*vn[0]+W[5]*vn[1]+W[6]*vn[2];

  schnaps_real p = (gam-1)*(W[2] - W[0]*(W[1]/W[0]*W[1]/W[0]
				 + W[3]/W[0]*W[3]/W[0]
				 + W[4]/W[0]*W[4]/W[0])/2
		    - (W[7]*W[7]+W[5]*W[5]+W[6]*W[6])/2);

  flux[0] = W[0]*un;
  flux[1] = W[0]*un*W[1]/W[0] + (p + (W[7]*W[7] + W[5]*W[5] + W[6]*W[6])/2)
    *vn[0] - bn*W[7];
  flux[2] = (W[2] + p + (W[7]*W[7] + W[5]*W[5] + W[6]*W[6])/2)*un
    - (W[7]*W[1]/W[0] + W[5]*W[3]/W[0] + W[6]*W[4]/W[0])*bn;
  flux[3] = W[0]*un*W[3]/W[0] + (p + (W[7]*W[7] + W[5]*W[5]
				      + W[6]*W[6])/2)*vn[1] - bn*W[5];
  flux[4] = W[0]*un*W[4]/W[0] + (p + (W[7]*W[7] + W[5]*W[5]
				      + W[6]*W[6])/2)*vn[2] - bn*W[6];

  flux[5] = -bn*W[3]/W[0] + un*W[5] + W[8]*vn[1];
  flux[6] = -bn*W[4]/W[0] + un*W[6] + W[8]*vn[2];
  flux[7] = -bn*W[1]/W[0] + un*W[7] + W[8]*vn[0];

  flux[8] = 6*6*bn;
}
#pragma end_opencl

#pragma start_opencl
void MHDNumFluxRusanov(schnaps_real *wL, schnaps_real *wR,schnaps_real *vnorm, schnaps_real *flux)
{
  schnaps_real fluxL[9];
  schnaps_real fluxR[9];
  
  fluxnum(wL, vnorm, fluxL);
  fluxnum(wR, vnorm, fluxR);

  for(int i = 0; i < 9; i++) {
    flux[i] = (fluxL[i] + fluxR[i]) / 2 - 6 * (wR[i] - wL[i]) / 2;
  }
}
#pragma end_opencl

#pragma start_opencl
void MHDNumFluxP2(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux)
{
  schnaps_real wmil[9];
  schnaps_real wRmwL[9];

  schnaps_real M[81];

  // Initialize matrix to 0
  for(int i = 0; i < 9 ; i++) {
    for(int j = 0; j < 9 ; j++) {
      M[i*9+j] = 0;
    }
  }

  // Compute the middle state
  for (int i = 0; i < 9; i++) {
    wmil[i] = (wL[i] + wR[i])/2;
  }

  // calcul de la matrice M
  jacobmhd(wmil, vn, M);

  schnaps_real coef[3] = {1./2, 0., 1./2};

  // calcul de (wR-wL)
  for (int i = 0; i < 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  schnaps_real dabs[9];

  for (int i = 0; i < 9 ; i++){
    dabs[i] = coef[2] * wRmwL[i];
  }

  schnaps_real Mw[9];
  for(int i = 1; i >= 0; i--){
    matrix_vector(M, dabs, Mw);
    for(int j = 0; j < 9; j++){
      dabs[j] = coef[i] * wRmwL[j] + Mw[j];
    }
  }

  schnaps_real fluxL[9];
  schnaps_real fluxR[9];

  fluxnum(wL, vn, fluxL);
  fluxnum(wR, vn, fluxR);

  for(int i = 0; i < 9; i++){
    flux[i] = (fluxL[i] + fluxR[i]) / 2 - 6 * dabs[i] / 2;
  }
}
#pragma end_opencl

#pragma start_opencl
void MHDNumFlux1D(schnaps_real *WL, schnaps_real *WR, schnaps_real *vn, schnaps_real *flux)
{
  schnaps_real gam = 1.6666666666;

  schnaps_real piL, piR, piyL, piyR, pizL, pizR;
  schnaps_real a, aR, aL, al0, ar0;
  schnaps_real cfL, cfR, cL, cR;
  schnaps_real Xl, Xr, pxl, pxr;
  schnaps_real alpha;
  schnaps_real us, pis, uys, uzs, piys, pizs;
  schnaps_real pi, piy, piz;
  schnaps_real sigma1, sigma2, sigma3;

  schnaps_real YL[9];
  schnaps_real YR[9];
  schnaps_real ymil[9];
  schnaps_real ystar[9];

  primitives(WL, YL);
  primitives(WR, YR);

  schnaps_real b = YL[7]; // En 1D BX est constant dans YL = YR

  // calcul des parametres issus de la relaxation
  piL = YL[2] + 0.5*fabs(YL[5]*YL[5] + YL[6]*YL[6]) - 0.5*b*b;
  piyL = -b*YL[5];
  pizL = -b*YL[6];

  piR = YR[2] + 0.5 * fabs(YR[5]*YR[5] + YR[6]*YR[6]) - 0.5*b*b;
  piyR = -b*YR[5];
  pizR = -b*YR[6];


  for(int i=0; i<9; i++){
    ymil[i] = 0.5*(YL[i]+YR[i]);
  }

  a = sqrt(gam*ymil[2]/ymil[0]);
  schnaps_real cf = sqrt(\
		 0.5*((b*b + ymil[5]*ymil[5] + ymil[6]*ymil[6])/(ymil[0])+a*a)
		 + sqrt(
			0.25*(pow((b*b + ymil[5]*ymil[5]
				   + ymil[6]*ymil[6])/(ymil[0])+a*a,2))
			-(a*a*b*b)/(ymil[0])
			)
		 );


  aL = sqrt(gam*YL[2]/YL[0]);
  cfL = sqrt(\
	     0.5*((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(YL[0])+aL*aL)
	     + sqrt(
		    0.25*(pow((b*b + YL[5]*YL[5]
			       + YL[6]*YL[6])/(YL[0])+aL*aL,2))
		    -(aL*aL*b*b)/(YL[0]) 
		    )
	     );

  aR = sqrt(gam*YR[2]/YR[0]);
  cfR = sqrt(\
	     0.5*((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(YR[0])+aR*aR)
	     + sqrt(
		    0.25*(pow((b*b + YR[5]*YR[5]
			       + YR[6]*YR[6])/(YR[0])+aR*aR,2))
		    -(aR*aR*b*b)/(YR[0])
		    )
	     );

  // calcul des vitesses relaxees a gauche et a droite
  alpha = (gam-1)/2.;

  Xl = (fmax(YL[1]-YR[1], 0)
	+ (fmax(piR-piL, 0))/(YL[0]*cfL+YR[0]*cfR))/cfL;
  Xr = (fmax(YL[1]-YR[1], 0)
	+ (fmax(piL-piR, 0))/(YL[0]*cfL+YR[0]*cfR))/cfR;

  pxl = 1 - Xl/(1+alpha*Xl);
  pxr = 1 - Xr/(1+alpha*Xr);

  al0 = sqrt(\
	     0.5*((b*b + YL[5]*YL[5] + YL[6]*YL[6])/(YL[0]*pxl)+aL*aL)\
	     + sqrt(
		    0.25*(pow((b*b + YL[5]*YL[5]
			       + YL[6]*YL[6])/(YL[0]*pxl)+aL*aL,2))
		    -(aL*aL*b*b)/(YL[0]*pxl) 
		    )
	     );

  ar0 = sqrt(\
	     0.5*((b*b + YR[5]*YR[5] + YR[6]*YR[6])/(YR[0]*pxr)+aR*aR)
	     + sqrt(
		    0.25*(pow((b*b + YR[5]*YR[5]
			       + YR[6]*YR[6])/(YR[0]*pxr)+aR*aR,2))
		    -(aR*aR*b*b)/(YR[0]*pxr)
		    )
	     );

  cL = al0*YL[0] + alpha*YL[0]*(fmax(YL[1]-YR[1],0)
				+ (fmax(piR-piL,0))/(YL[0]*cfL+YR[0]*cfR));
  cR = ar0*YR[0] + alpha*YR[0]*(fmax(YL[1]-YR[1],0)
				+ (fmax(piL-piR,0))/(YL[0]*cfL+YR[0]*cfR));

  // pour le 3-ondes ondes on prend des vitesses simples
  //real cA = cf;
  //real cB = cf;
  //real b2 = 0;

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
  if(sigma1 > 0){
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
    if(sigma2 > 0){
      ystar[0] = 1.0/(1.0/YL[0] + (piL-piR+cR*(YR[1]-YL[1]))/(cL*(cL+cR)));
      ystar[5] = ystar[0]*(YL[5]/YL[0] + b/(cL*cL)*piyL - b/(cL*cL)*piys);
      ystar[6] = ystar[0]*(YL[6]/YL[0] + b/(cL*cL)*pizL - b/(cL*cL)*pizs);
      ystar[1] = us;
      ystar[3] = uys;
      ystar[4] = uzs;
      pi = pis;
      piy = piys;
      piz = pizs;
      ystar[2] = ystar[0]*(gam-1)*(YL[2]/((gam-1)*YL[0])
				   + ((b*b + YL[5]*YL[5]
				       + YL[6]*YL[6])/(2.0*YL[0]))
				   - piL*piL/(2.0*cL*cL)
				   - (piyL*piyL+pizL*pizL)/(2.0*cL*cL)
				   - ((b*b+ystar[5]*ystar[5]
				       +ystar[6]*ystar[6])/(2.0*ystar[0]))
				   + pis*pis/(2.0*cL*cL)
				   + (piys*piys+pizs*pizs)/(2.0*cL*cL));
    }
    else
      if(sigma3 > 0){
	ystar[0] = 1.0/(1.0/YR[0] + (piR-piL+cL*(YR[1]-YL[1]))/(cR*(cL+cR)));
	ystar[5] = ystar[0]*(YR[5]/YR[0] + b/(cR*cR)*piyR - b/(cR*cR)*piys);
	ystar[6] = ystar[0]*(YR[6]/YR[0] + b/(cR*cR)*pizR - b/(cR*cR)*pizs);
	ystar[1] = us;
	ystar[3] = uys;
	ystar[4] = uzs;
	pi = pis;
	piy = piys;
	piz = pizs;
	ystar[2] = ystar[0]*(gam-1)*(YR[2]/((gam-1)*YR[0])
				     + ((b*b + YR[5]*YR[5]
					 + YR[6]*YR[6])/(2.0*YR[0]))
				     - piR*piR/(2.0*cR*cR)
				     - (piyR*piyR+pizR*pizR)/(2.0*cR*cR)
				     - ((b*b+ystar[5]*ystar[5]
					 +ystar[6]*ystar[6])/(2.0*ystar[0]))
				     + pis*pis/(2.0*cR*cR)
				     + (piys*piys+pizs*pizs)/(2.0*cR*cR));
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
  flux[2] = ystar[1]*(0.5*ystar[0]*(ystar[1]*ystar[1]
				    +ystar[3]*ystar[3]
				    +ystar[4]*ystar[4])
                      + ystar[2]/(gam-1)
		      + 0.5*(b*b+ystar[5]*ystar[5]+ystar[6]*ystar[6]) + pi)
    + piy*ystar[3] + piz*ystar[4];
  flux[3]  = ystar[0]*ystar[1]*ystar[3] + piy;
  flux[4]  = ystar[0]*ystar[1]*ystar[4] + piz;
  flux[5]  = ystar[1]*ystar[5] - b*ystar[3];
  flux[6]  = ystar[1]*ystar[6] - b*ystar[4];
  flux[7]  = 0;
  flux[8]  = 0;
}
#pragma end_opencl

#pragma start_opencl
void MHDBoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
		     schnaps_real *flux) {
  schnaps_real wR[9];

  if(vnorm[1] > 0.0001 || vnorm[1] < -0.0001){
    MHDImposedData(x,t,wR);
  }
  else if(vnorm[0] > 0.0001 || vnorm[0] < -0.0001){
    for(int i=0; i<9; i++){
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
void MHDInitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  MHDImposedData(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void MHDImposedData(const schnaps_real *x,const  schnaps_real t, schnaps_real *w) {
  schnaps_real gam = 1.6666666666;

  schnaps_real yL[9];
  //real yR[9];
  //real wL[9];
  //real wR[9];
  
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
  //    for(int i=0; i<9; i++){
  //      w[i] = wL[i];
  //    }
  //  else
  //    for(int i=0; i<9; i++){
  //      w[i] = wR[i];
  //    }

  yL[0] = gam * gam;
  yL[1] = -sin(x[1]);
  yL[2] = gam;
  yL[3] = sin(x[0]);
  yL[4] = 0;
  yL[5] = sin(2*(x[0]));
  yL[6] = 0;
  yL[7] = -sin(x[1]);
  yL[8] = 0;

  conservatives(yL, w);
}
#pragma end_opencl





////////////////////////////////////////////////////////////////////////////////
//                               Choc Fort                                    //
////////////////////////////////////////////////////////////////////////////////

#pragma start_opencl
void MHDBoundaryFluxChocFort(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
		     schnaps_real *flux) {
  schnaps_real wR[9];

  if(vnorm[1] > 0.0001 || vnorm[1] < -0.0001){
    MHDImposedDataChocFort(x,t,wR);
  }
  else{
    for(int i=0; i<9; i++){
      wR[i] = wL[i];
    }
  }

  MHDNumFluxRusanov(wL,wR,vnorm,flux);
}
#pragma end_opencl


#pragma start_opencl
void MHDInitDataChocFort(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  MHDImposedDataChocFort(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void MHDImposedDataChocFort(const schnaps_real *x,const  schnaps_real t, schnaps_real *w) {
  schnaps_real gam = 1.6666666666;

  schnaps_real yL[9];
  schnaps_real yR[9];
  schnaps_real wL[9];
  schnaps_real wR[9];


  yL[0] = 3;
  yL[1] = 1.3;
  yL[3] = 0;
  yL[4] = 0;
  yL[2] = 3;
  yL[5] = 1;
  yL[6] = 1;
  yL[7] = 1.5;
  yL[8] = 0;
  
  yR[0] = 1;
  yR[1] = 1.3;
  yR[3] = 0;
  yR[4] = 0;
  yR[2] = 1;
  yR[5] = 0.0707372016677029;
  yR[6] = 0.9974949866040544;
  yR[7] = 1.5;
  yR[8] = 0;
  
  conservatives(yL, wL);
  conservatives(yR, wR);
  
  if(x[0] < 5)
    for(int i=0; i<9; i++){
      w[i] = wL[i];
    }
  else
    for(int i=0; i<9; i++){
      w[i] = wR[i];
    }

}
#pragma end_opencl




////////////////////////////////////////////////////////////////////////////////
//                               Orszag-Tang                                  //
////////////////////////////////////////////////////////////////////////////////

#pragma start_opencl
void MHDBoundaryFluxOrszagTang(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
                                 schnaps_real *flux) {
  schnaps_real wR[9];

  MHDImposedDataOrszagTang(x,t,wR);
  MHDNumFluxRusanov(wL,wR,vnorm,flux);

}
#pragma end_opencl

#pragma start_opencl
void MHDInitDataOrszagTang(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  MHDImposedDataOrszagTang(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void MHDImposedDataOrszagTang(const schnaps_real *x,const  schnaps_real t, schnaps_real *w) {
  schnaps_real gam = 1.6666666666;
  schnaps_real yL[9];

  yL[0] = gam*gam;
  yL[1] = -sin(x[1]);
  yL[2] = gam;
  yL[3] = sin(x[0]);
  yL[4] = 0;
  yL[5] = sin(2*(x[0]));
  yL[6] = 0;
  yL[7] = -sin(x[1]);
  yL[8] = 0;

  conservatives(yL, w);
  
}
#pragma end_opencl


////////////////////////////////////////////////////////////////////////////////
//                               Reconnexion                                  //
////////////////////////////////////////////////////////////////////////////////

#pragma start_opencl
void MHDBoundaryFluxReconnexion(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
                                 schnaps_real *flux) {
  schnaps_real wR[9];
  MHDImposedDataReconnexion(x,t,wR);
  MHDNumFluxRusanov(wL,wR,vnorm,flux);
}
#pragma end_opencl

#pragma start_opencl
void MHDInitDataReconnexion(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  MHDImposedDataReconnexion(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void MHDImposedDataReconnexion(const schnaps_real *x,const  schnaps_real t, schnaps_real *w) {
  //real gam = 1.6666666666;
  schnaps_real yL[9];

  yL[0] = 1.0;
  yL[1] = 0.1*sin(3.14159265359*(x[1]));
  yL[2] = 0.1;
  yL[3] = 0;
  yL[4] = 0;
  if((x[0])<(2./3.)){
    yL[5] = 1.0;
  }
  else if(((x[0])>4./3.)){
    yL[5] = 1.0;
  }
  else{
    yL[5] = -1.0;
  }
  yL[6] = 0;
  yL[7] = 0;
  yL[8] = 0;

  conservatives(yL, w);
  
}
#pragma end_opencl


////////////////////////////////////////////////////////////////////////////////
//                            Kelvin-Helmotz                                  //
////////////////////////////////////////////////////////////////////////////////

#pragma start_opencl
void MHDBoundaryFluxKelvinHelmotz(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
                                 schnaps_real *flux) {

  schnaps_real wR[9];

  if(vnorm[1] > 0.0001 || vnorm[1] < -0.0001){
    for(int i=0; i<9; i++){
      wR[i] = wL[i];
    }
  }
  else{
    MHDImposedDataKelvinHelmotz(x,t,wR);
  }

  
  MHDNumFluxRusanov(wL,wR,vnorm,flux);
  
}
#pragma end_opencl

#pragma start_opencl
void MHDInitDataKelvinHelmotz(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  MHDImposedDataKelvinHelmotz(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void MHDImposedDataKelvinHelmotz(const schnaps_real *x,const  schnaps_real t, schnaps_real *w) {
  //real gam = 1.6666666666;
  schnaps_real yL[9];

  yL[0] = 1.0; //reprendre 1.29 a la place de 1.806 et prendre Bx = 1.29/3.33 = 0.39
  yL[1] = (1.29/2)*tanh((x[1])/0.05);
  yL[2] = 1.0;
  yL[3] = 0.01*sin(2*3.14159265359*(x[0]))*exp(-((x[1])/0.2)*((x[1])/0.2));
  yL[4] = 0;
  yL[7] = 0.39; // pour le Bx on prend le 1.29 et on divise pas le mac
  yL[5] = 0;
  yL[6] = 0;
  yL[8] = 0;
  
  conservatives(yL, w);
  
}
#pragma end_opencl


////////////////////////////////////////////////////////////////////////////////
//                            Double tearing                                  //
////////////////////////////////////////////////////////////////////////////////

#pragma start_opencl
void MHDBoundaryFluxDoubleTearing(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
                                 schnaps_real *flux) {
  schnaps_real wR[9];

  //if(vnorm[1] > 0.0001 || vnorm[1] < -0.0001){
  //  for(int i=0; i<9; i++){
  //    MHDImposedDataDoubleTearing(x,t,wR);
  //  }
  //}
  //else{
  //  MHDImposedDataDoubleTearing(x,t,wR);
  //}

  MHDImposedDataDoubleTearing(x,t,wR);
  MHDNumFluxRusanov(wL,wR,vnorm,flux);
}
#pragma end_opencl

#pragma start_opencl
void MHDInitDataDoubleTearing(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  MHDImposedDataDoubleTearing(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void MHDImposedDataDoubleTearing(const schnaps_real *x,const  schnaps_real t, schnaps_real *w) {
  //real gam = 1.6666666666;
  schnaps_real yL[9];

  yL[5] = 1 + tanh( (x[0]-0.5)/0.2 ) - tanh( (x[0]+0.5)/0.2 );     // By
  yL[2] = -yL[5]*yL[5]/2 + 0.6;                                  // P
  yL[0] = 2*yL[2]/0.2;                                          // Rho

  if((x[0]) < -1.5 && (x[0]) > 1.5)
    yL[7] = 0;
  else
    yL[7] = 0.003*sin(3.14159265359*(x[1])/2)*exp(-( (x[0])*(x[0])-0.5 )/0.1 );   // Bx essayer de forcer les bord a 0. si on est < -1.5 x>1.5 je met 0

  yL[1] = 0;                                                 // Ux
  yL[3] = 0;                                                 // Uy
  yL[4] = 0;                                                 // Uz
  yL[6] = 0;                                                 // Bz
  yL[8] = 0;                                                 // Psi

  
  conservatives(yL, w);
  
}
#pragma end_opencl

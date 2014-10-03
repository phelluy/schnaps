#include "interpolation.h"
#include "global.h"
/* gauss lobatto points */
#include<stdio.h>



const double gauss_lob_point[] = {
  0.5,
  0.,
  1.,
  0.,
  0.5,
  1.,
  0.,
  0.276393202250021030359082633127,
  0.723606797749978969640917366873,
  1.,
  0.,
  0.172673164646011428100853771877, 
  0.5, 
  0.827326835353988571899146228123,
  1.
};

const double gauss_lob_weight[] = {
  1.,
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

const int gauss_lob_offset[] = {
  0, 1, 3, 6, 10
};

const double gauss_lob_dpsi[] = {
  0.00000000000000000000000000000000000000000000000000000000000,
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

const int gauss_lob_dpsi_offset[] = {0, 1, 5, 14, 30};


void lagrange_polynomial(double* p,const double* subdiv,int deg,int ii,double x) {
  *p = 1;                                              
  for(int j=0;j<deg+1;j++){                 
    if (j != ii){                                 
      *p *= (x - subdiv[j]) /                  
	(subdiv[ii] - subdiv[j]);             
    }                                                 
  }
  
}




void dlagrange_polynomial(double* dp,const double* subdiv,int deg,int i,double x){
  double xj;
  *dp = 0;
  for(int k=0;k<((deg)+1);k++){
    if (k != (i)){
      double xk = (subdiv)[k];
      double dploc = ((double)1) / ((subdiv)[i] -  xk);
      for(int j=0;j<(deg)+1;j++){
	if (j != (i) && j != k){
	  xj = (subdiv)[j];
	  dploc *= ((x) - xj) / ((subdiv)[i] - xj);
	}
      }
      *dp += dploc;
    }
  }
}
                                                                





int NPG(int* param){
  return (param[0]+1)*(param[1]+1)*(param[2]+1) *
         (param[3])*(param[4])*(param[5]);
}


/// Number of interpolation points for each face
int NPGF(int* param, int ifa){
  // For each face, give the dimension index i
  int permut[6][4] = {
    {0, 2, 1, 0},
    {1, 2, 0, 1},
    {2, 0, 1, 1},
    {2, 1, 0, 0},
    {0, 1, 2, 1},
    {1, 0, 2, 0}
  };
  int i0 = permut[ifa][0];
  int i1 = permut[ifa][1];
  return (param[i0] + 1) * (param[i1] + 1) * param[i0+3] * param[i1+3]; 
}

// return the reference coordinates xpg[3] and weight wpg of the GLOP ipg
void ref_pg_vol(int* param,int ipg,double* xpg,double* wpg){
  int deg[3],offset[3],nraf[3];
  

  // approximation degree in each direction
  deg[0]=param[0];
  deg[1]=param[1];
  deg[2]=param[2];
  // number of subcells in each direction
  nraf[0]=param[3];
  nraf[1]=param[4];
  nraf[2]=param[5];
  
  
  int ix = ipg % (deg[0] + 1);  
  ipg/=(deg[0] + 1);

  int iy = ipg % (deg[1] + 1);
  ipg/=(deg[1] + 1);

  int iz = ipg % (deg[2] + 1);
  ipg/=(deg[2] + 1);

  int ncx= ipg % nraf[0];
  double hx=1/(double) nraf[0]; 
  ipg/=nraf[0];

  int ncy= ipg % nraf[1];
  double hy=1/(double) nraf[1]; 
  ipg/=nraf[1];

  int ncz= ipg;
  double hz=1/(double) nraf[2]; 


  offset[0]=gauss_lob_offset[deg[0]]+ix;
  offset[1]=gauss_lob_offset[deg[1]]+iy;
  offset[2]=gauss_lob_offset[deg[2]]+iz;


  xpg[0]=hx*(ncx+gauss_lob_point[offset[0]]);
  xpg[1]=hy*(ncy+gauss_lob_point[offset[1]]);
  xpg[2]=hz*(ncz+gauss_lob_point[offset[2]]);

  *wpg=hx*hy*hz*gauss_lob_weight[offset[0]]*
    gauss_lob_weight[offset[1]]*
    gauss_lob_weight[offset[2]];

  //  printf("%f %f %f\n", hx, hy, hz);//gauss_lob_weight[offset[0]],gauss_lob_weight[offset[1]],gauss_lob_weight[offset[2]]);
};
// same function for the face 
void ref_pg_face(int* param,int ifa,int ipg,double* xpg,double* wpg){
// For each face, give the dimension index i
  const int axis_permut[6][4] = {
    {0, 2, 1, 0},
    {1, 2, 0, 1},
    {2, 0, 1, 1},
    {2, 1, 0, 0},
    {0, 1, 2, 1},
    {1, 0, 2, 0}
  };
  int deg[3],offset[2],nraf[3];
  double h[3];
  int ipgxyz[3],ncpgxyz[3];
  int ipgf=ipg;

  // approximation degree in each direction
  deg[0]=param[axis_permut[ifa][0]];
  deg[1]=param[axis_permut[ifa][1]];
  deg[2]=param[axis_permut[ifa][2]];
  // number of subcells in each direction
  nraf[0]=param[3+axis_permut[ifa][0]];
  nraf[1]=param[3+axis_permut[ifa][1]];
  nraf[2]=param[3+axis_permut[ifa][2]];


  // Compute permuted indices
  int ix = ipg % (deg[0] + 1);  
  ipg/=(deg[0] + 1);

  int iy = ipg % (deg[1] + 1);
  ipg/=(deg[1] + 1);

  // Equals 0 or d depending on the face
  int iz = axis_permut[ifa][3]*deg[2];

  // Compute permuted indices of the subface
  int ncx= ipg % nraf[0];
  h[0]=1/(double) nraf[0]; 
  ipg/=nraf[0];

  int ncy= ipg;
  h[1]=1/(double) nraf[1]; 
 
  // Equals 0 or nraf-1 depending on the face
  int ncz = axis_permut[ifa][3]*(nraf[2]-1);
  h[2]=1/(double) nraf[2]; 

  // Compute non permuted indices for points and subfaces
  ipgxyz[axis_permut[ifa][0]]=ix;
  ipgxyz[axis_permut[ifa][1]]=iy;
  ipgxyz[axis_permut[ifa][2]]=iz;
  
  ncpgxyz[axis_permut[ifa][0]]=ncx;
  ncpgxyz[axis_permut[ifa][1]]=ncy;
  ncpgxyz[axis_permut[ifa][2]]=ncz;

  // Compute the global index of the 
  // Gauss-Lobatto point in the volume
  param[6]=ipgxyz[0]+(param[0]+1)*
    (ipgxyz[1]+(param[1]+1)*
     (ipgxyz[2]+(param[2]+1)*
      (ncpgxyz[0]+param[3]*
       (ncpgxyz[1]+param[4]*
	ncpgxyz[2]))));

  /* if (ifa==2 && ipgf == 3){ */
  /*   printf("ix=%d %d %d , ncx=%d %d %d ipgv=%d\n",ipgxyz[0],ipgxyz[1],ipgxyz[2],ncpgxyz[0],ncpgxyz[1],ncpgxyz[2], param[6]); */
  /* } */

  // Compute the reference coordinates of the 
  // Gauss-Lobatto point in the volume
  offset[0]=gauss_lob_offset[deg[0]]+ix;
  offset[1]=gauss_lob_offset[deg[1]]+iy;
  //printf("offset=%d\n",offset);

  xpg[axis_permut[ifa][0]] = h[0]*(ncx+gauss_lob_point[offset[0]]);
  xpg[axis_permut[ifa][1]] = h[1]*(ncy+gauss_lob_point[offset[1]]);
  xpg[axis_permut[ifa][2]] = axis_permut[ifa][3];

  *wpg=h[0]*h[1]*gauss_lob_weight[offset[0]]*
    gauss_lob_weight[offset[1]];

};
// return the value psi  and the gradient dpsi[3] of the basis 
// function ib at point xref[3]. Warning: the value of the gradient is
// not reliable if xref is on the boundary 
//of a subcell (because the gradient is discontinuous)
void psi_ref(int* param, int ib, double* xref, double* psi, double* dpsi){
  
  double dpsibx;
  double dpsiby;
  double dpsibz;

  int deg[3],offset[3],nraf[3];
  
  // approximation degree in each direction
  deg[0]=param[0];
  deg[1]=param[1];
  deg[2]=param[2];
  // number of subcells in each direction
  nraf[0]=param[3];
  nraf[1]=param[4];
  nraf[2]=param[5];
  // Starting Gauss-Lobatto point in each direction
  offset[0]=gauss_lob_offset[deg[0]];
  offset[1]=gauss_lob_offset[deg[1]];
  offset[2]=gauss_lob_offset[deg[2]];

  // basis functions indices
  int ibx = ib % (deg[0] + 1);  
  ib/=(deg[0] + 1);

  int iby = ib % (deg[1] + 1);
  ib/=(deg[1] + 1);

  int ibz = ib % (deg[2] + 1);
  ib/=(deg[2] + 1);

  int ncbx= ib % nraf[0];
  double hx=1/(double) nraf[0]; 
  ib/=nraf[0];

  int ncby= ib % nraf[1];
  double hy=1/(double) nraf[1]; 
  ib/=nraf[1];

  int ncbz= ib;
  double hz=1/(double) nraf[2]; 

  double psibx = 0;
  double psiby = 0;
  double psibz = 0;

  lagrange_polynomial(&psibx, gauss_lob_point + offset[0], deg[0], ibx, xref[0]/hx-ncbx);
  lagrange_polynomial(&psiby, gauss_lob_point + offset[1], deg[1], iby, xref[1]/hy-ncby);
  lagrange_polynomial(&psibz, gauss_lob_point + offset[2], deg[2], ibz, xref[2]/hz-ncbz);

  /* psibx *= (xref[0] <= (ncbx + 1) * hx)&&(xref[0] > ncbx * hx); */
  /* psiby *= (xref[1] <= (ncby + 1) * hy)&&(xref[1] > ncby * hy); */
  /* psibz *= (xref[2] <= (ncbz + 1) * hz)&&(xref[2] > ncbz * hz); */

  *psi = psibx * psiby * psibz;

  if (dpsi != NULL){

    dlagrange_polynomial(&dpsibx, gauss_lob_point + offset[0], deg[0], ibx, xref[0]);
    dlagrange_polynomial(&dpsiby, gauss_lob_point + offset[1], deg[1], iby, xref[1]);
    dlagrange_polynomial(&dpsibz, gauss_lob_point + offset[2], deg[2], ibz, xref[2]);

    dpsi[0] = dpsibx *  psiby *  psibz;
    dpsi[1] =  psibx * dpsiby *  psibz;
    dpsi[2] =  psibx *  psiby * dpsibz;
  }

};
// return the gradient dpsi[3] of the basis 
// function ib at GLOP ipg.
void grad_psi_pg(int* param,int ib,int ipg,double* dpsi){
  int deg[3],offset[3],nraf[3];
  
  // approximation degree in each direction
  deg[0]=param[0];
  deg[1]=param[1];
  deg[2]=param[2];
  // number of subcells in each direction
  nraf[0]=param[3];
  nraf[1]=param[4];
  nraf[2]=param[5];
  
  // indices of Gauss-Lobatto points in each subcell
  int ipgx = ipg % (deg[0] + 1);  
  ipg/=(deg[0] + 1);

  int ipgy = ipg % (deg[1] + 1);
  ipg/=(deg[1] + 1);

  int ipgz = ipg % (deg[2] + 1);
  ipg/=(deg[2] + 1);

  // indices of each subcell and space step in each direction
  int ncpgx= ipg % nraf[0];
  double hx=1/(double) nraf[0]; 
  ipg/=nraf[0];

  int ncpgy= ipg % nraf[1];
  double hy=1/(double) nraf[1]; 
  ipg/=nraf[1];

  int ncpgz= ipg;
  double hz=1/(double) nraf[2]; 

  // basis functions indices
  int ibx = ib % (deg[0] + 1);  
  ib/=(deg[0] + 1);

  int iby = ib % (deg[1] + 1);
  ib/=(deg[1] + 1);

  int ibz = ib % (deg[2] + 1);
  ib/=(deg[2] + 1);

  int ncbx= ib % nraf[0];
  ib/=nraf[0];

  int ncby= ib % nraf[1];
  ib/=nraf[1];

  int ncbz= ib;

  // Number of Gauss-Lobatto points in each direction
  offset[0]=gauss_lob_dpsi_offset[deg[0]];
  offset[1]=gauss_lob_dpsi_offset[deg[1]];
  offset[2]=gauss_lob_dpsi_offset[deg[2]];

  // Computation of the value of the interpollation polynomial gradient
  double psibx,psiby,psibz,dpsibx,dpsiby,dpsibz;

  psibx=(ipgx==ibx) * (ncpgx==ncbx);
  //dpsibx=dlagrange_polynomial(gauss_leg_point+offset,deg,ix,xref[0]);
  dpsibx = (ncpgx==ncbx) * gauss_lob_dpsi[offset[0]+ibx*(deg[0]+1)+ipgx] / hx;

  psiby=(ipgy==iby) * (ncpgy==ncby);
  dpsiby = (ncpgy==ncby) * gauss_lob_dpsi[offset[1]+iby*(deg[1]+1)+ipgy] / hy;

  psibz=(ipgz==ibz) * (ncpgz==ncbz);
  dpsibz = (ncpgz==ncbz) * gauss_lob_dpsi[offset[2]+ibz*(deg[2]+1)+ipgz] / hz;


  dpsi[0]=dpsibx*psiby*psibz;
  dpsi[1]=psibx*dpsiby*psibz;
  dpsi[2]=psibx*psiby*dpsibz;

  //printf("  ipgx=%d ibx=%d ncpgx=%d ncbx=%d psibx=%f dpsibx=%f dpsi=%f d=%d \n", ipgx, ibx, ncpgx, ncbx, psibx, dpsibx, dpsi[0], deg[0]);
};




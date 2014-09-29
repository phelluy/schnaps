#include "h20.h"
#include "global.h"
#include<stdio.h>
#include <math.h>
const int h20_refnormal[6][3]={{0,-1,0},{1,0,0},
                                           {0,1,0},{-1,0,0},
                                           {0,0,1},{0,0,-1}};


void Ref2Phy(double* physnode,
             double* xref,
             double* dphiref,
             int ifa,
             double* xphy,
             double* dtau,
             double* codtau,
             double* dphi,
             double* vnds){



  //HEXAN_REF2PHY(physnode, xref, 1, xphy, dtau);

  // compute the mapping and its jacobian
  double x,y,z;
  x=xref[0];
  y=xref[1];
  z=xref[2];

  //double phi[20];
  // gradient of the shape functions and value (4th component)
  // of the shape functions
  double gradphi[20][4];



  // this file fills the values of phi and gradphi
#include "h20phi.h"

  for(int ii=0;ii<3;ii++){
    xphy[ii]=0;
    for(int jj=0;jj<3;jj++){
      dtau[3*ii+jj]=0;
    }
    for(int i=0;i<20;i++){
      xphy[ii]+=physnode[3*i+ii]*gradphi[i][3];
      for(int jj=0;jj<3;jj++){
        dtau[3*ii+jj]+=physnode[3*i+ii]*gradphi[i][jj];;
      }
    }
  }
      

  if (codtau != NULL) {
    codtau[0] =  dtau[4] * dtau[8] - dtau[5] * dtau[7];
    codtau[1] = -dtau[3] * dtau[8] + dtau[5] * dtau[6];
    codtau[2] =  dtau[3] * dtau[7] - dtau[4] * dtau[6];
    codtau[3] = -dtau[1] * dtau[8] + dtau[2] * dtau[7];
    codtau[4] =  dtau[0] * dtau[8] - dtau[2] * dtau[6];
    codtau[5] = -dtau[0] * dtau[7] + dtau[1] * dtau[6];
    codtau[6] =  dtau[1] * dtau[5] - dtau[2] * dtau[4];
    codtau[7] = -dtau[0] * dtau[5] + dtau[2] * dtau[3];
    codtau[8] =  dtau[0] * dtau[4] - dtau[1] * dtau[3];
  }

  if (dphi != NULL){
    int ii;
    for(ii=0;ii<3;ii++){
      int jj;
      dphi[ii]=0;
      for(jj=0;jj<3;jj++){
        dphi[ii]+=codtau[3*ii+jj]*dphiref[jj];
      }
    }
  }

  if (vnds !=NULL) {

    int ii;
    for(ii=0;ii<3;ii++){
      int jj;
      vnds[ii]=0;
      for(jj=0;jj<3;jj++){
        vnds[ii]+=codtau[3*ii+jj]*h20_refnormal[ifa][jj];
      }
    }

  }

}


void Phy2Ref(double* physnode,double* xphy,double* xref){
#define ITERNEWTON 3

  double dtau[9], codtau[9];                                              
  double dxref[3], dxphy[3];                                              
  double det;                                                             
  int ifa =- 1;                                                         
  xref[0] = 0.5;                                                          
  xref[1] = 0.5;                                                          
  xref[2] = 0.5;                                                          
  for(int iter = 0;iter<ITERNEWTON;iter ++ ){                               
    Ref2Phy(physnode, xref, 0,ifa,                                      
            dxphy, dtau, codtau, 0,0);                                  
    dxphy[0] -= (xphy)[0];                                              
    dxphy[1] -= (xphy)[1];                                              
    dxphy[2] -= (xphy)[2];                                              
    det = dtau[0] * dtau[4] * dtau[8] - dtau[0] * dtau[5] * dtau[7]     
      - dtau[3] * dtau[1] * dtau[8]                                     
      + dtau[3] * dtau[2] * dtau[7] + dtau[6] * dtau[1] * dtau[5]       
      - dtau[6] * dtau[2] * dtau[4];     
    for(int ii = 0;ii<3;ii ++ ){                             
      dxref[ii] = 0;                                               
      for(int jj = 0;jj<3;jj ++ ){                                          
        dxref[ii] += codtau[3 * jj + ii] * dxphy[jj];         
      }      
      xref[ii] -= dxref[ii] / det;                            
    }
    printf("iter=%d erreur=%f\n",iter,sqrt(dxref[0]*dxref[0]+
				      dxref[1]*dxref[1]+
				      dxref[2]*dxref[2]));
  }

}

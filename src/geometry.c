#include "geometry.h"
#include "global.h"
#include<stdio.h>
#include <math.h>
#include <assert.h>
const int h20_refnormal[6][3]={{0,-1,0},{1,0,0},
                                           {0,1,0},{-1,0,0},
                                           {0,0,1},{0,0,-1}};
double Dist(double x1[3],double x2[3]){
  return sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+
              (x1[1]-x2[1])*(x1[1]-x2[1])+
              (x1[2]-x2[2])*(x1[2]-x2[2]));
}

void PrintPoint(double x[3]){
  printf("%f %f %f\n",x[0],x[1],x[2]);
}

void GeomRef2Phy(Geom* g){
  Ref2Phy(g->physnode,g->xref,g->dphiref,
          g->ifa,g->xphy,g->dtau,g->codtau,
          g->dphi,g->vnds);
  g->det =g->codtau[0][0]*g->dtau[0][0]+g->codtau[0][1]*g->dtau[0][1]+
	g->codtau[0][2]*g->dtau[0][2];
};

void Ref2Phy(double physnode[20][3],
             double xref[3],
             double dphiref[3],
             int ifa,
             double xphy[3],
             double dtau[3][3],
             double codtau[3][3],
             double dphi[3],
             double vnds[3]){

  // compute the mapping and its jacobian
  double x,y,z;
  x=xref[0];
  y=xref[1];
  z=xref[2];

  // gradient of the shape functions and value (4th component)
  // of the shape functions
  double gradphi[20][4];
#include "h20phi.h"  // this file fills the values of phi and gradphi

  if (xphy != NULL){
    for(int ii=0;ii<3;ii++){
      xphy[ii]=0;
      for(int i=0;i<20;i++){
	xphy[ii]+=physnode[i][ii]*gradphi[i][3];
      }
    }
  }

  if (dtau != NULL){
    for(int ii=0;ii<3;ii++){
      for(int jj=0;jj<3;jj++){
	dtau[ii][jj]=0;
      }
      for(int i=0;i<20;i++){
	for(int jj=0;jj<3;jj++){
	  dtau[ii][jj]+=physnode[i][ii]*gradphi[i][jj];;
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

  if (dphi != NULL){
    assert(codtau != NULL);
    int ii;
    for(ii=0;ii<3;ii++){
      int jj;
      dphi[ii]=0;
      for(jj=0;jj<3;jj++){
        dphi[ii]+=codtau[ii][jj]*dphiref[jj];
      }
    }
  }

  if (vnds !=NULL) {
    assert(codtau != NULL);
    assert(ifa >=0);
    int ii;
    for(ii=0;ii<3;ii++){
      int jj;
      vnds[ii]=0;
      for(jj=0;jj<3;jj++){
        vnds[ii]+=codtau[ii][jj]*h20_refnormal[ifa][jj];
      }
    }

  }

}

void GeomPhy2Ref(Geom* g){
  Phy2Ref(g->physnode,g->xphy,g->xref);
}
void Phy2Ref(double physnode[20][3],double xphy[3],double xref[3]){
#define ITERNEWTON 10

  double dtau[3][3], codtau[3][3];                                              
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
    det = dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
      dtau[0][2]*codtau[0][2];
    assert(det>0);
    for(int ii = 0;ii<3;ii ++ ){                             
      dxref[ii] = 0;                                               
      for(int jj = 0;jj<3;jj ++ ){                                          
        dxref[ii] += codtau[jj][ii] * dxphy[jj];         
      }      
      xref[ii] -= dxref[ii] / det;                            
    }
  }
  double eps=1e-2;  // may be to constraining...
  assert(xref[0]<1+eps && xref[0]>-eps);
  assert(xref[1]<1+eps && xref[1]>-eps);
  assert(xref[2]<1+eps && xref[2]>-eps);

}

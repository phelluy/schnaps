#include "h20.h"
#include "global.h"
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



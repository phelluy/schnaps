#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"



void VlasovP_Lagrangian_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  
  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    real vn = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j))*vnorm[0];
    
    real vnp = vn>0 ? vn : 0;
    real vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
  // do not change the potential !
  // and the electric field
  flux[_INDEX_PHI]=0;  // flux for phi
  flux[_INDEX_EX]=0; // flux for E
  flux[_INDEX_RHO]=0; // flux for rho
  flux[_INDEX_VELOCITY]=0; // flux for u
  flux[_INDEX_PRESSURE]=0; // flux for p
  flux[_INDEX_TEMP]=0; // flux for e ou T

};


//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void VlasovP_Lagrangian_Source(const real* x, const real t, const real* w, 
			       real* source) {

  real E=w[_INDEX_EX]; // electric field
  real Md[_INDEX_MAX_KIN+1];
  real db[_INDEX_MAX_KIN+1];
  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    Md[iv]=0;
    db[iv]=0;
  }
  
  
  for(int iv=0;iv<_MV+2;iv++){
    source[iv]=0;
  }
  // no source on the potential for the moment
  source[_INDEX_PHI]=0;
  source[_INDEX_EX]=0;
  source[_INDEX_RHO]=0; //rho init
  source[_INDEX_VELOCITY]=0; // u init
  source[_INDEX_PRESSURE]=0; // p init
  source[_INDEX_TEMP]=0; 
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<_DEG_V+1;kloc++){
      real omega=wglop(_DEG_V,kloc);
      int kpg=kloc+iel*_DEG_V;
      Md[kpg]+=omega*_DV;
      for(int iloc=0;iloc<_DEG_V+1;iloc++){
	int ipg=iloc+iel*_DEG_V;
	source[ipg]+=E*omega*w[kpg]*dlag(_DEG_V,iloc,kloc);
	if (iloc==kloc) db[ipg]+=E*omega*dlag(_DEG_V,iloc,kloc);
      }
    }
  }

  // upwinding
  if (E>0){
    source[_INDEX_MAX_KIN]-=E*w[_INDEX_MAX_KIN];
    db[_INDEX_MAX_KIN]-=E;
  }
  else {
    source[0]-=-E*w[0];
    db[0]-=-E;
  }

  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    source[iv]/=Md[iv];
    //printf("%f ",source[iv]);
  }
  //printf("\n");
  //assert(1==2);
  

};



#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"



void VlasovP_Lagrangian_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  KineticData * kd=&schnaps_kinetic_data;
  for(int i = 0;i < kd->index_max_kin+1;i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vn = (-kd->vmax+nel*kd->dv +
		 kd->dv* glop(kd->deg_v,j))*vnorm[0];
    
    schnaps_real vnp = vn>0 ? vn : 0;
    schnaps_real vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
  // do not change the potential !
  // and the electric field
  flux[kd->index_phi]=0;  // flux for phi
  flux[kd->index_ex]=0; // flux for E
  flux[kd->index_ey]=0;
  flux[kd->index_ez]=0;
  flux[kd->index_rho]=0; // flux for rho
  flux[kd->index_u]=0; // flux for u
  flux[kd->index_P]=0; // flux for p
  flux[kd->index_T]=0; // flux for e ou T

};


//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void VlasovP_Lagrangian_Source(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, 
			       schnaps_real* source) {
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real E=w[kd->index_ex]; // electric field
  schnaps_real Md[kd->index_max_kin +1];
  schnaps_real db[kd->index_max_kin +1];
  for(int iv=0;iv< kd->index_max_kin + 1;iv++){
    Md[iv]=0;
    db[iv]=0;
  }
  
  for(int iv=0;iv< kd->index_max_kin + 1;iv++){
    source[iv]=0;
  }
  // no source on the potential for the moment
  source[kd->index_phi]=0;
  source[kd->index_ex]=0;
  source[kd->index_ey]=0;
  source[kd->index_ez]=0;
  source[kd->index_rho]=0; //rho init
  source[kd->index_u]=0; // u init
  source[kd->index_P]=0; // p init
  source[kd->index_T]=0; 
  // loop on the finite emlements
  for(int iel=0;iel<kd->nb_elem_v;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<kd->deg_v+1;kloc++){
      schnaps_real omega=wglop(kd->deg_v,kloc);
      int kpg = kloc+iel*kd->deg_v;
      Md[kpg] += omega*kd->dv;
      for(int iloc=0;iloc<kd->deg_v+1;iloc++){
	int ipg=iloc+iel*kd->deg_v;
	source[ipg] += E*omega*w[kpg]*dlag(kd->deg_v,iloc,kloc);
	if (iloc==kloc) db[ipg]+=E*omega*dlag(kd->deg_v,iloc,kloc);
      }
    }
  }

  // upwinding
  if (E>0){
    source[kd->index_max_kin]-=E*w[kd->index_max_kin];
    db[kd->index_max_kin]-=E;
  }
  else {
    source[0]-=-E*w[0];
    db[0]-=-E;
  }

  for(int iv=0;iv<kd->index_max_kin+1;iv++){
    source[iv]/=Md[iv];
  }

  
};

void BGK_Source(const schnaps_real* w, schnaps_real* source) {

  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real Maxw=0;
  
  for(int i=0;i<kd->index_max;i++){
    source[i]=0;
  }
  
  for(int ielv=0;ielv<kd->nb_elem_v;ielv++){
    // loop on the local glops
    for(int iloc=0;iloc<kd->deg_v+1;iloc++){
      schnaps_real omega=wglop(kd->deg_v,iloc);
      schnaps_real vi=-kd->vmax+ielv*kd->dv+kd->dv*glop(kd->deg_v,iloc);
      int ipgv=iloc+ielv*kd->deg_v;
      Maxw=Computation_Maxwellian(w[kd->index_rho],w[kd->index_u],w[kd->index_T],vi);
      source[ipgv]=kd->knud*(Maxw-w[ipgv]);
    }
  }
  
  
};




//#include "collision.h"
#include "gyro.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"


void Gyro_Lagrangian_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  
  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    real vn = (nel*_DV +
		 _DV* glop(_DEG_V,j))*vnorm[0];
    
    real vnp = vn>0 ? vn : 0;
    real vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
  
};

//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
real GyroL2VelError(real* x,real t,real *w){


  real wex[_MV];
  real err2=0;
  GyroImposedData(x, t,wex);
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      real omega=wglop(_DEG_V,iloc);
      real vi=iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      err2+=omega*_DV*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
    }
  }
  return err2;
};



void Gyro_Lagrangian_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
				       real* flux){
  real wR[_MV];
  GyroImposedData(x,t,wR);
  Gyro_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};


void GyroInitData(real x[3],real w[]){

  real t=0;
  GyroImposedData(x,t,w);

};



void GyroImposedData(real x[3],real t,real w[]){

  for(int i=0;i<_MV;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    real vi = (nel*_DV +
		 _DV* glop(_DEG_V,j));
    w[i]=cos(x[0]-vi*t);
  }

};


real Gyro_ImposedKinetic_Data(real x[3],real t,real v){
  real f;
  f=cos(x[0]-v*t);
  return f;
};

real GyroL2_Kinetic_error(field* f){

  real error=0;
  real error_space=0;
  real moy=0; // mean value
  real moy_space=0;


  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      real xpgref[3],xphy[3],wpg;
      real dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      real w[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      error+=GyroL2VelError(xphy,f->tnow,w)*wpg*det;
    }
  }
  //moy=moy+weight*moy_space;

  return sqrt(error);
  //return sqrt(error)/sqrt(moy);
}

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void GyroSource(real* force,real* w, real* source){

  real E=force[0]; // electric field
  real Md[_MV];
  for(int iv=0;iv<_MV;iv++){
    Md[iv]=0;
    source[iv]=0;
  }
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<_DEG_V+1;kloc++){
      real omega=wglop(_DEG_V,kloc);
      int kpg=kloc+iel*_DEG_V;
      Md[kpg]+=omega*_DV;
      for(int iloc=0;iloc<_DEG_V+1;iloc++){
	int ipg=iloc+iel*_DEG_V;
	source[ipg]-=E*omega*w[kpg]*dlag(_DEG_V,iloc,kloc);
      }
    }
  }

  // upwinding
  if (E>0){
    source[_MV-1]+=E*w[_MV-1];
  }
  else {
    source[0]+=-E*w[0];
  }

  for(int iv=0;iv<_MV;iv++){
    source[iv]/=Md[iv];
  }
  

};



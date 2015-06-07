//#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "gyro.h"


/* void Gyro_Lagrangian_NumFlux(real wL[],real wR[],real* vnorm,real* flux){ */
  
/*   for(int i=0;i<_MV;i++){ */
/*     int j=i%_DEG_V; // local connectivity put in function */
/*     int nel=i/_DEG_V; // element num (TODO : function) */

/*     real vn = (nel*_DV + */
/* 		 _DV* glop(_DEG_V,j))*vnorm[0]; */
    
/*     real vnp = vn>0 ? vn : 0; */
/*     real vnm = vn-vnp; */

/*     flux[i] = vnp * wL[i] + vnm * wR[i]; */
/*   } */
  
/* } */

//flux num gyro
void Gyro_Lagrangian_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real eps =0; //if not equal 0 => decentered flux
  /* real E_x =0; //firstly consider the electric field is const */
  /* real E_y =1; */
  real E_x =wL[_INDEX_EX];
  real E_y =wL[_INDEX_EY];
  //real dv=_DV;
  //printf("test dv = %f \n",dv); 
  /* assert(1==2); */
  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)  
    real wm = (wL[i]+wR[i])/2;
    real flux1 = E_y*wm;
    real flux2 = -E_x*wm;
    real v =-_VMAX+ nel*_DV +
      _DV* glop(_DEG_V,j); // gauss_lob_point[j]
    real flux3 =v*wm;
  
    flux[i] = vnorm[0]*flux1+vnorm[1]*flux2+vnorm[2]*flux3-eps*(wR[i]-wL[i])/2;
  }
  flux[_INDEX_PHI] =0;
  flux[_INDEX_EX]=0;
  flux[_INDEX_EY]=0;
  flux[_INDEX_EZ]=0;
  
}

//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
real GyroL2VelError(real* x,real t,real *w){


  real wex[_INDEX_MAX];
  real err2=0;
  GyroImposedData(x, t,wex);
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      real omega=wglop(_DEG_V,iloc);
      real vi=-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      err2+=omega*_DV*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
    }
  }
  return err2;
}



void Gyro_Lagrangian_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
				       real* flux){
  real wR[_INDEX_MAX];
  GyroImposedData(x,t,wR);
  Gyro_Lagrangian_NumFlux(wL,wR,vnorm,flux);
}


void GyroInitData(real x[3],real w[]){

  real t=0;
  GyroImposedData(x,t,w);

}

void GyroImposedData(const real x[3], const real t,real w[]){

  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    real vi = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j));
    //w[i]=cos(x[0]-vi*t);
    real pi=4*atan(1.);
    real xi = x[0]-t;
    //pour transport 1d
    // w[i] = cos(2*pi*xi);//*exp(-0.5*pow(xi-0.5,2));
    //w[i] = cos(2*pi*xi)*exp(-1/(1-pow(2*xi-1,2)));
    //w[i] = exp(-4*pow(xi-0.5,2));
    real yi = x[1]+t;
    //pour transport 2D
    //w[i] = exp(-4*pow(yi-0.5,2))*exp(-4*pow(xi-0.5,2));
    //w[i] = cos(2*pi*yi);//*exp(-1/(1-pow(2*xi-1,2)));
    //w[i] = cos(2*pi*xi)*cos(2*pi*yi);//*exp(-0.5*pow(xi-0.5,2));
        //pour transport 3D
    //real zi=x[2]-vi*t;
    real zi=x[2]-vi*t;
     //w[i] = exp(-4*pow(zi-0.5,2));//exp(-4*pow(yi-0.5,2))*exp(-4*pow(xi-0.5,2))*exp(-4*pow(zi-0.5,2));

    w[i] = exp(-4*pow(zi-1.,2));// *Gyro_ImposedKinetic_Data(x,t,vi);
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI]=0;
  w[_INDEX_EX]=0;//1;
  w[_INDEX_EY]=0;//1;
  w[_INDEX_EZ]=0;

}


real Gyro_ImposedKinetic_Data(const real x[3], const real t, real v) {
  real f;
  f=exp(-pow((v-1.*t),2)/16.); //velocity transport, Ez=1
  //f=exp(-4*pow(xi-0.5,2))*exp(-pow((v-2.*t),2));
  return f;
}

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
void GyroSource(const real *x, const real t, const real *w, real *source) {
  real Ez=w[_INDEX_EZ]; // electric field
  real Md[_MV];
  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    Md[iv]=0;
    source[iv]=0;
  }
  // loop on the finite elements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<_DEG_V+1;kloc++){
      real omega=wglop(_DEG_V,kloc);
      int kpg=kloc+iel*_DEG_V;
      Md[kpg]+=omega*_DV;
      for(int iloc=0;iloc<_DEG_V+1;iloc++){
	int ipg=iloc+iel*_DEG_V;
	source[ipg]+=Ez*omega*w[kpg]*dlag(_DEG_V,iloc,kloc);
      }
    }
  }

  // upwinding !beta =1 (or 1/2?)
  if (Ez>0){
    source[_MV-1]+=-Ez*w[_MV-1]; 
    //source[0]+=Ez*w[0];
  }
  else {
    source[0]+=Ez*w[0];
    //source[_INDEX_MAX_KIN]+=-Ez*w[_INDEX_MAX_KIN];
  }

  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    source[iv]/=Md[iv];
  }
  source[_INDEX_PHI]=0;
  source[_INDEX_EX]=0;
  source[_INDEX_EY]=0;
  source[_INDEX_EZ]=0;

}



void Velocity_distribution_plot(real *w){
  
  FILE * ver;
  ver = fopen( "fvel.dat", "w" );

  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      real omega=wglop(_DEG_V,iloc);
      real dv=_DV;
      real vi=-_VMAX+iel*dv+dv*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      //err2+=omega*_DV*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
      fprintf(ver,"%f %f\n",vi,w[ipg]);
      //printf("dv is  %f \n", dv);
      //printf("vi est %d %d %d %f\n", iel,iloc,ipg, vi);
    }
  }
  fclose(ver);
}

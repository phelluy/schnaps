#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"


//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
//! \param[in] x : point of the mesh
//! \param[in] t : time
//! \param[in] t : type of L2norm. if type_norm=0 this is the numerical solution if type_norm=1 this is the error 
double L2VelError(field * f,double* x,double *w){

  double wex[_INDEX_MAX];
  double err2=0;
  double t=f->tnow;
  f->model.ImposedData(x,t,wex);
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      double omega=wglop(_DEG_V,iloc);
      double vi=-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      err2+=omega*_DV*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
    }
  }
  
  return err2;
};

double L2_Kinetic_error(field* f){

  double error=0;

  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],xphy[3],wpg;
      double dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      double w[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      error+=L2VelError(f,xphy,w)*wpg*det;
    }
  }
  return sqrt(error);
}







double local_kinetic_energy(field * f,double* x,double *w){

  double wex[_INDEX_MAX];
  double l_ke=0;
  double t=f->tnow;
  
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      double omega=wglop(_DEG_V,iloc);
      double vi=-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      l_ke+=omega * _DV * w[ipg] * vi*vi ;
     }
  }
  return l_ke;
};

void Energies(field* f,double k_energy, double e_energy, double t_energy){
  
  k_energy=0;
  e_energy=0;
  t_energy=0;


  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],xphy[3],wpg;
      double dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      double w[f->model.m];
      for(int iv=0;iv<_INDEX_MAX+1;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      k_energy+=local_kinetic_energy(f,xphy,w)*wpg*det;
      e_energy+=w[_MV+1]*w[_MV+1]*wpg*det;      
    }
  }   
  
  t_energy=0.5*(e_energy+k_energy);
  
  f->Diagnostics[f->iter_time]=0.5*k_energy;
  f->Diagnostics[f->iter_time+f->itermax]=0.5*e_energy;
  f->Diagnostics[f->iter_time+2*f->itermax]=t_energy;
}


void Plot_Energies(field* f){
  int nb_diag=0;
  double e_energy=0, k_energy=0, t_energy=0;
  FILE * Plot;
  Plot = fopen("Diagnostics.dat","w");

  for(int i=1;i<f->itermax+1;i++){
    f->tnow=i*f->dt;
    k_energy=f->Diagnostics[i];
    e_energy=f->Diagnostics[i+f->itermax];
    t_energy=f->Diagnostics[i+2*f->itermax];
    fprintf(Plot, "%lf %lf %lf %lf\n", f->tnow,k_energy,e_energy,t_energy);
  }
  fclose(Plot);
}

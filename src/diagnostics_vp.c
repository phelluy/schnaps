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
//! \param[in] t : type of L2norm. if type_norm=0 this is the
//! numerical solution if type_norm=1 this is the error
real L2VelError(field *f, real *x, real *w){

  real wex[_INDEX_MAX];
  real err2 = 0;
  real t = f->tnow;
  f->model.ImposedData(x, t, wex);
  // loop on the finite emlements
  for(int iel = 0; iel < _NB_ELEM_V; iel++){
    // loop on the local glops
    for(int iloc = 0; iloc < _DEG_V + 1; iloc++){
      real omega = wglop(_DEG_V, iloc);
      real vi = -_VMAX + iel*_DV + _DV * glop(_DEG_V, iloc);
      int ipg = iloc + iel * _DEG_V;
      err2 += omega * _DV * (w[ipg] - wex[ipg]) * (w[ipg] - wex[ipg]);
    }
  }
  
  return err2;
}

real L2_Kinetic_error(field* f){
  real error = 0;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++){
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++){
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++){
      real xpgref[3], xphy[3], wpg;
      real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL, -1, // dpsiref,ifa
	      xphy, dtau,  // xphy,dtau
	      codtau, NULL, NULL); // codtau,dpsi,vnds
      real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      real w[f->model.m];
      for(int iv = 0;iv < f->model.m; iv++){
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	w[iv] = f->wn[imem];
      }
      // get the exact value
      error += L2VelError(f, xphy, w) * wpg * det;
    }
  }
  return sqrt(error);
}

real local_kinetic_energy(field *f,real *x, real *w) {

  real wex[_INDEX_MAX];
  real l_ke=0;
  real t=f->tnow;
  
  // loop on the finite emlements
  for(int iel = 0; iel < _NB_ELEM_V; iel++){
    // loop on the local glops
    for(int iloc = 0; iloc < _DEG_V + 1; iloc++){
      real omega = wglop(_DEG_V, iloc);
      real vi = -_VMAX + iel * _DV + _DV * glop(_DEG_V, iloc);
      int ipg = iloc + iel * _DEG_V;
      l_ke += omega * _DV * w[ipg] * vi * vi ;
     }
  }
  return l_ke;
}


// TODO: do not store all diagnotics for all time, but instead just
// append to the output file.
void Energies(field *f, real *w, real k_energy, real e_energy, real t_energy) {
  
  k_energy = 0;
  e_energy = 0;
  t_energy = 0;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++){
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++){
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++){
      real xpgref[3], xphy[3], wpg;
      real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      real wn[f->model.m];
      for(int iv = 0; iv < _INDEX_MAX + 1; iv++){ 
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	wn[iv] = w[imem];
      }
      // get the exact value
      k_energy += local_kinetic_energy(f, xphy, wn) * wpg * det;
      e_energy += wn[_MV+1] * wn[_MV+1] * wpg * det;      
    }
  }   
  
  t_energy = 0.5 * (e_energy + k_energy);
  
  f->Diagnostics[f->iter_time] = 0.5 * k_energy;
  f->Diagnostics[f->iter_time + f->itermax] = 0.5 * e_energy;
  f->Diagnostics[f->iter_time + 2 * f->itermax] = t_energy;
}

void Plot_Energies(field *f, real dt) {
  int nb_diag = 0;
  real e_energy = 0, k_energy = 0, t_energy = 0;
  FILE *Plot;
  Plot = fopen("Diagnostics.dat","w");

  for(int i = 1; i < f->itermax + 1; i++){
    f->tnow = i * dt; // FIXME: this will break with adaptive time-stepping
    k_energy = f->Diagnostics[i];
    e_energy = f->Diagnostics[i + f->itermax];
    t_energy = f->Diagnostics[i + 2 * f->itermax];
    fprintf(Plot, "%lf %lf %lf %lf\n", f->tnow, k_energy, e_energy, t_energy);
  }
  fclose(Plot);
}

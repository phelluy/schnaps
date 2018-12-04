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

schnaps_real L2VelError(field *f, schnaps_real *x, schnaps_real *w){
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wex[kd->index_max];
  schnaps_real err2 = 0;
  schnaps_real t = f->tnow;
  f->model.ImposedData(x, t, wex);
  // loop on the finite emlements
  for(int iel = 0; iel < kd->nb_elem_v; iel++){
    // loop on the local glops
    for(int iloc = 0; iloc < kd->deg_v + 1; iloc++){
      schnaps_real omega = wglop(kd->deg_v, iloc);
      schnaps_real vi = -kd->vmax + iel*kd->dv + kd->dv * glop(kd->deg_v, iloc);
      int ipg = iloc + iel * kd->deg_v;
      err2 += omega * kd->dv * (w[ipg] - wex[ipg]) * (w[ipg] - wex[ipg]);
    }
  }
  
  return err2;
}

schnaps_real L2_Kinetic_error(Simulation* simu){
  schnaps_real error = 0;
  field *f =&simu->fd[0];

    // get the physical nodes of element ie

    // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){
      schnaps_real xpgref[3], xphy[3], wpg;
      schnaps_real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref,  // xref
	      NULL, -1, // dpsiref,ifa
	      xphy, dtau,  // xphy,dtau
	      codtau, NULL, NULL); // codtau,dpsi,vnds
      schnaps_real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      schnaps_real w[f->model.m];
      for(int iv = 0;iv < f->model.m; iv++){
	int imem = f->varindex(f->deg,f->raf, f->model.m, ipg, iv);
	w[iv] = f->wn[imem];
      }
      // get the exact value
      error += L2VelError(f, xphy, w) * wpg * det;
    }
  
  return sqrt(error);
}

schnaps_real local_kinetic_energy(field *f,schnaps_real *x, schnaps_real *w) {

  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wex[kd->index_max];
  schnaps_real l_ke=0;
  schnaps_real t=f->tnow;

  // loop on the finite emlements
  for(int iel = 0; iel < kd->nb_elem_v; iel++){
    // loop on the local glops
    for(int iloc = 0; iloc < kd->deg_v + 1; iloc++){
      schnaps_real omega = wglop(kd->deg_v, iloc);
      schnaps_real vi = -kd->vmax + iel * kd->dv + kd->dv * glop(kd->deg_v, iloc);
      int ipg = iloc + iel * kd->deg_v;
      l_ke += omega * kd->dv * w[ipg] * vi * vi ;
     }
  }
  return l_ke;
}


// TODO: do not store all diagnotics for all time, but instead just
// append to the output file.
void Energies(Simulation *simu, schnaps_real *w, schnaps_real k_energy,
	      schnaps_real e_energy, schnaps_real t_energy,int first_diag) {

  KineticData * kd = &schnaps_kinetic_data;
  k_energy = 0;
  e_energy = 0;
  t_energy = 0;

  field * f = &simu->fd[0];

  // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){

      schnaps_real xpgref[3], xphy[3], wpg;
      schnaps_real dtau[3][3], codtau[3][3];//,xpg[3];

      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds

      schnaps_real det	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2];

      schnaps_real wn[f->model.m];
      //printf("m=%d imax=%d\n",f->model.m,_INDEX_MAX + 1);
      for(int iv = 0; iv < kd->index_max + 1; iv++){
	int imem = f->varindex(f->deg,f->raf, f->model.m, ipg, iv);
	wn[iv] = w[imem];
      }
      // get the exact value
      k_energy += local_kinetic_energy(f, xphy, wn) * wpg * det;
      e_energy += wn[kd->index_ex] * wn[kd->index_ex] * wpg * det;
    }

  t_energy = 0.5 * (e_energy + k_energy);

  simu->Diagnostics[simu->iter_time_rk + (first_diag-1) * simu->itermax_rk]
    = 0.5 * k_energy;
  simu->Diagnostics[simu->iter_time_rk + (first_diag) * simu->itermax_rk]
    = 0.5 * e_energy;
  simu->Diagnostics[simu->iter_time_rk + (first_diag+1) * simu->itermax_rk]
    = t_energy;
}

void Charge_total(Simulation *simu, schnaps_real *w, schnaps_real t_charge,int first_diag) {
   KineticData * kd=&schnaps_kinetic_data;
  t_charge=0;
  field *f=&simu->fd[0];

    // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){
      schnaps_real xpgref[3], xphy[3], wpg;
      schnaps_real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      schnaps_real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      schnaps_real wn[f->model.m];
      for(int iv = 0; iv < kd->index_max + 1; iv++){ 
	int imem = f->varindex(f->deg,f->raf, f->model.m, ipg, iv);
	wn[iv] = w[imem];
      }
      t_charge += wn[kd->index_rho] * wpg * det;
    }

  simu->Diagnostics[simu->iter_time_rk + (first_diag-1) * simu->itermax_rk] = t_charge;
}

void Taux_instability(Simulation *simu, schnaps_real *w, schnaps_real m_1, schnaps_real taux_ins_m1,int first_diag) {

  KineticData * kd=&schnaps_kinetic_data;
  
  taux_ins_m1=0;
  schnaps_real taux_ins_r_m1 = 0;
  schnaps_real taux_ins_i_m1 = 0;
  
  schnaps_real pi = 4.0*atan(1.0);
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field * f = simu->fd + ie;

  // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){
      schnaps_real xpgref[3], xphy[3], wpg;
      schnaps_real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
  	      xpgref,  // xref
  	      NULL,-1, // dpsiref,ifa
  	      xphy,dtau,  // xphy,dtau
  	      codtau,NULL,NULL); // codtau,dpsi,vnds
      schnaps_real det
  	      = dtau[0][0] * codtau[0][0]
  	      + dtau[0][1] * codtau[0][1]
  	      + dtau[0][2] * codtau[0][2];
  	      
      schnaps_real r = sqrt(xphy[0]*xphy[0]+xphy[1]*xphy[1]);
      schnaps_real theta = atan2(xphy[1],xphy[0]);

      int imemphi = f->varindex(f->deg,f->raf, f->model.m, ipg, kd->index_phi);
      
      taux_ins_r_m1 += cos(m_1*theta) * f->wn[imemphi] * wpg * det;
      //taux_ins_r_m1 += wpg * det;
      taux_ins_i_m1 += sin(m_1*theta) * f->wn[imemphi] * wpg * det;

    }
  }
  
  taux_ins_r_m1 /= (pi*100);
  taux_ins_i_m1 /= (pi*100);
  taux_ins_m1 = sqrt(taux_ins_r_m1 * taux_ins_r_m1 + taux_ins_i_m1 * taux_ins_i_m1)/2;
  simu->Diagnostics[simu->iter_time_rk + (first_diag-1) * simu->itermax_rk] = taux_ins_m1;
/*  taux_ins_m1 =taux_ins_r_m1 ;
  simu->Diagnostics[simu->iter_time_rk + (first_diag-1) * simu->itermax_rk] = taux_ins_m1;
printf("integrale =  %.15e\n", taux_ins_m1);	
assert(1==2);*/
}
  


void Plot_Energies(Simulation *simu, schnaps_real dt) {

  schnaps_real e_energy = 0, k_energy = 0, t_energy = 0, t_charge=0;
  schnaps_real taux_ins_m1 = 0,  taux_ins_m2 = 0, val = 0;
  FILE *Plot;

  Plot = fopen("Diagnostics.dat","w");

  for(int i = 1; i < simu->itermax_rk ; i++){

    simu->tnow = i*dt;
    fprintf(Plot, "%.11e ", simu->tnow);

    for (int j=1; j<simu->nb_diags; j++) {
      val = simu->Diagnostics[i + simu->itermax_rk*(j-1)];
      fprintf(Plot, "%.11e ", simu->tnow);
    }
    fprintf(Plot, "\n");
  }
  fclose(Plot);
}

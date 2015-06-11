#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"





void Computation_charge_density(field *f, real * w){
  
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      int imemc=f->varindex(f->interp_param,ie,ipg,_MV+2);
      w[imemc]=0;
  
      for(int ielv=0;ielv<_NB_ELEM_V;ielv++){
	// loop on the local glops
	for(int iloc=0;iloc<_DEG_V+1;iloc++){
	  real omega=wglop(_DEG_V,iloc);
	  real vi=-_VMAX+ielv*_DV+_DV*glop(_DEG_V,iloc);
	  int ipgv=iloc+ielv*_DEG_V;
	  int imem=f->varindex(f->interp_param,ie,ipg,ipgv);
	  w[imemc]+=omega*_DV*w[imem];
	}
      }
    }
  }
  
}


real Computation_charge_average(field *f,real * w) {
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  real average = 0;
  real rho_imem = 0;
  real size_domain = 0;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->interp_param + 1);
    for(int ipg = 0; ipg < npg; ipg++) {
	int imem = f->varindex(f->interp_param, ie, ipg, _INDEX_RHO);
	rho_imem = f->wn[imem];
      
      real wpg, det;
      // Compute wpg, det, and the exact solution
      { 
	real xphy[3], xpgref[3];
	real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
	Ref2Phy(physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);
      }


        average += rho_imem * wpg * det;
	size_domain +=  wpg * det;

      
    }
  }
  return average/size_domain;
}


void ComputeElectricField(field* f){

  int nraf[3] = {f->interp_param[4], 
		 f->interp_param[5],
		 f->interp_param[6]};
  
  int npg[3] = {f->interp_param[1] + 1, 
		f->interp_param[2] + 1,
		f->interp_param[3] + 1};
    
  int nnodes = npg[0] * npg[1] * npg[2] ;
 
  int npgmacrocell = nnodes * nraf[0] * nraf[1] * nraf[2];


  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    for(int ipg = 0;ipg < npgmacrocell; ipg++){
      //real wpg;
      real xref[3];

      ref_pg_vol(f->interp_param+1,ipg,xref,NULL,NULL);
      int iex = f->varindex(f->interp_param,ie,
			    ipg,_INDEX_EX);
      f->wn[iex] = 0;
      
      for(int ib=0; ib < npgmacrocell; ib++){
	real dtau[3][3],codtau[3][3];
	real dphiref[3];
	real dphi[3];
	grad_psi_pg(f->interp_param+1,ib,ipg,dphiref);
	Ref2Phy(physnode,xref,dphiref,0,NULL,
		  dtau,codtau,dphi,NULL);
	real det = dot_product(dtau[0], codtau[0]);
	int ipot = f->varindex(f->interp_param,ie,
			   ib,_INDEX_PHI);
	f->wn[iex] -= f->wn[ipot] * dphi[0] / det;
      }
    }
 

  }
}

void Compute_electric_field(field* f, real * w){

  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    const int m = f->model.m;
    const int deg[3]={f->interp_param[1],
		      f->interp_param[2],
		      f->interp_param[3]};
    const int npg[3] = {deg[0]+1,
			deg[1]+1,
			deg[2]+1};
    const int nraf[3]={f->interp_param[4],
		       f->interp_param[5],
		       f->interp_param[6]};

    const unsigned int sc_npg=npg[0]*npg[1]*npg[2];


    int f_interp_param[8]= {f->interp_param[0],
			    f->interp_param[1],
			    f->interp_param[2],
			    f->interp_param[3],
			    f->interp_param[4],
			    f->interp_param[5],
			    f->interp_param[6],
			    f->interp_param[7]};

    // loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++){
      for(int icL1 = 0; icL1 < nraf[1]; icL1++){
	for(int icL2 = 0; icL2 < nraf[2]; icL2++){

	  int icL[3] = {icL0,icL1,icL2};
	  // get the L subcell id
	  int ncL=icL[0]+nraf[0]*(icL[1]+nraf[1]*icL[2]);
	  // first glop index in the subcell
	  int offsetL=npg[0]*npg[1]*npg[2]*ncL;

	  // compute all of the xref for the subcell
	  real *xref0 = malloc(sc_npg * sizeof(real));
	  real *xref1 = malloc(sc_npg * sizeof(real));
	  real *xref2 = malloc(sc_npg * sizeof(real));
	  int *imems = malloc(m * sc_npg * sizeof(int));
	  int pos=0;
	  for(unsigned int p=0; p < sc_npg; ++p) {
	    real xref[3];
	    real tomega;
	    ref_pg_vol(f->interp_param+1,offsetL+p,xref,&tomega,NULL);
	    xref0[p] = xref[0];
	    xref1[p] = xref[1];
	    xref2[p] = xref[2];
	    
	    for(int im=0; im < m; ++im) {
	      imems[pos++] = f->varindex(f_interp_param,ie,offsetL+p,im);
	    }
	  }

	  // loop in the "cross" in the three directions
	  // TO DO activate the 3D case !!!
	  //for(int dim0 = 0; dim0 < 3; dim0++){
	  for(int dim0 = 0; dim0 < 1; dim0++){
	    // point p at which we compute the flux
    
	    for(int p0 = 0; p0 < npg[0]; p0++){
	      for(int p1 = 0; p1 < npg[1]; p1++){
		for(int p2 = 0; p2 < npg[2]; p2++){
		  int p[3]={p0,p1,p2};
		  int ipgL=offsetL+p[0]+npg[0]*(p[1]+npg[1]*p[2]);
		  /* for(int iv=0; iv < m; iv++){ */
		  /*   ///int imemL=f->varindex(f_interp_param,ie,ipgL,iv); */
		  /*   wL[iv] = f->wn[imems[m*(ipgL-offsetL)+iv]];  */

		  /*   //wL[iv] = f->wn[imemL]; */
		  /* } */
		  //wn[imems[m*(ipgL-offsetL)+_MV+dim0]]=0;
		  real gradx[3]={0,0,0};

		  int q[3]={p[0],p[1],p[2]};
		  // loop on the direction dim0 on the "cross"
		  for(int iq = 0; iq < npg[dim0]; iq++){
		    q[dim0]=(p[dim0]+iq)%npg[dim0];
		    int ipgR=offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		    real phiq=w[imems[m*(ipgR-offsetL)+_INDEX_PHI]];
		    real dphiref[3]={0,0,0};
		    // compute grad phi_q at glop p
		    dphiref[dim0]=dlag(deg[dim0],q[dim0],p[dim0])*nraf[dim0];

		    real xrefL[3]={xref0[ipgL-offsetL],
		    		     xref1[ipgL-offsetL],
		    		     xref2[ipgL-offsetL]};
		    //real wpgL=omega[ipgL-offsetL];
		    /* real xrefL[3], wpgL; */
		    /* ref_pg_vol(f->interp_param+1,ipgL,xrefL,&wpgL,NULL); */

		    // mapping from the ref glop to the physical glop
		    real dtau[3][3],codtau[3][3],dphiL[3];
		    Ref2Phy(physnode,
			    xrefL,
			    dphiref,  // dphiref
			    -1,    // ifa                                 
			    NULL,  // xphy  
			    dtau,
			    codtau,
			    dphiL,  // dphi
			    NULL);  // vnds   
    
		    //f->model.NumFlux(wL,wL,dphiL,flux);

		    gradx[0]+=phiq*dphiL[0];
		  } // iq
		  w[imems[m*(ipgL-offsetL)+_INDEX_EX+dim0]]=gradx[0];		  
		  //printf("grad=%f\n",gradx[0]);
		} // p2
	      } // p1
	    } // p0
	    
	  } // dim loop
	  
	  free(xref0);
	  free(xref1);
	  free(xref2);
	  free(imems);

	} // icl2
      } //icl1
    } // icl0
  }
  

}


#include "implicit.h"
#include <stdlib.h>


//void InternalCoupling(Simulation *simu,  LinearSolver *solver, int itest);
//void FluxCoupling(Simulation *simu,  LinearSolver *solver,int itest);

/* void InternalAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void FluxAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void InterfaceAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void SourceAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void MassAssembly(Simulation *simu,  LinearSolver *solver); */

/* void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver,real theta, real dt); */


void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver){

  int neq = simu->wsize;

  MatrixStorage ms = SKYLINE;
  Solver st = LU; 
  InitLinearSolver(solver,neq,&ms,&st);

  int itest = 1;

  for (int isky=0 ; isky < itest; isky++){
  
    InternalCoupling(simu, solver, isky);
    FluxCoupling(simu, solver, isky);
    InterfaceCoupling(simu, solver, isky);

    if (isky == 0) AllocateLinearSolver(solver);
  }

  //DisplayLinearSolver(solver);

  

}

void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver,schnaps_real theta, schnaps_real dt){

  if(solver->mat_is_assembly == false){
    MassAssembly(simu, solver);
    InternalAssembly(simu, solver,theta,dt);
    FluxAssembly(simu, solver,theta,dt);
    InterfaceAssembly(simu, solver,theta,dt);
  }

  if(solver->rhs_is_assembly == false){
    for(int i=0;i<solver->neq;i++){
      solver->rhs[i]=0;
    }
      SourceAssembly(simu, solver,theta,dt);
      
  }
  //DisplayLinearSolver(solver);

}


void ThetaTimeScheme(Simulation *simu, schnaps_real tmax, schnaps_real dt){

  LinearSolver solver_implicit;
  LinearSolver solver_explicit;  

  schnaps_real theta=0.5;
  simu->dt=dt;
  
  int itermax=tmax/simu->dt+1;
  simu->itermax_rk=itermax;
  InitImplicitLinearSolver(simu, &solver_implicit);
  InitImplicitLinearSolver(simu, &solver_explicit);
  schnaps_real *res = calloc(simu->wsize, sizeof(schnaps_real));

  simu->tnow=0;
  for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
    simu->fd[ie].tnow=simu->tnow;
  } 

  for(int tstep=0;tstep<simu->itermax_rk;tstep++){
  

    if(tstep==0){ 
      solver_implicit.mat_is_assembly=false;
      solver_explicit.mat_is_assembly=false;
    } 
    else 
      { 
      solver_implicit.mat_is_assembly=true;
      solver_explicit.mat_is_assembly=true;
    } 

    solver_implicit.rhs_is_assembly=false;
    solver_explicit.rhs_is_assembly=false;

    
    AssemblyImplicitLinearSolver(simu, &solver_explicit,-(1.0-theta),simu->dt);
    simu->tnow=simu->tnow+simu->dt;
    for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
      simu->fd[ie].tnow=simu->tnow;
    } 
    AssemblyImplicitLinearSolver(simu, &solver_implicit,theta,simu->dt);
  

    MatVect(&solver_explicit, simu->w, res);

    for(int i=0;i<solver_implicit.neq;i++){
      solver_implicit.rhs[i]=-solver_explicit.rhs[i]+solver_implicit.rhs[i]+res[i];
    }
  
    SolveLinearSolver(&solver_implicit);

    for(int i=0;i<solver_implicit.neq;i++){
      simu->w[i]=solver_implicit.sol[i];
    }
    int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, tstep, simu->itermax_rk, dt);
  }
  
}

void InternalCoupling(Simulation *simu,  LinearSolver *solver, int isky){

  //for(int isky = 0; isky < itest; isky++){
  
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie; // &(simu->fd[ie])
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};

    //const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};
	  // get the L subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // first glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  // loop in the "cross" in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) {
	    //for(int dim0 = 0; dim0 < 2; dim0++) {  // TODO : return to 3d !

	    // point p at which we compute the flux
	    for(int p0 = 0; p0 < npg[0]; p0++) {
	      for(int p1 = 0; p1 < npg[1]; p1++) {
		for(int p2 = 0; p2 < npg[2]; p2++) {

		  int p[3] = {p0, p1, p2};
		  int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
		  for(int iv1 = 0; iv1 < m; iv1++) {
		    int imemL = f->varindex(deg, nraf, m, ipgL, iv1) + offsetw;

		    int q[3] = {p[0], p[1], p[2]};
		    // loop on the direction dim0 on the "cross"
		    for(int iq = 0; iq < npg[dim0]; iq++) {
		      q[dim0] = (p[dim0] + iq) % npg[dim0];

		      int ipgR = offsetL + q[0] + npg[0] * (q[1] + npg[1] * q[2]);
		      for(int iv2 = 0; iv2 < m; iv2++) {
			int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv2) + offsetw;
			if (isky ==0) IsNonZero(solver, imemL, imemR);
			if (isky ==1) AddLinearSolver(solver, imemL, imemR,1);
			// dtw[imems[temp]] += flux[iv] * wpgL;
		      }
		    } // iv1
		  } // iq
		} // p2
	      } // p1
	    } // p0

	    
	  } // dim loop


	} // icl2
      } //icl1
    } // icl0


  }



  
}





void FluxCoupling(Simulation *simu,  LinearSolver *solver, int isky){

  
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie; // &(simu->fd[ie])
    int offsetw = f->wsize * ie;
    const int m = f->model.m;


    const int nraf[3] = {f->raf[0],
			 f->raf[1],
			 f->raf[2]};
    const int deg[3] = {f->deg[0],
			f->deg[1],
			f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};

    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};

	  // Get the left subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	  // Sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	    // Compute the subface flux only if we do not touch the
	    // subcell boundary along the current direction dim0
	    if (icL[dim0] != nraf[dim0] - 1) {
	      int icR[3] = {icL[0], icL[1], icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
	      int offsetR = npg[0] * npg[1] * npg[2] * ncR;

	      // FIXME: write only write to L-values (and do both
	      // faces) to parallelise better.

	      const int altdim1[3] = {1, 0, 0};
	      const int altdim2[3] = {2, 2, 1};

	      // now loop on the left glops of the subface
	      //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	      int dim1 = altdim1[dim0];
	      int dim2 = altdim2[dim0];
	      int iL[3];
	      iL[dim0] = deg[dim0];
	      for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++) {
		for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++) {
		  // find the right and left glops volume indices

		  int iR[3] = {iL[0], iL[1], iL[2]};
		  iR[dim0] = 0;

		  int ipgL = offsetL 
		    + iL[0] + (deg[0] + 1) * (iL[1] + (deg[1] + 1) * iL[2]);
		  int ipgR = offsetR 
		    + iR[0] + (deg[0] + 1) * (iR[1] + (deg[1] + 1) * iR[2]);
		  //printf("ipgL=%d ipgR=%d\n", ipgL, ipgR);

		  for(int iv1 = 0; iv1 < m; iv1++) {
		    int imemL = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv1) + offsetw; 

		    // finally distribute the flux on the two sides
		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imemR = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2) + offsetw;
		      if (isky ==0) IsNonZero(solver, imemL, imemR);
		      if (isky ==1) {
			AddLinearSolver(solver, imemL, imemR,1);
			AddLinearSolver(solver, imemR, imemL,1);
		      }
		    }
		  }

		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop
	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

  }
  

}

void InternalAssembly(Simulation *simu,  LinearSolver *solver,schnaps_real theta, schnaps_real dt){

  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie;
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;

    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};

    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];


    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {
	  
	  int icL[3] = {icL0, icL1, icL2};
	  // get the L subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // first glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  
	  // compute all of the xref for the subcell
	  schnaps_real *xref0 = malloc(sc_npg * sizeof(schnaps_real));
	  schnaps_real *xref1 = malloc(sc_npg * sizeof(schnaps_real));
	  schnaps_real *xref2 = malloc(sc_npg * sizeof(schnaps_real));
	  schnaps_real *omega = malloc(sc_npg * sizeof(schnaps_real));
	  int *imems = malloc(m * sc_npg * sizeof(int));
	  int pos = 0;
	  for(unsigned int p = 0; p < sc_npg; ++p) {
	    schnaps_real xref[3];
	    schnaps_real tomega;
	    
	    ref_pg_vol(f->deg, f->raf, offsetL + p, xref, &tomega, NULL);
	    xref0[p] = xref[0];
	    xref1[p] = xref[1];
	    xref2[p] = xref[2];
	    omega[p] = tomega;
	    
	    for(int im = 0; im < m; ++im) {
	      imems[pos++] = f->varindex(f->deg,f->raf,f->model.m, offsetL + p, im) + offsetw;
	    }
	  }

	  // loop in the "cross" in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) {
	    //for(int dim0 = 0; dim0 < 2; dim0++) {  // TODO : return to 3d !
	    // point p at which we compute the flux
	    
	    for(int p0 = 0; p0 < npg[0]; p0++) {
	      for(int p1 = 0; p1 < npg[1]; p1++) {
		for(int p2 = 0; p2 < npg[2]; p2++) {
		  schnaps_real wL[m], flux[m];
		  int p[3] = {p0, p1, p2};
		  int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
		  int q[3] = {p[0], p[1], p[2]};
		  // loop on the direction dim0 on the "cross"
		  for(int iq = 0; iq < npg[dim0]; iq++) {
		    q[dim0] = (p[dim0] + iq) % npg[dim0];
		    schnaps_real dphiref[3] = {0, 0, 0};
		    // compute grad phi_q at glop p
		    dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) 
		      * nraf[dim0];
		    
		    schnaps_real xrefL[3] = {xref0[ipgL - offsetL],
				     xref1[ipgL - offsetL],
				     xref2[ipgL - offsetL]};
		    schnaps_real wpgL = omega[ipgL - offsetL];
		    /* real xrefL[3], wpgL; */
		    /* ref_pg_vol(f->interp_param+1,ipgL,xrefL, &wpgL, NULL); */

		    // mapping from the ref glop to the physical glop
		    schnaps_real dtau[3][3], codtau[3][3], dphiL[3];
		    schnaps_ref2phy(f->physnode,
			    xrefL,
			    dphiref, // dphiref
			    -1,  // ifa
			    NULL, // xphy
			    dtau,
			    codtau,
			    dphiL, // dphi
			    NULL);  // vnds


		    for(int iv1 = 0; iv1 < m; iv1++) {
		      int imemL = f->varindex(deg, nraf, m, ipgL, iv1) + offsetw;
		      for(int iv = 0; iv < m; iv++) {
			wL[iv] = (iv == iv1);
		      }

		      f->model.NumFlux(wL, wL, dphiL, flux);

		      int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		      for(int iv2 = 0; iv2 < m; iv2++) {
			schnaps_real val = theta * dt * flux[iv2] * wpgL;
			int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv2) + offsetw;
			AddLinearSolver(solver, imemR, imemL,-val);
		      }
		    }
		  } // iq
		} // p2
	      } // p1
	    } // p0

	  } // dim loop

	  free(omega);
	  free(xref0);
	  free(xref1);
	  free(xref2);
	  free(imems);

	} // icl2
      } //icl1
    } // icl0

      
  }

}


void FluxAssembly(Simulation *simu, LinearSolver *solver,schnaps_real theta, schnaps_real dt){


  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie; // &(simu->fd[ie])
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;

  

    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    int npg[3] = {deg[0] + 1,
		  deg[1] + 1,
		  deg[2] + 1};

 
    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {
	  
	  int icL[3] = {icL0, icL1, icL2};
	  
	  // Get the left subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  
	  // Sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	    // Compute the subface flux only if we do not touch the
	    // subcell boundary along the current direction dim0
	    if (icL[dim0] != nraf[dim0] - 1) {
	      int icR[3] = {icL[0], icL[1], icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
	      int offsetR = npg[0] * npg[1] * npg[2] * ncR;
	      
	      // FIXME: write only write to L-values (and do both
	      // faces) to parallelise better.
	      
	      const int altdim1[3] = {1, 0, 0};
	      const int altdim2[3] = {2, 2, 1};
	      
	      // now loop on the left glops of the subface
	      //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	      int dim1 = altdim1[dim0];
	      int dim2 = altdim2[dim0];
	      int iL[3];
	      iL[dim0] = deg[dim0];
	      for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++) {
		for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++) {
		  // find the right and left glops volume indices

		  int iR[3] = {iL[0], iL[1], iL[2]};
		  iR[dim0] = 0;
		  
		  int ipgL = offsetL 
		    + iL[0] + (deg[0] + 1) * (iL[1] + (deg[1] + 1) * iL[2]);
		  int ipgR = offsetR 
		    + iR[0] + (deg[0] + 1) * (iR[1] + (deg[1] + 1) * iR[2]);
		  //printf("ipgL=%d ipgR=%d\n", ipgL, ipgR);

		  // Compute the normal vector for integrating on the
		  // face
		  schnaps_real vnds[3];
		  {
		    schnaps_real xref[3], wpg3;
		    ref_pg_vol(f->deg, f->raf, ipgL, xref, &wpg3, NULL);
		    // mapping from the ref glop to the physical glop
		    schnaps_real dtau[3][3], codtau[3][3];
		    schnaps_ref2phy(f->physnode,
			    xref,
			    NULL, // dphiref
			    -1,  // ifa
			    NULL, // xphy
			    dtau,
			    codtau,
			    NULL, // dphi
			    NULL);  // vnds
		    // we compute ourself the normal vector because we
		    // have to take into account the subcell surface

		    schnaps_real h1h2 = 1. / nraf[dim1] / nraf[dim2];
		    vnds[0] = codtau[0][dim0] * h1h2;
		    vnds[1] = codtau[1][dim0] * h1h2;
		    vnds[2] = codtau[2][dim0] * h1h2;
		  }


		  schnaps_real wpg
		    = wglop(deg[dim1], iL[dim1])
		    * wglop(deg[dim2], iL[dim2]);
		  
		  
		  // numerical flux from the left and right state and
		  // normal vector
		  schnaps_real wL[m], wR[m], flux[m];

		  for (int iv1 = 0; iv1 < m; iv1++){

			    

		    for(int iv = 0; iv < m; iv++) {
		      wL[iv] = (iv == iv1);
		      wR[iv] = 0;
		    }
		    int imem1 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv1)+offsetw;

		    f->model.NumFlux(wL, wR, vnds, flux);	

		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv2)+offsetw;		  
		      schnaps_real val = theta * dt * flux[iv2] * wpg;		      
		      AddLinearSolver(solver, imem2, imem1, val);
		      
		      imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2)+offsetw;
		      val = theta * dt * flux[iv2] * wpg;		      
		      AddLinearSolver(solver, imem2, imem1, -val);
		    }
		  
		    for(int iv = 0; iv < m; iv++) {
		      wL[iv] = 0;
		      wR[iv] = (iv == iv1);
		    }
		    imem1 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv1)+offsetw;


		    f->model.NumFlux(wL, wR, vnds, flux);

		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv2)+offsetw;
		      schnaps_real val = theta * dt * flux[iv2] * wpg;
		      AddLinearSolver(solver, imem2, imem1, val);
		    
		      imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2)+offsetw;		    
		      val = theta *dt * flux[iv2] * wpg;		    
		      AddLinearSolver(solver, imem2, imem1, -val);
		    }
		  }

		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop
	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

  }
}



void MassAssembly(Simulation *simu,  LinearSolver *solver){

  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie;
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;

    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};

    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};


    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      schnaps_real det = dot_product(dtau[0], codtau[0]);
      for(int iv1 = 0; iv1 < m; iv1++) {
	int imem = f->varindex(deg, nraf, m, ipg, iv1)+offsetw;
	schnaps_real val = wpg * det;
	AddLinearSolver(solver, imem, imem,val);
      }
    }
 
      

  }

}

void SourceAssembly(Simulation *simu,  LinearSolver *solver, schnaps_real theta, schnaps_real dt){

  if(simu->fd[0].model.Source != NULL) {
    for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
      field *f = simu->fd + ie;
      int offsetw = f->wsize * ie;
    
    const int m = f->model.m;
    
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};

    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};


    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      schnaps_real det = dot_product(dtau[0], codtau[0]);
      schnaps_real wL[m], source[m];
      /* for(int iv = 0; iv < m; ++iv){ */
      /* 	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv); */
      /* 	wL[iv] = w[imem]; */
      /* } */
      f->model.Source(xphy, f->tnow, wL, source);

      for(int iv1 = 0; iv1 < m; iv1++) {
	int imem = f->varindex(deg, nraf, m, ipg, iv1)+offsetw;
	schnaps_real val = theta * dt * source[iv1] * wpg * det;
	solver->rhs[imem] += val;
      }
    }


 
     
  }
 }
  // assembly of the boundary terms

  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR < 0) {
 
      const unsigned int m = fL->model.m;
      
      // Loop over the points on a single macro cell interface.
      for(int ipgfL = 0; ipgfL < NPGF(fL->deg, fL->raf, locfaL); ipgfL++) {
	
	schnaps_real xpgref[3], xpgref_in[3], wpg;
	
	// Get the coordinates of the Gauss point and coordinates of a
	// point slightly inside the opposite element in xref_in
	int ipgL = ref_pg_face(fL->deg, fL->raf, locfaL, ipgfL, xpgref, &wpg, xpgref_in);
	
	// Normal vector at gauss point ipgL
	schnaps_real vnds[3], xpg[3];
	{
	  schnaps_real dtau[3][3], codtau[3][3];
	  schnaps_ref2phy(fL->physnode,
		  xpgref,
		  NULL, locfaL, // dpsiref, ifa
		  xpg, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
	}
	
 	// the boundary flux is an affine function
	schnaps_real flux0[m], wL[m];
	for(int iv = 0; iv < m; iv++) {
	  wL[iv] = 0;
	}

	fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux0);
	
	for(int iv2 = 0; iv2 < m; iv2++) {
	  int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2)+offsetL;
	  schnaps_real val = theta *dt * flux0[iv2] * wpg;
	  solver->rhs[imem2] -= val;
	}
      }
    } // if ier < 0
  } // macroface loop


 
} // SourceAssembly


/* void DGMacroCellInterface(int locfaL, */
/* 			  field *fL, int offsetL, field *fR, int offsetR, */
/* 			  real *w, real *dtw)  */
void InterfaceCoupling(Simulation *simu,  LinearSolver *solver, int itest)
{

  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      offsetR = fsize * ieR;
    }

  
    const unsigned int m = fL->model.m;


    // Loop over the points on a single macro cell interface.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ipgfL = 0; ipgfL < NPGF(fL->deg, fL->raf, locfaL); ipgfL++) {

      schnaps_real xpgref[3], xpgref_in[3], wpg;
    
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipgL = ref_pg_face(fL->deg, fL->raf, locfaL, ipgfL, xpgref, &wpg, xpgref_in);
    
      // Normal vector at gauss point ipgL
      schnaps_real vnds[3], xpg[3];
      {
	schnaps_real dtau[3][3], codtau[3][3];
	schnaps_ref2phy(fL->physnode,
		xpgref,
		NULL, locfaL, // dpsiref, ifa
		xpg, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
      }
    
      if (fR != NULL) {  // the right element exists
	schnaps_real xrefL[3];
	{
	  schnaps_real xpg_in[3];
	  schnaps_ref2phy(fL->physnode,
		  xpgref_in,
		  NULL, -1, // dpsiref, ifa
		  xpg_in, NULL,
		  NULL, NULL, NULL); // codtau, dpsi, vnds
	  PeriodicCorrection(xpg_in,fL->period);
	  schnaps_phy2ref(fR->physnode, xpg_in, xrefL);
	
	}
      
	int ipgR = ref_ipg(fR->deg,fR->raf, xrefL);



	for (int iv1 = 0; iv1 < m; iv1++){

			    
	  int imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv1) + offsetL;

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv2) + offsetL;		  
	    IsNonZero(solver, imem2, imem1);
		      
	    imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv2) + offsetR;
	    IsNonZero(solver, imem2, imem1);
	  }
		  
	  imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgR, iv1) + offsetR;


	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv2) + offsetL;
	    IsNonZero(solver, imem2, imem1);
		    
	    imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv2) + offsetR;		    
	    IsNonZero(solver, imem2, imem1);
	  }
	}

      } else { // The point is on the boundary.


	/* for(int iv2 = 0; iv2 < m; iv2++) { */
	/*   int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2); */
	/*   real val = theta *dt * flux0[iv2] * wpg; */
	/*   solver->rhs[imem2] -= val; */
	/* } */
      
	for(int iv1 = 0; iv1 < m; iv1++) {
	  int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv1) + offsetL;

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2) + offsetL;
	    IsNonZero(solver, imem2, imem1);
	  }
	} // iv1

      } // else

  
    } // ipgfl

  } // macroface loop

}

/* void DGMacroCellInterface(int locfaL, */
/* 			  field *fL, int offsetL, field *fR, int offsetR, */
/* 			  real *w, real *dtw)  */
void InterfaceAssembly(Simulation *simu,  LinearSolver *solver,schnaps_real theta, schnaps_real dt)
{

  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      offsetR = fsize * ieR;
    }

  
    const unsigned int m = fL->model.m;


    // Loop over the points on a single macro cell interface.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ipgfL = 0; ipgfL < NPGF(fL->deg, fL->raf, locfaL); ipgfL++) {

      schnaps_real xpgref[3], xpgref_in[3], wpg;
    
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipgL = ref_pg_face(fL->deg, fL->raf, locfaL, ipgfL, xpgref, &wpg, xpgref_in);
    
      // Normal vector at gauss point ipgL
      schnaps_real vnds[3], xpg[3];
      {
	schnaps_real dtau[3][3], codtau[3][3];
	schnaps_ref2phy(fL->physnode,
		xpgref,
		NULL, locfaL, // dpsiref, ifa
		xpg, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
      }
    
      if (fR != NULL) {  // the right element exists
	schnaps_real xrefL[3];
	{
	  schnaps_real xpg_in[3];
	  schnaps_ref2phy(fL->physnode,
		  xpgref_in,
		  NULL, -1, // dpsiref, ifa
		  xpg_in, NULL,
		  NULL, NULL, NULL); // codtau, dpsi, vnds
	  PeriodicCorrection(xpg_in,fL->period);
	  schnaps_phy2ref(fR->physnode, xpg_in, xrefL);
	
	}
      
	int ipgR = ref_ipg(fR->deg,fR->raf, xrefL);


	schnaps_real flux[m];
	schnaps_real wL[m];
	schnaps_real wR[m];

	for (int iv1 = 0; iv1 < m; iv1++){

			    

	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = (iv == iv1);
	    wR[iv] = 0;
	  }
	  int imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv1) + offsetL;

	  // int_dL F(wL, wR, grad phi_ib)

	  fL->model.NumFlux(wL, wR, vnds, flux);

	  // Add flux to both sides

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv2) + offsetL;		  
	    schnaps_real val = theta * dt * flux[iv2] * wpg;		      
	    AddLinearSolver(solver, imem2, imem1, val);
		      
	    imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv2) + offsetR;
	    val = theta * dt * flux[iv2] * wpg;		      
	    AddLinearSolver(solver, imem2, imem1, -val);
	  }
		  
	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = 0;
	    wR[iv] = (iv == iv1);
	  }
	  imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgR, iv1) + offsetR;


	  fL->model.NumFlux(wL, wR, vnds, flux);

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv2) + offsetL;
	    schnaps_real val = theta * dt * flux[iv2] * wpg;
	    AddLinearSolver(solver, imem2, imem1, val);
		    
	    imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv2) + offsetR;		    
	    val = theta *dt * flux[iv2] * wpg;		    
	    AddLinearSolver(solver, imem2, imem1, -val);
	  }
	}

      } else { // The point is on the boundary.

	// the boundary flux is an affine function
	schnaps_real flux0[m], wL[m];
	for(int iv = 0; iv < m; iv++) {
	  wL[iv] = 0;
	}
	fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux0);

	/* for(int iv2 = 0; iv2 < m; iv2++) { */
	/*   int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2); */
	/*   real val = theta *dt * flux0[iv2] * wpg; */
	/*   solver->rhs[imem2] -= val; */
	/* } */
      
	for(int iv1 = 0; iv1 < m; iv1++) {
	  int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv1) + offsetL;

	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = (iv == iv1);
	  }

	  schnaps_real flux[m];
	  fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux);

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    // The basis functions is also the gauss point index
	    int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2) + offsetL;
	    schnaps_real val = theta *dt * (flux[iv2]-flux0[iv2]) * wpg;		    
	    AddLinearSolver(solver, imem2, imem1, val);
	  }
	} // iv1

      } // else

  
    } // ipgfl

  } // macroface loop

}

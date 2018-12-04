#include "lbm_timeschemes.h"
#include "implicit.h"
/***********************************************************************************/
/*** small utility to recover multi-index from memory location**********************/
/**** WARNING need sanity checks to make varindex independent *********************/
/**** works only for varindex = iv + m * ipg for now ******************************/
void LBM_imem_to_glopmulti(Simulation * simu, int imem,
                           LBM_glop_multindex * glop)
{
  int fsize = simu->wsize / simu->macromesh.nbelems; // FIXME -> this should be precomputed in simu 
  // in an array fsize[nbelems], plus an array foffset[nb_elems+1]  
  // in the case where fsize is different for each macroelement. 
  // NB
  // -> the offset_el backward/forward computations should be modified accordingly (with a small bisection search
  // for eloffset)
  // -> the  interface routines should tackle the match between  macro elements with different RAF,DEG (and  even m for 
  // space dependent modelling) 
  glop->ie = imem / fsize;
  glop->ivarindex = imem % fsize;
  field *f = simu->fd + glop->ie;
  f->varindex_to_ipgiv(f->deg, f->raf, f->model.m, glop->ivarindex,
                       &(glop->ipg), &(glop->iv));
  ipg_to_xyz(f->raf, f->deg, glop->ic, glop->ix, &(glop->ipg));
}
// replace with void subroutine that compute all memoffsets from ie,ic,ix multiindex
int LBM_glopmuti_to_imem(Simulation * simu,LBM_glop_multindex *glop){
  int fsize = simu->wsize / simu->macromesh.nbelems; // FIXME -> this should be precomputed in simu 
  int offset_el=fsize * glop->ie;
  field *f=simu->fd+glop->ie;
  xyz_to_ipg(f->raf,f->deg,glop->ic,glop->ix,&(glop->ipg));
  glop->ivarindex= f->varindex(f->deg,f->raf,f->model.m,glop->ipg,glop->iv);
  int imem= offset_el+glop->ivarindex;
  return imem;
}
int CheckLBM_imem_multivarindex_constistency(Simulation *simu){
    // loop on macroelements
    LBM_glop_multindex iglop;
    for (iglop.ie = 0; iglop.ie < simu->macromesh.nbelems; iglop.ie++) {
      field *f =simu->fd+iglop.ie;
      for (iglop.ic[0] = 0; iglop.ic[0] < f->raf[0]; iglop.ic[0]++) {
      for (iglop.ic[1] = 0; iglop.ic[1] < f->raf[1]; iglop.ic[1]++) {
      for (iglop.ic[2] = 0; iglop.ic[2] < f->raf[2]; iglop.ic[2]++) {
        for (iglop.ix[0] = 0; iglop.ix[0] < f->deg[0]+1; iglop.ix[0]++) {
        for (iglop.ix[1] = 0; iglop.ix[1] < f->deg[1]+1; iglop.ix[1]++) {
        for (iglop.ix[2] = 0; iglop.ix[2] < f->deg[2]+1; iglop.ix[2]++) {
        for (iglop.iv=0;iglop.iv < f->model.m;iglop.iv++){
          int imem=LBM_glopmuti_to_imem(simu,&iglop);
          LBM_glop_multindex jglop;
          LBM_imem_to_glopmulti(simu,imem,&jglop);
          assert(iglop.ie==jglop.ie);
          assert(iglop.ivarindex==jglop.ivarindex);
          assert(iglop.ipg==jglop.ipg);
          for (int k=0;k<3;k++){
            assert(iglop.ic[k]==jglop.ic[k]);
            assert(iglop.ix[k]==jglop.ix[k]);
          }
        }
        }
        }
        }
      }
      }
      }
    }
    return 1;
}
/***********************************************************************************/
// Time schemes
/*************************************************************************************/
// *** adapted routines for the splitted lbm thetatimescheme
// (advection is done sequentially for each velocity node
// relaxation is applied separately  
// TODO give some name to that scheme to properly identify and avoid 
// confusion with the Nonlinear time scheme
void LBMThetaTimeScheme(LBMSimulation * lbsimu, schnaps_real theta,
                        schnaps_real tmax, schnaps_real dt)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int nb_nodes = lsd->lb_model->q;      // velocity nodes on the lattice(s)
  int ig_glob = 0, ig_node = 0;
  field *f_glob, *f_node;
  Simulation *micsimu = &(lbsimu->micro_simu);
  Simulation *macsimu = &(lbsimu->macro_simu);
  //
  MatrixStorage ms = SKYLINE;
  //
  int nraf[3] = { micsimu->fd[0].raf[0],
    micsimu->fd[0].raf[1],
    micsimu->fd[0].raf[2]
  };
  int deg[3] = { micsimu->fd[0].deg[0],
    micsimu->fd[0].deg[1],
    micsimu->fd[0].deg[2]
  };
  int nb_ipg = NPG(deg, nraf);
  //
  Simulation *simu_advec;
  simu_advec = malloc(sizeof(Simulation));
  EmptySimulation(simu_advec);
  InitSimulation(simu_advec, &(micsimu->macromesh), deg, nraf,
                 &(lbsimu->model_advec));
  //
  simu_advec->vmax = lbsimu->vmax;
  simu_advec->cfl = 0.0;
  simu_advec->nb_diags = 0;
  simu_advec->pre_dtfields = NULL;
  simu_advec->post_dtfields = NULL;
  simu_advec->update_after_rk = NULL;
  //
  LinearSolver solver_imp[nb_nodes];
  LinearSolver solver_exp[nb_nodes];
  schnaps_real *res_advec =
      calloc(simu_advec->wsize, sizeof(schnaps_real *));
  //
  lbsimu->dt = dt;
  int itermax = (int) (tmax / dt) + 1;
  printf("Called with tmax=%f dt=%f Nb iterations:%i \n", tmax, dt,
         itermax);
  lbsimu->itermax = itermax;
  // important sync micro simu params for diag compatibility with RK which uses only the micro simu);
  lbsimu->micro_simu.itermax_rk = itermax;
  lbsimu->micro_simu.dt = dt;
  int freq = (1 >= lbsimu->itermax / 10) ? 1 : lbsimu->itermax / 10;
  lbsimu->tnow = 0.0;
  simu_advec->dt = dt;
  simu_advec->itermax_rk = itermax;
  simu_advec->tnow = 0.0;
  for (int ie = 0; ie < micsimu->macromesh.nbelems; ++ie) {
    micsimu->fd[ie].tnow = lbsimu->tnow;
    macsimu->fd[ie].tnow = lbsimu->tnow;
    simu_advec->fd[ie].tnow = simu_advec->tnow;
  }
  // Diagnostics (this should be elsewhere, some timetraces  module ?
  int mac_size_diags = macsimu->nb_diags * itermax;
  int mic_size_diags = micsimu->nb_diags * itermax;
  if (macsimu->nb_diags != 0) {
    macsimu->Diagnostics = malloc(mac_size_diags * sizeof(schnaps_real));
  };
  if (micsimu->nb_diags != 0) {
    micsimu->Diagnostics = malloc(mic_size_diags * sizeof(schnaps_real));
  };
  //
  time_t t_start, t_end;        // time measurements for op factorization
  t_start = time(NULL);
  //  Solvers Init
  printf("Sparse Linear Solvers init");
  for (int isim = 0; isim < nb_nodes; isim++) {
    lsd->current_node_index = isim;
    //
    InitImplicitLinearSolver(simu_advec, &solver_imp[isim], ms);
    InitImplicitLinearSolver(simu_advec, &solver_exp[isim], ms);
    //
  }
  // End Operators init
  //
  // Time loop start
  for (int iter = 0; iter < itermax; iter++) {
    //
    if (iter % freq == 0) {
      printf(" iter %i/%i t=%f\n", iter, itermax, lbsimu->tnow);
    }
    //
    lbsimu->iter_time = iter;
    macsimu->iter_time_rk = iter;
    micsimu->iter_time_rk = iter;
    //
    if (iter == 0) {
      t_start = time(NULL);
    }
    //
    if (lbsimu->pre_advec != NULL) {
      lbsimu->pre_advec(lbsimu);
    };
    // now loop on velocity nodes
    for (int isim = 0; isim < nb_nodes; isim++) {
      lsd->current_node_index = isim;
      // dispatch main w to per node w's
      for (int ie = 0; ie < simu_advec->macromesh.nbelems; ie++) {
        f_glob = micsimu->fd + ie;
        for (int ipg = 0; ipg < nb_ipg; ipg++) {
          f_node = simu_advec->fd + ie;
          ig_glob =
              f_glob->varindex(f_glob->deg, f_glob->raf, f_glob->model.m,
                               ipg, isim);
          ig_node = f_node->varindex(f_node->deg, f_node->raf, 1, ipg, 0);
          f_node->wn[ig_node] = f_glob->wn[ig_glob];
          //
        };                      // ipg end loop glops
      };                        // ie end loop macroelements
      // end of data dispatch
      // Solvers assembly and factorization if necessary
      solver_imp[isim].rhs_is_assembly = false;
      solver_exp[isim].rhs_is_assembly = false;
      if (iter == 0) {
        solver_imp[isim].mat_is_assembly = false;
        solver_exp[isim].mat_is_assembly = false;
      } else {
        solver_imp[isim].mat_is_assembly = true;
        solver_exp[isim].mat_is_assembly = true;
      };
      //
      simu_advec->tnow = lbsimu->tnow;
      for (int ie = 0; ie < simu_advec->macromesh.nbelems; ++ie) {
        simu_advec->fd[ie].tnow = simu_advec->tnow;
      }
      LBM_AssemblyImplicitLinearSolver(simu_advec, &solver_exp[isim],
                                       -(1.0 - theta), simu_advec->dt);
      simu_advec->tnow = lbsimu->tnow + lbsimu->dt;
      for (int ie = 0; ie < simu_advec->macromesh.nbelems; ++ie) {
        simu_advec->fd[ie].tnow = simu_advec->tnow;
      }
      LBM_AssemblyImplicitLinearSolver(simu_advec, &solver_imp[isim],
                                       theta, simu_advec->dt);
      // compute residual
      MatVect(&solver_exp[isim], simu_advec->w, res_advec);
      //
      for (int i = 0; i < solver_imp[isim].neq; i++) {
        solver_imp[isim].rhs[i] =
            -solver_exp[isim].rhs[i] + solver_imp[isim].rhs[i] +
            res_advec[i];
      }
      //
      solver_imp[isim].solver_type = LU;
      Advanced_SolveLinearSolver(&solver_imp[isim], simu_advec);
      //
      for (int i = 0; i < solver_imp[isim].neq; i++) {
        simu_advec->w[i] = solver_imp[isim].sol[i];
      }
      // collect nodes ws to main w
      for (int ie = 0; ie < micsimu->macromesh.nbelems; ie++) {
        f_glob = micsimu->fd + ie;
        for (int ipg = 0; ipg < nb_ipg; ipg++) {
          f_node = simu_advec->fd + ie;
          ig_glob =
              f_glob->varindex(f_glob->deg, f_glob->raf, f_glob->model.m,
                               ipg, isim);
          ig_node = f_node->varindex(f_node->deg, f_node->raf, 1, ipg, 0);
          f_glob->wn[ig_glob] = f_node->wn[ig_node];
        };                      // ipg end loop glops
      };                        // ie end loop macroelements
      //
      if (lbsimu->post_advec_one_node) {
        lbsimu->post_advec_one_node(lbsimu);
      }
      //
    };                          // isim end loop on velocity node 
    // post advec ops
    if (lbsimu->post_advec != NULL) {
      lbsimu->post_advec(lbsimu);
    }
    if (lbsimu->post_tstep != NULL) {
      lbsimu->post_tstep(lbsimu, macsimu->w);
    }
    if (iter == 0) {
      t_end = time(NULL);
      printf("First step duration %ld\n", t_end - t_start);
    }
    lbsimu->tnow += lbsimu->dt;
    simu_advec->tnow = lbsimu->tnow;
    micsimu->tnow = lbsimu->tnow;
    macsimu->tnow = lbsimu->tnow;
    for (int ie = 0; ie < simu_advec->macromesh.nbelems; ++ie) {
      micsimu->fd[ie].tnow = lbsimu->tnow;
      macsimu->fd[ie].tnow = lbsimu->tnow;
      simu_advec->fd[ie].tnow = lbsimu->tnow;
    }
    //
  };                            // iter end of time loop 
  //
  if (res_advec != NULL) {
    free(res_advec);
  }
  if (simu_advec != NULL) {
    free(simu_advec);
  }
}
/***************************************************************************************************************/
void LBM_MassCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky){
  bool verbose = false;
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  //
  for (int i=0;i< simu->wsize;i++){
    IsNonZero(solver,i,i);
  }
}
/***************************************************************************************************************/
void LBM_InternalCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky){
  bool verbose = false;
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  //
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  //
  LBM_glop_multindex Iglop;
  LBM_glop_multindex Jglop;
  // restrict to advected variable iv_advec
  Iglop.iv = 0;
  Jglop.iv = 0;
  // loop on macro elements
  for(Iglop.ie = 0; Iglop.ie < simu->macromesh.nbelems; Iglop.ie++){
    field *f = simu->fd + Iglop.ie; // &(simu->fd[ie])
    Jglop.ie = Iglop.ie; // sync J element
    // Loop on the subcells
    for (Iglop.ic[0]=0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
      Jglop.ic[0]=Iglop.ic[0]; // sync cell ic0
    for (Iglop.ic[1]=0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
      Jglop.ic[1]=Iglop.ic[1]; // sync cell ic1
    for (Iglop.ic[2]=0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
      Jglop.ic[2]=Iglop.ic[2]; // sync cell ic2
    // Loop on glops in the subcell
    for (Iglop.ix[0] = 0; Iglop.ix[0] < f->deg[0]+1; Iglop.ix[0]++) {
    for (Iglop.ix[1] = 0; Iglop.ix[1] < f->deg[1]+1; Iglop.ix[1]++) {
    for (Iglop.ix[2] = 0; Iglop.ix[2] < f->deg[2]+1; Iglop.ix[2]++) {
    // recover memory index for Iglop;
    int imem = LBM_glopmuti_to_imem(simu,&Iglop);
    // loop on grapsi component
      for(int dim0 = 0; dim0 < 3; dim0++) {
          // sync Jglop.ix with Iglop.ix
          for (int k=0;k<3;k++){
            Jglop.ix[k]=Iglop.ix[k];
          }
          // add contribution only along dim0 due to separability in xyz of basis functions
          for (Jglop.ix[dim0] = 0; Jglop.ix[dim0] < f->deg[dim0]+1; Jglop.ix[dim0]++){
            int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
            if (verbose){
            printf("Int Switching on (%i,%i) in cell (%i,%i,%i) iv %i\n",
              imem,jmem,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],Iglop.iv);
            }
            IsNonZero(solver,imem,jmem);
          } // end loop Jglop.ix[dim0]
          //
      } // dim0 loop
    } // end Iglop.ix[0]
    } // end Iglop.ix[1]
    } // end Iglop.ix[2]
    // loop along grapsi components (x,y,z)
    } // end Iglop.ic[2]
    } // end Iglop.ic[1]
    } // end Iglop.ic[0]
  } // end Iglop.ie
} 
/***************************************************************************************************************/
void LBM_FluxCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky){
  //
  const bool verbose = false;
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  if (verbose){
    printf(" velocity:(%f,%f,%f) \n:",vel[0],vel[1],vel[2]);
  }
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  LBM_glop_multindex Iglop;
  LBM_glop_multindex Jglop;
  // restrict iv to the (unique) field id 
  Iglop.iv = 0;
  Jglop.iv = 0;
  // loop on macro elements
  for (Iglop.ie=0; Iglop.ie < simu->macromesh.nbelems;Iglop.ie++){
    field *f= simu->fd+Iglop.ie;
  // loop on subcells
  for (Iglop.ic[0]=0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
  for (Iglop.ic[1]=0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
  for (Iglop.ic[2]=0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
    if (verbose){
    printf(" ic:(%i;%i,%i) \n",Iglop.ic[0],Iglop.ic[1],Iglop.ic[2]);
    }
  // Sweeping subcell faces pairs in the three directions (ie two subfaces per dim)  
  for (int dim0 = 0; dim0 < 3; dim0++) {
    // do not treat last cell as a "left" cell
    if (Iglop.ic[dim0] == f->raf[dim0] - 1) {
      continue;
    }
    //
    const int altdim1[3] = { 1, 0, 0 };
    const int altdim2[3] = { 2, 2, 1 };
    int dim1= altdim1[dim0];
    int dim2= altdim2[dim0];
    // select right subcell along dim0
    Jglop.ie=Iglop.ie;
    Jglop.ic[dim0]=Iglop.ic[dim0]+1;
    Jglop.ic[dim1]=Iglop.ic[dim1];
    Jglop.ic[dim2]=Iglop.ic[dim2];
    // select rightmost glop in dim0 direction from the left subcell
    Iglop.ix[dim0] = f->deg[dim0];
    // select leftmost glop in dim0 direction from the right subcell
    Jglop.ix[dim0] = 0;
    // loop on left glop on the face
    for (Iglop.ix[dim1] = 0; Iglop.ix[dim1] < f->deg[dim1]+1; Iglop.ix[dim1]++) {
      // matching glop ix[idim1] in the right subcell
      Jglop.ix[dim1] = Iglop.ix[dim1];
      for (Iglop.ix[dim2] = 0; Iglop.ix[dim2] < f->deg[dim2]+1; Iglop.ix[dim2]++) {
        // matching glop ix[idim2] in the right subcell
        Jglop.ix[dim2] = Iglop.ix[dim2];
        // get memory index for the advected variable iv_advec
        int imem = LBM_glopmuti_to_imem(simu,&Iglop);
        int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
        // 
        // Compute the normal vector 
        schnaps_real vnds[3];
        {
          schnaps_real xref[3], wpg3;
          ref_pg_vol(f->deg, f->raf, Iglop.ipg, xref, &wpg3, NULL);
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
        }
        // compute the scalar product with advection velocity
        schnaps_real vidotn = vnds[0] * vel[0] + vnds[1] * vel[1] +vnds[2] * vel[2];
        //
        if (vidotn > _SMALL){
          if (verbose){
          printf("Flux Switching on (%i,%i) in cell (%i,%i,%i) iv %i\n",imem,jmem,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],Iglop.iv);
          }
          IsNonZero(solver, imem, jmem);
        }
        if (vidotn < -_SMALL){
          if (verbose){
          printf("Flux Switching on (%i,%i) in cell (%i,%i,%i) iv %i\n",jmem,imem,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],Iglop.iv);
          }
          IsNonZero(solver, jmem, imem);
        }
    } // end Iglop.ix[dim2]
    } // end Iglop.ix[dim1]
    //
  } // end dim0
  } // end ic2
  } // end ic1
  } // end ic0
  } // end ie
}
/***************************************************************************************************************/
void LBM_InterfaceCoupling_OneNode(Simulation *simu, LinearSolver *solver, int isky)
{
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    const unsigned int m = fL->model.m;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      offsetR = fsize * ieR;
    }
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
      //
      schnaps_real vidotn= vnds[0] * vel[0] + vnds[1] * vel[1] + vnds[2] * vel[2];
      //
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
	int imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, 0) + offsetL;
	int imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, 0) + offsetR;
	if (vidotn > _SMALL){
	  IsNonZero(solver, imem1, imem2);
	}
	else{
	if (vidotn < -_SMALL){
	  IsNonZero(solver, imem2, imem1);
	}
	}
	//
   } 
   else { // The point is on the boundary.
     // only the advected node varies
	   int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, 0) + offsetL;
	   IsNonZero(solver, imem1, imem1);
    } // end else end if
    } // ipgfl
  } // macroface loop
}
/***************************************************************************************************************/
void LBM_InitLinearSolver_OneNode(Simulation *simu, LinearSolver *solver, MatrixStorage ms){
  int neq = simu->wsize;
  //
  Solver st = LU;
  InitLinearSolver(solver,neq,&ms,&st);
  //solver->tol=1.e-8;
  solver->tol=_SMALL;
  int isky=0;
  LBM_MassCoupling_OneNode(simu,solver,isky);
  LBM_InternalCoupling_OneNode(simu,solver,isky);
  LBM_FluxCoupling_OneNode(simu,solver,isky);
  LBM_InterfaceCoupling_OneNode(simu,solver,isky);
  AllocateLinearSolver(solver);
  DisplayLinearSolver(solver);
}
/***************************************************************************************************************/
/* assembly routines for one node velocity advections */
void LBM_FluxAssembly_OneNode(Simulation * simu, LinearSolver * solver,schnaps_real theta, schnaps_real dt){
  const bool verbose = false;
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  LBM_glop_multindex Iglop;
  LBM_glop_multindex Jglop;
  // restrict iv to the (unique) field id 
  Iglop.iv = 0;
  Jglop.iv = 0;
  // loop on macro elements
  for (Iglop.ie=0; Iglop.ie < simu->macromesh.nbelems;Iglop.ie++){
    field *f= simu->fd+Iglop.ie;
  // loop on subcells
  for (Iglop.ic[0]=0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
  for (Iglop.ic[1]=0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
  for (Iglop.ic[2]=0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
  // Sweeping subcell faces pairs in the three directions (ie two subfaces per dim)  
  for (int dim0 = 0; dim0 < 3; dim0++) {
    // do not treat last cell as a "left" cell
    if (Iglop.ic[dim0] == f->raf[dim0] - 1) {
      continue;
    }
    //
    const int altdim1[3] = { 1, 0, 0 };
    const int altdim2[3] = { 2, 2, 1 };
    int dim1= altdim1[dim0];
    int dim2= altdim2[dim0];
    // select right subcell along dim0
    Jglop.ie=Iglop.ie;
    Jglop.ic[dim0]=Iglop.ic[dim0]+1;
    Jglop.ic[dim1]=Iglop.ic[dim1];
    Jglop.ic[dim2]=Iglop.ic[dim2];
    // select rightmost glop in dim0 direction from the left subcell
    Iglop.ix[dim0] = f->deg[dim0];
    // select leftmost glop in dim0 direction from the right subcell
    Jglop.ix[dim0] = 0;
    // loop on left glop on the face
    for (Iglop.ix[dim1] = 0; Iglop.ix[dim1] < f->deg[dim1]+1; Iglop.ix[dim1]++) {
      // matching glop ix[idim1] in the right subcell
      Jglop.ix[dim1] = Iglop.ix[dim1];
      for (Iglop.ix[dim2] = 0; Iglop.ix[dim2] < f->deg[dim2]+1; Iglop.ix[dim2]++) {
        // matching glop ix[idim2] in the right subcell
        Jglop.ix[dim2] = Iglop.ix[dim2];
        // get memory index for the advected variable iv_advec
        int imem = LBM_glopmuti_to_imem(simu,&Iglop);
        int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
        // 
        // Compute the normal vector 
        schnaps_real vnds[3];
        {
          schnaps_real xref[3], wpg3;
          ref_pg_vol(f->deg, f->raf, Iglop.ipg, xref, &wpg3, NULL);
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
          schnaps_real h1h2 = 1. / f->raf[dim1] / f->raf[dim2];
          vnds[0] = codtau[0][dim0] * h1h2;
          vnds[1] = codtau[1][dim0] * h1h2;
          vnds[2] = codtau[2][dim0] * h1h2;
        }
        //  WIP COMPUTE FLUX
        schnaps_real wpg
          = wglop(f->deg[dim1], Iglop.ix[dim1])
          * wglop(f->deg[dim2], Iglop.ix[dim2]);
		  schnaps_real wL[f->model.m], wR[f->model.m], flux[f->model.m];
		  //
		  //

    } // end Iglop.ix[dim2]
    } // end Iglop.ix[dim1]
    //
  } // end dim0
  } // end ic2
  } // end ic1
  } // end ic0
  } // end ie

}

/***************************************************************************************************************/
//!\brief Multistep time scheme based on splitted per-node advection and CN relaxation 
void LBMThetaTimeScheme_OneNode_Multistep(LBMSimulation * lbsimu, schnaps_real theta,
			schnaps_real tmax, schnaps_real dt){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int nb_nodes = lsd->lb_model->q;      // velocity nodes on the lattice(s)
  int ig_glob = 0, ig_node = 0;
  field *f_glob, *f_node;
  Simulation *micsimu = &(lbsimu->micro_simu);
  Simulation *macsimu = &(lbsimu->macro_simu);
  //
  MatrixStorage ms = KLU_CSR;
  //
  int nraf[3] = { micsimu->fd[0].raf[0],
    micsimu->fd[0].raf[1],
    micsimu->fd[0].raf[2]
  };
  int deg[3] = { micsimu->fd[0].deg[0],
    micsimu->fd[0].deg[1],
    micsimu->fd[0].deg[2]
  };
  int nb_ipg = NPG(deg, nraf);
  //
  Simulation *simu_advec;
  simu_advec = malloc(sizeof(Simulation));
  EmptySimulation(simu_advec);
  InitSimulation(simu_advec, &(micsimu->macromesh), deg, nraf,
                 &(lbsimu->model_advec));
  //
  simu_advec->vmax = lbsimu->vmax;
  simu_advec->cfl = 0.0;
  simu_advec->nb_diags = 0;
  simu_advec->pre_dtfields = NULL;
  simu_advec->post_dtfields = NULL;
  simu_advec->update_after_rk = NULL;
  //
  LinearSolver solver_imp[nb_nodes];
  LinearSolver solver_exp[nb_nodes];
  schnaps_real *res_advec =
      calloc(simu_advec->wsize, sizeof(schnaps_real *));
  //
  lbsimu->dt = dt;
  int itermax = (int) (tmax / dt) + 1;
  printf("Called with tmax=%f dt=%f Nb iterations:%i \n", tmax, dt,
         itermax);
  lbsimu->itermax = itermax;
  // important sync micro simu params for diag compatibility with RK which uses only the micro simu);
  lbsimu->micro_simu.itermax_rk = itermax;
  lbsimu->micro_simu.dt = dt;
  int freq = (1 >= lbsimu->itermax / 10) ? 1 : lbsimu->itermax / 10;
  lbsimu->tnow = 0.0;
  simu_advec->dt = dt;
  simu_advec->itermax_rk = itermax;
  simu_advec->tnow = 0.0;
  for (int ie = 0; ie < micsimu->macromesh.nbelems; ++ie) {
    micsimu->fd[ie].tnow = lbsimu->tnow;
    macsimu->fd[ie].tnow = lbsimu->tnow;
    simu_advec->fd[ie].tnow = simu_advec->tnow;
  }
  // Diagnostics (this should be elsewhere, some timetraces  module ?
  int mac_size_diags = macsimu->nb_diags * itermax;
  int mic_size_diags = micsimu->nb_diags * itermax;
  if (macsimu->nb_diags != 0) {
    macsimu->Diagnostics = malloc(mac_size_diags * sizeof(schnaps_real));
  };
  if (micsimu->nb_diags != 0) {
    micsimu->Diagnostics = malloc(mic_size_diags * sizeof(schnaps_real));
  };
  //
  time_t t_start, t_end;        // time measurements for op factorization
  t_start = time(NULL);
  //  Solvers Init
  printf("Sparse Linear Solvers init");
  for (int isim = 0; isim < nb_nodes; isim++) {
    lsd->current_node_index = isim;
    //
    LBM_InitLinearSolver_OneNode(simu_advec, &solver_imp[isim], ms);
    LBM_InitLinearSolver_OneNode(simu_advec, &solver_exp[isim], ms);
    //
  }
  //
  // DEBUG
  printf(" DEBUG : no time iterations");
  itermax =0;
  // End Operators init
  // Time loop start
  for (int iter = 0; iter < itermax; iter++) {
    //
    if (iter % freq == 0) {
      printf(" iter %i/%i t=%f\n", iter, itermax, lbsimu->tnow);
    }
    //
    lbsimu->iter_time = iter;
    macsimu->iter_time_rk = iter;
    micsimu->iter_time_rk = iter;
    //
    if (iter == 0) {
      t_start = time(NULL);
    }
    //
    if (lbsimu->pre_advec != NULL) {
      lbsimu->pre_advec(lbsimu);
    };
    // ADVECTION BLOCK START 
    // TODO export in function
    // now loop on velocity nodes
    for (int isim = 0; isim < nb_nodes; isim++) {
      lsd->current_node_index = isim;
      // dispatch main w to per node w's
      for (int ie = 0; ie < simu_advec->macromesh.nbelems; ie++) {
        f_glob = micsimu->fd + ie;
        for (int ipg = 0; ipg < nb_ipg; ipg++) {
          f_node = simu_advec->fd + ie;
          ig_glob =
              f_glob->varindex(f_glob->deg, f_glob->raf, f_glob->model.m,
                               ipg, isim);
          ig_node = f_node->varindex(f_node->deg, f_node->raf, 1, ipg, 0);
          f_node->wn[ig_node] = f_glob->wn[ig_glob];
          //
        };                      // ipg end loop glops
      };                        // ie end loop macroelements
      // end of data dispatch
      // Solvers assembly and factorization if necessary
      solver_imp[isim].rhs_is_assembly = false;
      solver_exp[isim].rhs_is_assembly = false;
      if (iter == 0) {
        solver_imp[isim].mat_is_assembly = false;
        solver_exp[isim].mat_is_assembly = false;
      } else {
        solver_imp[isim].mat_is_assembly = true;
        solver_exp[isim].mat_is_assembly = true;
      };
      //
      simu_advec->tnow = lbsimu->tnow;
      for (int ie = 0; ie < simu_advec->macromesh.nbelems; ++ie) {
        simu_advec->fd[ie].tnow = simu_advec->tnow;
      }
      LBM_AssemblyImplicitLinearSolver(simu_advec, &solver_exp[isim],
                                       -(1.0 - theta), simu_advec->dt);
      simu_advec->tnow = lbsimu->tnow + lbsimu->dt;
      for (int ie = 0; ie < simu_advec->macromesh.nbelems; ++ie) {
        simu_advec->fd[ie].tnow = simu_advec->tnow;
      }
      LBM_AssemblyImplicitLinearSolver(simu_advec, &solver_imp[isim],
                                       theta, simu_advec->dt);
      // compute residual
      MatVect(&solver_exp[isim], simu_advec->w, res_advec);
      //
      for (int i = 0; i < solver_imp[isim].neq; i++) {
        solver_imp[isim].rhs[i] =
            -solver_exp[isim].rhs[i] + solver_imp[isim].rhs[i] +
            res_advec[i];
      }
      //
      solver_imp[isim].solver_type = LU;
      Advanced_SolveLinearSolver(&solver_imp[isim], simu_advec);
      //
      for (int i = 0; i < solver_imp[isim].neq; i++) {
        simu_advec->w[i] = solver_imp[isim].sol[i];
      }
      // collect nodes ws to main w
      for (int ie = 0; ie < micsimu->macromesh.nbelems; ie++) {
        f_glob = micsimu->fd + ie;
        for (int ipg = 0; ipg < nb_ipg; ipg++) {
          f_node = simu_advec->fd + ie;
          ig_glob =
              f_glob->varindex(f_glob->deg, f_glob->raf, f_glob->model.m,
                               ipg, isim);
          ig_node = f_node->varindex(f_node->deg, f_node->raf, 1, ipg, 0);
          f_glob->wn[ig_glob] = f_node->wn[ig_node];
        };                      // ipg end loop glops
      };                        // ie end loop macroelements
      //
      if (lbsimu->post_advec_one_node) {
        lbsimu->post_advec_one_node(lbsimu);
      }
      //
    };                          // isim end loop on velocity node 
    // ADVECTION BLOCK END
    // post advec ops
    if (lbsimu->post_advec != NULL) {
      lbsimu->post_advec(lbsimu);
    }
    if (lbsimu->post_tstep != NULL) {
      lbsimu->post_tstep(lbsimu, macsimu->w);
    }
    if (iter == 0) {
      t_end = time(NULL);
      printf("First step duration %ld\n", t_end - t_start);
    }
    lbsimu->tnow += lbsimu->dt;
    simu_advec->tnow = lbsimu->tnow;
    micsimu->tnow = lbsimu->tnow;
    macsimu->tnow = lbsimu->tnow;
    for (int ie = 0; ie < simu_advec->macromesh.nbelems; ++ie) {
      micsimu->fd[ie].tnow = lbsimu->tnow;
      macsimu->fd[ie].tnow = lbsimu->tnow;
      simu_advec->fd[ie].tnow = lbsimu->tnow;
    }
    //
  };                            // iter end of time loop 
  //
  if (res_advec != NULL) {
    free(res_advec);
  }
  if (simu_advec != NULL) {
    free(simu_advec);
  }

  
}

/***************************************************************************************************************/
/***************************************************************************************************************/
void LBM_AssemblyImplicitLinearSolver(Simulation * simu,
                                      LinearSolver * solver,
                                      schnaps_real theta, schnaps_real dt)
{
  if (solver->mat_is_assembly == false) {
    MassAssembly(simu, solver);
    InternalAssembly(simu, solver, theta, dt);
    FluxAssembly(simu, solver, theta, dt);
    LBM_InterfaceAssembly(simu, solver, theta, dt);
  }
  if (solver->rhs_is_assembly == false) {
    for (int i = 0; i < solver->neq; i++) {
      solver->rhs[i] = 0;
    }
    LBM_SourceAssembly(simu, solver, theta, dt);
  }
  //DisplayLinearSolver(solver);
}

void LBM_SourceAssembly(Simulation * simu, LinearSolver * solver,
                        schnaps_real theta, schnaps_real dt)
{

  if (simu->fd[0].model.Source != NULL) {
    for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
      field *f = simu->fd + ie;
      int offsetw = f->wsize * ie;

      const int m = f->model.m;

      int deg[3] = { f->deg[0],
        f->deg[1],
        f->deg[2]
      };

      int nraf[3] = { f->raf[0],
        f->raf[1],
        f->raf[2]
      };

      for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
        schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
        ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
        schnaps_ref2phy(f->physnode,    // phys. nodes
                        xpgref, // xref
                        NULL, -1,       // dpsiref, ifa
                        xphy, dtau,     // xphy, dtau
                        codtau, NULL, NULL);    // codtau, dpsi, vnds
        schnaps_real det = dot_product(dtau[0], codtau[0]);
        schnaps_real wL[m], source[m];
        /* for(int iv = 0; iv < m; ++iv){ */
        /*      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv); */
        /*      wL[iv] = w[imem]; */
        /* } */
        f->model.Source(xphy, f->tnow, wL, source);

        for (int iv1 = 0; iv1 < m; iv1++) {
          int imem = f->varindex(deg, nraf, m, ipg, iv1) + offsetw;
          schnaps_real val = theta * dt * source[iv1] * wpg * det;
          solver->rhs[imem] += val;
        }
      }                         //ipg
    }
  }
  // assembly of the boundary terms

  int fsize = simu->wsize / simu->macromesh.nbelems;

  for (int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++) {
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR < 0) {

      const unsigned int m = fL->model.m;

      // Loop over the points on a single macro cell interface.
      for (int ipgfL = 0; ipgfL < NPGF(fL->deg, fL->raf, locfaL); ipgfL++) {

        schnaps_real xpgref[3], xpgref_in[3], wpg;

        // Get the coordinates of the Gauss point and coordinates of a
        // point slightly inside the opposite element in xref_in
        int ipgL =
            ref_pg_face(fL->deg, fL->raf, locfaL, ipgfL, xpgref, &wpg,
                        xpgref_in);

        // Normal vector at gauss point ipgL
        schnaps_real vnds[3], xpg[3];
        {
          schnaps_real dtau[3][3], codtau[3][3];
          schnaps_ref2phy(fL->physnode, xpgref, NULL, locfaL,   // dpsiref, ifa
                          xpg, dtau, codtau, NULL, vnds);       // codtau, dpsi, vnds
        }

        // the boundary flux is an affine function
        schnaps_real flux0[m], wL[m];
        for (int iv = 0; iv < m; iv++) {
          wL[iv] = 0;
        }
        // store the state of micro/macro simulation in global struct to allow access across nodes from 1node boundary flux
        LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
        LBMSimulation *lbsimu = lsd->current_lb_sim;
        for (int i = 0; i < lsd->lb_model->q; i++) {
          field *fmic = lbsimu->micro_simu.fd + ieL;
          int ivar =
              fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipgL, i);
          lbsimu->wmic_buffer[i] = fmic->wn[ivar];
        }
        for (int i = 0; i < lsd->lb_model->nb_macro; i++) {
          field *fmac = lbsimu->macro_simu.fd + ieL;
          int ivar =
              fmac->varindex(fmac->deg, fmac->raf, fmac->model.m, ipgL, i);
          lbsimu->wmac_buffer[i] = fmac->wn[ivar];
        }
        //
        fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux0);

        for (int iv2 = 0; iv2 < m; iv2++) {
          int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL,
                                   iv2) + offsetL;
          schnaps_real val = theta * dt * flux0[iv2] * wpg;
          solver->rhs[imem2] -= val;
        }
      }
    }                           // if ier < 0
  }                             // macroface loop
}                               // SourceAssembly

//
void LBM_InterfaceAssembly(Simulation * simu, LinearSolver * solver,
                           schnaps_real theta, schnaps_real dt)
{

  int fsize = simu->wsize / simu->macromesh.nbelems;

  for (int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++) {
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
    for (int ipgfL = 0; ipgfL < NPGF(fL->deg, fL->raf, locfaL); ipgfL++) {

      schnaps_real xpgref[3], xpgref_in[3], wpg;

      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipgL = ref_pg_face(fL->deg, fL->raf, locfaL, ipgfL, xpgref, &wpg,
                             xpgref_in);

      // Normal vector at gauss point ipgL
      schnaps_real vnds[3], xpg[3];
      {
        schnaps_real dtau[3][3], codtau[3][3];
        schnaps_ref2phy(fL->physnode, xpgref, NULL, locfaL,     // dpsiref, ifa
                        xpg, dtau, codtau, NULL, vnds); // codtau, dpsi, vnds
      }

      if (fR != NULL) {         // the right element exists
        schnaps_real xrefL[3];
        {
          schnaps_real xpg_in[3];
          schnaps_ref2phy(fL->physnode, xpgref_in, NULL, -1,    // dpsiref, ifa
                          xpg_in, NULL, NULL, NULL, NULL);      // codtau, dpsi, vnds
          PeriodicCorrection(xpg_in, fL->period);
          schnaps_phy2ref(fR->physnode, xpg_in, xrefL);

        }

        int ipgR = ref_ipg(fR->deg, fR->raf, xrefL);
        schnaps_real flux[m];
        schnaps_real wL[m];
        schnaps_real wR[m];

        for (int iv1 = 0; iv1 < m; iv1++) {



          for (int iv = 0; iv < m; iv++) {
            wL[iv] = (iv == iv1);
            wR[iv] = 0;
          }
          int imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL,
                                   iv1) + offsetL;

          // int_dL F(wL, wR, grad phi_ib)

          fL->model.NumFlux(wL, wR, vnds, flux);

          // Add flux to both sides

          for (int iv2 = 0; iv2 < m; iv2++) {
            int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL,
                                     iv2) + offsetL;
            schnaps_real val = theta * dt * flux[iv2] * wpg;
            AddLinearSolver(solver, imem2, imem1, val);

            imem2 =
                fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR,
                             iv2) + offsetR;
            val = theta * dt * flux[iv2] * wpg;
            AddLinearSolver(solver, imem2, imem1, -val);
          }

          for (int iv = 0; iv < m; iv++) {
            wL[iv] = 0;
            wR[iv] = (iv == iv1);
          }
          imem1 =
              fL->varindex(fL->deg, fL->raf, fL->model.m, ipgR,
                           iv1) + offsetR;


          fL->model.NumFlux(wL, wR, vnds, flux);

          for (int iv2 = 0; iv2 < m; iv2++) {
            int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL,
                                     iv2) + offsetL;
            schnaps_real val = theta * dt * flux[iv2] * wpg;
            AddLinearSolver(solver, imem2, imem1, val);

            imem2 =
                fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR,
                             iv2) + offsetR;
            val = theta * dt * flux[iv2] * wpg;
            AddLinearSolver(solver, imem2, imem1, -val);
          }
        }

      } else {                  // The point is on the boundary.

        // the boundary flux is an affine function
        schnaps_real flux0[m], wL[m];
        for (int iv = 0; iv < m; iv++) {
          wL[iv] = 0;
        }
        // store the state of micro simulation in global struct to allow access across nodes from 1node boundary flux
        LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
        LBMSimulation *lbsimu = lsd->current_lb_sim;
        for (int i = 0; i < lsd->lb_model->q; i++) {
          field *fmic = lbsimu->micro_simu.fd + ieL;
          int ivar =
              fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipgL, i);
          lbsimu->wmic_buffer[i] = fmic->wn[ivar];
        }
        for (int i = 0; i < lsd->lb_model->nb_macro; i++) {
          field *fmac = lbsimu->macro_simu.fd + ieL;
          int ivar =
              fmac->varindex(fmac->deg, fmac->raf, fmac->model.m, ipgL, i);
          lbsimu->wmac_buffer[i] = fmac->wn[ivar];
        }
        //
        fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux0);

        /* for(int iv2 = 0; iv2 < m; iv2++) { */
        /*   int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2); */
        /*   real val = theta *dt * flux0[iv2] * wpg; */
        /*   solver->rhs[imem2] -= val; */
        /* } */

        for (int iv1 = 0; iv1 < m; iv1++) {
          int imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL,
                                   iv1) + offsetL;

          for (int iv = 0; iv < m; iv++) {
            wL[iv] = (iv == iv1);
          }
          //
          //
          schnaps_real flux[m];
          fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux);

          for (int iv2 = 0; iv2 < m; iv2++) {
            // The basis functions is also the gauss point index
            int imem2 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL,
                                     iv2) + offsetL;
            schnaps_real val = theta * dt * (flux[iv2] - flux0[iv2]) * wpg;
            AddLinearSolver(solver, imem2, imem1, val);
          }
        }                       // iv1

      }                         // else


    }                           // ipgfl

  }                             // macroface loop

}
//**********************************************************************************//
//*** coupling routines for nonlinear per_velocity avection with full relaxation
void LBM_Test_Coupling_BlockCell(Simulation *simu, LinearSolver *solver,int isky){
  LBM_glop_multindex Iglop;
  LBM_glop_multindex Jglop;
  for(Iglop.ie = 0; Iglop.ie < simu->macromesh.nbelems; Iglop.ie++){
    field *f = simu->fd + Iglop.ie; // &(simu->fd[ie])
    for (Iglop.ic[0] = 0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
    for (Iglop.ic[1] = 0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
    for (Iglop.ic[2] = 0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
    for (Iglop.ix[0] = 0; Iglop.ix[0] < f->deg[0]+1; Iglop.ix[0]++) {
    for (Iglop.ix[1] = 0; Iglop.ix[1] < f->deg[1]+1; Iglop.ix[1]++) {
    for (Iglop.ix[2] = 0; Iglop.ix[2] < f->deg[2]+1; Iglop.ix[2]++) {
    for (Iglop.iv = 0;Iglop.iv < f->model.m;Iglop.iv++){
        int imem= LBM_glopmuti_to_imem(simu,&Iglop);
        //
        Jglop.ie=Iglop.ie;
        for (int k=0;k<3;k++){
          Jglop.ic[k]=Iglop.ic[k];
        }
        for (Jglop.ix[0] = 0; Jglop.ix[0] < f->deg[0]+1; Jglop.ix[0]++) {
        for (Jglop.ix[1] = 0; Jglop.ix[1] < f->deg[1]+1; Jglop.ix[1]++) {
        for (Jglop.ix[2] = 0; Jglop.ix[2] < f->deg[2]+1; Jglop.ix[2]++) {
        for (Jglop.iv = 0;Jglop.iv < f->model.m;Jglop.iv++){
          int jmem=LBM_glopmuti_to_imem(simu,&Jglop);
          IsNonZero(solver,imem,jmem);
        } // Jiv
        } // end Jix2
        } // end Jix1
        } // end Jix0
        //
    } // end iv
    } // end ix2
    } // end ix1
    } // end ix0
    } // end ic2
    } // end ic1
    } // end ic0
  }
}
//* 
void LBM_Mass_SourceJac_Coupling(Simulation *simu, LinearSolver *solver, int isky){
  LBM_glop_multindex Iglop;
  for(Iglop.ie = 0; Iglop.ie < simu->macromesh.nbelems; Iglop.ie++){
    field *f = simu->fd + Iglop.ie; // &(simu->fd[ie])
    for (Iglop.ic[0] = 0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
    for (Iglop.ic[1] = 0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
    for (Iglop.ic[2] = 0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
    for (Iglop.ix[0] = 0; Iglop.ix[0] < f->deg[0]+1; Iglop.ix[0]++) {
    for (Iglop.ix[1] = 0; Iglop.ix[1] < f->deg[1]+1; Iglop.ix[1]++) {
    for (Iglop.ix[2] = 0; Iglop.ix[2] < f->deg[2]+1; Iglop.ix[2]++) {
      for (int iv1 = 0;iv1 < f->model.m;iv1++){
        Iglop.iv=iv1;
        int imem= LBM_glopmuti_to_imem(simu,&Iglop);
        for (int iv2 = 0;iv2 < f->model.m;iv2++){
        Iglop.iv=iv2;
        int jmem= LBM_glopmuti_to_imem(simu,&Iglop);
        IsNonZero(solver,imem,jmem);
        }
      }
    } // end ix2
    } // end ix1
    } // end ix0
    } // end ic2
    } // end ic1
    } // end ic0
  }
}
void LBM_InternalCoupling_SelectVel(Simulation *simu,  LinearSolver *solver){
  //
  bool verbose = false;
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_advec_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  //
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  //
  LBM_glop_multindex Iglop;
  LBM_glop_multindex Jglop;
  // restrict to advected variable iv_advec
  Iglop.iv = iv_advec;
  Jglop.iv = iv_advec;
  // loop on macro elements
  for(Iglop.ie = 0; Iglop.ie < simu->macromesh.nbelems; Iglop.ie++){
    field *f = simu->fd + Iglop.ie; // &(simu->fd[ie])
    Jglop.ie = Iglop.ie; // sync J element
    // Loop on the subcells
    for (Iglop.ic[0]=0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
      Jglop.ic[0]=Iglop.ic[0]; // sync cell ic0
    for (Iglop.ic[1]=0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
      Jglop.ic[1]=Iglop.ic[1]; // sync cell ic1
    for (Iglop.ic[2]=0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
      Jglop.ic[2]=Iglop.ic[2]; // sync cell ic2
    // Loop on glops in the subcell
    for (Iglop.ix[0] = 0; Iglop.ix[0] < f->deg[0]+1; Iglop.ix[0]++) {
    for (Iglop.ix[1] = 0; Iglop.ix[1] < f->deg[1]+1; Iglop.ix[1]++) {
    for (Iglop.ix[2] = 0; Iglop.ix[2] < f->deg[2]+1; Iglop.ix[2]++) {
    // recover memory index for Iglop;
    int imem = LBM_glopmuti_to_imem(simu,&Iglop);
    // loop on grapsi component
      for(int dim0 = 0; dim0 < 3; dim0++) {
          // sync Jglop.ix with Iglop.ix
          for (int k=0;k<3;k++){
            Jglop.ix[k]=Iglop.ix[k];
          }
          // add contribution only along dim0 due to separability in xyz of basis functions
          for (Jglop.ix[dim0] = 0; Jglop.ix[dim0] < f->deg[dim0]+1; Jglop.ix[dim0]++){
            int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
            if (verbose){
            printf("Int Switching on (%i,%i) in cell (%i,%i,%i) iv %i\n",
              imem,jmem,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],Iglop.iv);
            }
            IsNonZero(solver,imem,jmem);
          } // end loop Jglop.ix[dim0]
          //
      } // dim0 loop
    } // end Iglop.ix[0]
    } // end Iglop.ix[1]
    } // end Iglop.ix[2]
    // loop along grapsi components (x,y,z)
    } // end Iglop.ic[2]
    } // end Iglop.ic[1]
    } // end Iglop.ic[0]
  } // end Iglop.ie
} 

//**********************************************************************************//
//*********************************************************************************//
void LBM_FluxCoupling_SelectVel(Simulation * simu, LinearSolver * solver){
  //
  const bool verbose = false;
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_advec_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  if (verbose){
    printf(" velocity:(%f,%f,%f) \n:",vel[0],vel[1],vel[2]);
  }
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  LBM_glop_multindex Iglop;
  LBM_glop_multindex Jglop;
  // restrict to advected variable iv
  Iglop.iv = iv_advec;
  Jglop.iv = iv_advec;
  // loop on macro elements
  for (Iglop.ie=0; Iglop.ie < simu->macromesh.nbelems;Iglop.ie++){
    field *f= simu->fd+Iglop.ie;
  // loop on subcells
  for (Iglop.ic[0]=0;Iglop.ic[0] < f->raf[0];Iglop.ic[0]++){
  for (Iglop.ic[1]=0;Iglop.ic[1] < f->raf[1];Iglop.ic[1]++){
  for (Iglop.ic[2]=0;Iglop.ic[2] < f->raf[2];Iglop.ic[2]++){
    if (verbose){
    printf(" ic:(%i;%i,%i) \n",Iglop.ic[0],Iglop.ic[1],Iglop.ic[2]);
    }
  // Sweeping subcell faces pairs in the three directions (ie two subfaces per dim)  
  for (int dim0 = 0; dim0 < 3; dim0++) {
    // do not treat last cell as a "left" cell
    if (Iglop.ic[dim0] == f->raf[dim0] - 1) {
      continue;
    }
    //
    const int altdim1[3] = { 1, 0, 0 };
    const int altdim2[3] = { 2, 2, 1 };
    int dim1= altdim1[dim0];
    int dim2= altdim2[dim0];
    // select right subcell along dim0
    Jglop.ie=Iglop.ie;
    Jglop.ic[dim0]=Iglop.ic[dim0]+1;
    Jglop.ic[dim1]=Iglop.ic[dim1];
    Jglop.ic[dim2]=Iglop.ic[dim2];
    // select rightmost glop in dim0 direction from the left subcell
    Iglop.ix[dim0] = f->deg[dim0];
    // select leftmost glop in dim0 direction from the right subcell
    Jglop.ix[dim0] = 0;
    // loop on left glop on the face
    for (Iglop.ix[dim1] = 0; Iglop.ix[dim1] < f->deg[dim1]+1; Iglop.ix[dim1]++) {
      // matching glop ix[idim1] in the right subcell
      Jglop.ix[dim1] = Iglop.ix[dim1];
      for (Iglop.ix[dim2] = 0; Iglop.ix[dim2] < f->deg[dim2]+1; Iglop.ix[dim2]++) {
        // matching glop ix[idim2] in the right subcell
        Jglop.ix[dim2] = Iglop.ix[dim2];
        // get memory index for the advected variable iv_advec
        int imem = LBM_glopmuti_to_imem(simu,&Iglop);
        int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
        // 
        // Compute the normal vector 
        schnaps_real vnds[3];
        {
          schnaps_real xref[3], wpg3;
          ref_pg_vol(f->deg, f->raf, Iglop.ipg, xref, &wpg3, NULL);
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
        }
        // compute the scalar product with advection velocity
        schnaps_real vidotn = vnds[0] * vel[0] + vnds[1] * vel[1] +vnds[2] * vel[2];
        //
        if (vidotn > _SMALL){
          if (verbose){
          printf("Flux Switching on (%i,%i) in cell (%i,%i,%i) iv %i\n",imem,jmem,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],Iglop.iv);
          }
          IsNonZero(solver, imem, jmem);
        }
        if (vidotn < -_SMALL){
          if (verbose){
          printf("Flux Switching on (%i,%i) in cell (%i,%i,%i) iv %i\n",jmem,imem,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],Iglop.iv);
          }
          IsNonZero(solver, jmem, imem);
        }
    } // end Iglop.ix[dim2]
    } // end Iglop.ix[dim1]
    //
  } // end dim0
  } // end ic2
  } // end ic1
  } // end ic0
  } // end ie
}
//*********************************************************************************//
void LBM_InterfaceCoupling_SelectVel(Simulation *simu,  LinearSolver *solver)
{
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_advec_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  printf(" velocity:(%f,%f,%f) \n:",vel[0],vel[1],vel[2]);
  if (velnorm2 < _SMALL) {
    return;                     //-> stagnant v=0 node no coupling
  }
  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    const unsigned int m = fL->model.m;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      offsetR = fsize * ieR;
    }
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
      //
      schnaps_real vidotn= vnds[0] * vel[0] + vnds[1] * vel[1] + vnds[2] * vel[2];
      //
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
	int imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgL, iv_advec) + offsetL;
	int imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv_advec) + offsetR;
	if (vidotn > _SMALL){
	  IsNonZero(solver, imem1, imem2);
	}
	else{
	if (vidotn < -_SMALL){
	  IsNonZero(solver, imem2, imem1);
	}
	}
	//
   } 
   else { // The point is on the boundary.
     // only the advected node varies
	   int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv_advec) + offsetL;
	   // the boundary flux may have a dependency on the local field values
	   bool with_bc_jacobian= false;
	   if (with_bc_jacobian){
       for(int iv2 = 0; iv2 < m; iv2++) { 
	     int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2) + offsetL;
	     IsNonZero(solver, imem1, imem2);
	     } // iv2
	   } 
	   else
	   {
  	   IsNonZero(solver, imem1, imem1);
	   }
    } // end else end if
    } // ipgfl
  } // macroface loop
}
//

// NL non splitted implicit time scheme *******************************************//
void LBMThetaTimeScheme_NL_LocalNewton(LBMSimulation * lbsimu,
                                       schnaps_real theta,
                                       schnaps_real tmax, schnaps_real dt)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int nb_nodes = lsd->lb_model->q;      // velocity nodes on the lattice(s)
  Simulation *micsimu = &(lbsimu->micro_simu);
  Simulation *macsimu = &(lbsimu->macro_simu);
  //
  int nraf[3] = { micsimu->fd[0].raf[0],
    micsimu->fd[0].raf[1],
    micsimu->fd[0].raf[2]
  };
  int deg[3] = { micsimu->fd[0].deg[0],
    micsimu->fd[0].deg[1],
    micsimu->fd[0].deg[2]
  };
  int nb_ipg = NPG(deg, nraf);
  //
  lbsimu->dt = dt;
  micsimu->dt = dt;
  macsimu->dt = dt;
  //
  int itermax = (int) (tmax / dt) + 1;
  printf("Called with tmax=%f dt=%f Nb iterations:%i \n", tmax, dt,
         itermax);
  lbsimu->itermax = itermax;
  micsimu->itermax_rk = itermax;
  macsimu->itermax_rk = itermax;
  int freq = (1 >= lbsimu->itermax / 10) ? 1 : lbsimu->itermax / 10;
  lbsimu->tnow = 0.0;
  micsimu->tnow = 0.0;
  macsimu->tnow = 0.0;
  //
  MatrixStorage ms = KLU_CSR;
  Solver st = LU;
  LinearSolver solver_imp[nb_nodes]; 
  // NOTE we will need another set for the inverse operators 
  // theta -> 1 - theta  dt -> -dt , swap tn tn+1 , this affects only values
  // but the op structure is different (using the opposite node velocity won't work)
  //LinearSolver solver_imp_inv[nb_nodes]; 
  // allocate buffer for explicit rhs storage
  schnaps_real *exp_rhs;
  exp_rhs = (schnaps_real *) calloc(micsimu->wsize,sizeof(schnaps_real));
  //
  for (int i = 0; i < nb_nodes; i++) {
    // initialize klu linear solver for a given transport pattern.
    lsd->current_advec_node_index = i;
    InitLinearSolver(&(solver_imp[i]), micsimu->wsize, &ms, &st);
    // define operator profile
    LBM_Mass_SourceJac_Coupling(micsimu,&(solver_imp[i]),0);
    LBM_InternalCoupling_SelectVel(micsimu, &(solver_imp[i]));
    LBM_FluxCoupling_SelectVel(micsimu, &(solver_imp[i]));
    LBM_InterfaceCoupling_SelectVel(micsimu, &(solver_imp[i]));
    AllocateLinearSolver(&(solver_imp[i]));
    // Display Ops (debug only, for small op size !!)
/*    printf("*********************************************\n");*/
/*    printf("******************i=%i***********************\n", i);*/
/*    printf("*********************************************\n");*/
/*    //DisplayLinearSolver(&(solver_imp[i]));*/
    // 
      // check match between topological and klu btf structure
      CheckKLUBLockCellConstitency(micsimu, &(solver_imp[i]),false);
  }
  // time iterations 
  for (int iter = 0; iter < itermax; iter++) {
    //
    if (iter % freq == 0) {
      printf(" iter %i/%i t=%f\n", iter, itermax, lbsimu->tnow);
    }
    //
    lbsimu->iter_time = iter;
    macsimu->iter_time_rk = iter;
    micsimu->iter_time_rk = iter;
    //
    // Full NL Step for each velocity
    // fully symmetrized call
    for (int i=0;i< nb_nodes;i++){
      lsd->current_advec_node_index=i;
          // Call routine
          // DEBUG
          if (i==0){ 
          printf("DEBUG : only iv=0\n");
          LBMThetaTimeScheme_NL_LocalNewton_solve(micsimu,&(solver_imp[i]),exp_rhs,theta,dt);
          }
    }
    // 
    for (int i=nb_nodes-1;i>-1;i--){
      lsd->current_advec_node_index=i;
          // Call routine inverse // FIXME need to pass tstart tend to allow for easy swapping
          //LBMThetaTimeScheme_NL_LocalNewton_solve(micsimu,&(solver_imp[i]),exp_rhs,1.0-theta,-dt);
    }
    // sync time at end of iteration
    lbsimu->tnow = 0.0;
    micsimu->tnow = 0.0;
    macsimu->tnow = 0.0;
  } // end of time iterations
  // cleanup
  for (int i = 0; i < nb_nodes; i++) {
    FreeLinearSolver(&(solver_imp[i]));
  }
  free(exp_rhs);
  //
}

void LBMThetaTimeScheme_NL_LocalNewton_solve(Simulation * simu,
                                             LinearSolver * impsolv,
                                             schnaps_real * exp_rhs,
                                             schnaps_real theta,
                                             schnaps_real dt)
{
  //
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  //
  assert(impsolv->storage_type == KLU_CSR);
  KLU *klu = (KLU *) impsolv->matrix;
  int nblocks = klu->symbolic->nblocks;
  // 
  int maxiter_newton = 1;
  schnaps_real newton_tol = 1.0e-4;
  schnaps_real err = 1.0;
  bool exp_is_assembled[nblocks];
  //
  // copy state in solution vector
  for (int imem=0;imem < simu->wsize;imem++){
    impsolv->sol[imem]=simu->w[imem];
  }
  // clear explicit rhs buffer
  for (int imem=0;imem < simu->wsize;imem++){
    exp_rhs[imem]=0.0;
  }
  // loop on blocks
  for (int iblock = nblocks - 1; iblock > -1; iblock--) {
    //
    // Computation of constant explicit rhs part at time n
    // Local newton iteration
    err = 1.0;
    exp_is_assembled[iblock]=false;
    // copy of current state to sol vector
    //
    //
    for (int iter = 0; iter < maxiter_newton; iter++) {
      // Implicit op block assembly
      printf("*********************************************\n");
      printf("******************i=%i***********************\n", lsd->current_advec_node_index);
      printf("*********************************************\n");
      //DisplayLinearSolver(impsolv);
      printf("Block wise Assembly for block %i \n", iblock);
      LBM_NL_LocalNewton_Assembly_Blockwise(simu, impsolv, theta,
                                            exp_rhs, exp_is_assembled[iblock], iblock);
      exp_is_assembled[iblock]=true;
      //DisplayLinearSolver(impsolv);
      // Implicit op block factorization TODO 
      // Block Solve TODO 
      // error computation
      //
      if (err < newton_tol)
        break;
    }
    //
    //
  } // end loop on blocks
  //
}
//
int CheckKLUBLockCellConstitency(Simulation *simu,LinearSolver *solv,bool verbose){
  // recover advection velocity
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int inode_advec = lsd->current_advec_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[inode_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  // check for null advection velocity (spatially local case)
  bool vel_is_zero = (velnorm2 < _SMALL);
  //
  KLU *(klumat) = solv->matrix;
  int nblocks = klumat->symbolic->nblocks;
  //
  const bool display_value=false;
  const bool display_only=false;
  // Optionnal display of some index value to check ordering
  if (display_value){
    LBM_glop_multindex Iglop;
    LBM_glop_multindex Jglop;
    for (int iblock = 0; iblock < nblocks; iblock++){
    int ibl_start = klumat->symbolic->R[iblock];
    int ibl_end = klumat->symbolic->R[iblock + 1];
    for (int i=ibl_start;i < ibl_end;i++){
      LBM_imem_to_glopmulti(simu,klumat->symbolic->P[i],&Iglop);
      for (int jblock = 0; jblock < nblocks; jblock++){
        int jbl_start = klumat->symbolic->R[jblock];
        int jbl_end = klumat->symbolic->R[jblock + 1];
      for (int j=jbl_start;j < jbl_end;j++){
        LBM_imem_to_glopmulti(simu,klumat->symbolic->Q[j],&Jglop);
        if (fabs(GetKLU(klumat, klumat->symbolic->P[i], klumat->symbolic->Q[j])) < 1e-8) {
          printf(" ");
        } else {
          printf("%i",Jglop.ic[2]);
        }
      } // j
      printf("|");
      } // jblock
      printf("\n");
      } // i
      for (int j=0; j < simu->wsize;j++){
        printf("-");
      }
      printf("\n");
    } // iblock
  if( display_only){
    return 1;
  }
  } // end if display_value


  // loop on klu generated blocks
  for (int iblock = nblocks-1; iblock > -1; iblock--){
    int bl_start = klumat->symbolic->R[iblock];
    int bl_end = klumat->symbolic->R[iblock + 1];
    // get multiindex of first point in block
    int imem = klumat->symbolic->P[bl_start]; // natural memory index
    //assert(test_blstart == bl_start);
    LBM_glop_multindex i_glop;
    LBM_imem_to_glopmulti(simu, imem, &i_glop);
    if (verbose){
    printf("I iblock : %i imem %i Ref ie %i ivarindex %i ,ic (%i,%i,%i) ix (%i,%i,%i) igp: %i,iv: %i \n",
      iblock,imem,i_glop.ie, i_glop.ivarindex,
      i_glop.ic[0],i_glop.ic[1],i_glop.ic[2],
      i_glop.ix[0],i_glop.ix[1],i_glop.ix[2],i_glop.ipg,i_glop.iv);
    }
    // loop on memory index inside the block
    for (int ibl = bl_start; ibl < bl_end; ibl++) {
      int jmem= klumat->symbolic->P[ibl];
      LBM_glop_multindex j_glop;
      LBM_imem_to_glopmulti(simu, jmem, &j_glop);
      if (verbose){
        printf("ibl %i jmem %i J ie %i ivarindex %i ,ic (%i,%i,%i) ix (%i,%i,%i) igp: %i,iv: %i \n",
          ibl,jmem,j_glop.ie,j_glop.ivarindex,
           j_glop.ic[0],j_glop.ic[1],j_glop.ic[2],
          j_glop.ix[0],j_glop.ix[1],j_glop.ix[2],j_glop.ipg,j_glop.iv);
      }
      bool check=true;
      // all memory indexes map to the same cell
      for (int k=0;k<3;k++){
        check = check & (i_glop.ic[k] == j_glop.ic[k]);
      }
      // in the v=0, all memory indexes in the block are also related to the same glop
      if (vel_is_zero){
        for (int k=0;k<3;k++){
          check = check & (i_glop.ix[k]== j_glop.ix[k]);
        }
      }
      assert(check);
      //
      jmem= klumat->symbolic->Q[ibl];
      LBM_imem_to_glopmulti(simu, jmem, &j_glop);
      check=true;
      // all memory indexes map to the same cell
      for (int k=0;k<3;k++){
        check = check & (i_glop.ic[k]== j_glop.ic[k]);
      }
      // in the v=0, all memory indexes in the block are also related to the same glop
      if (vel_is_zero){
        for (int k=0;k<3;k++){
          check = check & (i_glop.ix[k]== j_glop.ix[k]);
        }
      }
      assert(check);
    } // end loop ibl
    if (verbose){
      printf(" -> OK\n");
    }
  } // end loop on block
  //
  //
  //
  return 1;
} 
// CAUTION this routine assumes Block/Subcell  bijection when the advection velocity is non-zero
// and block/glop bijection in the 0 velocity case.
// for now all memory buffers and operators are assumed in the initial natural ordering
// the Blockwise factor/solve will reorder on entry and exit if necessary 
void LBM_NL_LocalNewton_Assembly_Blockwise(Simulation * simu,
                                           LinearSolver * impsolv,
                                           const schnaps_real theta,
                                           schnaps_real * exp_rhs,
                                           bool exp_is_assembled,
                                           int iblock)
{
  // 
const int h20_refnormal[6][3]={{0,-1,0},
			       {1,0,0},
			       {0,1,0},
			       {-1,0,0},
			       {0,0,1},
			       {0,0,-1} };
   // Implementation Checklist 
  // -> mass contribution DONE
  // -> source contribution DONE
  // -> jac source contribution DONE
  // -> flux contribution DONE
  // -> interface flux + boundary + bounday jacobian TODO
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int iv_advec = lsd->current_advec_node_index;
  schnaps_real vel[3];
  schnaps_real velnorm2 = 0.0;
  for (int i = 0; i < lsd->lb_model->d; i++) {
    vel[i] = lsd->lb_model->vi[iv_advec][i];
    velnorm2 += vel[i] * vel[i];
  }
  // check for null advection velocity (spatially local case)
  bool vel_is_zero = (velnorm2 < _SMALL);
  //
  // time quantities
  schnaps_real t_exp=simu->tnow;
  schnaps_real t_imp=simu->tnow+simu->dt;
  schnaps_real thdt_imp= theta * simu->dt;
  schnaps_real thdt_exp= -(1.-theta) * simu->dt;
  //
  // thdt_exp-thdt_imp = dt *(-1+theta-theta)= -dt
  //  
  bool verbose = true;
  //
  KLU *(klumat) = impsolv->matrix;
  int ibl_start = klumat->symbolic->R[iblock];
  // get multiindex of first point in block to recover element and cell multiindex;
  // nb the btf is unchanged by numerical factorization
  int imem_ref=klumat->symbolic->P[ibl_start]; // nb the btf is unchanged by numerical factorization
  // so we can use the symb pre-ordering even is the matrix is lu factorized 
  LBM_glop_multindex Iglop; // line multiindex
  LBM_glop_multindex Jglop; // column multiindex
  LBM_imem_to_glopmulti(simu,imem_ref,&Iglop);
  LBM_imem_to_glopmulti(simu,imem_ref,&Jglop);
  if (verbose){
    if (vel_is_zero){
    printf(" block :%i element : %i subcell :(%i,%i,%i) glop(%i,%i,%i)\n",
      iblock,Iglop.ie,
      Iglop.ic[0],Iglop.ic[1],Iglop.ic[2],
      Iglop.ix[0],Iglop.ix[1],Iglop.ix[2]);
    }
    else{
    printf(" block :%i element : %i subcell :(%i,%i,%i)\n",iblock,Iglop.ie,Iglop.ic[0],Iglop.ic[1],Iglop.ic[2]);
    }
  }
  field *f = simu->fd+ Iglop.ie;
  if (!vel_is_zero){
  // loop on points in subcell //
  for (Iglop.ix[0] = 0; Iglop.ix[0] < f->deg[0]+1; Iglop.ix[0]++) {
    for (Iglop.ix[1] = 0; Iglop.ix[1] < f->deg[1]+1; Iglop.ix[1]++) {
      for (Iglop.ix[2] = 0; Iglop.ix[2] < f->deg[2]+1; Iglop.ix[2]++) {
        // copy I glop base to Iglop for local couplings
        for (int k=0;k<3;k++){
          Jglop.ix[k]=Iglop.ix[k];
        }
        // glop volume index in subcell
        xyz_to_ipg(f->raf,f->deg,Iglop.ic,Iglop.ix,&(Iglop.ipg));
        // recover local geometric quantities 
        schnaps_real Idtau[3][3], Icodtau[3][3], Ixpgref[3], Ixphy[3], Iwpg;
        ref_pg_vol(f->deg, f->raf, Iglop.ipg, Ixpgref, &Iwpg, NULL);
        schnaps_ref2phy(f->physnode, // phys. nodes
        Ixpgref, // xref
        NULL, -1, // dpsiref, ifa
        Ixphy, Idtau, // xphy, dtau
        Icodtau, NULL, NULL); // codtau, dpsi, vnds
        schnaps_real Idet = dot_product(Idtau[0], Icodtau[0]);
        // Mass matrix element
        schnaps_real I_Mass= Iwpg * Idet;
        // Recover local values
        schnaps_real I_wloc[f->model.m];
        for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
          int imem= LBM_glopmuti_to_imem(simu,&Iglop);
          I_wloc[Iglop.iv]=impsolv->sol[imem]; 
        }
        // Source spatially diagonal term - Direct
        schnaps_real I_source_imp[f->model.m];
        f->model.Source(Ixphy, t_imp, I_wloc, I_source_imp);
        // Source spatially diagonal term jacobian
        schnaps_real I_Jacsource_imp[f->model.m * f->model.m];
        // recover jacobian of the source term
        f->model.SourceJac(Ixphy, t_imp, I_wloc, I_Jacsource_imp);
        // build spatially local part of rhs 
        for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
          int imem= LBM_glopmuti_to_imem(simu,&Iglop);
          schnaps_real val = I_source_imp[Iglop.iv]; 
          // add source jacobian contribution
          for (int jv=0; jv < f->model.m;jv++){
            val -= I_Jacsource_imp[jv * f->model.m+Iglop.iv] * I_wloc[jv];
          }
          impsolv->rhs[imem] = val * I_Mass * thdt_imp ;
        }
        // build spatially local par of explicit rhs if needed
        if (!exp_is_assembled){
          schnaps_real I_source_exp[f->model.m];
          f->model.Source(Ixphy, t_exp, I_wloc, I_source_imp);
          for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
            int imem= LBM_glopmuti_to_imem(simu,&Iglop);
            exp_rhs[imem] = -I_source_exp[Iglop.iv]; 
            exp_rhs[imem] *=  thdt_exp ;
            exp_rhs[imem] += I_wloc[Iglop.iv];
            exp_rhs[imem] *= I_Mass;
          } // end iv
         } // end if !exp_is_assembled
        // add spatially local contributions to implicit operator
        for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
          int imem= LBM_glopmuti_to_imem(simu,&Iglop);
          AddLinearSolver(impsolv,imem,imem,I_Mass);
          // add source jacobian contribution
          for (Jglop.iv=0; Jglop.iv < f->model.m;Jglop.iv++){
            schnaps_real tmp_val = I_Jacsource_imp[Jglop.iv * f->model.m +Iglop.iv];
            tmp_val *= thdt_imp;
            tmp_val *=I_Mass;
            int jmem= LBM_glopmuti_to_imem(simu,&Jglop);
            AddLinearSolver(impsolv,imem,jmem,tmp_val);
          } // end Jglop.iv
        } // end Iglop.iv
        // now add internal (volumic) transport contributions
        // transport is restricted to one velocity node
        schnaps_real wL[f->model.m];
        for (int iv=0;iv<f->model.m;iv++){
          wL[iv]=0.0;
        }
        wL[iv_advec]=1.0;
        // loop on gradient component direction
        for (int gdim=0;gdim < 3;gdim++){
          // loop on basis functions in the subcell
          for (int k=0;k<3;k++){
            Jglop.ix[k]=Iglop.ix[k];
          }
          // restrict to glops along the gradient component
          // due to separability in xyz of basis functions.
          for (Jglop.ix[gdim]=0;Jglop.ix[gdim] < f->raf[gdim];Jglop.ix[gdim]++){
            schnaps_real dphiref[3]={0.0,0.0,0.0}; // grad psi reference 
            schnaps_real dphiL[3]={0.0,0.0,0.0}; // grad psi / physical
            dphiref[gdim] = dlag(f->deg[gdim], Jglop.ix[gdim], Iglop.ix[gdim])* f->raf[gdim];
            schnaps_ref2phy(f->physnode,
              Ixpgref,
              dphiref, // dphiref
              -1,  // ifa
              NULL, // xphy
              Idtau,
              Icodtau,
              dphiL, // dphi
              NULL);  // vnds
            //
            schnaps_real flux[f->model.m];
            f->model.NumFlux(wL, wL, dphiL, flux);
            schnaps_real val= - thdt_imp * flux[iv_advec] * Iwpg;
            Iglop.iv = iv_advec;
            Jglop.iv = iv_advec;
            int imem = LBM_glopmuti_to_imem(simu,&Iglop);
            int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
            AddLinearSolver(impsolv,imem,jmem,val);
            // contribution to explicit rhs
            if (!exp_is_assembled){
              val = -thdt_exp * flux[iv_advec] * Iwpg;
              exp_rhs[imem] += val * impsolv->sol[jmem];
            } // end if !exp_is_assembled
          } // end Jglop.ix[gdim]
        } // end gdim
        // end of internal tranport contribution
      } // end iglop.ix[2]
    } // iglop.end ix[1]
  } // iglop.end ix[0]
  // now compute flux contribution of neighbour cells
  // loop on axis dimension (formal normal to faces)
  const int altdim1[3] = {1, 0, 0};
  const int altdim2[3] = {2, 2, 1};
  schnaps_real wL[f->model.m], wR[f->model.m];
  for (int iv=0;iv < f->model.m;iv++){
    wL[iv]=0.0;
    wR[iv]=0.0;
  }
  Iglop.iv=iv_advec; 
  Jglop.iv=iv_advec;
  //
  for (int dim0=0;dim0 < 3;dim0++){
    if (Iglop.ic[dim0] == f->raf[dim0]) continue; // 
    Jglop.ic[dim0]=Iglop.ic[dim0]+1; // right neighbour cell
    // face dimensions
    int dim1=altdim1[dim0];
    int dim2=altdim2[dim0];
    Iglop.ix[dim0]=f->deg[dim0]; // last point of left cell along dim0
    Jglop.ix[dim0]=0; // first points of right cell alond dim0
    for (Iglop.ix[dim1]=0;Iglop.ix[dim1] < f->deg[dim1]+1;Iglop.ix[dim1]++){
        Jglop.ix[dim1]=Iglop.ix[dim1];
      for (Iglop.ix[dim2]=0;Iglop.ix[dim2] < f->deg[dim2]+1;Iglop.ix[dim2]++){
          Jglop.ix[dim2]=Iglop.ix[dim2];
          // recover left/right values
          int imem=LBM_glopmuti_to_imem(simu,&Iglop); // nb this call updates Iglop.ipg
          int jmem=LBM_glopmuti_to_imem(simu,&Jglop); // nb this call updates Jglop.ipg
          wL[iv_advec]=impsolv->sol[imem];
          wR[iv_advec]=impsolv->sol[jmem];
          // Compute the normal vector
          schnaps_real vnds[3];
          {
            schnaps_real xref[3], wpg3;
            ref_pg_vol(f->deg, f->raf, Iglop.ipg, xref, &wpg3, NULL);
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
            schnaps_real h1h2 = 1. / f->raf[dim1] / f->raf[dim2];
            vnds[0] = codtau[0][dim0] * h1h2;
            vnds[1] = codtau[1][dim0] * h1h2;
            vnds[2] = codtau[2][dim0] * h1h2;
          }
          // weight
          schnaps_real wpg
            = wglop(f->deg[dim1], Iglop.ix[dim1])
            * wglop(f->deg[dim2], Iglop.ix[dim2]);
          // flux computation
          schnaps_real flux[f->model.m];
          f->model.NumFlux(wL,wR,vnds,flux);
          //
          schnaps_real val=  flux[iv_advec] * wpg;
          // contribution to implicit residual
          impsolv->rhs[imem] -= thdt_imp * val;
          // contribution to explicit rhs
          if (!exp_is_assembled){
            exp_rhs[imem] += thdt_exp * val;
          }
      }  // end Iglop.ix[dim2]
    } // end Iglop.ix[dim1]
    //
  } //end dim0
  // Compute Interface contribution.
  Iglop.iv= iv_advec; // restrict to the advected quantity 
  // 
  // we assume DEG,RAF constant across macroelements
  printf(" ie :%i \n",Iglop.ie);
  for (int ifa = 0 ; ifa < 6; ifa++){
    int dim0 =-1;
    int sgn = 0;
    printf(" ifa %i \n",ifa);
    // recover dimension and orientation
    for (int idim=0; idim<3;idim++){
      if (h20_refnormal[ifa][idim] ==1){
        dim0=idim;
        sgn = 1;
      }
      if (h20_refnormal[ifa][idim] ==-1){
        dim0=idim;
        sgn = -1;
      }
    } // end idim
    int dim1 = altdim1[dim0];
    int dim2 = altdim2[dim0];
    // initialize Jglop in same element/subcell
    Jglop.ie = Iglop.ie;
    Jglop.ic[dim1] = Iglop.ic[dim1];
    Jglop.ic[dim2] = Iglop.ic[dim2];
    Jglop.ic[dim0] = Iglop.ic[dim0];
    //
    bool is_boundary = false;
    bool is_border_cell = (Iglop.ic[dim0] == (f->raf[dim0]-1) * (1+sgn)/2);
    if (is_border_cell){
      // set glop index along dim0
      Iglop.ix[dim0] = (f->deg[dim0]) * (1-sgn)/2;
      // get neighbour element
      int ie_test = simu->macromesh.elem2elem[6 * Iglop.ie  + ifa];
      if (verbose){
        printf(" dim0 : %i sgn %i ifa %i  Neighbour element %i \n",
          dim0,sgn,ifa,simu->macromesh.elem2elem[6 * Iglop.ie  + ifa]);
        }
      // neighbour element field test
      if (ie_test > -1){
        is_boundary = false;
        Jglop.ie = ie_test;
        field *f_n = simu->fd + Jglop.ie; 
        // get neighbourg subcell 
        Jglop.ic[dim0] = (f_n->raf[dim0]-1) * (1-sgn)/2;
        // set glops index along dim0
        Jglop.ix[dim0] = (f_n->deg[dim0]) * (1-sgn)/2;
      }
      // boundary
      else{
        is_boundary = true;
      } // end if ie_test
      // now loop on glops on the subcell face
      for (Iglop.ix[dim1] = 0;Iglop.ix[dim1] < f->deg[dim1]+1; Iglop.ix[dim1]++){
        Jglop.ix[dim1] = Iglop.ix[dim1];
        for (Iglop.ix[dim2] = 0;Iglop.ix[dim2] < f->deg[dim2]+1; Iglop.ix[dim2]++){
          Jglop.ix[dim2] = Iglop.ix[dim2];
          // recover memory indexes
          int imem = LBM_glopmuti_to_imem(simu,&Iglop);
          int jmem = LBM_glopmuti_to_imem(simu,&Jglop);
          // compute normal vector
          schnaps_real vnds[3];
          {
            schnaps_real xref[3], wpg3;
            ref_pg_vol(f->deg, f->raf, Iglop.ipg, xref, &wpg3, NULL);
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
            schnaps_real h1h2 = 1. / f->raf[dim1] / f->raf[dim2];
            vnds[0] = codtau[0][dim0] * h1h2;
            vnds[1] = codtau[1][dim0] * h1h2;
            vnds[2] = codtau[2][dim0] * h1h2;
          }
          // weight
          schnaps_real wpg
            = wglop(f->deg[dim1], Iglop.ix[dim1])
            * wglop(f->deg[dim2], Iglop.ix[dim2]);
          // flux computation
          if (is_boundary){
            wL[iv_advec]=impsolv->sol[imem];
          }{
            wL[iv_advec]=impsolv->sol[imem];
            wR[iv_advec]=impsolv->sol[jmem];
          }
          //
        } // end Iglop.ix[dim2]
      } // end Iglop.ix[dim1]
    } // end if is_border_cell
  } // end ifa
  printf(" \n");
  //
  //
  //
  //
  } // end if !vel_is_zero 
  //
  // case of zero velocity 
  // the block is reduced to 1 glop 
  if (vel_is_zero){
    LBM_imem_to_glopmulti(simu,imem_ref,&Iglop);
    LBM_imem_to_glopmulti(simu,imem_ref,&Jglop);
    //
    schnaps_real Idtau[3][3], Icodtau[3][3], Ixpgref[3], Ixphy[3], Iwpg;
    ref_pg_vol(f->deg, f->raf, Iglop.ipg, Ixpgref, &Iwpg, NULL);
    schnaps_ref2phy(f->physnode, // phys. nodes
    Ixpgref, // xref
    NULL, -1, // dpsiref, ifa
    Ixphy, Idtau, // xphy, dtau
    Icodtau, NULL, NULL); // codtau, dpsi, vnds
    schnaps_real Idet = dot_product(Idtau[0], Icodtau[0]);
    // Mass matrix element
    schnaps_real I_Mass= Iwpg * Idet;
    // Recover local values
    schnaps_real I_wloc[f->model.m];
    for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
      int imem= LBM_glopmuti_to_imem(simu,&Iglop);
      I_wloc[Iglop.iv]=impsolv->sol[imem]; 
    }
    // Source spatially diagonal term - Direct
    schnaps_real I_source_imp[f->model.m];
    f->model.Source(Ixphy, t_imp, I_wloc, I_source_imp);
    // Source spatially diagonal term jacobian
    schnaps_real I_Jacsource_imp[f->model.m * f->model.m];
    // recover jacobian of the source term
    for (int k=0;k< f->model.m * f->model.m;k++){
      I_Jacsource_imp[k]=0.0;
    }
    // build spatially local part of rhs 
    for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
      int imem= LBM_glopmuti_to_imem(simu,&Iglop);
      impsolv->rhs[imem] = -I_source_imp[Iglop.iv]; 
      // add source jacobian contribution
      for (int jv=0; jv < f->model.m;jv++){
        impsolv->rhs[imem] += I_Jacsource_imp[jv * f->model.m+Iglop.iv] * I_wloc[jv];
      }
      impsolv->rhs[imem] *= I_Mass * thdt_imp ;
    }
    // build spatially local par of explicit rhs if needed
    if (!exp_is_assembled){
      schnaps_real I_source_exp[f->model.m];
      f->model.Source(Ixphy, t_exp, I_wloc, I_source_imp);
      for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
        int imem= LBM_glopmuti_to_imem(simu,&Iglop);
        exp_rhs[imem] = I_Mass * (I_wloc[Iglop.iv] + thdt_exp * I_source_exp[Iglop.iv]);
      } // end iv
     } // end if !exp_is_assembled
    // add spatially local contributions to implicit operator
    for (Iglop.iv=0;Iglop.iv < f->model.m;Iglop.iv++){
      int imem= LBM_glopmuti_to_imem(simu,&Iglop);
      AddLinearSolver(impsolv,imem,imem,I_Mass);
      // add source jacobian contribution
      for (Jglop.iv=0; Jglop.iv < f->model.m;Jglop.iv++){
        schnaps_real tmp_val = I_Jacsource_imp[Jglop.iv * f->model.m +Iglop.iv] * I_wloc[Jglop.iv];
      tmp_val *= thdt_imp;
      tmp_val *=I_Mass;
      int jmem= LBM_glopmuti_to_imem(simu,&Jglop);
      AddLinearSolver(impsolv,imem,jmem,tmp_val);
      } // end Jglop.iv
    } // end Iglop.iv
  }
  //
}

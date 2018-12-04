#include "lbm_generic.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "implicit.h"
//
/*****************************************************************************************/
/*****************************************************************************************/
// basic moments routines
/*****************************************************************************************/
/*****************************************************************************************/
void NewMPolyDescriptor(MPolyDescriptor * mpd, int ndim, int nc,
			char *name)
{
  DestroyMPolyDescriptor(mpd);
  mpd->ndim = ndim;
  mpd->nc = nc;
  mpd->c = (schnaps_real *) calloc(nc, sizeof(schnaps_real));
  mpd->e = (int **) calloc(mpd->nc, sizeof(int *));
  for (int i = 0; i < mpd->nc; i++) {
    mpd->e[i] = (int *) calloc(ndim, sizeof(int));
  };
  for (int i = 0; i < mpd->nc; i++) {
    for (int k = 0; k < mpd->ndim; k++) {
      mpd->e[i][k] = 1;
    };
  };
  if (name != NULL) {
    mpd->name = strdup(name);
  }
}

/*****************************************************************************************/
void SetMPolyDescriptor(MPolyDescriptor * mpd, schnaps_real c_in[mpd->nc],
			int e_in[mpd->nc][mpd->ndim])
{
  for (int i = 0; i < mpd->nc; i++) {
    mpd->c[i] = c_in[i];
    for (int k = 0; k < mpd->ndim; k++) {
      mpd->e[i][k] = e_in[i][k];
    };
  };
  //
  char tmp_poly[LBM_MAX_POLYSTRSIZE];
  tmp_poly[0] = '\0';
  for (int i = 0; i < mpd->nc; i++) {
    char tmp_str[LBM_MAX_POLYSTRSIZE];
    tmp_str[0] = '\0';
    sprintf(tmp_str, "%f", mpd->c[i]);
    bool has_nonnull_exp = false;
    for (int k = 0; k < mpd->ndim; k++) {
      if (mpd->e[i][k] > 0) {
	has_nonnull_exp = true;
      }				//endif
    }				// k
    if (has_nonnull_exp) {
      strcat(tmp_str, " * ");
      for (int k = 0; k < mpd->ndim; k++) {
	if (mpd->e[i][k] > 0) {
	  char varname[2] = { LBM_POLYVARNAMES[k], '\0' };
	  strcat(tmp_str, varname);
	  char tmp_exp[LBM_MAX_POLYSTRSIZE];
	  sprintf(tmp_exp, "^%i", mpd->e[i][k]);
	  strcat(tmp_str, tmp_exp);
	}			//endif
      }				//k
    }				// endif
    int l = strlen(tmp_poly) + strlen(tmp_str) + 1;
    if (l < LBM_MAX_POLYSTRSIZE) {
      strcat(tmp_poly, tmp_str);
    }
    l = strlen(tmp_poly) + 4;
    if (l < LBM_MAX_POLYSTRSIZE) {
      if (i < mpd->nc - 1) {
	strcat(tmp_poly, " + ");
      }
    }
  }				//i
  int len = strlen(tmp_poly);
  mpd->poly_string = (char *) malloc((len + 1));
  strcpy(mpd->poly_string, tmp_poly);
}

/*****************************************************************************************/
void NewMPolyDescriptorFromModel(MPolyDescriptor * mpd,
				 const MPolyData * model, char *name)
{
  if (name == NULL) {
    NewMPolyDescriptor(mpd, model->ndim, model->nc, model->s);
  } else {
    NewMPolyDescriptor(mpd, model->ndim, model->nc, name);
  }
  SetMPolyDescriptor(mpd, model->c, model->e);
}

/*****************************************************************************************/
void DisplayMPolyDescriptor(MPolyDescriptor * mpd)
{
  if (mpd) {
    if (mpd->poly_string) {
      char *temp_str =
	  calloc(LBM_MAX_POLYSTRSIZE + LBM_MAX_MACROVARNAMESIZE + 4,
		 sizeof(char));
      strcpy(temp_str, mpd->name);
      strcat(temp_str, " : ");
      strcat(temp_str, mpd->poly_string);
      puts(temp_str);
      free(temp_str);
    }
  }
}

/*****************************************************************************************/
schnaps_real EvaluatePoly(MPolyDescriptor * mpd, schnaps_real * V)
{
  schnaps_real result = 0.0;
  for (int i = 0; i < mpd->nc; i++) {
    schnaps_real temp_val = mpd->c[i];
    for (int k = 0; k < mpd->ndim; k++) {
      for (int q = 0; q < mpd->e[i][k]; q++) {
	temp_val = temp_val * V[k];
      }				//q
    }				//k
    result = result + temp_val;
  }				//i
  return result;
}

/*****************************************************************************************/
void DestroyMPolyDescriptor(MPolyDescriptor * mpd)
{
  if (mpd->c)
    free(mpd->c);
  if (mpd->e) {
    for (int i = 0; i < mpd->nc; i++) {
      free(mpd->e[i]);
    };
    free(mpd->e);
  }
  if (mpd->name)
    free(mpd->name);
  if (mpd->poly_string)
    free(mpd->poly_string);
  mpd->ndim = 0;
  mpd->nc = 0;
}

//
/*****************************************************************************************/
/*****************************************************************************************/
// LBModelDescriptor routines
/*****************************************************************************************/
/*****************************************************************************************/
void NewLBModelDescriptor(LBModelDescriptor * lb, int d, int nb_macro,
			  int q)
{
  DestroyLBModelDescriptor(lb);
  assert(d > 0);
  lb->d = d;
  assert(q > 0);
  lb->q = q;
  assert(nb_macro > 0);
  lb->nb_macro = nb_macro;
  //
  lb->macro_names = (char **) calloc(lb->nb_macro, sizeof(char *));
  for (int i = 0; i < lb->nb_macro; i++) {
    lb->macro_names[i] =
	(char *) calloc(LBM_MAX_MACROVARNAMESIZE, sizeof(char));
  }
  //
  lb->vi = (schnaps_real **) calloc(lb->q, sizeof(schnaps_real *));
  //
  for (int i = 0; i < lb->q; i++) {
    lb->vi[i] = (schnaps_real *) calloc(lb->d, sizeof(schnaps_real));
  };
  //
  lb->iopposite = (int *) calloc(lb->q, sizeof(int));
  //
  lb->Moments = (MPolyDescriptor *) calloc(lb->q, sizeof(MPolyDescriptor));
  for (int i = 0; i < lb->q; i++) {
    lb->Moments[i] = MPolyDescriptor_NULL;
  };
  //
  lb->inode_min = (int *) calloc(lb->q, sizeof(int));
  lb->inode_max = (int *) calloc(lb->q, sizeof(int));
  //
  lb->M = (schnaps_real **) calloc(lb->q, sizeof(schnaps_real *));
  for (int i = 0; i < lb->q; i++) {
    lb->M[i] = (schnaps_real *) calloc(lb->q, sizeof(schnaps_real));
  };
  //
  MatrixStorage ms = KLU_CSR;
  Solver sv = LU;
  lb->Msolv = calloc(1,sizeof(LinearSolver));
  InitLinearSolver(lb->Msolv,lb->q,&ms,&sv);
  //
  //lb->nu = (schnaps_real *) calloc(lb->q, sizeof(schnaps_real));
  lb->s = (schnaps_real *) calloc(lb->q, sizeof(schnaps_real));
  //
/*  lb->is_relaxed= (bool*) calloc(lb->q,sizeof(schnaps_real));*/
/*  //*/
/*  lb->is_conserved= (bool*) calloc(lb->q,sizeof(schnaps_real));*/
  //

}

/*****************************************************************************************/
/*****************************************************************************************/
void DestroyLBModelDescriptor(LBModelDescriptor * lb)
{
  if (lb->macro_names) {
    for (int i = 0; i < lb->nb_macro; i++) {
      free(lb->macro_names[i]);
    }
    free(lb->macro_names);
    lb->macro_names = NULL;
  }
  if (lb->vi) {
    for (int i = 0; i < lb->q; i++) {
      free(lb->vi[i]);
    }
    free(lb->vi);
    lb->vi = NULL;
  }
  if (lb->iopposite) {
    free(lb->iopposite);
    lb->iopposite=NULL;
  }
  if (lb->Moments) {
    for (int i = 0; i < lb->q; i++) {
      DestroyMPolyDescriptor(&(lb->Moments[i]));
    }
    free(lb->Moments);
    lb->Moments=NULL;
  }
  if (lb->M) {
    for (int i = 0; i < lb->q; i++) {
      free(lb->M[i]);
    };
    free(lb->M);
    lb->M=NULL;
  }
  if (lb->Msolv){
    FreeLinearSolver(lb->Msolv);
    free(lb->Msolv);
    lb->Msolv = NULL;
  }
  if (lb->inode_min){
    free(lb->inode_min);
    lb->inode_min=NULL;
    }
  if (lb->inode_max){
    free(lb->inode_max);
    lb->inode_max=NULL;
  }
  if (lb->s){
    free(lb->s);
    lb->s=NULL;
  }
  if (lb->model_spec_params){
    free(lb->model_spec_params);
  }
  lb->model_spec_params=NULL;
  lb->d = 0;
  lb->q = 0;
  lb->nb_macro = 0;
  lb->f_to_macro = NULL;
  lb->feq = NULL;
}

/******************************************************************************************/
void ComputeLBModelDescriptorMomentMatrix(LBModelDescriptor * lb)
{
  for (int i = 0; i < lb->q; i++) {
    for (int j = lb->inode_min[i]; j <= lb->inode_max[i]; j++) {
      lb->M[i][j] = EvaluatePoly(&(lb->Moments[i]), &(lb->vi[j][0]));
      if (fabs(lb->M[i][j]) > _SMALL){
        IsNonZero(lb->Msolv,i,j);
      }
    }
  }
  //
  AllocateLinearSolver(lb->Msolv);
  for (int i = 0; i < lb->q; i++) {
    for (int j = lb->inode_min[i]; j <= lb->inode_max[i]; j++) {
      if (fabs(lb->M[i][j]) > _SMALL){
        SetLinearSolver(lb->Msolv,i,j,lb->M[i][j]);
      }
    }
  }
  LUDecompLinearSolver(lb->Msolv);
}

/*******************************************************************************************/
void DisplayLBModelDescriptorMomentPoly(LBModelDescriptor * lb)
{
  for (int i = 0; i < lb->q; i++) {
    DisplayMPolyDescriptor(&(lb->Moments[i]));
  }
}

/******************************************************************************************/
void DisplayLBModelDescriptorMomentMatrix(LBModelDescriptor * lb)
{
  for (int i = 0; i < lb->q; i++) {
    for (int j = 0; j < lb->q; j++) {
      printf("%f\t", lb->M[i][j]);
    }
    printf("\n");
  }
}

/******************************************************************************************/
void ComputeLBModelDescriptorVmax(LBModelDescriptor * lb)
{
  schnaps_real vm = 0.0;
  for (int i = 0; i < lb->q; i++) {
    schnaps_real tmpv = 0.0;
    for (int k = 0; k < lb->d; k++) {
      tmpv += lb->vi[i][k] * lb->vi[i][k];
    }
    if (tmpv > vm)
      vm = tmpv;
  }
  lb->vmax = sqrt(vm);
}

/******************************************************************************************/
void CheckLBModelDescriptorMacroConservation(LBModelDescriptor * lb,
					     bool verbose)
{
  //
  schnaps_real Min[lb->nb_macro];
  schnaps_real Mout[lb->nb_macro];
  schnaps_real Merr[lb->nb_macro];
  schnaps_real feq[lb->q];
  //
  for (int l = 0; l < lb->nb_macro; l++) {
    Min[l] = 1.0 + ((schnaps_real) rand()) / ((schnaps_real) RAND_MAX);
  }
  for (int i = 0; i < lb->q; i++) {
    feq[i] = lb->feq_i(i, lb->nb_macro, Min);
  }
  lb->f_to_macro(feq, Mout);
  schnaps_real max_error = 0.0;
  int imaxerror = -1;
  for (int l = 0; l < lb->nb_macro; l++) {
    Merr[l] = Mout[l] / Min[l] - 1.0;
    if (Merr[l] > max_error) {
      max_error = Merr[l];
      imaxerror = l;
    }
  }
  if (verbose) {
    printf("Checking lb model conservation\n");
    for (int l = 0; l < lb->nb_macro; l++) {
      printf("%i: %s : %f\n", l, lb->macro_names[l], Merr[l]);
    }
  }
  if (max_error < _SMALL) {
    printf("LB model conservation check OK \n");
  } else {
    printf
	("LB model conservation check Warning  error %f on macro quatity %s \n",
	 max_error, lb->macro_names[imaxerror]);
  }
  //
}
/*****************************************************************************************/
void CheckLBMMomentMatrixInversion(LBModelDescriptor *lb, bool verbose){
  schnaps_real fin[lb->q];
  schnaps_real ferr[lb->q];
  //
  for (int l = 0; l < lb->q; l++) {
    fin[l] = 1.0 + ((schnaps_real) rand()) / ((schnaps_real) RAND_MAX);
  }
  //
  MatVect(lb->Msolv,fin,lb->Msolv->rhs);
  //
  if (verbose){
    DisplayLinearSolver(lb->Msolv);
  }
  SolveLinearSolver(lb->Msolv);
  //
  schnaps_real max_error = 0.0;
  int imaxerror = -1;
  for (int i=0;i< lb->q;i++){
    ferr[i] = lb->Msolv->sol[i]/fin[i]-1.0;
    if (ferr[i] > max_error){
      max_error= ferr[i];
      imaxerror = i;
    }
  }
  if (verbose){
    for (int i=0;i<lb->q;i++){
      printf("f%i orig :\t %f \t comp:\t %f \t error:\t %f\n",i,fin[i],lb->Msolv->sol[i],ferr[i]);
    }
  }
  if (max_error < _SMALL) {
    printf("LB model moment matrix inversion check OK \n");
  } else {
    printf
	("LB model matrix inversion check Warning  max error %f on components f%i \n",
	 max_error, imaxerror);
  }
  
}
/*****************************************************************************************/
void LBM_Dummy_InitMacroData(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->nb_macro; i++) {
    w[i] = 0.0;
  }
}

//
void LBM_Dummy_InitData(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->q; i++) {
    w[i] = 0.0;
  }
}

/*****************************************************************************************/
void LBM_Dummy_InitData_OneNode(schnaps_real x[3], schnaps_real w[])
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  w[0] = 0.0;
}

/*****************************************************************************************/
/*****************************************************************************************/
// Lattice Boltmann Simulation objects routines
/*****************************************************************************************/
/*****************************************************************************************/
void InitLBMSimulation(LBMSimulation * lbsimu,
		       LatticeBoltzmannSimData * lsd, MacroMesh * mesh,
		       int deg[3], int raf[3])
{
  lbsimu->d = lsd->lb_model->d;
  lbsimu->nb_macro_fields = lsd->lb_model->nb_macro;
  lbsimu->q = lsd->lb_model->q;
  lbsimu->vmax=lsd->lb_model->vmax;
  //
  lbsimu->macro_model.m = lsd->lb_model->nb_macro;
  lbsimu->macro_model.NumFlux = NULL;
  if (lbsimu->macro_model.InitData == NULL) {
    lbsimu->macro_model.InitData = LBM_Dummy_InitMacroData;
  }
  //
  EmptySimulation(&(lbsimu->macro_simu));
  InitSimulation(&(lbsimu->macro_simu), mesh, deg, raf,
		 &(lbsimu->macro_model));
  lbsimu->macro_simu.pre_dtfields = NULL;
  lbsimu->macro_simu.post_dtfields = NULL;
  lbsimu->macro_simu.update_after_rk = NULL;
  //
  lbsimu->micro_model.m = lsd->lb_model->q;
  //lbsimu->micro_model.NumFlux = NULL;
  //lbsimu->micro_model.InitData = LBM_Dummy_InitData;
  //lbsimu->micro_model.ImposedData = NULL;
  //lbsimu->micro_model.BoundaryFlux = NULL;
  //lbsimu->micro_model.Source = NULL;
  EmptySimulation(&(lbsimu->micro_simu));
  InitSimulation(&(lbsimu->micro_simu), mesh, deg, raf,
		 &(lbsimu->micro_model));
  lbsimu->micro_simu.pre_dtfields = LBM_pre_dtfields_wrapper;
  lbsimu->micro_simu.post_dtfields = LBM_post_dtfields_wrapper;
  lbsimu->micro_simu.update_after_rk = LBM_update_after_rk_wrapper;
  // set initial data to equilibrium one computed with initial macro data
  LB_Relaxation_bgk_f_full(lbsimu);
  lsd->current_lb_sim = lbsimu;
  //
  lbsimu->wmic_buffer =
      (schnaps_real *) calloc(lsd->lb_model->q, sizeof(schnaps_real));
  lbsimu->wmac_buffer =
      (schnaps_real *) calloc(lsd->lb_model->nb_macro, sizeof(schnaps_real));
  //
  lbsimu->pre_advec = NULL;
  lbsimu->post_advec = NULL;
  lbsimu->post_tstep = NULL;
  //
  lbsimu->diag_2d_period = 0.0;
  lbsimu->collect_diags = NULL;
}

void FreeBMSimulation(LBMSimulation * lbsimu)
{
  freeSimulation(&(lbsimu->macro_simu));
  freeSimulation(&(lbsimu->micro_simu));
  free(lbsimu->wmic_buffer);
  free(lbsimu->wmac_buffer);
}

//******************************************************************************//

void LB_Relaxation_bgk_f(void *lbs)
{
  LBMSimulation *lbsimu = lbs;
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  Simulation *macsimu = &(lbsimu->macro_simu);
  Simulation *micsimu = &(lbsimu->micro_simu);
  assert(lsd->lb_model->feq);
  schnaps_real rate = lsd->lb_model->s[0];	// TODO use something simpler ?
  //
  for (int ie = 0; ie < macsimu->macromesh.nbelems; ++ie) {
    field *f = macsimu->fd + ie;
    field *fmic = micsimu->fd + ie;
    for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      // recover macro quantities;
      schnaps_real wmac[lsd->lb_model->nb_macro];
      for (int iv = 0; iv < lsd->lb_model->nb_macro; iv++) {
	int imac = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	wmac[iv] = f->wn[imac];
      }				//iv macro quantities
      // now compute equilibria and relax for all micro quantities
      for (int inode = 0; inode < lsd->lb_model->q; inode++) {
	schnaps_real feq =
	    lsd->lb_model->feq_i(inode, lsd->lb_model->nb_macro, wmac);
	int imic =
	    fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipg,
			   inode);
	fmic->wn[imic] = fmic->wn[imic] - rate * (fmic->wn[imic] - feq);
      }
    }				//ipg
    //
  }				//ie
}

void LB_Relaxation_bgk_f_full(void *lbs)
{
  LBMSimulation *lbsimu = lbs;
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  assert(lsd->lb_model->feq);
  Simulation *macsimu = &(lbsimu->macro_simu);
  Simulation *micsimu = &(lbsimu->micro_simu);
  //
  for (int ie = 0; ie < macsimu->macromesh.nbelems; ++ie) {
    field *f = macsimu->fd + ie;
    field *fmic = micsimu->fd + ie;
    for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      // recover macro quantities;
      schnaps_real wmac[lsd->lb_model->nb_macro];
      for (int iv = 0; iv < lsd->lb_model->nb_macro; iv++) {
	int imac = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	wmac[iv] = f->wn[imac];
      }				//iv macro quantities
      // now compute equilibria and relax for all micro quantities
      for (int inode = 0; inode < lsd->lb_model->q; inode++) {
	int imic =
	    fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipg,
			   inode);
	fmic->wn[imic] =
	    lsd->lb_model->feq_i(inode, lsd->lb_model->nb_macro, wmac);
      }
    }				//ipg
    //
  }				//ie
}


void LB_ComputeMacroFromMicro(void *lbs)
{
  LBMSimulation *lbsimu = lbs;
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  Simulation *macsimu = &(lbsimu->macro_simu);
  Simulation *micsimu = &(lbsimu->micro_simu);
  for (int ie = 0; ie < macsimu->macromesh.nbelems; ++ie) {
    field *f = macsimu->fd + ie;
    field *fmic = micsimu->fd + ie;
    for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real wmac[lsd->lb_model->nb_macro];
      schnaps_real wmic[lsd->lb_model->q];
      for (int inode = 0; inode < lsd->lb_model->q; inode++) {
	int imic =
	    fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipg,
			   inode);
	wmic[inode] = fmic->wn[imic];
      }				//inode
      lsd->lb_model->f_to_macro(wmic, wmac);
      for (int iv = 0; iv < lsd->lb_model->nb_macro; iv++) {
	int imac = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	f->wn[imac] = wmac[iv];
      };
    }				//ipg
  }				//ie
}

/***********************************************************************************/
// wrappers for compatibility with RK schemes
void LBM_pre_dtfields_wrapper(void *simu)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  LBMSimulation *lbsimu = lsd->current_lb_sim;
  printf("pre dt wrapper");
  // update time data 
  lbsimu->tnow = lbsimu->micro_simu.tnow;
  lbsimu->iter_time = lbsimu->micro_simu.iter_time_rk;
  lbsimu->itermax = lbsimu->micro_simu.itermax_rk;
  // redirect call to global function
  lbsimu->pre_advec(lbsimu);
}

void LBM_post_dtfields_wrapper(void *simu)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  LBMSimulation *lbsimu = lsd->current_lb_sim;
  printf("post dt wrapper");
  // update time data 
  lbsimu->tnow = lbsimu->micro_simu.tnow;
  lbsimu->iter_time = lbsimu->micro_simu.iter_time_rk;
  lbsimu->itermax = lbsimu->micro_simu.itermax_rk;
  // redirect call to global function
  lbsimu->post_advec(lbsimu);
}

void LBM_update_after_rk_wrapper(void *simu, schnaps_real * w)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  LBMSimulation *lbsimu = lsd->current_lb_sim;
  printf("update_after_rk dt wrapper");
  // update time data
  lbsimu->tnow = lbsimu->micro_simu.tnow;
  lbsimu->iter_time = lbsimu->micro_simu.iter_time_rk;
  lbsimu->itermax = lbsimu->micro_simu.itermax_rk;
  // redirect call to global function
  lbsimu->post_tstep(lbsimu, w);
}

//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
// ******************flux functions for LB models ******************************//
//******************************************************************************//
//******************************************************************************//

void LBM_OneLatticeNumFlux(schnaps_real * wL, schnaps_real * wR,
			   schnaps_real * vnorm, schnaps_real * flux)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->q; i++) {
    schnaps_real vn = 0;
    for (int dim = 0; dim < lsd->lb_model->d; dim++) {
      vn += lsd->lb_model->vi[i][dim] * vnorm[dim];
    }
    schnaps_real vnp = vn > 0 ? vn : 0;
    schnaps_real vnm = vn - vnp;
    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
}
void LBM_OneLatticeNumFlux_SelectVel(schnaps_real * wL, schnaps_real * wR,
			   schnaps_real * vnorm, schnaps_real * flux)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->q; i++) {
    if (i==lsd->current_advec_node_index){
      schnaps_real vn = 0;
      for (int dim = 0; dim < lsd->lb_model->d; dim++) {
        vn += lsd->lb_model->vi[i][dim] * vnorm[dim];
      }
      schnaps_real vnp = vn > 0 ? vn : 0;
      schnaps_real vnm = vn - vnp;
      flux[i] = vnp * wL[i] + vnm * wR[i];
    }
    else
    {
      flux[i]=0.0;
    }
  }
  //printf("flux %f \t %f \t %f \n",flux[0],flux[1],flux[2]);
}

void LBM_OneNodeNumFlux(schnaps_real * wL, schnaps_real * wR,
			schnaps_real * vnorm, schnaps_real * flux)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  int inode = lsd->current_node_index;
  schnaps_real vn = 0;
  for (int dim = 0; dim < lsd->lb_model->d; dim++) {
    vn += lsd->lb_model->vi[inode][dim] * vnorm[dim];
  }
  schnaps_real vnp = vn > 0 ? vn : 0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
}
/***********************************************************************************/
// Basic sources
/*************************************************************************************/
void LBM_BGK_Relaxation_Source(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *s){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real srate = lsd->lb_model->s[0];
  //
  int nb_macro= lsd->lb_model->nb_macro;
  int q = lsd->lb_model->q;
  schnaps_real wmac[nb_macro];
  schnaps_real feq[q];
  lsd->lb_model->f_to_macro(w,wmac);
  lsd->lb_model->feq(q,nb_macro,wmac,feq);
  for (int i=0;i<q;i++){
    s[i] = srate * (feq[i]-w[i]);
  }
}
void LBM_BGK_Relaxation_Source_Jacobian(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *J){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real srate = lsd->lb_model->s[0];
  //
  int nb_macro= lsd->lb_model->nb_macro;
  int q = lsd->lb_model->q;
  schnaps_real wmac[nb_macro];
  schnaps_real Jfeq[q*q];
  lsd->lb_model->f_to_macro(w,wmac);
  lsd->lb_model->Jfeq(q,nb_macro,wmac,Jfeq);
  
  for (int k=0; k < q * q;k++){
      J[k] = srate * Jfeq[k];
  }
  // remove diagonal (linear part)
  for (int i=0;i< q ; i++){
    J[i*(q+1)] -= srate;  
  }
}
//*********************************************************************************//
void LBM_BGK_Relaxation_Source_TestDecay(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *s){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real srate = lsd->lb_model->s[0];
  //
  int q = lsd->lb_model->q;
  for (int i=0;i<q;i++){
    s[i] = -srate * w[i];
  }
}
void LBM_BGK_Relaxation_Source_Jacobian_TestDecay(const schnaps_real x[3], const schnaps_real t, const schnaps_real *w,schnaps_real *J){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real srate = lsd->lb_model->s[0];
  int q = lsd->lb_model->q;
  for (int k=0; k < q * q;k++){
      J[k] = 0.0;
  }
  // remove diagonal (linear part)
  for (int i=0;i< q ; i++){
    J[i*(q+1)] = -srate;  
  }
}
//*********************************************************************************//
//*********************************************************************************//
//********************************************************************************//
//********************************************************************************//
// per model specific routines and functions 
// TODO export to different file
//********************************************************************************//
//********************************************************************************//
//********************************************************************************//
schnaps_real LBM_dummy_zero_feq_i(int inode, int nb_macro, schnaps_real * w)
{
  return 0.0;
}
void LBM_dummy_zero_feq(int q,int nb_macro,schnaps_real *w, schnaps_real *f){
  for (int i=0;i< q;i++){
    f[i]=0.0;
  }
}
schnaps_real LBM_dummy_zero_Jfeq_ij(int inode, int jnode, int nb_macro,schnaps_real *w){
  return 0.0;
}
void LBM_dummy_zero_Jfeq(int q,int nb_macro, schnaps_real *w, schnaps_real *J){
  for (int i=0;i < q * q; i++){
    J[i] = 0.0;
  }
}
// D2Q3 TEST artificial model for debugging 
void LBM_Set_D2Q3_TEST_model(LBModelDescriptor *lb,schnaps_real alpha,schnaps_real cref){
  assert(lb->d=2);
  assert(lb->q=3);
  assert(lb->nb_macro=2);
  // 
  schnaps_real ca = cos(alpha * M_PI /180.0);
  schnaps_real sa = sin(alpha * M_PI /180.0);
  //
  NewMPolyDescriptorFromModel(&(lb->Moments[0]), &model_rho2D, NULL);
  NewMPolyDescriptor(&(lb->Moments[1]),lb->d,2,"jpar");
  schnaps_real c_in[2]= {ca,sa};
  int e_in[2][2];
  e_in[0][0]= 1;
  e_in[0][1]= 0;
  e_in[1][0]= 0;
  e_in[1][1]= 1;
  SetMPolyDescriptor(&(lb->Moments[1]),c_in,e_in);
  //
  NewMPolyDescriptorFromModel(&(lb->Moments[2]), &model_ekin2D, NULL);   
  sprintf(lb->macro_names[0], "rho");
  sprintf(lb->macro_names[0], "uparallel");
  // set velocity nodes
  lb->cref = cref;
  for (int i = 0; i < lb->q; i++) {
      lb->vi[i][0] = lb->cref * ca * LBM_D1Q3_nodes[i][0];
      lb->vi[i][1] = lb->cref * sa * LBM_D1Q3_nodes[i][0];
  }
  for (int i = 0; i < lb->q; i++) {
    lb->iopposite[i] = LBM_D1Q3_iopposite[i];
  }
  ComputeLBModelDescriptorVmax(lb);
  //
  for (int i = 0; i < lb->q; i++) {
    lb->inode_min[i] = 0;
    lb->inode_max[i] = lb->q - 1;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq_i = &LBM_feq_i_D2Q3_TEST;
  lb->feq = &LBM_feq_D2Q3_TEST;
  lb->Jfeq_ij= &LBM_dummy_zero_Jfeq_ij; // <-TODO actual implementation
  lb->Jfeq = &LBM_dummy_zero_Jfeq;  // <-TODO actual implementation
  lb->f_to_macro = &LBM_f_to_macro_D2Q3_TEST;
  //
  lb->model_spec_params = (params_D2Q3_TEST*) calloc(1,sizeof(params_D2Q3_TEST));
  params_D2Q3_TEST *p=lb->model_spec_params;
  p->alpha_deg=alpha;
  p->cosa = ca;
  p->sina = sa;
  p->theta=cref * cref /3.0;
  p->invtheta = 1.0 /p->theta;
}
void LBM_f_to_macro_D2Q3_TEST(const schnaps_real *f,schnaps_real *w){
  w[0]=0.0;
  w[1]=0.0;
  for (int i=0;i<3;i++){
   w[0]+=f[i];
   w[1]+=f[i] * LBM_D1Q3_nodes[i][0];
  }  
  // normalize
  if (w[0] != 0.0){
  w[1]=w[1]/w[0];
  }
}
schnaps_real LBM_feq_i_D2Q3_TEST(int inode,int nb_macro, schnaps_real *w){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real uparallel = w[1];
  params_D2Q3_TEST *param=lsd->lb_model->model_spec_params;
  schnaps_real invtemp = param->invtheta;
  schnaps_real ux = uparallel * param->cosa;
  schnaps_real uy = uparallel * param->sina;
  schnaps_real uv =
      (ux * lsd->lb_model->vi[inode][0] +
       uy * lsd->lb_model->vi[inode][1]) * invtemp;
  //
  schnaps_real feq= LBM_WEIGHTS_D1Q3_ISOTH[inode] * rho * (1.+ uv);
  return feq;
}
void LBM_feq_D2Q3_TEST(int q,int nb_macro,schnaps_real *w,schnaps_real *f){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  params_D2Q3_TEST *param=lsd->lb_model->model_spec_params;
  schnaps_real invtemp = param->invtheta;
  schnaps_real rho = w[0];
  schnaps_real uparallel = w[1];
  schnaps_real ux = uparallel * param->cosa;
  schnaps_real uy = uparallel * param->sina;
  for (int inode=0;inode < q;inode++){
    schnaps_real uv =
      (ux * lsd->lb_model->vi[inode][0] +
       uy * lsd->lb_model->vi[inode][1]) * invtemp;
    f[inode]=LBM_WEIGHTS_D1Q3_ISOTH[inode] * rho * (1.+ uv);
  }
}
// D2Q9 isothermal
//********************************************************************************//
void LBM_Set_D2Q9_ISOTH_model(LBModelDescriptor * lb, schnaps_real cref)
{
  assert(lb->d == 2);
  assert(lb->q == 9);
  assert(lb->nb_macro == 3);
  NewMPolyDescriptorFromModel(&(lb->Moments[0]), &model_rho2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[1]), &model_jx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[2]), &model_jy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[3]), &model_trMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[4]), &model_DifMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[5]), &model_Mxy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[6]), &model_Mabcx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[7]), &model_Mabcy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[8]), &model_V42D, NULL);
  //
  sprintf(lb->macro_names[0], "rho");
  sprintf(lb->macro_names[1], "ux");
  sprintf(lb->macro_names[2], "uy");
  // set velocity nodes
  lb->cref = cref;
  for (int i = 0; i < lb->q; i++) {
    for (int j = 0; j < lb->d; j++) {
      lb->vi[i][j] = lb->cref * LBM_D2Q9_nodes[i][j];
    }
  }
  for (int i = 0; i < lb->q; i++) {
    lb->iopposite[i] = LBM_D2Q9_iopposite[i];
  }
  ComputeLBModelDescriptorVmax(lb);
  //
  for (int i = 0; i < lb->q; i++) {
    lb->inode_min[i] = 0;
    lb->inode_max[i] = lb->q - 1;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq_i = &LBM_feq_i_D2Q9_ISOTH;
  lb->feq = &LBM_feq_D2Q9_ISOTH;
  lb->Jfeq_ij= &LBM_dummy_zero_Jfeq_ij; // <-TODO actual implementation
  lb->Jfeq = &LBM_dummy_zero_Jfeq;  // <-TODO actual implementation
  lb->f_to_macro = &LBM_f_to_macro_D2Q9_ISOTH;
  //
  lb->model_spec_params = (params_D2Q9_ISOTH*) calloc(1,sizeof(params_D2Q9_ISOTH));
  params_D2Q9_ISOTH *p=lb->model_spec_params;
  p->theta=cref * cref /3.0;
  p->invtheta = 1.0 /p->theta;
}

//
void LBM_f_to_macro_D2Q9_ISOTH(const schnaps_real * f, schnaps_real * w)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->nb_macro; i++) {
    w[i] = 0.0;
    for (int inode = 0; inode < lsd->lb_model->q; inode++) {
      w[i] += f[inode] * lsd->lb_model->M[i][inode];
    }
  }
  // normalize
  if (w[0] != 0.0) {
    w[1] = w[1] / w[0];
    w[2] = w[2] / w[0];
  }
}

//
schnaps_real LBM_feq_i_D2Q9_ISOTH(int i_node, int nb_macro, schnaps_real * w)
{
  //
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real ux = w[1];
  schnaps_real uy = w[2];
  params_D2Q9_ISOTH *param=lsd->lb_model->model_spec_params;
  schnaps_real invtemp = param->invtheta;
  //
  schnaps_real u2 = (ux * ux + uy * uy) * invtemp;
  schnaps_real uv =
      (ux * lsd->lb_model->vi[i_node][0] +
       uy * lsd->lb_model->vi[i_node][1]) * invtemp;
  schnaps_real feq =
      LBM_WEIGHTS_D2Q9_ISOTH[i_node] * rho * (1.0 + uv +
					      0.5 * (uv * uv - u2));
  return feq;
};
void LBM_feq_D2Q9_ISOTH(int q,int nb_macro, schnaps_real * w,schnaps_real *f)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real ux = w[1];
  schnaps_real uy = w[2];
  params_D2Q9_ISOTH *param=lsd->lb_model->model_spec_params;
  schnaps_real invtemp = param->invtheta;
  //
  schnaps_real u2 = (ux * ux + uy * uy) * invtemp;
  for (int i_node=0;i_node < q;i_node++){
  schnaps_real uv =
      (ux * lsd->lb_model->vi[i_node][0] +
       uy * lsd->lb_model->vi[i_node][1]) * invtemp;
  f[i_node] = LBM_WEIGHTS_D2Q9_ISOTH[i_node] * rho * (1.0 + uv + 0.5 * (uv * uv - u2));
  }
};

// D2Q9 isothermal with incompressibility mod , defined by parameter rho0
//********************************************************************************//
void LBM_Set_D2Q9_ISOTH_INC_model(LBModelDescriptor * lb, schnaps_real cref)
{
  assert(lb->d == 2);
  assert(lb->q == 9);
  assert(lb->nb_macro == 3);
  NewMPolyDescriptorFromModel(&(lb->Moments[0]), &model_rho2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[1]), &model_jx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[2]), &model_jy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[3]), &model_trMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[4]), &model_DifMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[5]), &model_Mxy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[6]), &model_Mabcx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[7]), &model_Mabcy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[8]), &model_V42D, NULL);
  //
  sprintf(lb->macro_names[0], "rho");
  sprintf(lb->macro_names[1], "ux");
  sprintf(lb->macro_names[2], "uy");
  // set velocity nodes
  lb->cref = cref;
  for (int i = 0; i < lb->q; i++) {
    for (int j = 0; j < lb->d; j++) {
      lb->vi[i][j] = lb->cref * LBM_D2Q9_nodes[i][j];
    }
  }
  for (int i = 0; i < lb->q; i++) {
    lb->iopposite[i] = LBM_D2Q9_iopposite[i];
  }
  ComputeLBModelDescriptorVmax(lb);
  //
  for (int i = 0; i < lb->q; i++) {
    lb->inode_min[i] = 0;
    lb->inode_max[i] = lb->q - 1;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq_i = &LBM_feq_i_D2Q9_ISOTH_INC;
  lb->feq = &LBM_feq_D2Q9_ISOTH_INC;
  lb->Jfeq_ij= &LBM_dummy_zero_Jfeq_ij; // <-TODO actual implementation
  lb->Jfeq = &LBM_dummy_zero_Jfeq;  // <-TODO actual implementation
  lb->f_to_macro = &LBM_f_to_macro_D2Q9_ISOTH_INC;
  //
  lb->model_spec_params =  (params_D2Q9_ISOTH_INC*) calloc(1,sizeof(params_D2Q9_ISOTH_INC));
  params_D2Q9_ISOTH_INC *p=lb->model_spec_params;
  p->theta=cref * cref /3.0;
  p->invtheta = 1.0 / p->theta;
  p->rho0 =1.0;
  //
}

//
void LBM_f_to_macro_D2Q9_ISOTH_INC(const schnaps_real * f, schnaps_real * w)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->nb_macro; i++) {
    w[i] = 0.0;
    for (int inode = 0; inode < lsd->lb_model->q; inode++) {
      w[i] += f[inode] * lsd->lb_model->M[i][inode];
    }
  }
  params_D2Q9_ISOTH_INC *param = lsd->lb_model->model_spec_params;
  schnaps_real rho0=param->rho0; //  
  // normalize
  if (w[0] != 0.0) {
    w[1] = w[1] / rho0;
    w[2] = w[2] / rho0;
  }
}
//
schnaps_real LBM_feq_i_D2Q9_ISOTH_INC(int i_node, int nb_macro, schnaps_real * w)
{
  //
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real ux = w[1];
  schnaps_real uy = w[2];
  params_D2Q9_ISOTH_INC *param = lsd->lb_model->model_spec_params;
  schnaps_real invtemp = param->invtheta;
  schnaps_real rho0 = param->rho0;
  //
  schnaps_real u2 = (ux * ux + uy * uy) * invtemp;
  schnaps_real uv =
      (ux * lsd->lb_model->vi[i_node][0] +
       uy * lsd->lb_model->vi[i_node][1]) * invtemp;
  schnaps_real feq =
      LBM_WEIGHTS_D2Q9_ISOTH[i_node] * (rho + rho0 * ( uv +
					      0.5 * (uv * uv - u2)));
  return feq;
};
void LBM_feq_D2Q9_ISOTH_INC(int q, int nb_macro, schnaps_real * w,schnaps_real *f){
  //
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real ux = w[1];
  schnaps_real uy = w[2];
  params_D2Q9_ISOTH_INC *param = lsd->lb_model->model_spec_params;
  schnaps_real invtemp = param->invtheta;
  schnaps_real rho0 = param->rho0;
  //
  schnaps_real u2 = (ux * ux + uy * uy) * invtemp;
  for (int i_node=0;i_node < q;i_node++){
    schnaps_real uv =
        (ux * lsd->lb_model->vi[i_node][0] +
         uy * lsd->lb_model->vi[i_node][1]) * invtemp;
    f[i_node] =
        LBM_WEIGHTS_D2Q9_ISOTH[i_node] * (rho + rho0 * ( uv +
					        0.5 * (uv * uv - u2)));
  }
};
//

// D2Q9 isothermal linearized (2D wave equation)
void LBM_Set_D2Q9_ISOTH_LINEARIZED_model(LBModelDescriptor * lb,
					 schnaps_real cref)
{
  assert(lb->d == 2);
  assert(lb->q == 9);
  assert(lb->nb_macro == 3);
  NewMPolyDescriptorFromModel(&(lb->Moments[0]), &model_rho2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[1]), &model_jx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[2]), &model_jy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[3]), &model_trMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[4]), &model_DifMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[5]), &model_Mxy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[6]), &model_Mabcx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[7]), &model_Mabcy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[8]), &model_V42D, NULL);
  //
  sprintf(lb->macro_names[0], "rho");
  sprintf(lb->macro_names[1], "jx");
  sprintf(lb->macro_names[2], "jy");
  // set velocity nodes
  lb->cref = cref;
  for (int i = 0; i < lb->q; i++) {
    for (int j = 0; j < lb->d; j++) {
      lb->vi[i][j] = lb->cref * LBM_D2Q9_nodes[i][j];
    }
  }
  for (int i = 0; i < lb->q; i++) {
    lb->iopposite[i] = LBM_D2Q9_iopposite[i];
  }
  ComputeLBModelDescriptorVmax(lb);
  //
  for (int i = 0; i < lb->q; i++) {
    lb->inode_min[i] = 0;
    lb->inode_max[i] = lb->q - 1;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq_i = &LBM_feq_i_D2Q9_ISOTH_LINEARIZED;
  lb->feq = &LBM_feq_D2Q9_ISOTH_LINEARIZED;
  lb->Jfeq_ij= &LBM_dummy_zero_Jfeq_ij; // <-TODO actual implementation
  lb->Jfeq = &LBM_dummy_zero_Jfeq;  // <-TODO actual implementation
  lb->f_to_macro = &LBM_f_to_macro_D2Q9_ISOTH_LINEARIZED;
  //
  lb->model_spec_params = (params_D2Q9_ISOTH_LINEARIZED*) calloc(1,sizeof(params_D2Q9_ISOTH_LINEARIZED));
  params_D2Q9_ISOTH_LINEARIZED *p=lb->model_spec_params;
  p->theta= cref * cref / 3.0; 
  p->invtheta = 1.0/p->theta;
}

schnaps_real LBM_feq_i_D2Q9_ISOTH_LINEARIZED(int i_node, int nb_macro,
					   schnaps_real * w)
{
  //
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real jx = w[1];
  schnaps_real jy = w[2];
  params_D2Q9_ISOTH_LINEARIZED *p=lsd->lb_model->model_spec_params;
  schnaps_real invtemp = p->invtheta;
  //
  schnaps_real jv =
      (jx * lsd->lb_model->vi[i_node][0] +
       jy * lsd->lb_model->vi[i_node][1]) * invtemp;
  schnaps_real feq = LBM_WEIGHTS_D2Q9_ISOTH[i_node] * (rho + jv);
  return feq;
};
void LBM_feq_D2Q9_ISOTH_LINEARIZED(int q, int nb_macro,
					   schnaps_real * w,schnaps_real *f)
{
  //
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  schnaps_real rho = w[0];
  schnaps_real jx = w[1];
  schnaps_real jy = w[2];
  params_D2Q9_ISOTH_LINEARIZED *p=lsd->lb_model->model_spec_params;
  schnaps_real invtemp = p->invtheta;
  //
  for (int i_node=0;i_node < q;i_node ++){
  schnaps_real jv =
      (jx * lsd->lb_model->vi[i_node][0] +
       jy * lsd->lb_model->vi[i_node][1]) * invtemp;
  f[i_node] = LBM_WEIGHTS_D2Q9_ISOTH[i_node] * (rho + jv);
  }
};

void LBM_f_to_macro_D2Q9_ISOTH_LINEARIZED(const schnaps_real * f,
					  schnaps_real * w)
{
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < lsd->lb_model->nb_macro; i++) {
    w[i] = 0.0;
    for (int inode = 0; inode < lsd->lb_model->q; inode++) {
      w[i] += f[inode] * lsd->lb_model->M[i][inode];
    }
  }
}
/******************************************************************************/
// MHD models
void LBM_Set_MHD_D2Q9_2D2Q5_model(LBModelDescriptor * lb, schnaps_real cref){
  assert(lb->d==2);
  assert(lb->q==19);
  assert(lb->nb_macro==5);
  // hydrodynamic part
  NewMPolyDescriptorFromModel(&(lb->Moments[0]), &model_rho2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[1]), &model_jx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[2]), &model_jy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[3]), &model_trMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[4]), &model_DifMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[5]), &model_Mxy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[6]), &model_Mabcx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[7]), &model_Mabcy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[8]), &model_V42D, NULL);
  // magnetic part
  // bx
  NewMPolyDescriptorFromModel(&(lb->Moments[9]), &model_rho2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[10]), &model_jx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[11]), &model_jy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[12]), &model_trMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[13]), &model_DifMab2D, NULL);
  // by
  NewMPolyDescriptorFromModel(&(lb->Moments[14]), &model_rho2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[15]), &model_jx2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[16]), &model_jy2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[17]), &model_trMab2D, NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[18]), &model_DifMab2D, NULL);
  //
  sprintf(lb->macro_names[0], "rho");
  sprintf(lb->macro_names[1], "jx");
  sprintf(lb->macro_names[2], "jy");
  sprintf(lb->macro_names[3], "Bx");
  sprintf(lb->macro_names[4], "By");
  //
  lb->cref = cref;
  // hydrodynamic nodes
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < lb->d; j++) {
      lb->vi[i][j] = lb->cref * LBM_D2Q9_nodes[i][j];
    }
    lb->iopposite[i] = LBM_D2Q9_iopposite[i];
  }
  // magnetic nodes
  for (int ibcomp=0;ibcomp<2;ibcomp++){
    for (int i=0;i<5;i++){
      int inode=9+ibcomp * 5 + i;
      for (int j=0;j<lb->d;j++){
      lb->vi[inode][j] = lb->cref * LBM_D2Q5_nodes[i][j];
      }
      lb->iopposite[inode]= 9+ibcomp * 5 + LBM_D2Q5_iopposite[i];
    }
  }
  //
  ComputeLBModelDescriptorVmax(lb);
  for (int i = 0; i < 9; i++) {
    lb->inode_min[i] = 0;
    lb->inode_max[i] = 8;
  }
  // Bx
  for (int i = 9; i < 14; i++) {
    lb->inode_min[i] = 9;
    lb->inode_max[i] = 13;
  }
  // By
  for (int i = 14; i < 19; i++) {
    lb->inode_min[i] = 14;
    lb->inode_max[i] = 18;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq_i = &LBM_feq_i_MHD_D2Q9_2D2Q5;
  lb->feq = &LBM_feq_MHD_D2Q9_2D2Q5;
  lb->Jfeq_ij= &LBM_dummy_zero_Jfeq_ij; // <-TODO actual implementation
  lb->Jfeq = &LBM_dummy_zero_Jfeq;  // <-TODO actual implementation
  lb->f_to_macro = &LBM_f_to_macro_MHD_D2Q9_2D2Q5;
  //
  lb->model_spec_params= (params_MHD_D2Q9_2D2Q5*) calloc(1,sizeof(params_MHD_D2Q9_2D2Q5));
  params_MHD_D2Q9_2D2Q5 *p=lb->model_spec_params;
  p->theta= cref * cref/3.0;
  p->invtheta = 1.0 / p->theta;
  p->theta_mag = cref * cref /3.0;
  p->invtheta_mag = 1.0 / p->theta_mag;
}
//
void LBM_f_to_macro_MHD_D2Q9_2D2Q5(const schnaps_real *f, schnaps_real *w){
  // hydrodynamic quatities
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  for (int i = 0; i < 3; i++) {
    w[i] = 0.0;
    for (int inode = 0; inode < lsd->lb_model->q; inode++) {
      w[i] += f[inode] * lsd->lb_model->M[i][inode];
    }
  }
  // normalize
  if (w[0] != 0.0) {
    w[1] = w[1] / w[0];
    w[2] = w[2] / w[0];
  }
  // magnetic quantities
  w[3]=0.0;
  for (int j=9;j<14;j++){
    w[3] += f[j] * lsd->lb_model->M[9][j]; 
  }
  w[4]=0.0;
  for (int j=14;j<19;j++){
    w[4] += f[j] * lsd->lb_model->M[14][j]; 
  }
}
schnaps_real LBM_feq_i_MHD_D2Q9_2D2Q5(int inode, int nb_macro, schnaps_real *w){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  params_MHD_D2Q9_2D2Q5 *p=lsd->lb_model->model_spec_params;
  schnaps_real rho = w[0];
  schnaps_real ux = w[1];
  schnaps_real uy = w[2];
  schnaps_real Bx = w[3];
  schnaps_real By = w[4];
  schnaps_real feq = -1;
  if (inode < 9){
    schnaps_real invtemp = 1.0/p->theta;
    schnaps_real vx = lsd->lb_model->vi[inode][0];
    schnaps_real vy = lsd->lb_model->vi[inode][1];
    //
    schnaps_real u2 = (ux * ux + uy * uy) * invtemp;
    schnaps_real v2 = (vx * vx + vy * vy) * invtemp;
    schnaps_real uv = (ux * vx + uy * vy) * invtemp;
    schnaps_real vB=  (vx * Bx + vy * By) * invtemp;
    schnaps_real B2= (Bx * Bx + By * By) * invtemp;
    feq = LBM_WEIGHTS_D2Q9_ISOTH[inode] *
      (rho * (1.0 + uv + 0.5 * (uv * uv - u2))
       + 0.5 * (0.5 * B2 * v2 - vB * vB));
  }
  else if (inode < 14){
    int inodeloc=inode-9;
    schnaps_real invtemp = 1.0/p->theta_mag;
    schnaps_real vy=lsd->lb_model->vi[inode][1];
    feq = LBM_WEIGHTS_D2Q5_ISOTH[inodeloc] *
      (Bx + invtemp * vy *(uy * Bx - ux * By));
  }
  else {
    // inode > 13
    int inodeloc=inode-14;
    schnaps_real invtemp = 1.0/p->theta_mag;
    schnaps_real vx=lsd->lb_model->vi[inode][0];
    feq = LBM_WEIGHTS_D2Q5_ISOTH[inodeloc] *
      (By + invtemp * vx *(ux * By - uy * Bx));
  }
  return feq;
}

void LBM_feq_MHD_D2Q9_2D2Q5(int q, int nb_macro, schnaps_real *w,schnaps_real *f){
  LatticeBoltzmannSimData *lsd = &schnaps_lbm_simdata;
  params_MHD_D2Q9_2D2Q5 *p=lsd->lb_model->model_spec_params;
  schnaps_real rho = w[0];
  schnaps_real ux = w[1];
  schnaps_real uy = w[2];
  schnaps_real Bx = w[3];
  schnaps_real By = w[4];
  schnaps_real invtemp = p->invtheta;
  schnaps_real u2 = (ux * ux + uy * uy) * invtemp;
  schnaps_real B2= (Bx * Bx + By * By) * invtemp;
  for (int inode=0;inode < 9;inode++){
    schnaps_real vx = lsd->lb_model->vi[inode][0]; 
    schnaps_real vy = lsd->lb_model->vi[inode][1]; 
    //
    schnaps_real v2 = (vx * vx + vy * vy) * invtemp;
    schnaps_real uv = (ux * vx + uy * vy) * invtemp;
    schnaps_real vB=  (vx * Bx + vy * By) * invtemp;
    f[inode] =
        LBM_WEIGHTS_D2Q9_ISOTH[inode] * 
        (rho * (1.0 + uv + 0.5 * (uv * uv - u2))
        + 0.5 * (0.5 * B2 * v2 - vB * vB));
  }
 schnaps_real invtemp_mag = p->invtheta_mag;
  for (int inode=9;inode < 14;inode ++){
    int inodeloc=inode-9;
    schnaps_real vy=lsd->lb_model->vi[inode][1];
    f[inode] = LBM_WEIGHTS_D2Q5_ISOTH[inodeloc] * (Bx + invtemp_mag * vy *(uy * Bx - ux * By));
  }
  for (int inode=14;inode < q;inode ++){
    int inodeloc=inode-14;
    schnaps_real vx=lsd->lb_model->vi[inode][0];
    f[inode] = LBM_WEIGHTS_D2Q5_ISOTH[inodeloc] * (By + invtemp_mag * vx *(ux * By - uy * Bx));
  }
}

#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "macromesh.h"


int CompareFatNode(const void* a,const void* b){

  FatNode* fna = (FatNode*) a;
  FatNode* fnb = (FatNode*) b;

  int r = fna->x_int[0]-fnb->x_int[0];
  if (r==0 || r==-1 || r==1)
    r = fna->x_int[1]-fnb->x_int[1];
  if (r==0 || r==-1 || r==1)
    r = fna->x_int[2]-fnb->x_int[2];
  if (r==0 || r==-1 || r==1)
    r=0;
  return r;

}

void slice_projection(schnaps_real* xin,schnaps_real* xout){
  // flat tokamak
  xout[0] = xin[2];
  xout[1] = xin[0];
  xout[2] = xin[1];
}

int BuildFatNodeList(Simulation *simu,FatNode* fn_list){


  int big_int = 1 << 28; // 2**28 = 268 435 456

  int ino=0;
  schnaps_real* xmin=simu->macromesh.xmin;
  schnaps_real* xmax=simu->macromesh.xmax;
  int nb_dg_nodes=0;

  schnaps_real smin[3], smax[3];

  smin[0] = 1e10;
  smin[1] = 1e10;
  smin[2] = 1e10;
  smax[0] = -1e10;
  smax[1] = -1e10;
  smax[2] = -1e10;
  
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    
    field *f = &simu->fd[ie];

    nb_dg_nodes =  NPG(f->deg, f->raf) * simu->macromesh.nbelems;

    
    
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real xpg[3];
      schnaps_real xref[3];
      ref_pg_vol(f->deg, f->raf, ipg, xref, NULL, NULL);
      schnaps_ref2phy(f->physnode,
		      xref,
		      0, -1, // dphiref, ifa
		      xpg, NULL,
		      NULL, NULL, NULL); // codtau, dphi, vnds


      slice_projection(xpg, fn_list[ino].x);

      for(int dim = 0; dim < 3; dim++){
	smin[dim] = smin[dim] >  fn_list[ino].x[dim] ? fn_list[ino].x[dim] : smin[dim];
	smax[dim] = smax[dim] <  fn_list[ino].x[dim] ? fn_list[ino].x[dim] : smax[dim];
      }
      
      fn_list[ino].dg_index = ino;

      ino++;
    }
  }

  assert(ino == nb_dg_nodes);

  for(int ino = 0; ino < nb_dg_nodes; ino++){
    // convert points to "pixels" and round to nearest integer
    schnaps_real xloc[3]={fn_list[ino].x[0],
			  fn_list[ino].x[1],
			  fn_list[ino].x[2]};
    
    fn_list[ino].x_int[0]=(int) ((xloc[0]-smin[0])/(smax[0]-smin[0]) * big_int)  ;
    fn_list[ino].x_int[1]=(int) ((xloc[1]-smin[1])/(smax[1]-smin[1]) * big_int)  ;
    fn_list[ino].x_int[2]=(int) ((xloc[2]-smin[2])/(smax[2]-smin[2]) * big_int)  ;
  }

  qsort(fn_list, nb_dg_nodes, sizeof(FatNode),CompareFatNode);

  fn_list[0].fe_index=0;
  int fe_index=0;
  
  for(int ino=1;ino<nb_dg_nodes;ino++){
    if (CompareFatNode(fn_list+ino,fn_list+ino-1)!=0) fe_index++;
    fn_list[ino].fe_index=fe_index;
  }

  return fe_index+1;

}

void InitContinuousSolver(void * cs, Simulation* simu,int type_bc,int nb_phy_vars, int * listvar){
  ContinuousSolver * ps=cs;
  ps->simu = simu;

  ps->type_bc=type_bc;
  ps->nb_phy_vars=nb_phy_vars;
  ps->list_of_var = calloc(nb_phy_vars,sizeof(int));
  assert(sizeof(ps->list_of_var)==sizeof(listvar));
  for (int i=0; i<nb_phy_vars;i++){
    ps->list_of_var[i]=listvar[i];
  }

  ps->postcomputation_assembly=NULL;
  ps->matrix_assembly=NULL;
  ps->rhs_assembly=NULL;
 
  field *f0 = &simu->fd[0];

  ps->nb_dg_nodes = NPG(f0->deg, f0->raf)
    * simu->macromesh.nbelems;

  ps->nb_dg_dof= ps->nb_dg_nodes * ps->nb_phy_vars;
  
  ps->fn_list = malloc(ps->nb_dg_nodes * sizeof(FatNode));
  assert(ps->fn_list);
  // paste the nodes of the DG mesh
  ps->nb_fe_nodes=BuildFatNodeList(simu,ps->fn_list);
  
  ps->nb_fe_dof= ps->nb_fe_nodes * ps->nb_phy_vars;
  ps->dg_to_fe_index = malloc(ps->nb_dg_nodes * sizeof(int));
  assert(ps->dg_to_fe_index);

  for(int ino=0;ino<ps->nb_dg_nodes;ino++){
    ps->dg_to_fe_index[ps->fn_list[ino].dg_index]=ps->fn_list[ino].fe_index;
  }

  // count the number of nodes in each slice
  // first build a list of integer points of
  // finite element nodes
  int *xfe = malloc(3 * ps->nb_fe_nodes * sizeof(int));
  
  for(int ino = 0; ino < ps->nb_dg_nodes; ino++){
    int ife = ps->dg_to_fe_index[ps->fn_list[ino].dg_index];
    xfe[3 * ife + 0] = ps->fn_list[ino].x_int[0];
    xfe[3 * ife + 1] = ps->fn_list[ino].x_int[1];
    xfe[3 * ife + 2] = ps->fn_list[ino].x_int[2];
  }
  
  ps->slice_size = 0;
  ps->nb_slices=0;
  
  if (schnaps_kinetic_data.solve_quasineutrality){
    
    for(int ino = 0; ino < ps->nb_fe_nodes; ino++){
      int i1 = xfe[3 * ino + 0];
      if (i1 == 0 || i1 == -1 || i1 == 1) {
	ps->slice_size++;
      }
      
    }
    
    
    if(ps->slice_size != 0){
      ps->nb_slices = ps->nb_fe_nodes / ps->slice_size;
      if (ps->nb_fe_nodes % ps->slice_size != 0){
	printf("Incorrect slice size or slice useless:%d %d\n",
	       ps->nb_fe_nodes, ps->slice_size);
	assert(ps->nb_fe_nodes % ps->slice_size == 0);
      }
    }
    else{
      ps->nb_slices=0;
    }

  }
  
  free(xfe);


  // now construct the list of boundary nodes
  ps->is_boundary_node = malloc(ps->nb_fe_nodes * sizeof(int));
  assert(ps->is_boundary_node);
  for(int ino = 0; ino < ps->nb_fe_nodes; ino++){
    ps->is_boundary_node[ino] = 0;
  }

  int nraf[3] = {f0->raf[0],f0->raf[1],f0->raf[2]};
  
  int deg[3] = {f0->deg[0],f0->deg[1],f0->deg[2]};

  int npg[3] = {f0->npg[0],f0->npg[1],f0->npg[2]};

  ps->nnodes = npg[0] * npg[1] * npg[2];
  ps->npgmacrocell = ps->nnodes *  nraf[0] * nraf[1] * nraf[2];

  
  ps->nbel = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2];

  for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    int nbfa = 6;
    if (simu->macromesh.is2d)  nbfa = 4; 
    if (simu->macromesh.is1d) nbfa = 2;
    for(int ii = 0; ii < nbfa; ii++) {
      int ifa = ii;
      if (simu->macromesh.is1d) ifa = 2 * ii + 1;
      int ieR = simu->macromesh.elem2elem[6*ie+ifa];
      if (ieR < 0) {
	for(int ipgf = 0; ipgf < NPGF(deg,nraf, ifa); ipgf++) {
	  int ipg = ref_pg_face(deg,nraf, ifa, ipgf,
				NULL, NULL, NULL);
	  int ino_dg = ipg + ie * ps->npgmacrocell;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  ps->is_boundary_node[ino_fe] = 1;

	}
      }
    }
  }

  int count_boundary = 0;

  for(int ino = 0; ino < ps->nb_fe_nodes; ino++){
    count_boundary += ps->is_boundary_node[ino];
  }

  InitLinearSolver(&ps->lsol,ps->nb_fe_dof,NULL,NULL); //InitSkyline(&sky, neq);

  ps->lsol.rhs = malloc(ps->nb_fe_dof * sizeof(schnaps_real));
  assert(ps->lsol.rhs);
  ps->lsol.sol = calloc(ps->nb_fe_dof, sizeof(schnaps_real));
  assert(ps->lsol.sol);

  AllocateContinuousMatrix(ps);

}


void SolveContinuous2D(void* cs){

  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];

  //printf("Init...\n");
  
  int nraf[3] = {f0->raf[0],f0->raf[1],f0->raf[2]};
  int deg[3] = {f0->deg[0],f0->deg[1],f0->deg[2]};

  
  // number of equation of the Poisson solver
  int neq=ps->nb_fe_dof; //(nodes)
 
  //printf("Matrix assembly nb_eqs=%d.....\n",ps->nb_fe_dof);
  if(ps->matrix_assembly != NULL){
    ps->matrix_assembly(ps);
  }
  else
    {
      printf("no matrix assembly ");
      exit(1);
    }

  
  //printf("RHS assembly.....\n");
  if(ps->rhs_assembly != NULL){
    ps->rhs_assembly(ps);
  }
  else
    {
      printf("no rhs assembly ");
      exit(1);
    }
  
  //printf("BC assembly.....\n");
  if(ps->bc_assembly != NULL){
    ps->bc_assembly(ps);
  }
  else
    {
      printf("no bc assembly ");
      exit(1);
    }
  
  for(int i=0;i<ps->nb_fe_dof;i++){
    //printf(" rhs %d %f\n",i,ps->lsol.rhs[i]);
  }
  //assert(1==2);
   
  Advanced_SolveLinearSolver(&ps->lsol,ps->simu);

  for(int i=0;i<ps->nb_fe_dof;i++){
    //  printf(" sol %d %f\n",i,ps->lsol.sol[i]);
  }
  //assert(1==2);  
  //printf("post computation assembly.....\n");
  if(ps->postcomputation_assembly != NULL){
    ps->postcomputation_assembly(ps);
  }


  //
  // printf("End Solve2D.\n");

}


void ContinuousToDiscontinuous_Copy(ContinuousSolver * cs){
  
  field* f0 = &cs->simu->fd[0];

  //printf("Copy...\n");

  
  // copy the potential at the right place
  for(int var =0; var < cs->nb_phy_vars; var++){ 
    for(int ie = 0; ie < cs->simu->macromesh.nbelems; ie++){  
      for(int ipg = 0;ipg < cs->npgmacrocell; ipg++){
	int ino_dg = ipg + ie * cs->npgmacrocell;
	int ino_fe = cs->dg_to_fe_index[ino_dg];
	int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
				ipg,cs->list_of_var[var]);
	cs->simu->fd[ie].wn[ipot]=cs->lsol.sol[ino_fe];
      }
    }
  }

}

void ExactDirichletContinuousMatrix(void * cs){
  ContinuousSolver * ps=cs;

  field* f0 = &ps->simu->fd[0];

  for(int ino=0; ino<ps->nb_fe_nodes; ino++){
    if (ps->is_boundary_node[ino]){
      for (int iv=0; iv<ps->nb_phy_vars;iv++){
        int iBord = ps->nb_phy_vars*ino+iv;
        for(int i=0; i<ps->nb_fe_dof; i++){
          SetLinearSolver(&ps->lsol,iBord,i,0.0);
        }
	SetLinearSolver(&ps->lsol,iBord,iBord,1.0);
      }
    }
  }
  ps->lsol.mat_is_assembly=true;
  
  for(int var =0; var < ps->nb_phy_vars; var++){ 
    //for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){  
    for (int i=0; i<ps->simu->macromesh.nboundaryfaces;i++){
      int ifa = ps->simu->macromesh.boundaryface[i];
      int locfaL = ps->simu->macromesh.face2elem[4 * ifa + 1];
      int ie = ps->simu->macromesh.face2elem[4 * ifa ];
      int ieR = ps->simu->macromesh.face2elem[4 * ifa + 2];
      assert(ieR < 0);
      if (ieR<0){
        field *f = &ps->simu->fd[ie];

        for(int ipglf = 0;ipglf < NPGF(f->deg,f->raf,locfaL); ipglf++){
          schnaps_real xpgref[3], xpgref_in[3], wpg;
          
          // Get the coordinates of the Gauss point and coordinates of a
          // point slightly inside the opposite element in xref_in
          int ipg = ref_pg_face(f->deg, f->raf, locfaL, ipglf, xpgref, &wpg, xpgref_in);
          int ino_dg = ipg + ie * ps->npgmacrocell;
          int ino_fe = ps->dg_to_fe_index[ino_dg];
          int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
				  ipg,ps->list_of_var[var]);
          int ipot_fe = ino_fe*ps->nb_phy_vars + var;
          // Normal vector at gauss point ipgL
          schnaps_real vnds[3], xpg[3];
          {
            schnaps_real dtau[3][3], codtau[3][3];
            schnaps_ref2phy(f->physnode,
			    xpgref,
			    NULL, locfaL, // dpsiref, ifa
			    xpg, dtau,
			    codtau, NULL, vnds); // codtau, dpsi, vnds
          }
          
          // the boundary flux is an affine function
          schnaps_real flux0[f->model.m];
          f->model.ImposedData(xpg, f->tnow, flux0);
          ps->lsol.rhs[ipot_fe] = flux0[ps->list_of_var[var]];
        }
      }
    }
  }
}




void PenalizedDirichletContinuousMatrix(void * cs){
  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];
  schnaps_real bigval = 1.e20;//.e16;
  for(int ino=0; ino<ps->nb_fe_nodes; ino++){
   
    if (ps->is_boundary_node[ino]){
      for (int iv=0; iv<ps->nb_phy_vars;iv++){
        SetLinearSolver(&ps->lsol,ps->nb_phy_vars*ino+iv,ps->nb_phy_vars*ino+iv,bigval);
      }
    }
  }

  ps->lsol.mat_is_assembly=true;
  
  for(int var =0; var < ps->nb_phy_vars; var++){ 
    //for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){  
    for (int i=0; i<ps->simu->macromesh.nboundaryfaces;i++){
      int ifa = ps->simu->macromesh.boundaryface[i];
      int locfaL = ps->simu->macromesh.face2elem[4 * ifa + 1];
      int ie = ps->simu->macromesh.face2elem[4 * ifa ];
      int ieR = ps->simu->macromesh.face2elem[4 * ifa + 2];
      if (ieR<0){
        field *f = &ps->simu->fd[ie];

        for(int ipglf = 0;ipglf < NPGF(f->deg,f->raf,locfaL); ipglf++){
          schnaps_real xpgref[3], xpgref_in[3], wpg;
          
          // Get the coordinates of the Gauss point and coordinates of a
          // point slightly inside the opposite element in xref_in
          int ipg = ref_pg_face(f->deg, f->raf, locfaL, ipglf, xpgref, &wpg, xpgref_in);
          int ino_dg = ipg + ie * ps->npgmacrocell;
          int ino_fe = ps->dg_to_fe_index[ino_dg];
          int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
				  ipg,ps->list_of_var[var]);
          int ipot_fe = ino_fe*ps->nb_phy_vars + var;
          // Normal vector at gauss point ipgL
          schnaps_real vnds[3], xpg[3];
          {
            schnaps_real dtau[3][3], codtau[3][3];
            schnaps_ref2phy(f->physnode,
			    xpgref,
			    NULL, locfaL, // dpsiref, ifa
			    xpg, dtau,
			    codtau, NULL, vnds); // codtau, dpsi, vnds
          }
          
          // the boundary flux is an affine function
          schnaps_real flux0[f->model.m];
          f->model.ImposedData(xpg, f->tnow, flux0);
          ps->lsol.rhs[ipot_fe] = flux0[ps->list_of_var[var]] * bigval;
        }
      }
    }
  }
}



void AllocateContinuousMatrix(void *cs){
  ContinuousSolver * ps=cs;

  //static bool is_lu = false;

  field* f0 = &ps->simu->fd[0];
  //printf("Init...\n");

  // number of equation of the Poisson solver
  int neq=ps->nb_fe_dof; //(nodes)
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m = ps->nb_phy_vars;

  int nraf[3] = {f0->raf[0],f0->raf[1],f0->raf[2]};
  int deg[3] = {f0->deg[0],f0->deg[1],f0->deg[2]};
  
  if (ps->simu->macromesh.is2d) {
    assert( nraf[0] == nraf[1]);
    assert( nraf[2] == 1);
    assert(deg[2] == 0);
  }
  
  if (ps->simu->macromesh.is1d) {
    assert( nraf[1] == 1);
    assert(deg[1] == 0);
    assert( nraf[2] == 1);
    assert(deg[2] == 0);
  }
 
  //printf("Allocation...\n");
  if(!ps->lsol.is_alloc){

    // compute the profile of the matrix
    for(int ie = 0; ie < ps->nbel; ie++){
      for(int iloc = 0; iloc < ps->nnodes; iloc++){
	for(int jloc = 0; jloc < ps->nnodes; jloc++){
	  int ino_dg = iloc + ie * ps->nnodes;
	  int jno_dg = jloc + ie * ps->nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  int ino_slice = 0;
	  int jno_slice = 0;
	  if (schnaps_kinetic_data.solve_quasineutrality){
	    ino_slice = ino_fe / ps->slice_size;
	    jno_slice = jno_fe / ps->slice_size;
	  }
	  for (int iv1=0;iv1<ps->nb_phy_vars;iv1++){
	    for (int iv2=0;iv2<ps->nb_phy_vars;iv2++){
	      if (ino_slice == jno_slice)
		IsNonZero(&ps->lsol, ino_fe*ps->nb_phy_vars+iv1, jno_fe*ps->nb_phy_vars+iv2);
	    }
	  }
	}
      }
    }
    
    
      AllocateLinearSolver(&ps->lsol);
  } 
}



void GenericOperator_Continuous(void * cs){

  ContinuousSolver * ps=cs;

  field* f0 = &ps->simu->fd[0];

  if(!ps->lsol.mat_is_assembly){
    for(int ie = 0; ie < ps->nbel; ie++){  

      // local matrix 
      schnaps_real aloc[ps->nnodes*ps->nb_phy_vars][ps->nnodes*ps->nb_phy_vars];
      for(int iloc = 0; iloc < ps->nnodes*ps->nb_phy_vars; iloc++){
	for(int jloc = 0; jloc < ps->nnodes*ps->nb_phy_vars; jloc++){
	  aloc[iloc][jloc] = 0.0;
	}
      }

      int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
      int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);

      for(int ipg = 0;ipg < ps->nnodes; ipg++){
	schnaps_real wpg;
	schnaps_real xref[3];
	int ipgmacro= ipg + isubcell * ps->nnodes;

	ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

	for(int iloc = 0; iloc < ps->nnodes; iloc++){
	  schnaps_real dtau[3][3],codtau[3][3];
	  schnaps_real dphiref_i[3],dphiref_j[3];
	  schnaps_real dphi_i[3],dphi_j[3];
	  schnaps_real basisPhi_i[4], basisPhi_j[4];
	  int ilocmacro = iloc + isubcell * ps->nnodes;
	  grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
	  schnaps_ref2phy(ps->simu->fd[iemacro].physnode,
			  xref,dphiref_i,0,NULL,
			  dtau,codtau,dphi_i,NULL);
        
	  schnaps_real det = dot_product(dtau[0], codtau[0]);
	  if (ilocmacro==ipgmacro){
	    basisPhi_i[0]=1;
	  }
	  else
	    {
	      basisPhi_i[0]=0;
	    }
	  basisPhi_i[1]=dphi_i[0]/det;
	  basisPhi_i[2]=dphi_i[1]/det;
	  basisPhi_i[3]=dphi_i[2]/det;
	  for(int jloc = 0; jloc < ps->nnodes; jloc++){
	    int jlocmacro = jloc + isubcell * ps->nnodes;
	    grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
	    schnaps_ref2phy(ps->simu->fd[iemacro].physnode,
			    xref,dphiref_j,0,NULL,
			    dtau,codtau,dphi_j,NULL);
	    if (jlocmacro==ipgmacro){
	      basisPhi_j[0]=1;
	    }
	    else
	      {
		basisPhi_j[0]=0;
	      }
	    basisPhi_j[1]=dphi_j[0]/det;
	    basisPhi_j[2]=dphi_j[1]/det;
	    basisPhi_j[3]=dphi_j[2]/det;
	    for (int iv1=0; iv1<ps->nb_phy_vars; iv1++){
	      for (int iv2=0; iv2<ps->nb_phy_vars; iv2++){
		schnaps_real res[4] = {0, 0, 0, 0};
		for (int i=0; i<4; i++){
		  for (int j=0; j<4; j++){
		    res[i]+=basisPhi_j[j]*ps->diff_op[ps->nb_phy_vars*iv1+iv2].DO[i][j];
		  }
		}
		aloc[iv1+iloc*ps->nb_phy_vars][iv2+jloc*ps->nb_phy_vars] += dot_product(basisPhi_i, res) * wpg * det  ;
	      }
	    }
	  }
	}
      }

      for(int iloc = 0; iloc < ps->nnodes; iloc++){
	for(int jloc = 0; jloc < ps->nnodes; jloc++){
	  int ino_dg = iloc + ie * ps->nnodes;
	  int jno_dg = jloc + ie * ps->nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  int ino_slice = 0;
	  int jno_slice = 0;
	  if (schnaps_kinetic_data.solve_quasineutrality){
	    ino_slice = ino_fe / ps->slice_size;
	    jno_slice = jno_fe / ps->slice_size;
	  }
	  //printf("ino=%d islice=%d nbslices=%d\n",ino_fe, ino_fe / ps->slice_size, ps->nb_slices);
	  for (int iv1=0;iv1<ps->nb_phy_vars;iv1++){
	    for (int iv2=0;iv2<ps->nb_phy_vars;iv2++){
	      schnaps_real val = aloc[iloc*ps->nb_phy_vars+iv1][jloc*ps->nb_phy_vars+iv2];
	      if (ino_slice == jno_slice)
		AddLinearSolver(&ps->lsol,ino_fe*ps->nb_phy_vars+iv1,jno_fe*ps->nb_phy_vars+iv2,val);
	    }
	  }
	}
      }
    }
  }
  ps->lsol.mat_is_assembly = true;
}



void freeContinuousSolver(ContinuousSolver* cs){

  LinearSolver* lsol = &(cs->lsol);
  FreeLinearSolver(lsol);
  if (cs->fn_list!=NULL){
    free(cs->fn_list);
  }
  if (cs->is_boundary_node!=NULL){
    free(cs->is_boundary_node);
  }
  if (cs->list_of_var!=NULL){
    free(cs->list_of_var);
  }
  if (cs->dg_to_fe_index!=NULL){
    free(cs->dg_to_fe_index);
  }
  cs->bc_assembly=NULL;
  cs->matrix_assembly=NULL;
  cs->rhs_assembly=NULL;
  cs->postcomputation_assembly=NULL;
}





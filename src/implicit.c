#include "implicit.h"
#include <stdlib.h>
#include <string.h>
//#include "csparse.h"


//void InternalCoupling(Simulation *simu,  LinearSolver *solver, int itest);
//void FluxCoupling(Simulation *simu,  LinearSolver *solver,int itest);

/* void InternalAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void FluxAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void InterfaceAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void SourceAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt); */
/* void MassAssembly(Simulation *simu,  LinearSolver *solver); */

/* void AssemblyImplilcitLinearSolver(Simulation *simu, LinearSolver *solver,real theta, real dt); */


void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver, MatrixStorage ms){

  int neq = simu->wsize;

  //MatrixStorage ms = SKYLINE;
  //MatrixStorage ms = KLU_CSR;

  Solver st = LU;
  InitLinearSolver(solver,neq,&ms,&st);
  //solver->tol=1.e-8;
  solver->tol=_SMALL;

  int itest = 1;
  for (int isky=0 ; isky < itest; isky++){

    InternalCoupling(simu, solver, isky);
    FluxCoupling(simu, solver, isky);
    InterfaceCoupling(simu, solver, isky);

    if (isky == 0) AllocateLinearSolver(solver);
  }

  //DisplayLinearSolver(solver);


}

void InitFieldImplicitSolver(field *fd, MatrixStorage ms){

  int neq = fd->wsize;

  //MatrixStorage ms = SKYLINE;
  Solver st = LU;
  if (fd->solver == NULL) fd->solver = malloc(sizeof(LinearSolver));
  if (fd->rmat == NULL) fd->rmat = malloc(sizeof(LinearSolver));
  InitLinearSolver(fd->solver,neq,&ms,&st);
  InitLinearSolver(fd->rmat,neq,&ms,&st);

  // int itest = 2; // for testing
  int itest = 1;

  for (int isky=0 ; isky < itest; isky++){

    InternalLocalCoupling(fd, isky);
    FluxLocalCoupling(fd, isky);

    if (isky == 0) {
      AllocateLinearSolver(fd->solver);
      AllocateLinearSolver(fd->rmat);
    }
  }

  //DisplayLinearSolver(fd->solver);
  //assert(1==2);


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


void AssemblyFieldImplicitSolver(field *fd,schnaps_real theta, schnaps_real dt)
{

  if(fd->solver->mat_is_assembly == false){
    assert(fd->rmat->mat_is_assembly == false);
    MassLocalAssembly(fd);
    InternalLocalAssembly(fd, theta, dt);
    FluxLocalAssembly(fd, theta, dt);


    fd->solver->mat_is_assembly = true;
    fd->rmat->mat_is_assembly = true;
  }



  if(fd->solver->rhs_is_assembly == false){
    for(int i=0;i<fd->solver->neq;i++){
      fd->solver->rhs[i]=0;
    }
  }
  /* DisplayLinearSolver(fd->solver); */
  /* assert(1==2); */
}

void LocalThetaTimeScheme(Simulation *simu, schnaps_real tmax, schnaps_real dt)
{

  schnaps_real theta=0.5;
  simu->dt=dt;
  int itermax=tmax / simu->dt;
  simu->itermax_rk=itermax;
  simu->tnow=0;
  simu->tmax = tmax;

  // assembly of the volume part of the matrix
  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field *f = simu->fd + ie;
    InitFieldImplicitSolver(f, SKYLINE);
    //performs: InternalLocalCoupling(fd, isky);
    //and FluxLocalCoupling(fd, isky);

    AssemblyFieldImplicitSolver(f, theta, dt);
    // performs: MassLocalAssembly(fd);
    //InternalLocalAssembly(fd, theta, dt);
    //FluxLocalAssembly(fd, theta, dt);
  }

  // assembly of the boundary condition part of the matrix
  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    assert(inter->fL);
    InterfaceLocalAssembly(inter, theta, dt);
  }

  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  freq = 1;
  int iter = 0;

  time_t start;
  while(simu->tnow < tmax) {


    if (iter == 1) start = time(NULL);

    for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      f->dt = simu->dt;
      MatVect(f->rmat, f->wn, f->solver->rhs);
      //for(int i=0;i<f->solver->neq;i++) f->solver->rhs[i]=0;
    }



    simu->tnow += theta * dt;

    for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
      field *f = simu->fd + ie;
      /* DisplayLinearSolver(f->solver); */
      /* DisplayLinearSolver(f->rmat); */
      /* assert(1==3); */
      f->tnow = simu->tnow;
      f->dt = simu->dt;
      SourceLocalAssembly(f, 1. , dt);
    }

    for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
      Interface* inter = simu->interface + ifa;
      // left = 0  right = 1
      ExtractInterface(inter, 0);
      ExtractInterface(inter, 1);
      InterfaceExplicitFlux(inter, 0);
      InterfaceExplicitFlux(inter, 1);
    }

    /* for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){ */
    /*   field *f = simu->fd + ie; */
    /*   DisplayLinearSolver(f->solver); */
    /*   DisplayLinearSolver(f->rmat); */
    /* } */



    simu->tnow += (1 - theta) * dt;

    for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      SolveLinearSolver(f->solver);
      for(int i=0;i<f->solver->neq;i++){
	f->wn[i] = f->solver->sol[i];
	//printf("i=%d sol=%f\n",i,f->solver->sol[i]);
      }
    }

    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n",
	     simu->tnow, iter+1, simu->itermax_rk, dt);
    iter++;

  }

  printf("Elapsed time=%f\n", (double) (time(NULL) -start));

}

void LocalThetaTimeScheme_SPU(Simulation *simu, schnaps_real tmax, schnaps_real dt){

  schnaps_real theta=0.5;
  simu->dt=dt;

  int itermax=tmax / simu->dt;
  simu->itermax_rk=itermax;
  simu->tnow=0;
  simu->tmax = tmax;

  // assembly of the volume part of the matrix
  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field *f = simu->fd + ie;
    InitFieldImplicitSolver(f, SKYLINE_SPU);
    AssemblyFieldImplicitSolver(f, theta, dt);
  }

  // assembly of the boundary condition part of the matrix
  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    assert(inter->fL);
    InterfaceLocalAssembly(inter, theta, dt);
  }

  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  //freq = 10;
  int iter = 0;


  // clean the allocated arrays in each field
  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field *f = simu->fd + ie;
    f->local_source_cl_init = false;
    Skyline_SPU* solver = f->solver->matrix;
    FactoLU_SPU(solver);
    free(solver->sol);
    solver->sol = f->wn;
    Skyline_SPU* rmat = f->rmat->matrix;
    free(rmat->sol);
    free(rmat->rhs);

  }

  time_t start;
  while(simu->tnow < tmax) {


    if (iter == 1) start = time(NULL);


    for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      f->dt = simu->dt;
      Skyline_SPU* sky_spu = f->solver->matrix;
      RegisterSkyline_SPU(sky_spu);
      MatVect_SPU(f->rmat,
		  sky_spu->sol_handle,
		  sky_spu->rhs_handle); // note: matvect implies a register
    }



    simu->tnow += theta * dt;

    for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      f->dt = simu->dt;
      SourceLocalAssembly_SPU(f, 1. , dt);
    }

    for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
      Interface* inter = simu->interface + ifa;
      // left = 0  right = 1
      ExtractInterface_SPU(inter, 0);
      ExtractInterface_SPU(inter, 1);
      if (inter->fR != NULL) {
      	InterfaceExplicitFlux_SPU(inter, 0);
      	InterfaceExplicitFlux_SPU(inter, 1);
      }
      else{
      	InterfaceBoundaryFlux_SPU(inter);
      }
    }

    simu->tnow += (1 - theta) * dt;

    for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      SolveLinearSolver(f->solver);
      for(int i=0;i<f->solver->neq;i++){
      }
    }

    if (iter % freq == 0) {
      printf("t=%f iter=%d/%d dt=%f\n",
	     simu->tnow, iter+1, simu->itermax_rk, dt);
    }
    iter++;
  }

  starpu_task_wait_for_all();
  printf("Elapsed time=%f\n", (double) (time(NULL) -start));

  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field* f = simu->fd + ie;
    Skyline_SPU* sky_spu = f->solver->matrix;
    UnRegisterSkyline_SPU(sky_spu);
  }

  destroy_global_arbiter();
  starpu_shutdown();



};


void GraphThetaTimeScheme_SPU(Simulation *simu, schnaps_real tmax, schnaps_real dt){

  schnaps_real theta=0.5;
  //schnaps_real theta=1;
  simu->dt=dt;

  int itermax=tmax / simu->dt;
  simu->itermax_rk=itermax;
  simu->tnow=0;
  simu->tmax = tmax;

  // assembly of the volume part of the matrix
  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field *f = simu->fd + ie;
    InitFieldImplicitSolver(f, SKYLINE_SPU);
    AssemblyFieldImplicitSolver(f, theta, dt);
  }

  // assembly of the boundary condition part of the matrix
  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    assert(inter->fL);
    InterfaceLocalAssembly(inter, theta, dt);
  }

  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  //freq = 10;
  int iter = 0;


  // clean the allocated arrays in each field
  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field *f = simu->fd + ie;
    f->local_source_cl_init = false;
    Skyline_SPU* solver = f->solver->matrix;
    FactoLU_SPU(solver);
    free(solver->sol);
    solver->sol = f->wn;
    free(solver->rhs);
    solver->rhs = f->res;
    RegisterSkyline_SPU(solver);
    Skyline_SPU* rmat = f->rmat->matrix;
    free(rmat->sol);
    free(rmat->rhs);

  }

  time_t start;
  start = time(NULL);
  while(simu->tnow < tmax) {

    if (iter <= 1) start = time(NULL);

    simu->tnow += theta * dt;

    for(int ies=0; ies <  simu->macromesh.nbelems; ++ies){
      // -1 because the ghost upwind vertex is always the first
      // and the ghost downwind vertex is always the last
      int ie = simu->macromesh.topo_order[ies+1];
      assert(ie >=0 && ie < simu->macromesh.nbelems);
      //printf("ies=%d ie=%d\n",ies,ie);
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      f->dt = simu->dt;
      Skyline_SPU* sky_spu = f->solver->matrix;
      MatVect_SPU(f->rmat,
		  sky_spu->sol_handle,
		  sky_spu->rhs_handle);

      SourceLocalAssembly_SPU(f, 1. , dt);

      // search upwind boundaries
      igraph_t* graph = &(simu->macromesh.connect_graph);
      igraph_vector_t found_edge;
      igraph_vector_init(&found_edge, 6);
      igraph_incident(graph, &found_edge, ie, IGRAPH_IN);
      int nup = igraph_vector_size(&found_edge);
      for(int ii=0; ii < nup; ii++) {
	int eid = (int)VECTOR(found_edge)[ii];
	int iefrom, ieto;
	igraph_edge(graph, eid, &iefrom, &ieto);
	assert(ie == ieto);
	int side = simu->macromesh.edge_dir[eid];
	int interid = simu->macromesh.edge2face[eid];
	Interface *inter = simu->interface + interid;
	if (iefrom >= simu->macromesh.nbelems) {
	  InterfaceBoundaryFlux_SPU(inter);
	}
	else {
	  InterfaceExplicitFlux_SPU(inter, 1 - side);
	}
      }

      SolveLinearSolver(f->solver);

      igraph_incident(graph, &found_edge, ie, IGRAPH_OUT);
      int ndown = igraph_vector_size(&found_edge);
      for(int ii=0; ii < ndown; ii++) {
	int eid = (int)VECTOR(found_edge)[ii];
	int iefrom, ieto;
	igraph_edge(graph, eid, &iefrom, &ieto);
	assert(ie == iefrom);
	int side = simu->macromesh.edge_dir[eid];
	int interid = simu->macromesh.edge2face[eid];
	Interface *inter = simu->interface + interid;
	if (ieto < simu->macromesh.nbelems) {
	  ExtractInterface_SPU(inter, side);
	}
      }
    }

    simu->tnow += (1 - theta) * dt;



    if (iter % freq == 0) {
      printf("t=%f iter=%d/%d dt=%f\n",
	     simu->tnow, iter+1, simu->itermax_rk, dt);
    }
    iter++;
    if (iter == 1) {
      starpu_task_wait_for_all();
      printf("Elapsed time first iter=%f\n", (double) (time(NULL) -start));
    }
  }

  starpu_task_wait_for_all();
  //starpu_task_wait_for_all();
  printf("Elapsed time=%f\n", (double) (time(NULL) -start));

  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field* f = simu->fd + ie;
    Skyline_SPU* sky_spu = f->solver->matrix;
    UnRegisterSkyline_SPU(sky_spu);
  }


  destroy_global_arbiter();
  starpu_shutdown();



};

void GraphThetaTimeSchemeSubCell_SPU(Simulation *simu, schnaps_real tmax, schnaps_real dt){

  schnaps_real theta=0.5;
  simu->dt=dt;

  int itermax=tmax / simu->dt;
  simu->itermax_rk=itermax;
  simu->tnow=0;
  simu->tmax = tmax;


  for(int ie=0; ie <  simu->macromesh.nbelems; ++ie){
    field *f = simu->fd + ie;
    f->local_source_cl_init = false;
  }


  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  //freq = 10;
  int iter = 0;

  time_t start;
  start = time(NULL);
  while(simu->tnow < tmax) {

    if (iter <= 1) start = time(NULL);

    simu->tnow += theta * dt;

    for(int ies=0; ies <  simu->macromesh.nbelems; ++ies){
      // -1 because the ghost upwind vertex is always the first
      // and the ghost downwind vertex is always the last
      int ie = simu->macromesh.topo_order[ies+1];
      assert(ie >=0 && ie < simu->macromesh.nbelems);
      //printf("ies=%d ie=%d\n",ies,ie);
      field *f = simu->fd + ie;
      f->tnow = simu->tnow;
      f->dt = simu->dt;

      // explicit part of the Crank-Nicolson algo
      // All the macrocell internal terms are computed:
      // fluxes, sources, derivatives
      FieldResidual(f, 1. - theta, dt);

      // search upwind boundaries
      igraph_t* graph = &(simu->macromesh.connect_graph);
      igraph_vector_t found_edge;
      igraph_vector_init(&found_edge, 6);
      igraph_incident(graph, &found_edge, ie, IGRAPH_IN);
      int nup = igraph_vector_size(&found_edge);
      for(int ii=0; ii < nup; ii++) {
	int eid = (int)VECTOR(found_edge)[ii];
	int iefrom, ieto;
	igraph_edge(graph, eid, &iefrom, &ieto);
	assert(ie == ieto);
	int side = simu->macromesh.edge_dir[eid];
	int interid = simu->macromesh.edge2face[eid];
	Interface *inter = simu->interface + interid;
	printf("iefrom=%d ieto=%d side=%d\n",iefrom,ieto,side);
	if (iefrom >= simu->macromesh.nbelems) {
	  //InterfaceBoundaryFlux_SPU(inter);
	  InterfaceExplicitFlux(inter, 0);	}
	else {
	  InterfaceExplicitFlux(inter, 1 - side);
	}
      }

      FieldDownwindSolve(f, theta, dt);

      igraph_incident(graph, &found_edge, ie, IGRAPH_OUT);
      int ndown = igraph_vector_size(&found_edge);
      for(int ii=0; ii < ndown; ii++) {
	int eid = (int)VECTOR(found_edge)[ii];
	int iefrom, ieto;
	igraph_edge(graph, eid, &iefrom, &ieto);
	assert(ie == iefrom);
	int side = simu->macromesh.edge_dir[eid];
	int interid = simu->macromesh.edge2face[eid];
	Interface *inter = simu->interface + interid;
	printf("extract iefrom=%d ieto=%d side=%d\n",iefrom,ieto,side);
	if (ieto < simu->macromesh.nbelems) {
	  ExtractInterface(inter, side);
	}
      }
    }

    simu->tnow += (1 - theta) * dt;



    if (iter % freq == 0) {
      printf("t=%f iter=%d/%d dt=%f\n",
	     simu->tnow, iter+1, simu->itermax_rk, dt);
    }
    iter++;
    if (iter == 1) {
      //starpu_task_wait_for_all();
      printf("Elapsed time first iter=%f\n", (double) (time(NULL) -start));
    }
  }

  //starpu_task_wait_for_all();
  //starpu_task_wait_for_all();
  printf("Elapsed time=%f\n", (double) (time(NULL) -start));


  //destroy_global_arbiter();
  //starpu_shutdown();



};


void  build_subcells_downwind_graph(field *fd, schnaps_real vit[3], int *);


void  build_subcells_downwind_graph(field *f, schnaps_real vit[3], int* subcell_order){


  int nraf[3] = {f->raf[0],
		 f->raf[1],
		 f->raf[2]};
  int deg[3] = {f->deg[0],
		f->deg[1],
		f->deg[2]};
  int npg[3] = {deg[0] + 1,
		deg[1] + 1,
		deg[2] + 1};


  igraph_t dwgraph;

  int nbverts = nraf[0] * nraf[1] * nraf[2] + 2;

  igraph_empty(&dwgraph, nbverts, IGRAPH_DIRECTED);





  int count_edge = 0;

  // Loop on the edges
  // Sweeping subcell faces in the three directions
  for(int dim0 = 0; dim0 < 3; dim0++) {
    int dim_offset[3] = {0, 0, 0};
    // in direction dim0 we have an aditional face
    dim_offset[dim0] = -1;
    for(int icL0 = dim_offset[0]; icL0 < nraf[0]; icL0++) {
      for(int icL1 = dim_offset[1]; icL1 < nraf[1]; icL1++) {
	for(int icL2 = dim_offset[2]; icL2 < nraf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};


	  bool is_first = (icL[dim0] == -1);
	  bool is_last = (icL[dim0] == nraf[dim0] - 1);

	  // Get the left subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;


	  int icR[3] = {icL[0], icL[1], icL[2]};
	  // The right cell index corresponds to an increment in
	  // the dim0 direction
	  icR[dim0]++;
	  int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
	  int offsetR = npg[0] * npg[1] * npg[2] * ncR;

	  const int altdim1[3] = {1, 0, 0};
	  const int altdim2[3] = {2, 2, 1};

	  // now loop on the left glops of the subface
	  //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	  int dim1 = altdim1[dim0];
	  int dim2 = altdim2[dim0];
	  int iL[3];
	  iL[dim0] = deg[dim0];

	  bool edgeLtoR = false;
	  bool edgeRtoL = false;
	  bool upwind = false;
	  bool downwind = false;


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
		// change the glop index if we are on the first cell
		if (is_first) ipgL = ipgR;
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
		vnds[0] = - codtau[0][dim0] * h1h2;
		vnds[1] = - codtau[1][dim0] * h1h2;
		vnds[2] = - codtau[2][dim0] * h1h2;
	      }


	      schnaps_real v_dot_n = vnds[0] * vit[0] + vnds[1] * vit[1] + vnds[2] * vit[2];



	      if (v_dot_n > _SMALL && is_first) upwind = true;
	      if (v_dot_n > _SMALL && is_last) downwind = true;
	      if (v_dot_n < -_SMALL && is_first) downwind = true;
	      if (v_dot_n < -_SMALL && is_last) upwind = true;
	      if (v_dot_n > _SMALL && !(is_first || is_last)) edgeLtoR = true;
	      if (v_dot_n < -_SMALL && !(is_first || is_last)) edgeRtoL = true;

	    }  // face yhat loop
	  } // face xhat loop

	  if (edgeLtoR) {
	    igraph_add_edge(&dwgraph, ncL, ncR);
	    //m->edge2face[count_edge] = ifa;
	    //m->edge_dir[count_edge++] = 0;
	  }
	  if (edgeRtoL) {
	    igraph_add_edge(&dwgraph, ncR, ncL);
	    //m->edge2face[count_edge] = ifa;
	    //m->edge_dir[count_edge++] = 1;
	  }
	  if (upwind && is_first) {
	    igraph_add_edge(&dwgraph, nbverts - 2, ncR);
	    //m->edge2face[count_edge] = ifa;
	    //m->edge_dir[count_edge++] = 1;
	  }
	  if (upwind && is_last) {
	    igraph_add_edge(&dwgraph, nbverts - 2, ncL);
	    //m->edge2face[count_edge] = ifa;
	    //m->edge_dir[count_edge++] = 1;
	  }
	  if (downwind && is_first) {
	    //printf("ncR=%d nbv-1=%d \n",ncR,nbverts - 1);
	    igraph_add_edge(&dwgraph, ncR, nbverts - 1);
	    //m->edge2face[count_edge] = ifa;
	    //m->edge_dir[count_edge++] = 0;
	  }
	  if (downwind && is_last) {
	    igraph_add_edge(&dwgraph, ncL, nbverts - 1);
	    //m->edge2face[count_edge] = ifa;
	    //m->edge_dir[count_edge++] = 0;
	  }


	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop
  } // dim loop

  // drawing
  FILE *ff = fopen("subcell_graph.dot","w");;
  igraph_write_graph_dot(&dwgraph, ff);
  fclose(ff);

  printf("For drawing the graph:\n");
  printf("dot -Tpng subcell_graph.dot -o subcell_graph.png\n");

  igraph_bool_t is_dag;
  igraph_is_dag(&dwgraph, &is_dag);
  igraph_vector_t sorting;
  if (is_dag) {
    printf("The graph is a DAG. Topological sorting...\n");
    igraph_vector_init(&sorting, nbverts);
    igraph_topological_sorting(&dwgraph, &sorting,
			       IGRAPH_OUT);
    //print_vector(&sorting, stdout);
  }

  for(int nid = 0; nid < nbverts; nid++){
    subcell_order[nid] = (int)VECTOR(sorting)[nid];
  }

  igraph_vector_destroy(&sorting);




}


// accumulate the fluxes of the neighbour subcells
// into the residual of cell isub
void get_neighbor_fluxes(field* fd, int isub);


void get_neighbor_fluxes(field* fd, int isub){

  int nraf[3] = {fd->raf[0],
		 fd->raf[1],
		 fd->raf[2]};

  int deg[3] = {fd->deg[0],
		fd->deg[1],
		fd->deg[2]};
  int npg[3] = {deg[0] + 1,
		deg[1] + 1,
		deg[2] + 1};

  const int m = fd->model.m;


  //int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
  int ncL = isub;
  // First glop index in the subcell
  int offsetL = npg[0] * npg[1] * npg[2] * ncL;

  int icL[3];
  icL[0] = isub % nraf[0];
  isub /= nraf[0];
  icL[1] = isub % nraf[1];
  isub /= nraf[1];
  icL[2] = isub;
  // dim loop
  for(int dim0 = 0; dim0 < 3; dim0++) {
    // left or right
    for(int s = -1; s <= 1; s += 2){
      int icR[3] = {icL[0], icL[1], icL[2]};
      icR[dim0] += s;
      // the boundary fluxes have been previously computed
      if ((icR[dim0] != -1) && (icR[dim0] != nraf[dim0])){

	int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
	int offsetR = npg[0] * npg[1] * npg[2] * ncR;

	const int altdim1[3] = {1, 0, 0};
	const int altdim2[3] = {2, 2, 1};

	// now loop on the left glops of the subface
	//int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	int dim1 = altdim1[dim0];
	int dim2 = altdim2[dim0];
	int iL[3];
	// we want to perform
	//if (s == 1) {
	//  iL[dim0] = deg[dim0];
	//} else {
	//  iL[dim0] = 0;
	//}
	// faster way:
	iL[dim0] = (s+1)/2 * deg[dim0];

	for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++) {
	  for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++) {
	    // find the right and left glops volume indices

	    int iR[3] = {iL[0], iL[1], iL[2]};
	    iR[dim0] = (1-s)/2 * deg[dim0];

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
	      // change the glop index if we are on the first cell
	      ref_pg_vol(deg, nraf, ipgL, xref, &wpg3, NULL);
	      // mapping from the ref glop to the physical glop
	      schnaps_real dtau[3][3], codtau[3][3];
	      schnaps_ref2phy(fd->physnode,
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
	    // numerical flux from the left and right state and
	    // normal vector
	    schnaps_real wL[m], wR[m], flux[m];
	    for(int iv = 0; iv < m; iv++) {
	      int imemL = fd->varindex(deg, nraf, m, ipgL, iv);
	      int imemR = fd->varindex(deg, nraf, m, ipgR, iv);
	      //wL[iv] = fd->wn[imemL];
	      // the left flux is implicit and computed elsewhere
	      wL[iv] = 0;
	      wR[iv] = fd->wn[imemR];
	    }
	    fd->model.NumFlux(wL, wR, vnds, flux);

	    // subcell ref surface glop weight
	    schnaps_real wpg
	      = wglop(deg[dim1], iL[dim1])
	      * wglop(deg[dim2], iL[dim2]);

	    /* printf("vnds %f %f %f flux %f wpg %f\n", */
	    /* 	 vnds[0], vnds[1], vnds[2], */
	    /* 	 flux[0], wpg); */

	    // finally distribute the flux on the L side
	    for(int iv = 0; iv < m; iv++) {
	      int imemL = fd->varindex(deg, nraf, fd->model.m, ipgL, iv);
	      // the s sign is here because vnds
	      // is not necessarily oriented from
	      // inside to outside of the subcell
	      fd->res[imemL] -= s * flux[iv] * wpg;
	    }
	  }
	}
      }
    }
  }

}



void subcell_solve(field* fd, int isub);

void subcell_solve(field* f, int isub){


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

  int sc_npg = npg[0] * npg[1] * npg[2];


  const int wsub_size = sc_npg * m;


  // local "rigidity" matrix in triplet form
  //cs *T;

  // assembly of the flux part of the matrix

  // assembly of the volume part of the matrix

  // local "rigidity" matrix in compressed-column form
  //cs *K;
  //K = cs_triplet(T);
  /* cs_lusol(K,  // matrix */
  /* 	   0, // order */
  /* 	   f->res + isub * wsub_size, // rhs and solution... */
  /* 	   0.001   // tolerance */
  /* 	   ); */
}

void FieldDownwindSolve(field *fd,schnaps_real theta, schnaps_real dt){


  // build the downwind graph
  schnaps_real vit[3] = {1, 1, 0};

  int nraf[3] = {fd->raf[0],
		 fd->raf[1],
		 fd->raf[2]};

  int nbverts = nraf[0] * nraf[1] * nraf[2] + 2;


  int* subcell_order = malloc(nbverts * 3 * sizeof(int));
  build_subcells_downwind_graph(fd, vit, subcell_order);

  // loop on the subcells in the good order
  for( int isubt = 0; isubt < nbverts - 2; isubt++){
    // +1 because the ghost upwind vertex is always the first
    // and the ghost downwind vertex is always the last
    int isub = subcell_order[isubt + 1];
    // accumulate the fluxes of the neighbour subcells
    // into the residual
    get_neighbor_fluxes(fd, isub);

    // assembly and solve the subcell linear system
    subcell_solve(fd, isub);

  } // end loop

    // free memory
  free(subcell_order);


}




void ThetaTimeScheme(Simulation *simu, schnaps_real tmax, schnaps_real dt){

  LinearSolver solver_implicit;
  LinearSolver solver_explicit;
  // TODO Matrix Storage should be passed either as additionnal arg or through a simu data member.
  MatrixStorage ms= SKYLINE;
  //MatrixStorage ms= KLU_CSR;
  schnaps_real theta=0.5;
  simu->dt=dt;
  int itermax=tmax/simu->dt;
  simu->itermax_rk=itermax;
  InitImplicitLinearSolver(simu, &solver_implicit,ms);
  InitImplicitLinearSolver(simu, &solver_explicit,ms);
  schnaps_real *res = calloc(simu->wsize, sizeof(schnaps_real));
  simu->tnow=0;
  for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
    simu->fd[ie].tnow=simu->tnow;
  }
  time_t start;
  int iter=0;
  for(int tstep=0;tstep<simu->itermax_rk;tstep++){

    if (iter <= 1) start = time(NULL);
    iter++;

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
    if (simu->pre_dtfields !=NULL){
      simu->pre_dtfields(simu);
    }
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

    solver_implicit.solver_type=LU;
    Advanced_SolveLinearSolver(&solver_implicit,simu);

    for(int i=0;i<solver_implicit.neq;i++){
      simu->w[i]=solver_implicit.sol[i];
    }
    //
    if (simu->post_dtfields !=NULL){
      simu->post_dtfields(simu);
    }
    if (simu->update_after_rk !=NULL){
      simu->update_after_rk(simu,simu->w);
    }
    //
    int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, tstep+1, simu->itermax_rk, dt);
    if (iter == 1)
      printf("Elapsed time first iter=%f\n", (double) (time(NULL) -start));
  }
  printf("Elapsed time=%f\n", (double) (time(NULL) -start));
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


void InternalLocalCoupling(field *f, int itest)
{

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
		  int imemL = f->varindex(deg, nraf, m, ipgL, iv1);

		  int q[3] = {p[0], p[1], p[2]};
		  // loop on the direction dim0 on the "cross"
		  for(int iq = 0; iq < npg[dim0]; iq++) {
		    q[dim0] = (p[dim0] + iq) % npg[dim0];

		    int ipgR = offsetL + q[0] + npg[0] * (q[1] + npg[1] * q[2]);
		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv2);
		      if (itest ==0) IsNonZero(f->solver, imemL, imemR);
		      if (itest ==0) IsNonZero(f->rmat, imemL, imemR);
		      if (itest ==1) AddLinearSolver(f->solver, imemL, imemR,1);
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



void FluxLocalCoupling(field *f,int itest)
{

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
		  int imemL = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv1);

		  // finally distribute the flux on the two sides
		  for(int iv2 = 0; iv2 < m; iv2++) {
		    int imemR = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2);
		    if (itest ==0) {
		      IsNonZero(f->solver, imemL, imemR);
		      IsNonZero(f->rmat, imemL, imemR);
		    }
		    if (itest ==1) {
		      AddLinearSolver(f->solver, imemL, imemR,1);
		      AddLinearSolver(f->solver, imemR, imemL,1);
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


void InternalLocalAssembly(field *f, schnaps_real theta, schnaps_real dt)
{

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
	    imems[pos++] = f->varindex(f->deg,f->raf,f->model.m, offsetL + p, im);
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
		    int imemL = f->varindex(deg, nraf, m, ipgL, iv1);
		    for(int iv = 0; iv < m; iv++) {
		      wL[iv] = (iv == iv1);
		    }

		    f->model.NumFlux(wL, wL, dphiL, flux);

		    int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		    for(int iv2 = 0; iv2 < m; iv2++) {
		      schnaps_real val =  - flux[iv2] * wpgL;
		      int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv2);
		      /* AddLinearSolver(f->solver, imemR, imemL, theta * dt * val); */
		      /* AddLinearSolver(f->rmat, imemR, imemL, - dt * val); */
		      AddLinearSolver(f->solver, imemR, imemL, theta * dt * val);
		      AddLinearSolver(f->rmat, imemR, imemL, -(1-theta) * dt * val);
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
		      } //iv2
		    }// iv1
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

void FluxLocalAssembly(field* f,schnaps_real theta, schnaps_real dt)
{


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
		  int imem1 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv1);

		  f->model.NumFlux(wL, wR, vnds, flux);

		  for(int iv2 = 0; iv2 < m; iv2++) {
		    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv2);
		    schnaps_real val = flux[iv2] * wpg;
		    /* AddLinearSolver(f->solver, imem2, imem1, theta * dt * val); */
		    /* AddLinearSolver(f->rmat, imem2, imem1, -dt * val); */
		    AddLinearSolver(f->solver, imem2, imem1, theta * dt * val);
		    AddLinearSolver(f->rmat, imem2, imem1, -(1-theta) * dt * val);

		    imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2);
		    val = flux[iv2] * wpg;
		    /* AddLinearSolver(f->solver, imem2, imem1, -theta * dt * val); */
		    /* AddLinearSolver(f->rmat, imem2, imem1, dt * val); */
		    AddLinearSolver(f->solver, imem2, imem1, theta * dt * (-val));
		    AddLinearSolver(f->rmat, imem2, imem1, -(1-theta) * dt * (-val));
		  }

		  for(int iv = 0; iv < m; iv++) {
		    wL[iv] = 0;
		    wR[iv] = (iv == iv1);
		  }
		  imem1 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv1);


		  f->model.NumFlux(wL, wR, vnds, flux);

		  for(int iv2 = 0; iv2 < m; iv2++) {
		    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv2);
		    schnaps_real val =  flux[iv2] * wpg;
		    /* AddLinearSolver(f->solver, imem2, imem1, theta * dt * val); */
		    /* AddLinearSolver(f->rmat, imem2, imem1, -dt * val); */
		    AddLinearSolver(f->solver, imem2, imem1, theta * dt * val);
		    AddLinearSolver(f->rmat, imem2, imem1, -(1-theta) * dt * val);

		    imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2);
		    val =  flux[iv2] * wpg;
		    /* AddLinearSolver(f->solver, imem2, imem1, -theta * dt * val); */
		    /* AddLinearSolver(f->rmat, imem2, imem1, dt *val); */
		    AddLinearSolver(f->solver, imem2, imem1, theta * dt * (-val));
		    AddLinearSolver(f->rmat, imem2, imem1, -(1-theta) * dt * (-val));
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

void MassLocalAssembly(field *f)
{


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
      int imem = f->varindex(deg, nraf, m, ipg, iv1);
      schnaps_real val = wpg * det;
      AddLinearSolver(f->solver, imem, imem,val);
      AddLinearSolver(f->rmat, imem, imem,  val);
      //printf("val local imp =%f\n",val);
    }
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
void SourceLocalAssembly_C(void *buffers[], void *cl_arg);

void SourceLocalAssembly_SPU(field *f, schnaps_real theta, schnaps_real dt){

  f->dt = dt;
  f->theta = theta;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet SourceLocalAssembly...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = SourceLocalAssembly_C;
    codelet.nbuffers = 1;
    codelet.modes[0] = STARPU_RW;
    codelet.name="SourceLocalAssembly";
  }
  /* if (!f->local_source_cl_init){ */
  /*   printf("register rhs %d %d...\n",f->wsize,f->solver->neq); */
  /*   f->local_source_cl_init = true; */
  /*   starpu_vector_data_register(&(f->rhs_handle), // mem handle */
  /* 				0, // location: CPU */
  /* 				(uintptr_t)(f->solver->rhs), // vector location */
  /* 				f->wsize,  // size */
  /* 				sizeof(real));  // type */
  /* register_data_arbiter(f->rhs_handle); */
  /*   printf("end register...\n"); */
  /* } */


  if(f->model.Source != NULL) {

    /* int m; */
    /* int deg[3]; */
    /* int raf[3]; */
    /* real physnode[20][3]; */
    /* varindexptr Varindex; */
    /* sourceptr Source; */
    /* real tnow,dt,theta; */

    void* arg_buffer;
    size_t arg_buffer_size;

    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			     STARPU_VALUE, &f->model.m, sizeof(int),
			     STARPU_VALUE, f->deg, 3 * sizeof(int),
			     STARPU_VALUE, f->raf, 3 * sizeof(int),
			     STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			     STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			     STARPU_VALUE, &f->model.Source, sizeof(sourceptr),
			     STARPU_VALUE, &f->tnow, sizeof(schnaps_real),
			     STARPU_VALUE, &dt, sizeof(schnaps_real),
			     STARPU_VALUE, &theta, sizeof(schnaps_real),
			     0);


    task = starpu_task_create();
    task->cl = &codelet;
    task->cl_arg = arg_buffer;
    task->cl_arg_size = arg_buffer_size;
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[0] = sky_spu->rhs_handle;


    int ret = starpu_task_submit(task);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
    //printf("task submitted\n");

    /* void* buffers[1]; */
    /* buffers[0] = f->solver->rhs; */
    //SourceLocalAssembly_C(buffers, f);


  }
}


void SourceLocalAssembly(field *f, schnaps_real theta, schnaps_real dt){

  f->dt = dt;
  f->theta = theta;


  if(f->model.Source != NULL) {

    const int m = f->model.m;

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
      // we should extract wL and add a buffer to the codelet TODO
      f->model.Source(xphy, f->tnow, wL, source);

      for(int iv1 = 0; iv1 < m; iv1++) {
	int imem = f->varindex(f->deg, f->raf, m, ipg, iv1);
	schnaps_real val = source[iv1] * wpg * det;
	f->solver->rhs[imem] += f->theta * f->dt * val;
      }
    }



  }
}


void SourceLocalAssembly_C(void *buffers[], void *cl_arg) {

  //field *f = cl_arg;
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  schnaps_real tnow,dt,theta;
  varindexptr Varindex;
  sourceptr Source;

  starpu_codelet_unpack_args(cl_arg,
			     &m, deg, raf, physnode, &Varindex, &Source,
			     &tnow, &dt, &theta);

  free(cl_arg);

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffers[0];
  schnaps_real* rhs = (schnaps_real *)STARPU_VECTOR_GET_PTR(rhs_v);

  for(int ipg = 0; ipg < NPG(deg, raf); ipg++) {
    schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(deg, raf, ipg, xpgref, &wpg, NULL);
    schnaps_ref2phy(physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    schnaps_real det = dot_product(dtau[0], codtau[0]);
    schnaps_real wL[m], source[m];
    Source(xphy, tnow, wL, source);

    for(int iv1 = 0; iv1 < m; iv1++) {
      int imem = Varindex(deg, raf, m, ipg, iv1);
      schnaps_real val = source[iv1] * wpg * det;
      rhs[imem] += theta * dt * val;
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


void InterfaceLocalAssembly(Interface *inter,  schnaps_real theta, schnaps_real dt)
{


  field* fL = inter->fL;
  field* fR = inter->fR;

  int locfaL = inter->locfaL;

  const unsigned int m = fL->model.m;


  for(int ipgf = 0; ipgf < NPGF(fL->deg, fL->raf, locfaL); ipgf++) {

    schnaps_real xpgref[3], xpgref_in[3], wpg;

    int ipgL = inter->vol_indexL[ipgf];
    schnaps_real vnds[3] = {
      inter->vnds[3 * ipgf + 0] * inter->wpg[ipgf],
      inter->vnds[3 * ipgf + 1] * inter->wpg[ipgf],
      inter->vnds[3 * ipgf + 2] * inter->wpg[ipgf],
    };
    schnaps_real* xpg = inter->xpg + 3 * ipgf;

    int offsetL = 0;
    int offsetR = 0;

    if (fR != NULL) {

      int ipgR = inter->vol_indexR[ipgf];
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
	  schnaps_real val =  flux[iv2];
	  AddLinearSolver(fL->solver, imem2, imem1, theta * dt * val);
	  AddLinearSolver(fL->rmat, imem2, imem1, -(1-theta) * dt * val);

	}

	for(int iv = 0; iv < m; iv++) {
	  wL[iv] = 0;
	  wR[iv] = (iv == iv1);
	}
	imem1 = fL->varindex(fL->deg, fL->raf, fL->model.m, ipgR, iv1) + offsetR;


	fL->model.NumFlux(wL, wR, vnds, flux);

	for(int iv2 = 0; iv2 < m; iv2++) {

	  int imem2 = fR->varindex(fR->deg, fR->raf, fR->model.m, ipgR, iv2) + offsetR;
	  schnaps_real val =  flux[iv2];
	  AddLinearSolver(fR->solver, imem2, imem1, theta * dt * (-val));
	  AddLinearSolver(fR->rmat, imem2, imem1, -(1-theta) * dt * (-val));
	}
      }
    }
    else{ // case of a boundary condition


      // the boundary flux is an affine function
      schnaps_real flux0[m], wL[m];
      for(int iv = 0; iv < m; iv++) {
	wL[iv] = 0;
      }
      fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux0);

      for(int iv1 = 0; iv1 < m; iv1++) {
	int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv1);

	for(int iv = 0; iv < m; iv++) {
	  wL[iv] = (iv == iv1);
	}

	schnaps_real flux[m];
	fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux);

	for(int iv2 = 0; iv2 < m; iv2++) {
	  // The basis functions is also the gauss point index
	  int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2);
	  schnaps_real val =  (flux[iv2]-flux0[iv2]);
	  AddLinearSolver(fL->solver, imem2, imem1, theta * dt * val);
	  AddLinearSolver(fL->rmat, imem2, imem1,  -(1-theta) * dt * val);
	  //printf("val=%f",val);
	}
      } // iv1

    } // ipgf


  } // if (fR == NULL)

}


void FieldResidual(field *fd,schnaps_real theta, schnaps_real dt){

  DGSubCellInterface(fd, fd->wn, fd->res);
  DGVolume(fd, fd->wn, fd->res);
  DGSource(fd, fd->wn, fd->res);

  // TODO: think and multiply by dt and/or +/-(1-theta)

}

void InitImplicitJFLinearSolver(Simulation *simu, JFLinearSolver *solver){

  int neq = simu->wsize;

  Solver st = LU;
  InitJFLinearSolver(solver,neq,&st);
  solver->MatVecProduct=MatVecJacobianFree;
  solver->NonlinearVector_computation=ImplicitNonlinearVector_computation;
  solver->soln=simu->w;
  solver->iter_max=10000;

  solver->tol=_SMALL*10;
  solver->eps=0.00001;

}

void ImplicitNonlinearVector_computation(Simulation * simu,void* lsol,schnaps_real * solvector,schnaps_real *nlvector){


  NonlinearThetaVector_assembly(simu,solvector,nlvector,simu->theta,simu->dt);

}


void AssemblyImplicitJFLinearSolver(Simulation *simu, JFLinearSolver *solver, schnaps_real dt){
  schnaps_real * rhs_implicit;
  schnaps_real * rhs_explicit;

  rhs_implicit = calloc(simu->wsize, sizeof(schnaps_real));
  rhs_explicit = calloc(simu->wsize, sizeof(schnaps_real));

  NonlinearThetaVector_assembly(simu,simu->w,rhs_explicit,-(1.0-simu->theta),simu->dt);

  simu->tnow=simu->tnow+simu->dt;
  for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
    simu->fd[ie].tnow=simu->tnow;
  }

  NonlinearThetaVector_assembly(simu,simu->w,rhs_implicit,simu->theta,simu->dt);

  DtFields_old(simu, simu->w, simu->dtw);

  for(int i=0;i<solver->neq;i++){
    solver->rhs[i]=-rhs_implicit[i]+rhs_explicit[i];
  }

  free(rhs_implicit);
  free(rhs_explicit);

}

void ThetaTimeScheme_WithJF(Simulation *simu, schnaps_real tmax, schnaps_real dt){

  JFLinearSolver solver_implicit;

  schnaps_real theta=0.5;
  simu->dt=dt;
  simu->theta=theta;

  int itermax=tmax/simu->dt;
  simu->itermax_rk=itermax;

  InitImplicitJFLinearSolver(simu, &solver_implicit);
  simu->tnow=0;
  for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
    simu->fd[ie].tnow=simu->tnow;
  }

  for(int tstep=0;tstep<simu->itermax_rk;tstep++){

    AssemblyImplicitJFLinearSolver(simu,&solver_implicit,simu->dt);


    SolveJFLinearSolver(&solver_implicit,simu);

    for(int i=0;i<solver_implicit.neq;i++){
     simu->w[i]=simu->w[i]+solver_implicit.sol[i];
    }

    int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, tstep, simu->itermax_rk, dt);
  }


}

void NonlinearThetaVector_assembly(Simulation * simu,schnaps_real * solvector,schnaps_real *nlvector,schnaps_real theta, schnaps_real dt){
  const unsigned int m = simu->fd[0].model.m;
    int fsize =  simu->wsize / simu->macromesh.nbelems;

    ///// Compute the nonlinear vecteur M uL - coef*dt F(uL) +  coef*dt F(uL,uR) - coef*dt M B(uL)
    ///// which correspond to thediscretization of dt u + div(F(u))=B(u)

    /////////// We compute the fllux + F(uL,uR) between the macro cell interface//////
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

    schnaps_real *fwL = solvector + offsetL;
    schnaps_real *fwR = solvector + offsetR;

    schnaps_real *nlvL = nlvector + offsetL;
    schnaps_real *nlvR = nlvector + offsetR;

    for(int ipgfL = 0; ipgfL < NPGF(fL->deg, fL->raf, locfaL); ipgfL++) {

      schnaps_real xpgref[3], xpgref_in[3], wpg;

      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipgL = ref_pg_face(fL->deg, fL->raf, locfaL, ipgfL, xpgref, &wpg, xpgref_in);

      schnaps_real flux[m];
      schnaps_real wL[m];

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

	schnaps_real wR[m];
	for(int iv = 0; iv < m; iv++) {
	  int imemL = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv);
	  wL[iv] = fwL[imemL];
	  int imemR = fR->varindex(fR->deg, fR->raf,fR->model.m, ipgR, iv);
	  wR[iv] = fwR[imemR];
	}

	fL->model.NumFlux(wL, wR, vnds, flux);

	// Add flux to both sides
	for(int iv = 0; iv < m; iv++) {
	  // The basis functions is also the gauss point index
	  int imemL = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv);
	  int imemR = fR->varindex(fR->deg, fR->raf,fR->model.m, ipgR, iv);
	  nlvL[imemL] += theta * dt * flux[iv] * wpg;
	  nlvR[imemR] -= theta * dt * flux[iv] * wpg;
	  //nlvL[imemL] -= theta * dt * flux[iv] * wpg;
	  //nlvR[imemR] += theta * dt * flux[iv] * wpg;
	}

      }
      else
	{ // The point is on the boundary.
	  for(int iv = 0; iv < m; iv++) {
	    int imemL = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv);
	    wL[iv] = fwL[imemL];
	  }


	  fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux);

	  for(int iv = 0; iv < m; iv++) {
	// The basis functions is also the gauss point index
	    int imemL = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv);
	    nlvL[imemL] += theta * dt * flux[iv] * wpg;
	    //nlvL[imemL] -= theta * dt * flux[iv] * wpg;
	  }
	}

    }

 }

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd + ie;
    schnaps_real * w = solvector + ie * fsize;
    schnaps_real * nlv = nlvector + ie * fsize;

    const int nraf[3] = {f->raf[0],
		       f->raf[1],
		       f->raf[2]};
    const int deg[3] = {f->deg[0],
		      f->deg[1],
		      f->deg[2]};
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};

    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];



    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};

	  // Get the left subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;

        /////////// We compute the fllux + F(uL,uR) between the sub-cell interface//////
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

		  // numerical flux from the left and right state and
		  // normal vector
		  schnaps_real wL[m], wR[m], flux[m];
		  for(int iv = 0; iv < m; iv++) {
		    // TO DO change the varindex signature
		    int imemL = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv);
		    int imemR = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv);
		    // end TO DO
		    wL[iv] = w[imemL];
		    wR[iv] = w[imemR];
		  }
		  f->model.NumFlux(wL, wR, vnds, flux);

		  // subcell ref surface glop weight
		  schnaps_real wpg
		  = wglop(deg[dim1], iL[dim1])
		    * wglop(deg[dim2], iL[dim2]);

		  /* printf("vnds %f %f %f flux %f wpg %f\n", */
		/* 	 vnds[0], vnds[1], vnds[2], */
		/* 	 flux[0], wpg); */

		// finally distribute the flux on the two sides
		  for(int iv = 0; iv < m; iv++) {
		    // TO DO change the varindex signature
		    int imemL = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv);
		    int imemR = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv);
		    // end TO DO
		     nlv[imemL] += theta * dt * flux[iv] * wpg;
		    nlv[imemR] -= theta * dt * flux[iv] * wpg;
		    //nlv[imemL] -= theta * dt * flux[iv] * wpg;
		    //nlv[imemR] += theta * dt * flux[iv] * wpg;
		  }

		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop

          /////////// We compute the volume term - F(uL) in the subs-cell//////
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
	    imems[pos++] = f->varindex(f->deg,f->raf,f->model.m, offsetL + p, im);
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
		for(int iv = 0; iv < m; iv++) {
		  ///int imemL = f->varindex(f_interp_param, ie, ipgL, iv);
		  wL[iv] = w[imems[m * (ipgL - offsetL) + iv]];
		}
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

		  f->model.NumFlux(wL, wL, dphiL, flux);

		  int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		  for(int iv = 0; iv < m; iv++) {
		    int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv);
		    int temp = m * (ipgR - offsetL) + iv;
		    assert(imemR == imems[temp]);
		    nlv[imems[temp]] -= theta * dt * flux[iv] * wpgL;
		    //nlv[imems[temp]] += theta * dt * flux[iv] * wpgL;
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

	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

     /////////// We compute the source term - B(uL) in the subs-cell and the mass term + M uL//////
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      schnaps_real det = dot_product(dtau[0], codtau[0]);  //// temp !!!!!!!!!!!!!!!
      schnaps_real wL[m], source[m];
      for(int iv = 0; iv < m; ++iv){
	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	wL[iv] = w[imem];
      }

      f->model.Source(xphy, f->tnow, wL, source);
      // printf("tnow=%f\n",f->tnow);

      for(int iv = 0; iv < m; ++iv) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	nlv[imem] -= theta * dt *source[iv] * det * wpg; //nlv[imem] += theta * dt *source[iv] * det * wpg;
	nlv[imem] += wL[iv] * (wpg * det);
	//nlv[imem] /=(wpg * det); //+=wL[iv] * (wpg * det);
      }
    }


  }


}

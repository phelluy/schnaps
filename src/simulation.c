#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

#include "field_cl.h"

void EmptySimulation(Simulation *simu){

#ifdef _WITH_OPENCL
   sprintf(numflux_cl_name,"%s",""); // FIXME: move to field struct.
  sprintf(cl_buildoptions,"%s",""); // FIXME: move to field struct.
  simu->sourcename_cl = malloc(1024 * sizeof(char)); // TODO set to NULL
  sprintf(simu->sourcename_cl,"%s"," "); // FIXME: remove
#endif
  simu->w_handle = NULL;
  simu->dtw_handle = NULL;
  simu->res_handle = NULL;
}

void RegisterSimulation_SPU(Simulation *simu){

  if (simu->w_handle == NULL){
    simu->w_handle = malloc(simu->macromesh.nbelems
			    * sizeof( starpu_data_handle_t));
    simu->dtw_handle = malloc(simu->macromesh.nbelems
			    * sizeof( starpu_data_handle_t));
    simu->res_handle = malloc(simu->macromesh.nbelems
			    * sizeof( starpu_data_handle_t));
    for(int ie = 0; ie < simu->macromesh.nbelems; ie++) {
      Registerfield_SPU(simu->fd + ie);
      simu->w_handle[ie] = simu->fd[ie].wn_handle;
      simu->dtw_handle[ie] = simu->fd[ie].dtwn_handle;
      simu->res_handle[ie] = simu->fd[ie].res_handle;

    }
  }

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++) {
    RegisterInterface_SPU(simu->interface + ifa);      
  }

  
}

void UnregisterSimulation_SPU(Simulation *simu){

  if (simu->w_handle != NULL){
    for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
      Unregisterfield_SPU(simu->fd + ie);
    }
    free(simu->w_handle);
    free(simu->dtw_handle);
    free(simu->res_handle);
    simu->w_handle = NULL;
    simu->dtw_handle = NULL;
    simu->res_handle = NULL;
  }

  for (int ifa = 0; ifa < simu->macromesh.nbfaces; ++ifa) {
    Interface *inter = simu->interface + ifa;
    if (inter->starpu_registered == true)
      UnregisterInterface_SPU(inter);
  }
}



void InitInterfaces(Simulation *simu){


  int nb_interfaces = simu->macromesh.nbfaces;
  printf("Alloc %d interfaces of size %zd...\n",simu->macromesh.nbfaces, sizeof(Interface));
  simu->interface = malloc(sizeof(Interface) * nb_interfaces);
  Interface *inter = simu->interface;

  printf("Compute interface points and normals: ");
  fflush(stdout);
   for(int ifa = 0; ifa < nb_interfaces; ifa++){
    printf("\%d ",ifa);
    if (ifa%100 == 0) fflush(stdout);
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    field *fR = NULL;
    int locfaR = -1;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      locfaR = simu->macromesh.face2elem[4 * ifa + 3];
    }

    inter[ifa].fL = fL;
    inter[ifa].fR = fR;

    inter[ifa].ieL = ieL;
    inter[ifa].ieR = ieR;

    inter[ifa].locfaL = locfaL;
    inter[ifa].locfaR = locfaR;

    int npgfL = NPGF(fL->deg, fL->raf, locfaL);
    inter[ifa].npgL = npgfL;

    int npgfR = 0;
    if (ieR >= 0) {
      npgfR = NPGF(fR->deg, fR->raf, locfaR);
    }
    inter[ifa].npgR = npgfR;

    inter[ifa].wsizeL = npgfL * fL->model.m;
    inter[ifa].wsizeR = 0;
    if (fR != NULL) inter[ifa].wsizeR = npgfR * fR->model.m;

    inter[ifa].wL = calloc(inter[ifa].wsizeL, sizeof(schnaps_real));
    inter[ifa].wR = NULL;
    if (npgfR > 0)  inter[ifa].wR =
		      calloc(inter[ifa].wsizeR, sizeof(schnaps_real));

    inter[ifa].vol_indexL = malloc(sizeof(int) * npgfL);
    inter[ifa].vol_indexR = NULL;
    if (fR != NULL) inter[ifa].vol_indexR = malloc(sizeof(int) * npgfR);

    inter[ifa].vnds = malloc(3* sizeof(schnaps_real) * npgfL);
    inter[ifa].xpg = malloc(3* sizeof(schnaps_real) * npgfL);
    inter[ifa].wpg = malloc(sizeof(schnaps_real) * npgfL);

    schnaps_real xpgref[3], xpgref_in[3],xpg[3], xpg_in[3];
    for(int ipgf = 0; ipgf < npgfL; ipgf++){
      schnaps_real wpg;
      schnaps_real vnds[3];
      int ipgv =
	ref_pg_face(fL->deg, fL->raf, locfaL, ipgf,
		    xpgref, &wpg, xpgref_in);
      inter[ifa].vol_indexL[ipgf] = ipgv;
      {
	schnaps_real codtau[3][3], dtau[3][3];
	schnaps_ref2phy(fL->physnode,
			xpgref_in,
			NULL, locfaL, // dpsiref, ifa
			xpg_in, dtau,
			codtau, NULL, vnds); // codtau, dpsi, vnds
	schnaps_ref2phy(fL->physnode,
			xpgref,
			NULL, locfaL, // dpsiref, ifa
			xpg, dtau,
			codtau, NULL, vnds); // codtau, dpsi, vnds
	PeriodicCorrection(xpg_in,fL->period);
      }
      /* printf("ifa=%d ieL=%d ieR=%d ipgf=%d  */
      /*        vnds=%f %f %f xpg=%f %f %f\n", */
      /* 	      ifa,ieL,ieR,ipgf,vnds[0],vnds[1],vnds[2], */
      /* 	      xpg[0],xpg[1],xpg[2]); */

      inter[ifa].vnds[3 * ipgf + 0] = vnds[0];
      inter[ifa].vnds[3 * ipgf + 1] = vnds[1];
      inter[ifa].vnds[3 * ipgf + 2] = vnds[2];
      inter[ifa].xpg[3 * ipgf + 0] = xpg[0];
      inter[ifa].xpg[3 * ipgf + 1] = xpg[1];
      inter[ifa].xpg[3 * ipgf + 2] = xpg[2];
      inter[ifa].wpg[ipgf] = wpg;

      if (fR != NULL){
	schnaps_phy2ref(fR->physnode, xpg_in, xpgref_in);
	inter[ifa].vol_indexR[ipgf] =
	  ref_ipg(fR->deg, fR->raf, xpgref_in);
      }
    }

    /* if (fR != NULL){ */
    /*   for(int ipgf = 0; ipgf < npgfR; ipgf++){ */
    /* 	 int ipgv1 = */
    /* 	   ref_pg_face(fR->deg, fR->raf, locfaR, ipgf, */
    /* 		       xpgref, NULL, NULL); */
    /* 	 int ipgv2 = ref_ipg(fR->deg, fR->raf, xpgref); */
    /* 	 assert(ipgv1 == ipgv2); */
    /* 	 inter[ifa].vol_indexR[ipgf] = ipgv2; */
    /*   } */

    inter[ifa].starpu_registered = false;

    //RegisterInterface_SPU(inter + ifa);
  }
  printf("\n");
}


void InitSimulation(Simulation *simu, MacroMesh *mesh,
		    int *deg, int *raf, Model *model){


  simu->macromesh = *mesh;

  simu->tnow = 0;

  simu->iter_time_rk =0;

  simu->fd = malloc(mesh->nbelems * sizeof(field));

  field *fd = simu->fd;

  int field_size = NPG(deg,raf) * model->m;

  simu->wsize = field_size * mesh->nbelems;

  printf("field_size = %d reals simusize = %d reals\n",field_size,simu->wsize);
  printf("Memory for saving one time step: %f MB\n", simu->wsize * sizeof(schnaps_real) / 1e6);
  printf("Memory for RK2 algo: %f MB\n", 3 * simu->wsize * sizeof(schnaps_real) / 1e6);
  printf("NPG = %d m = %d\n",NPG(deg,raf),model->m);
  simu->w = calloc(simu->wsize, sizeof(schnaps_real));
  simu->dtw = calloc(simu->wsize, sizeof(schnaps_real));
  simu->res = calloc(simu->wsize, sizeof(schnaps_real));

  schnaps_real *w = simu->w;
  schnaps_real *dtw = simu->dtw;
  schnaps_real *res = simu->res;

  schnaps_real physnode[20][3];

  simu->hmin = FLT_MAX;
  simu->vcfl = 1;

#ifdef _WITH_OPENCL
  simu->use_source_cl = false;
  if (model->Source != NULL) simu->use_source_cl = true;
#endif

  printf("Init macro elem");
  fflush(stdout);
  for(int ie=0; ie < mesh->nbelems; ++ie){
    printf(" %d",ie);
    if (ie%100 == 0) fflush(stdout);
    for(int inoloc = 0; inoloc < 20; ++inoloc){
      int ino = mesh->elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = mesh->node[3 * ino + 0];
      physnode[inoloc][1] = mesh->node[3 * ino + 1];
      physnode[inoloc][2] = mesh->node[3 * ino + 2];
    }

    init_empty_field(fd + ie);

    fd[ie].period[0] = simu->macromesh.period[0];
    fd[ie].period[1] = simu->macromesh.period[1];
    fd[ie].period[2] = simu->macromesh.period[2];

    fd[ie].id = ie;

    fd[ie].store_det = false;

    Initfield(fd + ie, *model, physnode, deg, raf,
	      w + ie * field_size, dtw + ie * field_size
	      );

    fd[ie].res = res + ie * field_size;

    simu->hmin = simu->hmin > fd[ie].hmin ? fd[ie].hmin : simu->hmin;
  }
  printf("\n");

  
  /* if (!starpu_is_init && starpu_use){ */
  /*   printf("StarPU registration...\n"); */
  /*   int ret; */
  /*   ret = starpu_init(NULL); */
  /*   init_global_arbiter(); */
  /*   assert(ret != -ENODEV) ; */
  /*   starpu_is_init = true; */

  /*   for(int ie=0; ie < mesh->nbelems; ++ie){ */
  /*     // starpu registration (must be called after init_data) */
  /*     Initfield_SPU(f); */
  /*   } */

    
  /* } */



  simu->pre_dtfields = NULL;
  simu->post_dtfields = NULL;
  simu->update_after_rk = NULL;
  simu->Diagnostics = NULL;
  simu->pic =NULL;
  simu->nb_diags = 0;


  // to do remove the legacy interp_param array
  field *f = simu->fd;
  simu->interp_param[0] = f->model.m;
  simu->interp_param[1] = f->deg[0];
  simu->interp_param[2] = f->deg[1];
  simu->interp_param[3] = f->deg[2];
  simu->interp_param[4] = f->raf[0];
  simu->interp_param[5] = f->raf[1];
  simu->interp_param[6] = f->raf[2];


#ifdef _WITH_OPENCL
  // opencl inits
  /* if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) { */
  /*   printf("OpenCL device not acceptable; OpenCL initialization disabled.\n"); */
  /* } else { */
    
    //init_field_cl(simu);  /// to do deactivate legacy opencl...
    
    
  //}
#endif // _WITH_OPENCL
  init_field_cl(simu);
  printf("Init interfaces...\n");
  InitInterfaces(simu);

}


void DisplaySimulation(Simulation *simu){

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie){
    printf("Field %d\n",ie);

    Displayfield(simu->fd + ie);
  }

}

// Save the results in the gmsh format typplot: index of the plotted
// variable int compare == true -> compare with the exact value.  If
// fieldname is NULL, then the fieldname is typpplot.
void PlotFields(int typplot, int compare, Simulation* simu, char *fieldname,
	       char *filename) {

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;


  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );
  


  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(schnaps_real));

  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);

  schnaps_real *value = malloc(nb_plotnodes * sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;

  // Nodes
  int npgv = NPG(deg, nraf);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L

    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {

	  for(int ino = 0; ino < 64; ino++) {
	    schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
			     hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
			     hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };

	    for(int ii = 0; ii < 3; ii++) {
	      assert(Xr[ii] < 1 +  1e-10);
	      assert(Xr[ii] > -1e-10);
	    }

	    schnaps_real Xphy[3];
	    schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

	    schnaps_real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};

	    value[nodecount] = 0;
	    schnaps_real testpsi = 0;
	    for(int ib = 0; ib < npgv; ib++) {
	      schnaps_real psi;
	      psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
	      testpsi += psi;
	      int vi = f->varindex(f->deg, f->raf, f->model.m, ib, typplot);
	      value[nodecount] += psi * f->wn[vi];
	      //printf("a");
	    }
	    assert(fabs(testpsi-1) < _SMALL);

	    // Compare with an exact solution
	    if (compare) {
	      schnaps_real wex[f->model.m];
	      f->model.ImposedData(Xphy, f->tnow, wex);
	      value[nodecount] -= wex[typplot];

	    }
	    nodecount++;
	    fprintf(gmshfile, "%d %f %f %f\n", nodecount,
		    Xplot[0], Xplot[1], Xplot[2]);
	  }
	}
      }
    }
  }

  fprintf(gmshfile, "$EndNodes\n");

  // Elements
  fprintf(gmshfile, "$Elements\n");
  fprintf(gmshfile, "%d\n",
	  simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);

  int elm_type = 92;
  int num_tags = 0;

  // fwrite((char*) &elm_type, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_elm_follow, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_tags, sizeof(int), 1, gmshfile);

  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
	  // Get the subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);

	  // Global subcell id
	  int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;

	  //fwrite((char*) &numelem, sizeof(int), 1, gmshfile);
	  fprintf(gmshfile, "%d ", numelem);
	  fprintf(gmshfile, "%d ", elm_type);
	  fprintf(gmshfile, "%d ", num_tags);

	  for(int ii = 0; ii < 64; ii++) {
	    int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
	    //fwrite((char*) &numnoe, sizeof(int), 1, gmshfile);
	    fprintf(gmshfile, "%d ", numnoe);
	  }
	  fprintf(gmshfile, "\n");
	}
      }
    }
  }

  fprintf(gmshfile, "$EndElements\n");

  // Now display data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else
    fprintf(gmshfile, "\"field: %s\"\n", fieldname);

  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);

  for(int ino = 0; ino < nb_plotnodes; ino++) {
    //fwrite(const void *ptr, size_t size_of_elements,
    // size_t number_of_elements, FILE *a_file);
    //fwrite((char*) &nodenumber, sizeof(int), 1, gmshfile);
    //fwrite((char*) &value, sizeof(real), 1, gmshfile);
    //fprintf(gmshfile, "%d %f\n", nodenumber, value);
    fprintf(gmshfile, "%d %f\n", ino + 1, value[ino]);
  }

  /* for(int i = 0;i < f->macromesh.nbelems;i++) { */
  /*   for(int ino = 0;ino < 20;ino++) { */
  /* 	int numnoe = elem2nodes[nnodes*i+ino]; */
  /* 	for(int ii = 0;ii < 3;ii++) { */
  /* 	  physnode[ino][ii] = node[3 * numnoe+ii]; */
  /* 	} */
  /*   } */

  /*   // data at the eight nodes */
  /*   for(int ii = 0;ii < 64;ii++) { */
  /* 	int nodenumber = 64*i + ii  + 1; */

  /* 	Xr[0] = (real) (hexa64ref[3 * ii+0]) / 3; */
  /* 	Xr[1] = (real) (hexa64ref[3 * ii + 1]) / 3; */
  /* 	Xr[2] = (real) (hexa64ref[3 * ii+2]) / 3; */

  /* 	Ref2Phy(physnode, */
  /* 		Xr, */
  /* 		NULL, */
  /* 		-1, */
  /* 		Xphy, */
  /* 		NULL, */
  /* 		NULL, */
  /* 		NULL, */
  /* 		NULL); */


  /* 	real value = 0; */
  /* 	for(int ib = 0;ib < npgv;ib++) { */
  /* 	  real psi; */
  /* 	  psi_ref(f->deg, f->raf, ib, Xr, &psi, NULL); */

  /* 	  int vi  =  f->varindex(f->interp_param, i, ib, typplot); */
  /* 	  value += psi * f->wn[vi]; */
  /* 	} */

  /* 	// compare with an */
  /* 	// exact solution */
  /*     if (compare) { */
  /*       real wex[f->model.m]; */
  /*       f->model.ImposedData(Xphy, f->tnow, wex); */
  /*       value -= wex[typplot]; */
  /*     } */


  /* 	//fwrite(const void *ptr, size_t size_of_elements, */
  /* 	// size_t number_of_elements, FILE *a_file); */
  /* 	//fwrite((char*) &nodenumber, sizeof(int), 1, gmshfile); */
  /* 	//fwrite((char*) &value, sizeof(real), 1, gmshfile); */
  /* 	//fprintf(gmshfile, "%d %f\n", nodenumber, value); */
  /* 	fprintf(gmshfile, "%d %f\n", nodenumber, value); */
  /*   } */

  /* } */

  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
    
  free(value);
}

void freeSimulation(Simulation* simu){

  for(int i=0;i<simu->macromesh.nbelems;i++){
    simu->fd[i].wn=NULL;
    simu->fd[i].dtwn=NULL;
  }


  free(simu->w);   //Already freed in field.
  free(simu->dtw); //Already freed in field.
  free(simu->fd);

  if(simu->Diagnostics != NULL){
    free(simu->Diagnostics);
  }
  if(simu->pic != NULL){
    free(simu->pic);
  }



}


// Apply the Discontinuous Galerkin approximation to compute the
// time derivative of the field
void DtFields_old(Simulation *simu, schnaps_real *w, schnaps_real *dtw) {

#ifdef _OPENMP
#pragma omp parallel
#endif

  //real *w = simu->fd[0].wn;
  //real *dtw = simu->fd[0].dtwn;

  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int iw = 0; iw < simu->wsize; iw++)
    dtw[iw] = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }
  // the field pointers must be updated
  /* for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) { */
  /*   simu->fd[ie].wn = w + ie * fsize; */
  /*   simu->fd[ie].dtwn = dtw + ie * fsize; */
  /* } */

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
    DGMacroCellInterface(locfaL,
    			 fL, offsetL, fR, offsetR,
    			 w, dtw);
  }


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    DGSubCellInterface(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGVolume(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGSource(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
    DGMass(simu->fd + ie, w + ie * fsize, dtw + ie * fsize);
  }
}


// Apply the Discontinuous Galerkin approximation to compute the
// time derivative of the field
void DtFields(Simulation *simu) {

#ifdef _OPENMP
#pragma omp parallel
#endif

  //real *w = simu->fd[0].wn;
  //real *dtw = simu->fd[0].dtwn;

  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int iw = 0; iw < simu->wsize; iw++)
    simu->dtw[iw] = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }
  // the field pointers must be updated
  /* for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) { */
  /*   simu->fd[ie].wn = w + ie * fsize; */
  /*   simu->fd[ie].dtwn = dtw + ie * fsize; */
  /* } */

  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    int ieL = simu->macromesh.face2elem[4 * ifa + 0];
    int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
    int ieR = simu->macromesh.face2elem[4 * ifa + 2];
    field *fL = simu->fd + ieL;
    field *fR = NULL;
    int offsetR = -1;
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      fR = simu->fd + ieR;
      offsetR = fsize * ieR;
    }
    DGMacroCellInterface(locfaL,
    			 fL, offsetL, fR, offsetR,
    			 simu->w, simu->dtw);
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    DGSubCellInterface(simu->fd + ie, simu->w + ie * fsize, simu->dtw + ie * fsize);
    DGVolume(simu->fd + ie, simu->w + ie * fsize, simu->dtw + ie * fsize);
    DGSource(simu->fd + ie, simu->w + ie * fsize, simu->dtw + ie * fsize);
    DGMass(simu->fd + ie, simu->w + ie * fsize, simu->dtw + ie * fsize);
  }
}



schnaps_real L2error(Simulation *simu) {

  schnaps_real error = 0;
  schnaps_real mean = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {

    field *f = simu->fd + ie;

    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      schnaps_real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	w[iv] = f->wn[imem];
      }

      schnaps_real wex[f->model.m];
      schnaps_real wpg, det;
      // Compute wpg, det, and the exact solution
      {
	schnaps_real xphy[3], xpgref[3];
	schnaps_real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
	schnaps_ref2phy(f->physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);

	// Get the exact value
	f->model.ImposedData(xphy, simu->tnow, wex);
      }

      for(int iv = 0; iv < f->model.m; iv++) {
	//for(int iv = 0; iv < 4; iv++) {   ///////error here for coil2d
	schnaps_real diff = w[iv] - wex[iv];
       error += diff * diff * wpg * det;
        mean += wex[iv] * wex[iv] * wpg * det;
	//printf("ie=%d ipg=%d iv=%d err=%f \n",ie,ipg,iv,diff);
        }
    }
  }
  //printf("errl2=%f\n",sqrt(error) / (sqrt(mean)  + 1e-14));
  return sqrt(error) / (sqrt(mean)  + 1e-14);
}

// Time integration by a second-order Runge-Kutta algorithm
void RK1(Simulation *simu, schnaps_real tmax){

  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;
  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt + 1;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  simu->tnow = 0;
  
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

   if(simu->nb_diags != 0) {
     simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
   }

   while(simu->tnow < tmax) {
     if (iter % freq == 0)
       printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

     //------ Pre time step ---- //
     if(simu->pre_dtfields != NULL) {
       simu->pre_dtfields(simu);
     }
    
     //------ Step 1 ---- //
     DtFields(simu);
     RK_out(simu,simu->w);
     simu->tnow += dt;

     //------ Post time step ---- //
     if(simu->update_after_rk != NULL){
       simu->update_after_rk(simu,simu->w);
     }

     iter++;
     simu->iter_time_rk = iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
}


// Time integration by a second-order Runge-Kutta algorithm
void RK2(Simulation *simu, schnaps_real tmax){

  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;
  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt + 1;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  simu->tnow = 0;

  schnaps_real *wn = calloc(simu->wsize, sizeof(schnaps_real));
  assert(wn);

  // FIXME: remove
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

   if(simu->nb_diags != 0) {
     simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
   }

  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);


    //--------- Pre RK step --------//   
    simu->dt=0.5*dt;
    
    if(simu->pre_dtfields != NULL) {
      simu->pre_dtfields(simu);
    }

    //--------- Begin RK ---------//
     //----- Rk Step 1 ----- //
    RK_Copy(simu,simu->w,wn);
    DtFields(simu);
    RK_in(simu);

    simu->tnow += 0.5 * dt;

     //----- Step 2 ----- //
    DtFields(simu);
    simu->dt=dt;
    RK_out(simu,wn);

    simu->tnow += 0.5 * dt;
    //---------- End RK ----------//

    
    //-------- Post RK step -------//  
    simu->dt=0.5*dt;

    if(simu->post_dtfields != NULL) {
       simu->post_dtfields(simu);
    }

    //------ Final post processing -----//
    simu->dt=dt;

    if(simu->update_after_rk != NULL){
      simu->update_after_rk(simu, simu->w);
    }

    iter++;
    simu->iter_time_rk = iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
    free(wn);
}



// Time integration by a fourth-order Runge-Kutta algorithm
void RK4_old(Simulation *simu, schnaps_real tmax)
{

  simu->dt = Get_Dt_RK(simu);

  schnaps_real dt = simu->dt;

  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  // Allocate memory for RK time-stepping
  schnaps_real *l1, *l2, *l3;
  l1 = calloc(simu->wsize, sizeof(schnaps_real));
  l2 = calloc(simu->wsize, sizeof(schnaps_real));
  l3 = calloc(simu->wsize, sizeof(schnaps_real));

  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

    if(simu->nb_diags != 0)
    simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));

  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    // l_1 = w_n + 0.5dt * S(w_n, t_0)
    DtFields_old(simu, simu->w, simu->dtw);
    RK_out_old(l1, simu->w, simu->dtw, 0.5 * dt, simu->wsize);

    simu->tnow += 0.5 * dt;

    // l_2 = w_n + 0.5dt * S(l_1, t_0 + 0.5 * dt)
    DtFields_old(simu, l1, simu->dtw);
    RK_out_old(l2, simu->w, simu->dtw, 0.5 * dt, simu->wsize);

    // l_3 = w_n + dt * S(l_2, t_0 + 0.5 * dt)
    DtFields_old(simu, l2, simu->dtw);
    RK_out_old(l3, simu->w, simu->dtw, dt, simu->wsize);

    simu->tnow += 0.5 * dt;

    // Compute S(l_3, t_0 + dt)
    DtFields_old(simu, l3, simu->dtw);
    RK4_final_inplace_old(simu->w, l1, l2, l3, simu->dtw, dt, simu->wsize);


     if(simu->update_after_rk != NULL){
      simu->update_after_rk(simu, simu->w);
    }

    iter++;
     simu->iter_time_rk=iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

  free(l3);
  free(l2);
  free(l1);
}



// Time integration by a fourth-order Runge-Kutta algorithm
void RK4(Simulation *simu, schnaps_real tmax)
{
  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;
  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  // Allocate memory for RK time-stepping
  schnaps_real *l1, *l2, *l3, *l4;
  l1 = calloc(simu->wsize, sizeof(schnaps_real));
  l2 = calloc(simu->wsize, sizeof(schnaps_real));
  l3 = calloc(simu->wsize, sizeof(schnaps_real));
  l4 = calloc(simu->wsize, sizeof(schnaps_real));

  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

    if(simu->nb_diags != 0)
    simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));

  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    // l_1 = w_n 
    RK_Copy(simu,simu->w,l1);

    // l_2 = l_1 + 0.5dt * S(l_1)
    simu->dt=0.5*dt;
     if(simu->pre_dtfields != NULL) {
      simu->pre_dtfields(simu);
    }
     DtFields(simu);
    RK_out(simu,l1);

    simu->tnow += 0.5 * dt;

    // l_3 = w_n + 0.5dt * S(l_2)
    RK_Copy(simu,simu->w,l2);
    DtFields(simu);
    RK_out(simu,l1);

    // l_4 = w_n + dt * S(l_3)
    simu->dt=dt;
    RK_Copy(simu,simu->w,l3);
    DtFields(simu);
    RK_out(simu,l1);

    simu->tnow += 0.5 * dt;

    // Compute yn+1= yn + dt(l1,l2,l3,l4)
    RK_Copy(simu,simu->w,l4);
    DtFields(simu);
    RK4_final_inplace(simu, l1, l2, l3,l4);


     if(simu->update_after_rk != NULL){
      simu->update_after_rk(simu, simu->w);
    }

    iter++;
     simu->iter_time_rk=iter;
  }
 
  free(l4);
  free(l3);
  free(l2);
  free(l1);
}





// An out-of-place RK step
void RK_Copy(Simulation * simu,schnaps_real * w, schnaps_real * w_temp)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < simu->wsize; iw++) {
    w_temp[iw] = w[iw];
  }
}

void RK_out(Simulation *simu, schnaps_real * wn)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < simu->wsize; iw++) {
    simu->w[iw] = wn[iw]+ simu->dt * simu->dtw[iw];
  }
}

void RK_in(Simulation *simu)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < simu->wsize; iw++) {
    simu->w[iw] = simu->w[iw]+ simu->dt * simu->dtw[iw];
  }
}

void RK_out_old(schnaps_real *dest, schnaps_real *fwn, schnaps_real *fdtwn, const schnaps_real dt,
	    const int sizew)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < sizew; iw++) {
    dest[iw] = fwn[iw] + dt * fdtwn[iw];
  }
}


// An in-place RK step
void RK_in_old(schnaps_real *fwnp1, schnaps_real *fdtwn, const schnaps_real dt, const int sizew)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < sizew; iw++) {
    fwnp1[iw] += dt * fdtwn[iw];
  }
}

void RK4_final_inplace(Simulation * simu, schnaps_real *l1, schnaps_real *l2, schnaps_real *l3,schnaps_real *l4)
{
  const schnaps_real a[] = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, simu->dt / 6.0};
  const schnaps_real b=-1.0 / 3.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < simu->wsize; ++i) {
    simu->w[i] =
      b * l1[i] +
      a[0] * l2[i] +
      a[1] * l3[i] +
      a[2] * l4[i] +
      a[3] * simu->dtw[i];
  }
}

void RK4_final_inplace_old(schnaps_real *w, schnaps_real *l1, schnaps_real *l2, schnaps_real *l3,
		       schnaps_real *dtw, const schnaps_real dt, const int sizew)
{
  const schnaps_real b = -1.0 / 3.0;
  const schnaps_real a[] = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, dt / 6.0};
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < sizew; ++i) {
    w[i] =
      b * w[i] +
      a[0] * l1[i] +
      a[1] * l2[i] +
      a[2] * l3[i] +
      a[3] * dtw[i];
  }
}


schnaps_real Get_Dt_RK(Simulation *simu)
{
  
  printf("cfl=%f hmin=%f \n",simu->cfl,simu->hmin);
  //assert(1==2);
  return simu->cfl * simu->hmin / simu->vcfl;

}


void DisplayArray(schnaps_real* array,
                  size_t size,
                  const char* name) {
  for (int i = 0; i < size; ++i)
    printf("%s[%d]: %f\n", name, i , array[i]);
}


void Compute_derivative(Simulation *simu, schnaps_real * wd, int nbfield){
  int nb_dof=0;
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {

    field *f = simu->fd + ie;
    nb_dof=simu->wsize/f->model.m;

    const int nraf[3] = {f->raf[0],
		       f->raf[1],
		       f->raf[2]};
    const int deg[3] = {f->deg[0],
		      f->deg[1],
		      f->deg[2]};
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const int m = f->model.m;
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Loop on the subcells
  for(int ic0 = 0; ic0 < nraf[0]; ic0++) {
    for(int ic1 = 0; ic1 < nraf[1]; ic1++) {
      for(int ic2 = 0; ic2 < nraf[2]; ic2++) {

	int ic[3] = {ic0, ic1, ic2};

	// index cell
	int nc = ic[0] + nraf[0] * (ic[1] + nraf[1] * ic[2]);
	// first glop index in the subcell
	int first_gp_cell = npg[0] * npg[1] * npg[2] * nc;

	for(int ipg=0;ipg<sc_npg;ipg++){
	  int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
	  schnaps_real derivate[3] ={0,0,0};
    schnaps_real wpg;
    schnaps_real xref[3];
    ref_pg_vol(f->deg,f->raf,index_glob_igp,xref,&wpg,NULL);
	  for(int jpg=0;jpg<sc_npg;jpg++){
	    int index_glob_jgp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + jpg, nbfield);
	    schnaps_real w=f->wn[index_glob_jgp];
	    schnaps_real dtau[3][3],codtau[3][3];
	    schnaps_real dphiref_j[3];
	    schnaps_real dphi_j[3];
	    grad_psi_pg(f->deg,f->raf,jpg,ipg,dphiref_j);
	    //
	    //
	    schnaps_ref2phy(f->physnode,
			    xref,dphiref_j,0,NULL,
			    dtau,codtau,dphi_j,NULL);
	    schnaps_real det = dot_product(dtau[0], codtau[0]);
	    
	    derivate[0]+=w * (dphi_j[0]/det);
	    derivate[1]+=w * (dphi_j[1]/det);
	    derivate[2]+=w * (dphi_j[2]/det);
	    //printf(" index glob:%i dphiref_j[0], jpg %i  %f \n",index_glob_igp,jpg,dphiref_j[0]);
	  }

	  wd[index_glob_igp]=derivate[0];
	  wd[index_glob_igp+nb_dof]=derivate[1];
	  wd[index_glob_igp+2*nb_dof]=derivate[2];
	}

      } // icl2
    } //icl1
  } // icl0

    
  }
  
}
/////// alternatives to PlotFields with a restriction to neighbouring subcells when projecting on the fe mesh
/////// this saves quite a lot of time
////// Still no optimal as the projection stencil should be point-based  but far cheaper than going through
///// all macromesh gauss points.
void PlotFieldsBinSparse(int typplot, int compare, Simulation* simu, char *fieldname,char *filename){

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  int one = 1;
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
            schnaps_real psi;
            psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
            testpsi += psi;
            int vi = f->varindex(f->deg, f->raf, f->model.m, index_glob_igp, typplot);
            value[nodecount] += psi * f->wn[vi];
          };
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        // Compare with an exact solution
        if (compare) {
          schnaps_real wex[f->model.m];
          f->model.ImposedData(Xphy, f->tnow, wex);
          value[nodecount] -= wex[typplot];
        }
        nodecount++;
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);

  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
  free(value);
} 
///
void PlotFieldsAsciiSparse(int typplot, int compare, Simulation* simu, char *fieldname,
	       char *filename) {

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
    simu->fd[0].raf[1],
    simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
    simu->fd[0].deg[1],
    simu->fd[0].deg[2]};
  const int npg[3] = {deg[0] + 1,
    deg[1] + 1,
    deg[2] + 1};
  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];
  //
  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(schnaps_real));
  //
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  //int npgv = NPG(deg, nraf);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        }
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

        schnaps_real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};
        schnaps_real testpsi = 0;
        ////////////////////////////////////////
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projection to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
              schnaps_real psi;
              psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
              testpsi += psi;
              int vi = f->varindex(f->deg, f->raf, f->model.m, index_glob_igp, typplot);
              value[nodecount] += psi * f->wn[vi];
          }; // end loop subcell gauss points
        }; //end loop neighbour subcell 2
        };//end loop neighbour subcell 1
        }; //end loop neighbour subcell 0
        assert(fabs(testpsi-1) < _SMALL);
      ///////////////////////////////////////
      // Compare with an exact solution
      if (compare) {
        schnaps_real wex[f->model.m];
        f->model.ImposedData(Xphy, f->tnow, wex);
        value[nodecount] -= wex[typplot];
      }
      nodecount++;
      fprintf(gmshfile, "%d %f %f %f\n", nodecount,
        Xplot[0], Xplot[1], Xplot[2]);
      } // end loop Hex64 fe nodes
    } // end loop subcell 2
    } // end loop subcell 1
    } // end loop subcell 0
  } // end loop macrocells

  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  fprintf(gmshfile, "%d\n",
  simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);
  int elm_type = 92;
  int num_tags = 0;
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
  	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    // Get the subcell id
    int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
    // Global subcell id
    int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
    fprintf(gmshfile, "%d ", numelem);
    fprintf(gmshfile, "%d ", elm_type);
    fprintf(gmshfile, "%d ", num_tags);
    for(int ii = 0; ii < 64; ii++) {
      int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
      fprintf(gmshfile, "%d ", numnoe);
    }
    fprintf(gmshfile, "\n");
    }
    }
    }
  } 
  fprintf(gmshfile, "$EndElements\n");
  // Now display data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else
    fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);
  //
  for(int ino = 0; ino < nb_plotnodes; ino++) {
    fprintf(gmshfile, "%d %f\n", ino + 1, value[ino]);
  }
  fprintf(gmshfile, "\n$EndNodeData\n");
  //
  fclose(gmshfile);
  free(value);
}

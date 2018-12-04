#include "field.h"
#include "geometry.h"
#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "global.h"
#include <math.h>
#include <float.h>
#include <string.h>
#include "quantities_vp.h"
#include "solverpoisson.h"
#include "model.h"

#ifdef _WITH_OPENCL
#include "clutils.h"
#include "clinfo.h"
#endif


// param[0] = M
// param[1] = deg x
// param[2] = deg y
// param[3] = deg z
// param[4] = raf x
// param[5] = raf y
// param[6] = raf z

#pragma start_opencl
int GenericVarindex(__constant int *deg, __constant int *raf, int m,
		    int ipg, int iv) {

  return iv + m * ipg;
  /* int npg = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1) */
  /*  * raf[0] * raf[1] * raf[2]; */
  /* return ipg + npg * iv; */
}
#pragma end_opencl

void GenericVarindex_to_ipgiv(int *deg, int *raf, int m, int imem,int *ipg,int *iv){

  //*iv = imem%m;
  //*ipg = imem/m;
  int npg = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1)
    * raf[0] * raf[1] * raf[2];
  *iv = imem/npg;
  *ipg= imem%npg;
  //
}

int GenericVarindex_CG(__constant int *deg, __constant int *raf, int m,
		    int ipg, int iv) {
  return iv + m * ipg;

  /* int npg = NPG_CG(deg, raf); */
  /* return ipg + npg * iv; */
}



#pragma start_opencl
int GenericVarindex3d(int m, int *nx, int *nc,
		      int elem,
		      int iv, int *ix, int *ic)
{
  // NB: passing nx and nx separately may allow one to deal with faces better.
  // Also, keeping the values in registers is theoretically faster than
  // even __local memory.

  // int nx[3] = {param[1] + 1, param[2] + 1, param[3] + 1};
  // int nc[3] = {param[4], param[5], param[6]};

  // number of glops in subcell
  int npgc = nx[0] * nx[1] * nx[2];
  // number of glops in macrocell:
  int npg = nc[0] * nc[1] * nc[2] * npgc;

  // index in subcell:
  int ipgc = ix[0] + nx[0] * (ix[1] + nx[1] * ix[2]);
  // index of subcell in macrocell:
  int nsubcell = ic[0] + nc[0] * (ic[1] + nc[1] * ic[2]);

  // index of point in macrocell:
  int ipg = ipgc + npgc * nsubcell;
  return iv + m * (ipg + npg * elem);
}
#pragma end_opencl

// Given a the index ipg of a poing in a subcell, determine the three
// logical coordinates of that point in the subcell.
/* void ipg_to_xyz(int ipg, int *p, int *npg) */
/* { */
/*   p[0] = ipg % npg[0]; */
/*   p[1] = (ipg / npg[0]) % npg[1]; */
/*   p[2] = ipg / npg[0] / npg[1]; */
/* }  */

schnaps_real min_grid_spacing(field *f)
{
  schnaps_real hmin = FLT_MAX;

    schnaps_real vol = 0, surf = 0;

    // Loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real xpgref[3], wpg;
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_real codtau[3][3], dtau[3][3];
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      NULL, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      schnaps_real det = dot_product(dtau[0], codtau[0]);
      vol += wpg * det;
    }
    for(int ifa = 0; ifa < 6; ifa++) {
      // loop on the faces
      for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, ifa); ipgf++) {
	schnaps_real xpgref[3], wpg;
	// get the coordinates of the Gauss point
	ref_pg_face(f->deg, f->raf, ifa, ipgf, xpgref, &wpg, NULL);
	schnaps_real vnds[3];
	{
	  schnaps_real codtau[3][3], dtau[3][3];
	  schnaps_ref2phy(f->physnode,
		  xpgref,
		  NULL, ifa, // dpsiref, ifa
		  NULL, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
	}
	surf += norm(vnds) * wpg;
      }
    }
    hmin = hmin < vol/surf ? hmin : vol/surf;
  // Now take into account the polynomial degree and the refinement
  int maxd = f->deg[0];
  maxd = maxd > f->deg[1] ? maxd : f->deg[1];
  maxd = maxd > f->deg[2] ? maxd : f->deg[2];
  hmin /= ((maxd + 1) * f->raf[0]);
  return hmin;
}

void init_empty_field(field *f)
{
/* #ifdef _WITH_OPENCL */
/*   f->use_source_cl = false; */
/* #endif */

  f->period[0] = -1;
  f->period[1] = -1;
  f->period[2] = -1;

  f->jacobian = NULL;
  f->det = NULL;
  f->store_det = false;
  f->starpu_registered = false;
  f->wn = NULL;
  f->dtwn = NULL;
  f->res = NULL;
  f->solver = NULL;
  f->solver_spu = NULL;
  f->rmat = NULL;
  f->rmat_spu = NULL;

}

void init_data(field *f)
{
  int deg[3];
  deg[0] = f->deg[0];
  deg[1] = f->deg[1];
  deg[2] = f->deg[2];
  int raf[3];
  raf[0] = f->raf[0];
  raf[1] = f->raf[1];
  raf[2] = f->raf[2];
  int n_npg = NPG(deg, raf);

    for(int ipg = 0; ipg < n_npg; ipg++) {
      schnaps_real xpg[3];
      schnaps_real xref[3], omega;
      ref_pg_vol(deg, raf, ipg, xref, &omega, NULL);
      schnaps_real dtau[3][3];
      schnaps_ref2phy(f->physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds
      { // Check the reverse transform at all the GLOPS
 	/* schnaps_real xref2[3]; */
	/* schnaps_phy2ref(f->physnode, xpg, xref2); */
	/* if (Dist(xref, xref2) >= _SMALL) { */
	/*   printf("xref %f %f %f xref2 %f %f %f \n", */
	/* 	 xref[0],xref[1],xref[2], */
	/* 	 xref2[0],xref2[1],xref[2]); */
	/* } */
	//assert(Dist(xref, xref2) < _SMALL);
      }
      schnaps_real w[f->model.m];
      f->model.InitData(xpg, w);
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(deg, raf, f->model.m, ipg, iv);
	f->wn[imem] = w[iv];
	f->dtwn[imem] = 0;
      }
    }
}

#ifdef _WITH_OPENCL
void set_physnodes_cl(Simulation *simu)
{
  const int nmacro = simu->macromesh.nbelems;
  schnaps_real buf_size = sizeof(schnaps_real) * 60 * nmacro;
  schnaps_real *physnode = malloc(buf_size);

  for(int ie = 0; ie < nmacro; ++ie) {
    int ie20 = 20 * ie;
    for(int inoloc = 0; inoloc < 20; ++inoloc) {
      int ino = 3 * simu->macromesh.elem2node[ie20 + inoloc];
      schnaps_real *iphysnode = physnode + 3 * ie20 + 3 * inoloc;
      schnaps_real *nodeino = simu->macromesh.node + ino;
      iphysnode[0] = nodeino[0];
      iphysnode[1] = nodeino[1];
      iphysnode[2] = nodeino[2];
    }
  }

  cl_int status;
  status = clEnqueueWriteBuffer(simu->cli.commandqueue,
  				simu->physnodes_cl, // cl_mem buffer,
  				CL_TRUE,// cl_bool blocking_write,
  				0, // size_t offset
  				buf_size, // size_t cb
  				physnode, //  	void *ptr,
  				0, 0, 0);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  free(physnode);


}
#endif

void Initfield(field *f, Model model,
	       schnaps_real physnode[][3],
	       int *deg, int *raf, schnaps_real *w,
	       schnaps_real* dtw) {
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};

  //f->vmax = 1.0; // FIXME: make this variable ??????

  f->model = model;

  f->interp.m = model.m;

  f->interp.deg[0] = deg[0];
  f->interp.deg[1] = deg[1];
  f->interp.deg[2] = deg[2];

  f->interp.raf[0] = raf[0];
  f->interp.raf[1] = raf[1];
  f->interp.raf[2] = raf[2];

  f->deg[0] = deg[0];
  f->deg[1] = deg[1];
  f->deg[2] = deg[2];

  f->raf[0] = raf[0];
  f->raf[1] = raf[1];
  f->raf[2] = raf[2];

  f->npg[0] = f->deg[0] + 1;
  f->npg[1] = f->deg[1] + 1;
  f->npg[2] = f->deg[2] + 1;

  for(int inoloc = 0 ; inoloc < 20; ++inoloc){
    f->physnode[inoloc][0] = physnode[inoloc][0];
    f->physnode[inoloc][1] = physnode[inoloc][1];
    f->physnode[inoloc][2] = physnode[inoloc][2];
  }

  f->varindex = GenericVarindex;
  f->varindex_to_ipgiv= GenericVarindex_to_ipgiv;
  int nmem = f->model.m * NPG(f->deg, f->raf);
  f->wsize = nmem;
  //real g_memsize = nmem * sizeof(real) * 1e-9;
  /* if(sizeof(real) == sizeof(real)) */
  /*   printf("Allocating %d doubles per array (%f GB).\n", nmem, g_memsize); */
  /* else */
  /*   printf("Allocating %d floats per array (%f GB)\n", nmem, g_memsize); */

  if (w == NULL){
    f->wn = calloc(nmem, sizeof(schnaps_real));
  } else {
    f->wn = w;
  }
  assert(f->wn);

  if (dtw == NULL){
    f->dtwn = calloc(nmem, sizeof(schnaps_real));
  } else {
    f->dtwn = dtw;
  }
  assert(f->dtwn);

  f->res = f->dtwn;
  f->pre_dtfield = NULL;
  f->post_dtfield = NULL;
  f->update_after_rk = NULL;
  //f->model.Source = NULL;
  //f->pic = NULL;

  // TODO: move this to the integrator code
  f->tnow=0;

  if (f->store_det) {
    assert(1==2);
    int npg = NPG(f->deg, f->raf);
    f->jacobian = malloc(npg * 9 * sizeof(schnaps_real));
    f->det = malloc(npg * sizeof(schnaps_real));
    for(int ipg = 0; ipg < npg; ipg++) {
      schnaps_real xpg[3];
      schnaps_real xref[3], omega;
      ref_pg_vol(f->deg, f->raf, ipg, xref, &omega, NULL);
      schnaps_real codtau[9];
      schnaps_real* dtau = f->jacobian + 9 * ipg;
      schnaps_ref2phy(f->physnode,
	      xref,
	      NULL, -1, // dphiref, ifa
              NULL, (schnaps_real(*)[3])dtau, // xphy , dtau
	      (schnaps_real(*)[3])codtau, NULL, NULL); // codtau, dphi, vnds
      f->det[ipg] = dot_product(dtau, codtau);
    }
  }

  
  init_data(f);   // to be exchanged with spu init xxx

  // Compute cfl parameter min_i vol_i/surf_i
  f->hmin = min_grid_spacing(f);

  //printf("hmin=%f\n", f->hmin);

  //printf("field init done\n");
  // default values for the linear solvers: not used
  f->solver = NULL;
  f->rmat = NULL;
  f->solver_spu = NULL;
  f->rmat_spu = NULL;


}


void Registerfield_SPU(field *f){

  if (!starpu_is_init && starpu_use){
    int ret;
    ret = starpu_init(NULL);
    assert(ret != -ENODEV) ;
    starpu_is_init = true;
  }

  if (starpu_use && !f->starpu_registered){

    
    starpu_vector_data_register(&(f->wn_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(f->wn), // vector location
				f->wsize,  // size
				sizeof(schnaps_real));  // type

    
    starpu_vector_data_register(&(f->dtwn_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(f->dtwn), // vector location
				f->wsize,  // size
				sizeof(schnaps_real));  // type

    starpu_vector_data_register(&(f->res_handle), // mem handle
				0, // location: CPU
				(uintptr_t)(f->res), // vector location
				f->wsize,  // size
				sizeof(schnaps_real));  // type
    f->starpu_registered = true;
  }
    
}

void Unregisterfield_SPU(field *f){


  if (starpu_use && f->starpu_registered){

    starpu_data_unregister(f->wn_handle);
    starpu_data_unregister(f->dtwn_handle);
    starpu_data_unregister(f->res_handle);
    f->starpu_registered = false;
  }


}


void Initfield_SPU(field *f){

  if (!starpu_is_init && starpu_use){
    int ret;
    ret = starpu_init(NULL);
    assert(ret != -ENODEV) ;
    starpu_is_init = true;
  }

  /* if (starpu_use){ */

  /*   starpu_vector_data_register(&(f->wn_handle), // mem handle */
  /* 				0, // location: CPU */
  /* 				(uintptr_t)(f->wn), // vector location */
  /* 				f->wsize,  // size */
  /* 				sizeof(real));  // type */

  /*   for(int ii=0; ii < f->wsize; ii++) f->res[ii] = 0; */


  /*   starpu_vector_data_register(&(f->res_handle), // mem handle */
  /* 				0, // location: CPU */
  /* 				(uintptr_t)(f->res), // vector location */
  /* 				f->wsize,  // size */
  /* 				sizeof(real));  // type */
  /*   assert(f->res != NULL); */
  /*   assert(f->res == f->dtwn); */
  /* } */


}

// This is the destructor for a field
void free_field(field *f)
{

#ifdef _WITH_OPENCL
  cl_int status;

  //status = clReleaseMemObject(f->physnode_cl);
  //if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  //assert(status >= CL_SUCCESS);
  //free(f->physnode);
#endif
}

// Display the field on screen
void Displayfield(field *f) {
  printf("Display field...\n");
    printf("elem data \n");
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real xref[3], wpg;
      ref_pg_vol(f->deg, f->raf, ipg, xref, &wpg, NULL);

      schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3];
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      schnaps_ref2phy(f->physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      printf("Gauss point %d xref=%f %f %f \n", ipg, xref[0], xref[1], xref[2]);
      printf("Gauss point %d xphy=%f %f %f \n", ipg, xphy[0], xphy[1], xphy[2]);
      printf("dtw= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	printf("%f ", f->dtwn[imem]);
      }
      printf("\n");
      printf("w= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	printf("%f ", f->wn[imem]);
      }
      printf("\n");

  }
};

// Save the results in a text file
// in order plot it with Gnuplot
void Gnuplot(Simulation *simu,int dir, schnaps_real fixval, char* filename) {

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  printf("Save for Gnuplot...\n");

  for(int ie=0; ie < simu->macromesh.nbelems; ie++){
    field *f = &simu->fd[ie];
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {

      schnaps_real xref[3], xphy[3], wpg;
      schnaps_real dtau[3][3];
      ref_pg_vol(f->deg, f->raf, ipg, xref, &wpg, NULL);

      schnaps_ref2phy(f->physnode,
	      xref,
	      0, -1, // dphiref, ifa
	      xphy, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      if(xphy[dir] > -(fixval + 0.0001) && xphy[dir] < (fixval + 0.00001)){

	fprintf(gmshfile, "%f ",xphy[1-dir]);

	for(int iv = 0; iv < f->model.m; iv++) {
	  int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	  fprintf(gmshfile, "%f ",f->wn[imem]);
	}
	fprintf(gmshfile, "\n");

      }
    }
}

  fclose(gmshfile);
};



// Compute inter-subcell fluxes
void DGSubCellInterface(field *f, schnaps_real *w, schnaps_real *dtw)
{


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

  // Loop on the subcells
  for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
    for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
      for(int icL2 = 0; icL2 < nraf[2]; icL2++) {

	int icL[3] = {icL0, icL1, icL2};

	// Get the left subcell id
	int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	// First glop index in the subcell
	// TODO change this wrong approach !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
		  dtw[imemL] -= flux[iv] * wpg;
		  dtw[imemR] += flux[iv] * wpg;
		}

	      }  // face yhat loop
	    } // face xhat loop
	  } // endif internal face
	} // dim loop
      } // subcell icl2 loop
    } // subcell icl1 loop
  } // subcell icl0 loop

}


// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
void DGMacroCellInterface(int locfaL,
			  field *fL, int offsetL, field *fR, int offsetR,
			  schnaps_real *w, schnaps_real *dtw)
{
  const unsigned int m = fL->model.m;
  // Assembly of the surface terms loop on the macrocells faces
  //int ieL = msh->face2elem[4 * ifa + 0];
  //int locfaL = msh->face2elem[4 * ifa + 1];


  //int ieR = msh->face2elem[4 * ifa + 2];
  //int locfaR = msh->face2elem[4 * ifa + 3];

  schnaps_real *fwL = w + offsetL;
  schnaps_real *fwR = w + offsetR;

  schnaps_real *fdtwL = dtw + offsetL;
  schnaps_real *fdtwR = dtw + offsetR;

  // Loop over the points on a single macro cell interface.
/* #ifdef _OPENMP */
/* #pragma omp parallel for */
/* #endif */
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
      //printf("ipgL=%d ipgR=%d\n",ipgL,ipgR);

      //Uncomment to check that the neighbour-finding algorithm worked.
      /* { */
      /*   real xpgR[3], xrefR[3], wpgR; */
      /*   ref_pg_vol(iparam + 1, ipgR, xrefR, &wpgR, NULL); */
      /*   Ref2Phy(physnodeR, */
      /* 	  xrefR, */
      /* 	  NULL, -1, dphiref, ifa */
      /* 	  xpgR, NULL, */
      /* 	  NULL, NULL, NULL); codtau, dphi, vnds */
      /*  #ifdef _PERIOD */
      /*   assert(fabs(Dist(xpg,xpgR)-_PERIOD)<_SMALL); */
      /*   #else */
      /* assert(Dist(xpg,xpgR)<1e-11); */
      /* #endif */
      /* } */

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
	fdtwL[imemL] -= flux[iv] * wpg;
	fdtwR[imemR] += flux[iv] * wpg;
      }

    } else { // The point is on the boundary.
      for(int iv = 0; iv < m; iv++) {
	int imemL = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv);
	wL[iv] = fwL[imemL];
      }

      fL->model.BoundaryFlux(xpg, fL->tnow, wL, vnds, flux);

      for(int iv = 0; iv < m; iv++) {
	// The basis functions is also the gauss point index
	int imemL = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv);
	fdtwL[imemL] -= flux[iv] * wpg;
      }
    }

  }
}

// Apply division by the mass matrix
void DGMass(field *f, schnaps_real *w, schnaps_real *dtw)
{

  int m = f->model.m;

  for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
    schnaps_real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
    schnaps_ref2phy(f->physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    schnaps_real det = dot_product(dtau[0], codtau[0]);
    for(int iv = 0; iv < m; iv++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      dtw[imem] /= (wpg * det);
    }
  }

}

// Apply the source term
void DGSource(field *f, schnaps_real *w, schnaps_real *dtw)
{
  if (f->model.Source == NULL) {
    return;
  }

  const int m = f->model.m;

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
      dtw[imem] += source[iv] * det * wpg;

    }
  }

}

// Compute the Discontinuous Galerkin volume terms, fast version
void DGVolume(field *f, schnaps_real *w, schnaps_real *dtw)
{


  const int m = f->model.m;
  const int deg[3] = {f->deg[0],
		      f->deg[1],
		      f->deg[2]};
  const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
  const int nraf[3] = {f->raf[0],
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

	// TODO: wrong approach not compatible with xyz_to_ipg !!!!!!!!!!!!!!!!!!!!!!!!!
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
		// TODO rather call xyz_to_ipg
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
		    dtw[imems[temp]] += flux[iv] * wpgL;
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

//cmake . -DUSE_OPENCL:BOOL=OFF

// Compute the average in a subcell : mimic 0 order
void DGAverage(field *f, schnaps_real *w)
{


  const int m = f->model.m;
  const int deg[3] = {f->deg[0], f->deg[1], f->deg[2]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {f->raf[0], f->raf[1], f->raf[2]};

  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];


  /* Loop on the subcells */
  for (int icL0 = 0; icL0 < nraf[0]; icL0++) {
    for (int icL1 = 0; icL1 < nraf[1]; icL1++) {
      for (int icL2 = 0; icL2 < nraf[2]; icL2++) {

        int icL[3] = {icL0, icL1, icL2};
        /* Get the L subcell id */
        int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);

        /* Get the first glop index in the subcell */
        int offsetL = npg[0] * npg[1] * npg[2] * ncL;

        /* Compute all xref associated to the subcell */
        schnaps_real *xref0 = malloc(sc_npg * sizeof(schnaps_real));
        schnaps_real *xref1 = malloc(sc_npg * sizeof(schnaps_real));
        schnaps_real *xref2 = malloc(sc_npg * sizeof(schnaps_real));
        schnaps_real *omega = malloc(sc_npg * sizeof(schnaps_real));
        int *imems = malloc(m * sc_npg * sizeof(int));
        int pos = 0;
        
        
        /* Get imem */
        for(unsigned int p = 0; p < sc_npg; p++) {
          schnaps_real xref[3];
          schnaps_real tomega;

          ref_pg_vol(f->deg, f->raf, offsetL + p, xref, &tomega, NULL);
          xref0[p] = xref[0];
          xref1[p] = xref[1];
          xref2[p] = xref[2];
          omega[p] = tomega;

          for(int im = 0; im < m; im++) {
            
            imems[pos++] = f->varindex(f->deg,
                                       f->raf,
                                       f->model.m,
                                       offsetL + p,
                                       im);
          }
        }

        schnaps_real vol_subcell = 0;
        schnaps_real wmean[m];
        
        for (int iv = 0; iv < m; ++iv) wmean[iv] = 0;
        bool test = false;
        /* Loop over the subcell gauss points */
        for (int p0 = 0; p0 < npg[0]; p0++) {
          for (int p1 = 0; p1 < npg[1]; p1++) {
            for (int p2 = 0; p2 < npg[2]; p2++) {
              schnaps_real wL[m];
              int p[3] = {p0, p1, p2};
              int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
              for(int iv = 0; iv < m; iv++) {
                wL[iv] = w[imems[m * (ipgL - offsetL) + iv]];
              }
              
              schnaps_real xrefL[3] = {xref0[ipgL - offsetL],
                                       xref1[ipgL - offsetL],
                                       xref2[ipgL - offsetL]};
                                       
              schnaps_real wpgL = omega[ipgL - offsetL];

              schnaps_real dtau[3][3], codtau[3][3], dphiL[3];
              schnaps_ref2phy(f->physnode,
                              xrefL,
                              NULL,
                              -1, /* ifa */
                              NULL, /* xphy */
                              dtau,
                              codtau,
                              NULL,
                              NULL); /* vnds */

              schnaps_real det = dot_product(dtau[0], codtau[0]);
              
              

              if (wL[0] <= 10) {
                test = true;
              }
              
              if (wL[1] / wL[0] >= 1) { 
                test = true;
              }
              if (wL[2] / wL[0] >= 1) {
                test = true;
              }
              if ( wL[3] / wL[0] >= 1) {
                test = true;             
              }
              for (int iv = 0; iv < m; ++iv) {
                wmean[iv] += wL[iv] * wpgL * det;
               
              }
              vol_subcell += wpgL * det;
            }  /* p2 */
          }  /* p1 */
        }  /* p0 */
        // printf("vol cell = %f\n",vol_subcell);
        // printf("w_mean = %f\n",wmean[0]);
        /* Loop over the subcell gauss points */
        if ( test ){
          // printf("AVERAGING\n");
          for (int p0 = 0; p0 < npg[0]; p0++) {
            for (int p1 = 0; p1 < npg[1]; p1++) {
              for (int p2 = 0; p2 < npg[2]; p2++) {
                int p[3] = {p0, p1, p2};
                int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
                
                for (int iv = 0; iv < m; ++iv)
                  w[imems[m * (ipgL - offsetL) + iv]] = wmean[iv] / vol_subcell;
              } /* p2 */
            } /* p1 */
          } /* p0 */
        }
        free(omega);
        free(xref0);
        free(xref1);
        free(xref2);
        free(imems);

      } /* icL2 */
    } /* icL1 */
  } /* icL0 */
}  





// Compute the normalized L2 distance with the imposed data

// Compute the normalized L2 distance with the imposed data
schnaps_real L2error_onefield(Simulation *simu, int nbfield) {
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  schnaps_real error = 0;
  schnaps_real mean = 0;

  for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    // Get the physical nodes of element ie

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
	f->model.ImposedData(xphy, f->tnow, wex);
      }

      int iv = nbfield;
      schnaps_real diff = w[iv] - wex[iv];
      error += diff * diff * wpg * det;
      mean += w[iv] * w[iv] * wpg * det;
    }
  }
  return sqrt(error);
}


void InterpField(field *f, schnaps_real *xref, schnaps_real *w){
  int deg[3] = {f->deg[0],
		      f->deg[1],
		      f->deg[2]};

  int nraf[3] = {f->raf[0],
		       f->raf[1],
		       f->raf[2]};

  for(int iv = 0; iv < f->model.m; iv++)
    w[iv] = 0;

  int is[3];

  for(int ii = 0; ii < 3; ii++){
    is[ii] = xref[ii] * nraf[ii];
    assert(is[ii] < nraf[ii] && is[ii]>= 0);
  }

  int npgv = NPG(deg, nraf);
  // TODO: loop only on non zero basis function
  for(int ib = 0; ib < npgv; ib++) {
    schnaps_real psi;
    psi_ref_subcell(deg, nraf, is, ib, xref, &psi, NULL);

    for(int iv=0;iv<f->model.m;iv++){
      int imem = f->varindex(deg, nraf, f->model.m, ib, iv);
      w[iv] += psi * f->wn[imem];
    }
  }
}


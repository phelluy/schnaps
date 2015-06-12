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

#ifdef _WITH_OPENCL 
#include "clutils.h"
#include "clinfo.h"
#endif

#ifdef _WITH_PTHREAD
#include <pthread.h>
#endif

// param[0] = M
// param[1] = deg x
// param[2] = deg y
// param[3] = deg z
// param[4] = raf x
// param[5] = raf y
// param[6] = raf z

#pragma start_opencl
int GenericVarindex(__constant int *param, int elem, int ipg, int iv) {
  int npg = (param[1] + 1) * (param[2] + 1) * (param[3] + 1)
    * param[4] * param[5] * param[6];
  return iv + param[0] * (ipg + npg * elem);
}
#pragma end_opencl

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

real min_grid_spacing(field *f)
{
  real hmin = FLT_MAX;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    real vol = 0, surf = 0;
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // Loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real xpgref[3], wpg;
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
      real codtau[3][3], dtau[3][3];
      Ref2Phy(physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      NULL, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real det = dot_product(dtau[0], codtau[0]);
      vol += wpg * det;
    }
    for(int ifa = 0; ifa < 6; ifa++) {
      // loop on the faces
      for(int ipgf = 0; ipgf < NPGF(f->interp_param + 1, ifa); ipgf++) {
	real xpgref[3], wpg;
	// get the coordinates of the Gauss point
	ref_pg_face(f->interp_param + 1, ifa, ipgf, xpgref, &wpg, NULL);
	real vnds[3];
	{
	  real codtau[3][3], dtau[3][3];
	  Ref2Phy(physnode,
		  xpgref,
		  NULL, ifa, // dpsiref, ifa
		  NULL, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
	}
	surf += norm(vnds) * wpg;
      }
    }
    hmin = hmin < vol/surf ? hmin : vol/surf;

  }
 
  // Now take into account the polynomial degree and the refinement
  int maxd = f->interp_param[1];
  maxd = maxd > f->interp_param[2] ? maxd : f->interp_param[2];
  maxd = maxd > f->interp_param[3] ? maxd : f->interp_param[3];

  hmin /= ((maxd + 1) * f->interp_param[4]);

  return hmin;
}

void init_empty_field(field *f)
{
#ifdef _WITH_OPENCL
  f->use_source_cl = false;
#endif
}

void init_data(field *f)
{
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real xpg[3];
      real xref[3], omega;
      ref_pg_vol(f->interp_param + 1, ipg, xref, &omega, NULL);
      real dtau[3][3];
      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds
      { // Check the reverse transform at all the GLOPS
 	real xref2[3];
	Phy2Ref(physnode, xpg, xref2);
	assert(Dist(xref, xref2) < 1e-8);
      }

      real w[f->model.m];
      f->model.InitData(xpg, w);
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem;
	imem = f->varindex(f->interp_param, ie, ipg, iv);
	f->wn[imem] = w[iv];

      }
    }
  }
}

#ifdef _WITH_OPENCL
void init_field_cl(field *f)
{
  InitCLInfo(&f->cli, nplatform_cl, ndevice_cl);
  cl_int status;

  f->wn_cl = clCreateBuffer(f->cli.context,
			    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
			    sizeof(real) * f->wsize,
			    f->wn,
			    &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dtwn_cl = clCreateBuffer(f->cli.context,
			      CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
			      sizeof(real) * f->wsize,
			      f->dtwn,
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->param_cl = clCreateBuffer(f->cli.context,
			       CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			       sizeof(int) * 7,
			       f->interp_param,
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->physnode = calloc(60, sizeof(real));

  f->physnode_cl = clCreateBuffer(f->cli.context,
				  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				  sizeof(real) * 60,
				  f->physnode,
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->physnodeR = calloc(60, sizeof(real));

  f->physnodeR_cl = clCreateBuffer(f->cli.context,
				   CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				   sizeof(real) * 60,
				   f->physnodeR,
				   &status);


  // Program compilation
  char *strprog;
  GetOpenCLCode();
  ReadFile("schnaps.cl", &strprog);

  printf("\t%s\n", numflux_cl_name);
  //printf("\t%s\n", strprog);

  // If the source term is set (via set_source_CL) then add it to the
  // buildoptions and compile using the new buildoptions.
  if(f->use_source_cl) {
    char *temp;
    int len0 = strlen(cl_buildoptions);
    char *D_SOURCE_FUNC = " -D_SOURCE_FUNC=";
    int len1 = strlen(D_SOURCE_FUNC);
    int len2 = strlen(f->sourcename_cl);
    temp = calloc(sizeof(char), len0 + len1 + len2 + 2);
    strcat(temp, cl_buildoptions);
    strcat(temp, D_SOURCE_FUNC);
    strcat(temp, f->sourcename_cl);
    strcat(temp, " ");
    BuildKernels(&f->cli, strprog, temp);
  } else {
    printf("No source term\n");
    BuildKernels(&f->cli, strprog, cl_buildoptions);
  }

  f->dgmass = clCreateKernel(f->cli.program,
			     "DGMass",
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgflux = clCreateKernel(f->cli.program,
			     "DGFlux",
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgvolume = clCreateKernel(f->cli.program,
			       "DGVolume",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgsource = clCreateKernel(f->cli.program,
			       "DGSource",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dginterface = clCreateKernel(f->cli.program,
				  "DGMacroCellInterface",
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->dgboundary = clCreateKernel(f->cli.program,
				 "DGBoundary",
				 &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->RK_out_CL = clCreateKernel(f->cli.program,
				"RK_out_CL",
				&status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->RK4_final_stage = clCreateKernel(f->cli.program,
				      "RK4_final_stage",
				      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->RK_in_CL = clCreateKernel(f->cli.program,
			       "RK_in_CL",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->zero_buf = clCreateKernel(f->cli.program,
			       "set_buffer_to_zero",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Initialize events. // FIXME: free on exit
  f->clv_zbuf = clCreateUserEvent(f->cli.context, &status);
  
  f->clv_mapdone = clCreateUserEvent(f->cli.context, &status);

  f->clv_physnodeupdate = clCreateUserEvent(f->cli.context, &status);

  f->clv_mci = clCreateUserEvent(f->cli.context, &status);
  f->clv_interkernel = clCreateUserEvent(f->cli.context, &status);
  f->clv_interupdate = clCreateUserEvent(f->cli.context, &status);
  f->clv_interupdateR = clCreateUserEvent(f->cli.context, &status);

  f->clv_mass = clCreateUserEvent(f->cli.context, &status);

  f->clv_flux = calloc(3, sizeof(cl_event));
  f->clv_flux[0] = clCreateUserEvent(f->cli.context, &status);
  f->clv_flux[1] = clCreateUserEvent(f->cli.context, &status);
  f->clv_flux[2] = clCreateUserEvent(f->cli.context, &status);

  f->clv_volume = clCreateUserEvent(f->cli.context, &status);
  f->clv_source = clCreateUserEvent(f->cli.context, &status);

  // Set timers to zero
  f->zbuf_time = 0;
  f->mass_time = 0;
  f->vol_time = 0;
  f->flux_time = 0;
  f->minter_time = 0;
  f->boundary_time = 0;
  f->source_time = 0;
  f->rk_time = 0;

  // Set roofline counts to zero
  f->flops_vol = 0;
  f->flops_flux = 0;
  f->flops_mass = 0; 
  f->reads_vol = 0;
  f->reads_flux = 0;
  f->reads_mass = 0; 
}
#endif

void Initfield(field *f) {
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  
  //f->vmax = 1.0; // FIXME: make this variable ??????

  // a copy for avoiding too much "->"
  for(int ip = 0; ip < 8; ip++)
    f->interp_param[ip] = f->interp.interp_param[ip];

  int nmem = f->model.m * f->macromesh.nbelems * NPG(f->interp_param + 1);
  f->wsize = nmem;
  printf("allocate %d reals\n", nmem);
  f->wn = calloc(nmem, sizeof(real));
  assert(f->wn);
  f->dtwn = calloc(nmem, sizeof(real));
  assert(f->dtwn);
  f->Diagnostics = NULL;
  f->pre_dtfield = NULL;
  f->update_after_rk = NULL;
  f->model.Source = NULL;

  // TODO: move this to the integrator code
  f->tnow=0;
  f->itermax=0;
  f->iter_time=0;
  f->nb_diags = 0;

  f->tnow = 0;

  init_data(f);

  // Compute cfl parameter min_i vol_i/surf_i
  f->hmin = min_grid_spacing(f);

  printf("hmin=%f\n", f->hmin);

  // Allocate and set MacroFaces
  f->mface = calloc(f->macromesh.nbfaces, sizeof(MacroFace*));
  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++) {
    f->mface[ifa].first = ifa;
    f->mface[ifa].last_p1 = ifa + 1;
  }

  // Allocate and set MacroCells
  f->mcell = calloc(f->macromesh.nbelems, sizeof(MacroCell*));
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    f->mcell[ie].first = ie;
    f->mcell[ie].last_p1 = ie + 1;
  }


#ifdef _WITH_OPENCL
  // opencl inits
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable; OpenCL initialization disabled.\n");
  } else {

    init_field_cl(f);
  }
#endif // _WITH_OPENCL
  
  printf("field init done\n");
}

// This is the destructor for a field
void free_field(field *f) 
{
  free(f->mcell);
  free(f->mface);

#ifdef _WITH_OPENCL
  cl_int status;

  status = clReleaseMemObject(f->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  free(f->physnode);
#endif
}

// Display the field on screen
void Displayfield(field *f) {
  printf("Display field...\n");
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    printf("elem %d\n", ie);
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real xref[3], wpg;
      ref_pg_vol(f->interp_param + 1, ipg, xref, &wpg, NULL);

      printf("Gauss point %d %f %f %f \n", ipg, xref[0], xref[1], xref[2]);
      printf("dtw= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	printf("%f ", f->dtwn[imem]);
      }
      printf("\n");
      printf("w= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	printf("%f ", f->wn[imem]);
      }
      printf("\n");
    }
  }
};

// Save the results in a text file
// in order plot it with Gnuplot
void Gnuplot(field* f,int dir, real fixval, char* filename) {

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  printf("Save for Gnuplot...\n");
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {

    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {

      real xref[3], xphy[3], wpg;
      real dtau[3][3];
      ref_pg_vol(f->interp_param + 1, ipg, xref, &wpg, NULL);

      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
	      xphy, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      if(xphy[dir] > -(fixval + 0.0001) && xphy[dir] < (fixval + 0.00001)){

	fprintf(gmshfile, "%f ",xphy[1-dir]);

	for(int iv = 0; iv < f->model.m; iv++) {
	  int imem = f->varindex(f->interp_param, ie, ipg, iv);
	  fprintf(gmshfile, "%f ",f->wn[imem]);
	}
	fprintf(gmshfile, "\n");

      }
    }
  }
  fclose(gmshfile);
};


// Save the results in the gmsh format typplot: index of the plotted
// variable int compare == true -> compare with the exact value.  If
// fieldname is NULL, then the fieldname is typpplot.
void Plotfield(int typplot, int compare, field* f, char *fieldname,
	       char *filename) {

  real hexa64ref[3 * 64] = { 
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

  int *elem2nodes = f->macromesh.elem2node;
  real *node = f->macromesh.node;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {f->interp_param[4], f->interp_param[5], f->interp_param[6]};
  // Refinement size in each direction
  real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(real));

  int nb_plotnodes = f->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);

  real *value = malloc(nb_plotnodes * sizeof(real));
  assert(value);
  int nodecount = 0;

  // Nodes
  int npgv = NPG(f->interp_param + 1);
  for(int i = 0; i < f->macromesh.nbelems; i++) {
    // Get the nodes of element L
    int nnodes = 20;
    real physnode[nnodes][3];
    for(int ino = 0; ino < nnodes; ino++) {
      int numnoe = elem2nodes[nnodes * i + ino];
      for(int ii = 0; ii < 3; ii++) {
        physnode[ino][ii] = node[3 * numnoe + ii];
      }
    }

    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {

	  for(int ino = 0; ino < 64; ino++) {
	    real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
			     hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
			     hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
	    
	    for(int ii = 0; ii < 3; ii++) {
	      assert(Xr[ii] < 1 +  1e-10);
	      assert(Xr[ii] > -1e-10);
	    }

	    real Xphy[3];
	    Ref2Phy(physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

	    real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};

	    value[nodecount] = 0;
	    real testpsi = 0;
	    for(int ib = 0; ib < npgv; ib++) {
	      real psi;
	      psi_ref_subcell(f->interp_param + 1, icL, ib, Xr, &psi, NULL);
	      testpsi += psi;
	      int vi = f->varindex(f->interp_param, i, ib, typplot);
	      value[nodecount] += psi * f->wn[vi];
	    }
	    assert(fabs(testpsi-1) < 1e-10);

	    // Compare with an exact solution
	    if (compare) {
	      real wex[f->model.m];
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
	  f->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);

  int elm_type = 92;
  int num_tags = 0;

  // fwrite((char*) &elm_type, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_elm_follow, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_tags, sizeof(int), 1, gmshfile);

  for(int i = 0; i < f->macromesh.nbelems; i++) {
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

  real t = 0;
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
  /* 	  psi_ref(f->interp_param + 1, ib, Xr, &psi, NULL); */

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

// Compute inter-subcell fluxes
void DGSubCellInterface(void *mc, field *f, real *w, real *dtw) 
{
  MacroCell *mcell = (MacroCell*) mc;

  // Loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    const int nraf[3] = {f->interp_param[4],
			 f->interp_param[5],
			 f->interp_param[6]};
    const int deg[3] = {f->interp_param[1],
			f->interp_param[2],
			f->interp_param[3]};
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
		  real vnds[3];
		  {
		    real xref[3], wpg3;
		    ref_pg_vol(f->interp_param + 1, ipgL, xref, &wpg3, NULL);
		    // mapping from the ref glop to the physical glop
		    real dtau[3][3], codtau[3][3];
		    Ref2Phy(physnode,
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

		    real h1h2 = 1. / nraf[dim1] / nraf[dim2];
		    vnds[0] = codtau[0][dim0] * h1h2;
		    vnds[1] = codtau[1][dim0] * h1h2;
		    vnds[2] = codtau[2][dim0] * h1h2;
		  }

		  // numerical flux from the left and right state and
		  // normal vector
		  real wL[m], wR[m], flux[m];
		  for(int iv = 0; iv < m; iv++) {
		    // TO DO change the varindex signature
		    int imemL = f->varindex(f->interp_param, ie, ipgL, iv); 
		    int imemR = f->varindex(f->interp_param, ie, ipgR, iv);
		    // end TO DO
		    wL[iv] = w[imemL];
		    wR[iv] = w[imemR];
		  }
		  f->model.NumFlux(wL, wR, vnds, flux);

		  // subcell ref surface glop weight
		  real wpg
		    = wglop(deg[dim1], iL[dim1])
		    * wglop(deg[dim2], iL[dim2]);

		  /* printf("vnds %f %f %f flux %f wpg %f\n", */
		  /* 	 vnds[0], vnds[1], vnds[2], */
		  /* 	 flux[0], wpg); */

		  // finally distribute the flux on the two sides
		  for(int iv = 0; iv < m; iv++) {
		    // TO DO change the varindex signature
		    int imemL = f->varindex(f->interp_param, ie, ipgL, iv);
		    int imemR = f->varindex(f->interp_param, ie, ipgR, iv);
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

  } // macro elem loop
}

// compute the Discontinuous Galerkin inter-macrocells boundary terms
void DGMacroCellInterfaceSlow(void *mc, field *f, real *w, real *dtw) {
  MacroCell *mcell = (MacroCell*) mc;

  // Local copy of the interpretation parameters
  int iparam[8];
  for(int ip = 0; ip < 8; ip++)
    iparam[ip] = f->interp_param[ip];

  // Init to zero the time derivative
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    for(int ipg = 0; ipg < NPG(iparam + 1); ipg++) {
      for(int iv = 0; iv < f->model.m; iv++) {
	//int imem = f->varindex(iparam, ie, ipg, iv);
	//f->dtwn[imem] = 0;
      }
    }
  }
  //assert(sizew==f->macromesh.nbelems * f->model.m * NPG(iparam + 1));

  // Assembly of the surface terms
  // Loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // loop on the 6 faces
    // or four faces for 2d computations
    int nbfa = 6;
    if (f->macromesh.is2d) nbfa = 4;
    if (f->macromesh.is1d) nbfa = 2;
    for(int nifa = 0; nifa < nbfa; nifa++) {
      // get the right elem or the boundary id
      int ifa= f->macromesh.is1d ?  2 * nifa + 1 : nifa;
      int ieR = f->macromesh.elem2elem[6*ie+ifa];
      real physnodeR[20][3];
      if (ieR >= 0) {
      	for(int inoloc = 0; inoloc < 20; inoloc++) {
      	  int ino = f->macromesh.elem2node[20 * ieR + inoloc];
      	  physnodeR[inoloc][0] = f->macromesh.node[3 * ino + 0];
      	  physnodeR[inoloc][1] = f->macromesh.node[3 * ino + 1];
      	  physnodeR[inoloc][2] = f->macromesh.node[3 * ino + 2];
      	}
      }

      // loop on the glops (numerical integration)
      // of the face ifa
      for(int ipgf = 0; ipgf < NPGF(f->interp_param + 1, ifa); ipgf++) {
  	real xpgref[3], xpgref_in[3], wpg;
  	//real xpgref2[3], wpg2;
  	// get the coordinates of the Gauss point
	// and coordinates of a point slightly inside the
	// opposite element in xref_in
  	ref_pg_face(iparam + 1, ifa, ipgf, xpgref, &wpg, xpgref_in);

/* #ifdef _PERIOD */
/* 	assert(f->is1d); // TODO: generalize to 2d */
/* 	if (xpgref_in[0] > _PERIOD) { */
/* 	  //printf("à droite ifa= %d x=%f ipgf=%d ieL=%d ieR=%d\n", */
/* 	  //	 ifa,xpgref_in[0],ipgf,ie,ieR); */
/* 	  xpgref_in[0] -= _PERIOD; */
/* 	} */
/* 	else if (xpgref_in[0] < 0) {  */
/* 	  //printf("à gauche ifa= %d  x=%f ieL=%d ieR=%d \n", */
/* 	  //ifa,xpgref_in[0],ie,ieR); */
/* 	  xpgref_in[0] += _PERIOD; */
/* 	} */
/* #endif */

  	// recover the volume gauss point from
  	// the face index
  	int ipg = iparam[7];
  	// get the left value of w at the gauss point
  	real wL[f->model.m], wR[f->model.m];
  	for(int iv = 0; iv < f->model.m; iv++) {
  	  int imem = f->varindex(iparam, ie, ipg, iv);
  	  wL[iv] = f->wn[imem];
  	}
  	// the basis functions is also the gauss point index
  	int ib = ipg;
        int ipgL = ipg;
  	// normal vector at gauss point ipg
  	real dtau[3][3], codtau[3][3], xpg[3];
  	real vnds[3];
  	Ref2Phy(physnode,
  		xpgref,
  		NULL, ifa, // dpsiref, ifa
  		xpg, dtau,
  		codtau, NULL, vnds); // codtau, dpsi, vnds
  	real flux[f->model.m];
  	if (ieR >=0) {  // the right element exists
  	  // find the corresponding point in the right elem
	  real xpg_in[3];
	  Ref2Phy(physnode,
		  xpgref_in,
		  NULL, ifa, // dpsiref, ifa
		  xpg_in, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
  	  real xref[3];
	  PeriodicCorrection(xpg_in,f->macromesh.period);
	  Phy2Ref(physnodeR, xpg_in, xref);
  	  int ipgR = ref_ipg(iparam + 1, xref);
	  real xpgR[3], xrefR[3], wpgR;
	  ref_pg_vol(iparam + 1, ipgR, xrefR, &wpgR, NULL);
	  Ref2Phy(physnodeR,
		  xrefR,
		  NULL, -1, // dphiref, ifa
		  xpgR, NULL,
		  NULL, NULL, NULL); // codtau, dphi, vnds

/* #ifdef _PERIOD */
/* 	  assert(fabs(Dist(xpg,xpgR)-_PERIOD)<1e-11); */
/* #else */
/*           assert(Dist(xpg,xpgR)<1e-11); */
/* #endif */
	  
  	  for(int iv = 0; iv < f->model.m; iv++) {
  	    int imem = f->varindex(iparam, ieR, ipgR, iv);
  	    wR[iv] = f->wn[imem];
  	  }
  	  // int_dL F(wL, wR, grad phi_ib )
  	  f->model.NumFlux(wL, wR, vnds, flux);

  	}
  	else { //the right element does not exist
  	  f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
  	}
  	for(int iv = 0; iv < f->model.m; iv++) {
  	  int imem = f->varindex(iparam, ie, ib, iv);
  	  f->dtwn[imem] -= flux[iv] * wpg;
  	}

      }

    }
  }
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
void DGMacroCellInterface(void *mc, field *f, real *w, real *dtw) 
{
  MacroFace *mface = (MacroFace*) mc;
  MacroMesh *msh = &f->macromesh;
  const unsigned int m = f->model.m;

  int iparam[8];
  for(int ip = 0; ip < 8; ip++)
    iparam[ip] = f->interp_param[ip];

  // Assembly of the surface terms loop on the macrocells faces
  for (int ifa = mface->first; ifa < mface->last_p1; ifa++) {
    int ieL = msh->face2elem[4 * ifa + 0];
    int locfaL = msh->face2elem[4 * ifa + 1];

    // Get the physical nodes of element ieL
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = msh->elem2node[20 * ieL + inoloc];
      physnode[inoloc][0] = msh->node[3 * ino + 0];
      physnode[inoloc][1] = msh->node[3 * ino + 1];
      physnode[inoloc][2] = msh->node[3 * ino + 2];
    }

    int ieR = msh->face2elem[4 * ifa + 2];
    int locfaR = msh->face2elem[4 * ifa + 3];
    real physnodeR[20][3];
    if (ieR >= 0) {
      for(int inoloc = 0; inoloc < 20; inoloc++) {
        int ino = msh->elem2node[20 * ieR + inoloc];
        physnodeR[inoloc][0] = msh->node[3 * ino + 0];
        physnodeR[inoloc][1] = msh->node[3 * ino + 1];
        physnodeR[inoloc][2] = msh->node[3 * ino + 2];
      }
    }

    // Loop over the points on a single macro cell interface.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ipgfL = 0; ipgfL < NPGF(f->interp_param + 1, locfaL); ipgfL++) {

      int iparam[8];
      for(int ip = 0; ip < 8; ip++)
	iparam[ip] = f->interp_param[ip];

      real xpgref[3], xpgref_in[3], wpg;
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      ref_pg_face(iparam + 1, locfaL, ipgfL, xpgref, &wpg, xpgref_in);

/* #ifdef _PERIOD */
/*       assert(f->is1d); // TODO: generalize to 2d */
/*       if (xpgref_in[0] > _PERIOD) { */
/* 	//printf("à droite ifa= %d x=%f ipgf=%d ieL=%d ieR=%d\n", */
/* 	//	 ifa,xpgref_in[0],ipgf,ie,ieR); */
/* 	xpgref_in[0] -= _PERIOD; */
/*       } */
/*       else if (xpgref_in[0] < 0) {  */
/* 	//printf("à gauche ifa= %d  x=%f ieL=%d ieR=%d \n", */
/* 	//ifa,xpgref_in[0],ie,ieR); */
/* 	xpgref_in[0] += _PERIOD; */
/*       } */
/* #endif */

      // Recover the volume gauss point from the face index
      int ipgL = iparam[7];

      real flux[m];
      real wL[m];

      // Normal vector at gauss point ipgL
      real vnds[3], xpg[3];
      {
	real dtau[3][3], codtau[3][3];
	Ref2Phy(physnode,
		xpgref,
		NULL, locfaL, // dpsiref, ifa
		xpg, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
      }

      if (ieR >= 0) {  // the right element exists
        real xrefL[3];
	{
	  real xpg_in[3];
	  Ref2Phy(physnode,
		  xpgref_in,
		  NULL, -1, // dpsiref, ifa
		  xpg_in, NULL,
		  NULL, NULL, NULL); // codtau, dpsi, vnds
	  PeriodicCorrection(xpg_in,f->macromesh.period);
	  Phy2Ref(physnodeR, xpg_in, xrefL);

	}

        int ipgR = ref_ipg(iparam + 1, xrefL);

	//printf("ipgL=%d ipgR=%d\n",ipgL,ipgR);

	// Uncomment to check that the neighbour-finding algorithm worked.
	/* { */
	/*   real xpgR[3], xrefR[3], wpgR; */
	/*   ref_pg_vol(iparam + 1, ipgR, xrefR, &wpgR, NULL); */
	/*   Ref2Phy(physnodeR, */
	/* 	  xrefR, */
	/* 	  NULL, -1, // dphiref, ifa */
	/* 	  xpgR, NULL, */
	/* 	  NULL, NULL, NULL); // codtau, dphi, vnds */
	/*  #ifdef _PERIOD */
	/*   assert(fabs(Dist(xpg,xpgR)-_PERIOD)<1e-11); */
	/*   #else */
	/* assert(Dist(xpg,xpgR)<1e-11); */
	/* #endif */
	/* }	 */

	real wR[m];
        for(int iv = 0; iv < m; iv++) {
	  int imemL = f->varindex(iparam, ieL, ipgL, iv);
	  wL[iv] = w[imemL];
          int imemR = f->varindex(iparam, ieR, ipgR, iv);
          wR[iv] = w[imemR];
        }

        // int_dL F(wL, wR, grad phi_ib)

        f->model.NumFlux(wL, wR, vnds, flux);

	// Add flux to both sides
	for(int iv = 0; iv < m; iv++) {
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(iparam, ieL, ipgL, iv);
          int imemR = f->varindex(iparam, ieR, ipgR, iv);
	  dtw[imemL] -= flux[iv] * wpg;
          dtw[imemR] += flux[iv] * wpg;
	}

      } else { // The point is on the boundary.
	for(int iv = 0; iv < m; iv++) {
	  int imemL = f->varindex(iparam, ieL, ipgL, iv);
	  wL[iv] = w[imemL];
	}

        f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);

	for(int iv = 0; iv < m; iv++) {
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(iparam, ieL, ipgL, iv);
	  dtw[imemL] -= flux[iv] * wpg;
	}
      }

    }

  }
}

// Apply division by the mass matrix
void DGMass(void *mc, field *f, real *dtw) 
{
  MacroCell *mcell = (MacroCell*)mc;

  int m = f->model.m;

  // loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
      ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real det = dot_product(dtau[0], codtau[0]);
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	dtw[imem] /= (wpg * det);
      }
    }
  }
}

// Apply the source term
void DGSource(void *mc, field *f, real *w, real *dtw) 
{
  if (f->model.Source == NULL) {
    return;
  }

  MacroCell *mcell = (MacroCell*)mc;

  const int m = f->model.m;

  // Loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
      ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real wL[m], source[m];
      for(int iv = 0; iv < m; ++iv){
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	wL[iv] = w[imem];
      }
      
      f->model.Source(xphy, f->tnow, wL, source);
      
      for(int iv = 0; iv < m; ++iv) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	dtw[imem] += source[iv];
      }
    }
  }
}

// Compute the Discontinuous Galerkin volume terms, fast version
void DGVolume(void *mc, field *f, real *w, real *dtw) 
{
  MacroCell *mcell = (MacroCell*) mc;

  // loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    const int m = f->model.m;
    const int deg[3] = {f->interp_param[1],
			f->interp_param[2],
			f->interp_param[3]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    const int nraf[3] = {f->interp_param[4],
			 f->interp_param[5],
			 f->interp_param[6]};

    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

    int f_interp_param[8] = {f->interp_param[0],
			     f->interp_param[1],
			     f->interp_param[2],
			     f->interp_param[3],
			     f->interp_param[4],
			     f->interp_param[5],
			     f->interp_param[6],
			     f->interp_param[7]};

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
	  real *xref0 = malloc(sc_npg * sizeof(real));
	  real *xref1 = malloc(sc_npg * sizeof(real));
	  real *xref2 = malloc(sc_npg * sizeof(real));
	  real *omega = malloc(sc_npg * sizeof(real));
	  int *imems = malloc(m * sc_npg * sizeof(int));
	  int pos = 0;
	  for(unsigned int p = 0; p < sc_npg; ++p) {
	    real xref[3];
	    real tomega;

	    ref_pg_vol(f->interp_param + 1, offsetL + p, xref, &tomega, NULL);
	    xref0[p] = xref[0];
	    xref1[p] = xref[1];
	    xref2[p] = xref[2];
	    omega[p] = tomega;

	    for(int im = 0; im < m; ++im) {
	      imems[pos++] = f->varindex(f_interp_param, ie, offsetL + p, im);
	    }
	  }

	  // loop in the "cross" in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) {
	    // for(int dim0 = 0; dim0 < 2; dim0++) {  // TODO : return to 3d !
	    // point p at which we compute the flux

	    for(int p0 = 0; p0 < npg[0]; p0++) {
	      for(int p1 = 0; p1 < npg[1]; p1++) {
		for(int p2 = 0; p2 < npg[2]; p2++) {
		  real wL[m], flux[m];
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
		    real dphiref[3] = {0, 0, 0};
		    // compute grad phi_q at glop p
		    dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) 
		      * nraf[dim0];

		    real xrefL[3] = {xref0[ipgL - offsetL],
				       xref1[ipgL - offsetL],
				       xref2[ipgL - offsetL]};
		    real wpgL = omega[ipgL - offsetL];
		    /* real xrefL[3], wpgL; */
		    /* ref_pg_vol(f->interp_param+1,ipgL,xrefL, &wpgL, NULL); */

		    // mapping from the ref glop to the physical glop
		    real dtau[3][3], codtau[3][3], dphiL[3];
		    Ref2Phy(physnode,
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
		      int imemR = f->varindex(f_interp_param, ie, ipgR, iv);
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
}


void dtfield_pthread(field *f) 
{
#ifdef _WITH_PTHREAD
  bool facealgo = true;
  //facealgo=false;
  if(facealgo) {
    for(int iw = 0; iw < f->wsize; iw++)
      f->dtwn[iw] = 0;
    for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++)
      DGMacroCellInterface((void*) (f->mface + ifa), f, f->wn, f->dtwn);
  }


  // One flying thread per macrocell
  pthread_t tmcell[f->macromesh.nbelems];
  int status;

  // Launch a thread for each macro cell computation of the inter
  // subcell fluxes
  if (!facealgo) {
    for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
      status = pthread_create (&tmcell[ie], // thread
			       NULL,                 // default attributes
			       DGMacroCellInterfaceSlow,  // called function
			       (void*) (mcell+ie));  // function params
      assert(status==0);
    }
    // Wait the end of the threads before next step
    for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
      pthread_join(tmcell[ie], NULL);
    }
  }

  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    status = pthread_create (&tmcell[ie], // thread
			     NULL,                 // default attributes
			     DGSubCellInterface,  // called function
			     (void*) (mcell+ie));  // function params
    assert(status==0);
  }
  // wait the end of the threads before next step
  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    pthread_join(tmcell[ie], NULL);
  }

  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    status = pthread_create (&tmcell[ie], // thread
			     NULL,                 // default attributes
			     DGVolume,  // called function
			     (void*) (mcell+ie));  // function params
    assert(status==0);
  }
  // Wait the end of the threads before next step
  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    pthread_join(tmcell[ie], NULL);
  }

  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    status = pthread_create (&tmcell[ie], // thread
			     NULL,                 // default attributes
			     DGMass,  // called function
			     (void*) (mcell+ie));  // function params
    assert(status==0);
  }
  // Wait the end of the threads before next step
  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    pthread_join(tmcell[ie], NULL);
  }
#endif
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field
void dtfield(field *f, real *w, real *dtw) {
  if(f->pre_dtfield != NULL) // FIXME: rename to before dtfield
      f->pre_dtfield(f, w);

#ifdef _WITH_PTHREAD
  dtfield_pthread(f, w, dtw);
#else

#ifdef _OPENMP
#pragma omp parallel
#endif
  for(int iw = 0; iw < f->wsize; iw++)
    dtw[iw] = 0;

  bool facealgo = true;
  if(facealgo)
    for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++){
      DGMacroCellInterface((void*) (f->mface + ifa), f, w, dtw);
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    MacroCell *mcelli = f->mcell + ie;
    if(!facealgo) DGMacroCellInterfaceSlow(mcelli, f, w, dtw);
    DGSubCellInterface(mcelli, f, w, dtw);
    DGVolume(mcelli, f, w, dtw);
    DGMass(mcelli, f, dtw);
    DGSource(mcelli, f, w, dtw);
  }
#endif
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field
void dtfieldSlow(field* f) 
{
  // Interpolation params
  // Warning: this is ugly, but the last parameter is used for
  // computing the volume GLOP index from the face GLOP index...
  // Ugly too: the first parameter is not used by all
  // utilities. we have sometimes to jump over : pass param + 1
  // instead of param...

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};

  // init to zero the time derivative
  int sizew = 0;
  //#pragma omp parallel for
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	f->dtwn[imem] = 0;
        sizew++;
      }
    }
  }
  assert(sizew == f->macromesh.nbelems * f->model.m * NPG(f->interp_param + 1));

  // assembly of the subrface terms loop on the elements
  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // loop on the 6 faces
    // or four faces for 2d computations
    int nbfa = 6;
    if (f->macromesh.is2d) nbfa = 4;
    for(int ifa = 0; ifa < nbfa; ifa++) {
      // get the right elem or the boundary id
      int ieR = f->macromesh.elem2elem[6*ie+ifa];
      real physnodeR[20][3];
      if (ieR >= 0) {
      	for(int inoloc = 0; inoloc < 20; inoloc++) {
      	  int ino = f->macromesh.elem2node[20 * ieR + inoloc];
      	  physnodeR[inoloc][0] = f->macromesh.node[3 * ino + 0];
      	  physnodeR[inoloc][1] = f->macromesh.node[3 * ino + 1];
      	  physnodeR[inoloc][2] = f->macromesh.node[3 * ino + 2];
      	}
      }

      // Loop on the glops (numerical integration) of the face ifa
      for(int ipgf = 0; ipgf < NPGF(f->interp_param + 1, ifa); ipgf++) {
  	real xpgref[3], xpgref_in[3], wpg;
  	//real xpgref2[3], wpg2;
  	// get the coordinates of the Gauss point
	// and coordinates of a point slightly inside the
	// opposite element in xref_in
  	ref_pg_face(f->interp_param + 1, ifa, ipgf, xpgref, &wpg, xpgref_in);

  	// recover the volume gauss point from
  	// the face index
  	int ipg = f->interp_param[7];
  	// get the left value of w at the gauss point
  	real wL[f->model.m], wR[f->model.m];
  	for(int iv = 0; iv < f->model.m; iv++) {
  	  int imem = f->varindex(f->interp_param, ie, ipg, iv);
  	  wL[iv] = f->wn[imem];
  	}
  	// the basis functions is also the gauss point index
  	int ib = ipg;
        int ipgL = ipg;
  	// normal vector at gauss point ipg
  	real dtau[3][3], codtau[3][3], xpg[3];
  	real vnds[3];
  	Ref2Phy(physnode,
  		xpgref,
  		NULL, ifa, // dpsiref, ifa
  		xpg, dtau,
  		codtau, NULL, vnds); // codtau, dpsi, vnds
  	real flux[f->model.m];
  	if (ieR >=0) {  // the right element exists
  	  // find the corresponding point in the right elem
	  real xpg_in[3];
	  Ref2Phy(physnode,
		  xpgref_in,
		  NULL, ifa, // dpsiref, ifa
		  xpg_in, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
  	  real xref[3];
	  Phy2Ref(physnodeR, xpg_in, xref);
  	  int ipgR = ref_ipg(f->interp_param + 1, xref);
	  real xpgR[3], xrefR[3], wpgR;
	  ref_pg_vol(f->interp_param + 1, ipgR, xrefR, &wpgR, NULL);
	  Ref2Phy(physnodeR,
		  xrefR,
		  NULL, -1, // dphiref, ifa
		  xpgR, NULL,
		  NULL, NULL, NULL); // codtau, dphi, vnds
	  assert(Dist(xpgR, xpg) < 1e-10);
  	  for(int iv = 0; iv < f->model.m; iv++) {
  	    int imem = f->varindex(f->interp_param, ieR, ipgR, iv);
  	    wR[iv] = f->wn[imem];
  	  }
  	  // int_dL F(wL, wR, grad phi_ib )
  	  f->model.NumFlux(wL, wR, vnds, flux);

  	}
  	else { //the right element does not exist
  	  f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
  	}
  	for(int iv = 0; iv < f->model.m; iv++) {
  	  int imem = f->varindex(f->interp_param, ie, ib, iv);
  	  f->dtwn[imem] -= flux[iv]*wpg;
  	}

      }

    }
  }

  // Assembly of the volume terms loop on the elements
  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // Mass matrix
    real masspg[NPG(f->interp_param + 1)];
    // Loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real xpgref[3], wpg;
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);

      // Get the value of w at the gauss point
      real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	w[iv] = f->wn[imem];
      }

      // Loop on the basis functions
      for(int ib = 0; ib < NPG(f->interp_param + 1); ib++) {
	// gradient of psi_ib at gauss point ipg
	real dpsiref[3], dpsi[3];
	real dtau[3][3], codtau[3][3];//, xpg[3];
	grad_psi_pg(f->interp_param + 1, ib, ipg, dpsiref);
	Ref2Phy(physnode, // phys. nodes
		xpgref, // xref
		dpsiref, -1, // dpsiref, ifa
		NULL, dtau, // xphy, dtau
		codtau, dpsi, NULL); // codtau, dpsi, vnds
	// remember the diagonal mass term
	if (ib == ipg) {
	  real det = dot_product(dtau[0], codtau[0]);
	  masspg[ipg] = wpg * det;
	}
	// int_L F(w, w, grad phi_ib )
	real flux[f->model.m];
	f->model.NumFlux(w, w, dpsi, flux);

	for(int iv = 0; iv < f->model.m; iv++) {
	  int imem = f->varindex(f->interp_param, ie, ib, iv);
	  f->dtwn[imem] += flux[iv] * wpg;
	}
      }
    }

    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      // apply the inverse of the diagonal mass matrix
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	f->dtwn[imem] /= masspg[ipg];
      }
    }

  }
}

// An out-of-place RK step
void RK_out(real *dest, real *fwn, real *fdtwn, const real dt, 
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
void RK_in(real *fwnp1, real *fdtwn, const real dt, const int sizew)
{
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int iw = 0; iw < sizew; iw++) {
    fwnp1[iw] += dt * fdtwn[iw];
  }
}

real set_dt(field *f)
{
  return f->model.cfl * f->hmin / f->vmax; 
}

// Time integration by a second-order Runge-Kutta algorithm
void RK2(field *f, real tmax, real dt) 
{
  if(dt <= 0)
    dt = set_dt(f);

  f->itermax = tmax / dt;
  int size_diags;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(f->interp_param + 1);
  int iter = 0;

  real *wnp1 = calloc(f->wsize, sizeof(real));
  assert(wnp1);

  // FIXME: remove
  size_diags = f->nb_diags * f->itermax;
  f->iter_time = iter;
  
  if(f->nb_diags != 0)
    f->Diagnostics = malloc(size_diags * sizeof(real));

  while(f->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    dtfield(f, f->wn, f->dtwn);
    RK_out(wnp1, f->wn, f->dtwn, 0.5 * dt, sizew);

    f->tnow += 0.5 * dt;

    dtfield(f, wnp1, f->dtwn);
    RK_in(f->wn, f->dtwn, dt, sizew);

    f->tnow += 0.5 * dt;

    if(f->update_after_rk != NULL)
      f->update_after_rk(f, f->wn);

    iter++;
    f->iter_time=iter;
  }
  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
  free(wnp1);
}

void RK4_final_inplace(real *w, real *l1, real *l2, real *l3, 
		       real *dtw, const real dt, const int sizew)
{
  const real b = -1.0 / 3.0;
  const real a[] = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, dt / 6.0};
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

// Time integration by a fourth-order Runge-Kutta algorithm
void RK4(field *f, real tmax, real dt) 
{
  if(dt <= 0)
    dt = set_dt(f);

  f->itermax = tmax / dt;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(f->interp_param + 1);
  int iter = 0;

  // Allocate memory for RK time-stepping
  real *l1, *l2, *l3;
  l1 = calloc(sizew, sizeof(real));
  l2 = calloc(sizew, sizeof(real));
  l3 = calloc(sizew, sizeof(real));

  while(f->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    // l_1 = w_n + 0.5dt * S(w_n, t_0)
    dtfield(f, f->wn, f->dtwn);
    RK_out(l1, f->wn, f->dtwn, 0.5 * dt, sizew);

    f->tnow += 0.5 * dt;

    // l_2 = w_n + 0.5dt * S(l_1, t_0 + 0.5 * dt)
    dtfield(f, l1, f->dtwn);
    RK_out(l2, f->wn, f->dtwn, 0.5 * dt, sizew);

    // l_3 = w_n + dt * S(l_2, t_0 + 0.5 * dt)
    dtfield(f, l2, f->dtwn);
    RK_out(l3, f->wn, f->dtwn, dt, sizew);

    f->tnow += 0.5 * dt;

    // Compute S(l_3, t_0 + dt)
    dtfield(f, l3, f->dtwn);
    RK4_final_inplace(f->wn, l1, l2, l3, f->dtwn, dt, sizew);

    iter++;
  }
  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

  free(l3);
  free(l2);
  free(l1);
}

// Compute the normalized L2 distance with the imposed data
real L2error(field *f) {
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  real error = 0;
  real mean = 0;

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
      real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ie, ipg, iv);
	w[iv] = f->wn[imem];
      }

      real wex[f->model.m];
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

	// Get the exact value
	f->model.ImposedData(xphy, f->tnow, wex);
      }

      for(int iv = 0; iv < f->model.m; iv++) {
	real diff = w[iv] - wex[iv];
        error += diff * diff * wpg * det;
        mean += wex[iv] * wex[iv] * wpg * det;
      }
    }
  }
  return sqrt(error) / (sqrt(mean)  + 1e-16);
}


void InterpField(field *f, int ie, real *xref, real *w){
  const int nraf[3] = {f->interp_param[4],
		       f->interp_param[5],
		       f->interp_param[6]};
  const int deg[3] = {f->interp_param[1],
		      f->interp_param[2],
		      f->interp_param[3]};

  for(int iv = 0; iv < f->model.m; iv++)
    w[iv] = 0;

  int is[3];

  for(int ii = 0; ii < 3; ii++){
    is[ii] = xref[ii] * nraf[ii];
    assert(is[ii] < nraf[ii] && is[ii]>= 0);
  }
  
  int npgv = NPG(f->interp_param + 1);
  // TODO: loop only on non zero basis function
  for(int ib = 0; ib < npgv; ib++) { 
    real psi;
    psi_ref_subcell(f->interp_param + 1, is, ib, xref, &psi, NULL);
    
    for(int iv=0;iv<f->model.m;iv++){
      int imem = f->varindex(f->interp_param, ie, ib, iv);
      w[iv] += psi * f->wn[imem];
    }
  }
}

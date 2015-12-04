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
int GenericVarindex(__constant int *param, int ipg, int iv)
{
  // param[0] = m
  return iv + param[0] * ipg;
}
#pragma end_opencl

real min_grid_spacing(field *f)
{
  real hmin = FLT_MAX;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    MacroCell *mcell = f->mcell + ie;

    real vol = 0, surf = 0;

    // Loop on the glops (for numerical integration)
    for(int ipg = 0; ipg < NPG(mcell->raf, mcell->deg); ipg++) {
      real xpgref[3], wpg;
      // Get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param + 4, f->interp_param + 1,
		 ipg, xpgref, &wpg, NULL);
      real codtau[3][3], dtau[3][3];
      Ref2Phy(mcell->physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      NULL, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real det = dot_product(dtau[0], codtau[0]);
      vol += wpg * det;
    }
    
    for(int ifa = 0; ifa < 6; ifa++) {
      // loop on the faces
      for(int ipgf = 0; ipgf < NPGF(mcell->raf, mcell->deg, ifa); ipgf++) {
	real xpgref[3], wpg;
	// get the coordinates of the Gauss point
	int *deg = f->interp_param + 1;
	int *raf = f->interp_param + 4;
	ref_pg_face(raf, deg, ifa, ipgf, xpgref, &wpg, NULL);
	real vnds[3];
	{
	  real codtau[3][3], dtau[3][3];
	  Ref2Phy(mcell->physnode,
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
    MacroCell *mcell = f->mcell + ie;
    
    for(int ipg = 0; ipg < mcell->npg; ipg++) {
      real xpg[3];
      real xref[3], omega;
      ref_pg_vol(f->interp_param + 4, f->interp_param + 1,
		 ipg, xref, &omega, NULL);
      real dtau[3][3];
      Ref2Phy(mcell->physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds
      { // Check the reverse transform at all the GLOPS
 	real xref2[3];
	Phy2Ref(mcell->physnode, xpg, xref2);

	real tolerance = sizeof(real) == sizeof(double) ? 1e-8 : 1e-5;
	assert(Dist(xref, xref2) < tolerance);
      }

      real w[f->model.m];
      f->model.InitData(xpg, w);
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ipg, iv) + mcell->woffset;
	f->wn[imem] = w[iv];
      }
    }
  }
}

#ifdef _WITH_OPENCL
void set_physnodes_cl(field *f) 
{
  const int nmacro = f->macromesh.nbelems;
  real *physnode = malloc(60 * sizeof(real));

  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;

    cl_int status;
    const size_t buf_size = sizeof(real) * 60;
    mcell->physnode_cl = clCreateBuffer(f->cli.context,
					CL_MEM_READ_ONLY,
					buf_size,
					NULL,
					&status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    int ie20 = 20 * ie;
    for(int inoloc = 0; inoloc < 20; ++inoloc) {
      int ino = 3 * f->macromesh.elem2node[ie20 + inoloc];
      real *iphysnode = physnode + 3 * inoloc;
      real *nodeino = f->macromesh.node + ino;
      iphysnode[0] = nodeino[0];
      iphysnode[1] = nodeino[1];
      iphysnode[2] = nodeino[2];
    }
   
    status = clEnqueueWriteBuffer(f->cli.commandqueue,
				  mcell->physnode_cl, // cl_mem buffer,
				  CL_TRUE,// cl_bool blocking_read,
				  0, // size_t offset
				  buf_size, // size_t cb
				  physnode, //  	void *ptr,
				  0, 0, 0);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

  }

  free(physnode);
}

void init_field_buffers_cl(field *f)
{
  cl_int status;

  const int nmacro = f->macromesh.nbelems;

  f->wn_cl = calloc(nmacro, sizeof(cl_mem));
  f->dtwn_cl = calloc(nmacro, sizeof(cl_mem));

  f->param_cl = clCreateBuffer(f->cli.context,
			       CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			       sizeof(int) * 7,
			       f->interp_param,
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  set_physnodes_cl(f);
}

void init_field_kernels_cl(field *f)
{
  cl_int status;
  
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

  f->ExtractInterface = clCreateKernel(f->cli.program,
				       "ExtractInterface",
				       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->InsertInterface = clCreateKernel(f->cli.program,
				      "InsertInterface",
				      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->ExtractedDGInterfaceFlux = clCreateKernel(f->cli.program,
					       "ExtractedDGInterfaceFlux",
					       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  f->ExtractedDGBoundaryFlux = clCreateKernel(f->cli.program,
					      "ExtractedDGBoundaryFlux",
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

  f->RK4_first_stages = clCreateKernel(f->cli.program,
				       "RK4_first_stages",
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

  f->empty_kernel = clCreateKernel(f->cli.program,
				   "empty_kernel",
				   &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void init_field_events_cl(field *f)
{
  cl_int status;
  
  const int nmacro = f->macromesh.nbelems;
  
  f->clv_zbuf = calloc(nmacro, sizeof(cl_event));
  
  const int ninterfaces = f->macromesh.nmacrointerfaces;
  if(ninterfaces > 0)
    f->clv_mci = calloc(ninterfaces, sizeof(cl_event));
    
  const int nbound = f->macromesh.nboundaryfaces;
  if(nbound > 0)
    f->clv_boundary = calloc(nbound, sizeof(cl_event));

  f->clv_flux = calloc(3, sizeof(cl_event*));
  for(int dim = 0; dim < 3; ++dim) {
    f->clv_flux[dim] = calloc(nmacro, sizeof(cl_event));
  }
  
  f->clv_volume = calloc(nmacro, sizeof(cl_event));
  f->clv_source = calloc(nmacro, sizeof(cl_event));
  f->clv_mass = calloc(nmacro, sizeof(cl_event));
}

void init_field_macrocells_cl(field *f)
{
  cl_int status;
  const int m = f->model.m;
  
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };
  
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    MacroCell *mcell = f->mcell + ie; 

    mcell->wn_cl = clCreateBuffer(f->cli.context,
				  CL_MEM_READ_WRITE,
				  sizeof(real) * mcell->nreal,
				  NULL,
				  &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    f->wn_cl[ie] = mcell->wn_cl;
    
    mcell->dtwn_cl = clCreateBuffer(f->cli.context,
				    CL_MEM_READ_WRITE,
				    sizeof(real) * mcell->nreal,
				    NULL,
				    &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    f->dtwn_cl[ie] = mcell->dtwn_cl;

    const size_t buf_size = sizeof(real) * mcell->npg;
    mcell->mass_cl = clCreateBuffer(f->cli.context,
				    CL_MEM_READ_ONLY,
				    buf_size,
				    NULL,
				    &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    status = clEnqueueWriteBuffer(f->cli.commandqueue,
				  mcell->mass_cl,
				  CL_TRUE,
				  0,
				  mcell->npg * sizeof(real),
				  mcell->mass,
				  0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    status = clFinish(f->cli.commandqueue);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    mcell->interface_cl = malloc(6 * sizeof(cl_mem));
    for(int ifa = 0; ifa < 6; ++ifa) {
      int d0 = axis_permut[ifa][0];
      int d1 = axis_permut[ifa][1];
      int npgf = mcell->raf[d0] * (mcell->deg[d0] + 1)
	* mcell->raf[d1] * (mcell->deg[d1] + 1); 
      size_t bufsize = sizeof(real) * m * npgf;
      assert(bufsize > 0);
      
      mcell->interface_cl[ifa] = clCreateBuffer(f->cli.context,
						CL_MEM_READ_WRITE,
						bufsize,
						NULL,
						&status);
      if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
      assert(status >= CL_SUCCESS);
    }
    
  }
}
  
void init_field_MacroFaces_cl(field *f)
{
  cl_int status;

  int *param = f->interp_param;
  int m = param[0];
  
  // Set up interfaces
  const int ninterfaces = f->macromesh.nmacrointerfaces;
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f->macromesh.macrointerface[i];
    MacroFace *mface = f->mface + ifa;
  
    MacroCell *mcellL = f->mcell + mface->ieL;
    mface->wL_cl = mcellL->interface_cl[mface->locfaL];
    
    MacroCell *mcellR = f->mcell + mface->ieR;
    mface->wR_cl = mcellR->interface_cl[mface->locfaR];
  }
}
  
void init_field_cl(field *f)
{
  InitCLInfo(&f->cli, nplatform_cl, ndevice_cl);

  init_field_buffers_cl(f);
  // Program compilation
  char *strprog;
  GetOpenCLCode();
  ReadFile("schnaps.cl", &strprog);

  printf("\t%s\n", numflux_cl_name);

  // If the source term is set (via set_source_CL) then add it to the
  // buildoptions and compile using the new buildoptions.
  if(f->use_source_cl) {
    char *temp;
    int len0 = strlen(cl_buildoptions) + 2;
    char *D_SOURCE_FUNC = " -D_SOURCE_FUNC=";
    int len1 = strlen(D_SOURCE_FUNC) + 2;
    int len2 = strlen(f->sourcename_cl) + 2;
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

  init_field_kernels_cl(f);

  init_field_events_cl(f);

  init_field_macrocells_cl(f);

  init_field_MacroFaces_cl(f);
  
  // Set timers to zero
  f->zbuf_time = 0;
  f->mass_time = 0;
  f->vol_time = 0;
  f->flux_time = 0;
  f->minter_time = 0;
  f->boundary_time = 0;
  f->source_time = 0;
  f->rk_time = 0;
}
#endif

void ipgf_to_xphy(MacroCell *mcell, int locfa, int ipgf, real *xphy)
{
  real xref[3];
  int ipgvL = ref_pg_face(mcell->raf, mcell->deg, locfa, ipgf, 
			  xref,  NULL, NULL);
  Ref2Phy(mcell->physnode, xref, NULL, -1, xphy, NULL, NULL, NULL, NULL);
}

void icix_to_xphy(MacroCell *mcell, int *ic, int *ix, real *xphy)
{
  int ipg = xyz_to_ipg(mcell->raf, mcell->deg, ic, ix);
  real xref[3];
  ref_pg_vol(mcell->raf, mcell->deg, ipg, xref, NULL, NULL);
  Ref2Phy(mcell->physnode, xref, NULL, -1, xphy, NULL, NULL, NULL, NULL);
}

void test_MacroFace_orientation(field *f, MacroFace *mface)
{
  real tol = 1e-8;
  
  MacroCell *mcellL = f->mcell + mface->ieL;
  MacroCell *mcellR = f->mcell + mface->ieR;
  
  //printf("mface->Rcorner: %d\n", mface->Rcorner);
  
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };

  int paxisL[4] = {axis_permut[mface->locfaL][0],
		   axis_permut[mface->locfaL][1],
		   axis_permut[mface->locfaL][2],
		   axis_permut[mface->locfaL][3] } ;
    
  int paxisR[4] = {axis_permut[mface->locfaR][0],
		   axis_permut[mface->locfaR][1],
		   axis_permut[mface->locfaR][2],
		   axis_permut[mface->locfaR][3] } ;

  int d0L = paxisL[0];
  int d1L = paxisL[1];
  int d2L = paxisL[2];
  int signL = paxisL[3] == 0 ? -1 : 1;

  int d0R = paxisR[0];
  int d1R = paxisR[1];
  int d2R = paxisR[2];
  int signR = paxisR[3] == 0 ? -1 : 1;
  
  // Test that each left-surface point actually matches the expected
  // right-surface point.
  for(unsigned int picL0 = 0; picL0 < mcellL->raf[d0L]; ++picL0) {
    for(unsigned int picL1 = 0; picL1 < mcellL->raf[d1L]; ++picL1) {
      for(unsigned int pixL0 = 0; pixL0 <= mcellL->deg[d0L]; ++pixL0) {
	for(unsigned int pixL1 = 0; pixL1 <= mcellL->deg[d1L]; ++pixL1) {

	  int icL[3];
	  icL[d0L] = picL0;
	  icL[d1L] = picL1;
	  icL[d2L] = signL == -1 ? 0 : mcellL->raf[d2L] - 1;

	  int ixL[3];
	  ixL[d0L] = pixL0;
	  ixL[d1L] = pixL1;
	  ixL[d2L] = signL == -1 ? 0 : mcellL->deg[d2L];

	  // printf("icL: %d %d %d\n", icL[0], icL[1], icL[2]);
	  // printf("ixL: %d %d %d\n", ixL[0], ixL[1], ixL[2]);
	    
	  int icR[3];
	  icR[d2R] = signR == -1 ?  0 : mcellR->raf[d2R] - 1;

	  int ixR[3];
	  ixR[d2R] = signR == -1 ?  0 : mcellR->dnpg[d2R] - 1;

	  switch(mface->Rcorner) {
	  case 0:
	    // Rcorner = 0: R0 = L1, R1 = L0
	    icR[d0R] = icL[d1L];
	    icR[d1R] = icL[d0L];

	    ixR[d0R] = ixL[d1L];
	    ixR[d1R] = ixL[d0L];
	    break;

	  case 1:
	    // Rcorner = 1: R0 = L0, R1 = -L1
	    icR[d0R] = icL[d0L];
	    icR[d1R] = mcellL->raf[d1L] - icL[d1L] - 1;
	      
	    ixR[d0R] = ixL[d0L];
	    ixR[d1R] = mcellL->dnpg[d1L] - ixL[d1L] - 1;
	    break;

	  case 2:
	    // Rcorner = 2: R0 = -L1, R1 = -L0
	    icR[d0R] = mcellR->raf[d0L] - icL[d1L] - 1;
	    icR[d1R] = mcellR->raf[d1L] - icL[d0L] - 1;
	      
	    ixR[d0R] = mcellR->dnpg[d0L] - ixL[d1L] - 1;
	    ixR[d1R] = mcellR->dnpg[d1L] - ixL[d0L] - 1;
	    break;
	      
	  case 3:
	    // Rcorner = 3: R0 = -L0, R1 = -L1
	    icR[d0R] = mcellR->raf[d0L] - icL[d0L] - 1;
	    icR[d1R] = icL[d1L];
	      
	    ixR[d0R] = mcellR->dnpg[d0L] - ixL[d0L] - 1;
	    ixR[d1R] = ixL[d1L];
	    break;
	      
	  default:
	    assert(0);
	  }

	  //printf("icR: %d %d %d\n", icR[0], icR[1], icR[2]);
	  //printf("ixR: %d %d %d\n", ixR[0], ixR[1], ixR[2]);
	    
	  real xphyR[3];
	  icix_to_xphy(mcellR, icR, ixR, xphyR);

	  real xphyL[3];
	  icix_to_xphy(mcellL, icL, ixL, xphyL);

	  if(DistPeriodic(xphyR, xphyL, f->macromesh.period, 2 * tol) > tol) {
	    printf("icL: %d %d %d\t", icL[0], icL[1], icL[2]);
	    printf("ixL: %d %d %d\n", ixL[0], ixL[1], ixL[2]);
	    printf("icR: %d %d %d\t", icR[0], icR[1], icR[2]);
	    printf("ixR: %d %d %d\n", ixR[0], ixR[1], ixR[2]);
    
	    printf("\tdist: %f\n",
		   DistPeriodic(xphyR, xphyL, f->macromesh.period, 2 * tol));
	  }
	      
	  assert(DistPeriodic(xphyR, xphyL, f->macromesh.period, 2 * tol)
		 < tol);
	}
      }
    }
  }
}

void init_field_macrointerfaces(field *f)
{
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };

  const int ninterfaces = f->macromesh.nmacrointerfaces;
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f->macromesh.macrointerface[i];
    MacroFace *mface = f->mface + ifa;

    assert(mface->ieR != -1);
    
    // Determine the relative orientation of the faces.

    MacroCell *mcellL = f->mcell + mface->ieL;
    MacroCell *mcellR = f->mcell + mface->ieR;
    
    int paxisL[4] = {axis_permut[mface->locfaL][0],
		     axis_permut[mface->locfaL][1],
		     axis_permut[mface->locfaL][2],
		     axis_permut[mface->locfaL][3] } ;
    
    int paxisR[4] = {axis_permut[mface->locfaR][0],
		     axis_permut[mface->locfaR][1],
		     axis_permut[mface->locfaR][2],
		     axis_permut[mface->locfaR][3] } ;

    int d0L = paxisL[0];
    int d1L = paxisL[1];
    int d2L = paxisL[2];
    int signL = paxisL[3] == 0 ? -1 : 1;
    
    int icL[3] = {0, 0, 0};
    icL[d2L] = signL == -1 ? 0 : mcellL->raf[d2L] - 1;
    int ixL[3] = {0, 0, 0};
    ixL[d2L] = signL == -1 ? 0 : mcellL->deg[d2L];

    real tol = 1e-8;
    
    /*
    printf("period: %f %f %f\t",
	   f->macromesh.period[0],
	   f->macromesh.period[1],
	   f->macromesh.period[2]);
    */

    // printf("icL: %d %d %d\t", icL[0], icL[1], icL[2]);
    // printf("ixL: %d %d %d\n", ixL[0], ixL[1], ixL[2]);

    real xphyL0[3];
    icix_to_xphy(mcellL, icL, ixL, xphyL0);
    
    int d0R = paxisR[0];
    int d1R = paxisR[1];
    int d2R = paxisR[2];
    int signR = paxisR[3] == 0 ? -1 : 1;
        
    // Sanity check that both MacroCells think they share the same
    // number of points
    int npgfR = NPGF(mcellR->raf, mcellR->deg, mface->locfaR);
    assert(mface->npgf == npgfR);

    assert(mface->npgf == mcellR->npgdir[d0R] * mcellR->npgdir[d1R]);

    int icR[3];
    icR[d2R] = signR == -1 ? 0 : mcellR->raf[d2R] - 1;
    int ixR[3];
    ixR[d2R] = signR == -1 ? 0 : mcellR->deg[d2R];
    
    // Distances from the four corners of the R face to the origin of
    // the L face.
    real dist[4];

    real xphyR[3];

    // corner 0
    icR[d0R] = 0;
    ixR[d0R] = 0;
    icR[d1R] = 0;
    ixR[d1R] = 0;
    // printf("icR: %d %d %d\t", icR[0], icR[1], icR[2]);
    // printf("ixR: %d %d %d\n", ixR[0], ixR[1], ixR[2]);
    
    icix_to_xphy(mcellR, icR, ixR, xphyR);
    dist[0] = DistPeriodic(xphyL0, xphyR, f->macromesh.period, 2 * tol);
    
    // corner 1
    icR[d0R] = 0;
    ixR[d0R] = 0;
    icR[d1R] = mcellR->raf[d1R] - 1;
    ixR[d1R] = mcellR->deg[d1R];
    // printf("icR: %d %d %d\t", icR[0], icR[1], icR[2]);
    // printf("ixR: %d %d %d\n", ixR[0], ixR[1], ixR[2]);
    icix_to_xphy(mcellR, icR, ixR, xphyR);
    dist[1] = DistPeriodic(xphyL0, xphyR, f->macromesh.period, 2 * tol);
    
    // corner 2 
    icR[d0R] = mcellR->raf[d0R] - 1;
    ixR[d0R] = mcellR->deg[d0R];
    icR[d1R] = mcellR->raf[d1R] - 1;
    ixR[d1R] = mcellR->deg[d1R];
    // printf("icR: %d %d %d\t", icR[0], icR[1], icR[2]);
    // printf("ixR: %d %d %d\n", ixR[0], ixR[1], ixR[2]);
    icix_to_xphy(mcellR, icR, ixR, xphyR);
    dist[2] = DistPeriodic(xphyL0, xphyR, f->macromesh.period, 2 * tol);

    // corner 3
    icR[d0R] = mcellR->raf[d0R] - 1;
    ixR[d0R] = mcellR->deg[d0R];
    icR[d1R] = 0;
    ixR[d1R] = 0;
    //printf("icR: %d %d %d\t", icR[0], icR[1], icR[2]);
    //printf("ixR: %d %d %d\n", ixR[0], ixR[1], ixR[2]);
    icix_to_xphy(mcellR, icR, ixR, xphyR);
    dist[3] = DistPeriodic(xphyL0, xphyR, f->macromesh.period, 2 * tol);

    /*
    printf("d0: %f, \td1: %f, \td2: %f, \td3: %f\n",
	   dist[0], dist[1], dist[2], dist[3]);
    printf("locfa L/R: %d %d\n", mface->locfaL, mface->locfaR);
    printf("mcell L/R: %d %d\n", mface->ieL, mface->ieR);
    printf("paxisL: %d %d %d %d\n", paxisL[0], paxisL[1], paxisL[2], paxisL[3]);
    printf("paxisR: %d %d %d %d\n", paxisR[0], paxisR[1], paxisR[2], paxisR[3]);
    printf("rafL: %d %d %d \n", mcellL->raf[0], mcellL->raf[1], mcellL->raf[2]);
    printf("degL: %d %d %d \n", mcellL->deg[0], mcellL->deg[1], mcellL->deg[2]);
    printf("rafR: %d %d %d \n", mcellR->raf[0], mcellR->raf[1], mcellR->raf[2]);
    printf("degR: %d %d %d \n", mcellR->deg[0], mcellR->deg[1], mcellR->deg[2]);
    */

    real mindist = MIN(dist[0], dist[1]);
    mindist = MIN(mindist, dist[2]);
    mindist = MIN(mindist, dist[3]);
      
    // verify that the corner points actually find the closest point.
    real fmindist = FLT_MAX;
    for(int ipgfR = 0; ipgfR < npgfR; ++ipgfR) {
      real xphyR[3];
      ipgf_to_xphy(mcellR, mface->locfaR, ipgfR, xphyR);
      real d = DistPeriodic(xphyL0, xphyR, f->macromesh.period, 2 * tol);
      fmindist = MIN(d, fmindist);
    }
    assert(ABS(mindist -fmindist) < 1e-8);

    int dimension = 3;
    int nzero = 0;

    for(int i = 0; i < 4; ++i) {
      if(ABS(dist[i] - fmindist) < tol)
	nzero++;
    }

    mface->Rcorner = -1;
    
    int ndimensions;
    switch(nzero) {
      case 1:
	{
	  ndimensions = 3;
	  assert(!f->macromesh.is2d);
	  assert(!f->macromesh.is1d);
	  
	  // Find the unique corner which is close.  This is enough
	  // information to determine the loop variables.
	  for(int i = 0; i < 4; ++i) {
	    if(ABS(dist[i] - fmindist) < tol) {
	      mface->Rcorner = i;
	    }
	  }
	}
	break;
      case 2:
	{
	  /*
	  printf("npgdirL: %d, %d\n", mcellL->npgdir[d0L], mcellL->npgdir[d1L]);
	  printf("npgdirR: %d, %d\n", mcellR->npgdir[d0R], mcellR->npgdir[d1R]);
	  */
	  mface->Rcorner = -1;

	  if(mcellL->npgdir[d0L] != 1) {
	    if(ABS(dist[0] - fmindist) < tol && ABS(dist[1] - fmindist) < tol)
	      mface->Rcorner = 1;
	    if(ABS(dist[0] - fmindist) < tol && ABS(dist[3] - fmindist) < tol)
	      mface->Rcorner = 0;
	    if(ABS(dist[1] - fmindist) < tol && ABS(dist[2] - fmindist) < tol)
	      mface->Rcorner = 2;
	    if(ABS(dist[2] - fmindist) < tol && ABS(dist[3] - fmindist) < tol)
	      mface->Rcorner = 3;
	  } else {
	    if(ABS(dist[0] - fmindist) < tol && ABS(dist[1] - fmindist) < tol)
	      mface->Rcorner = 0;
	    if(ABS(dist[0] - fmindist) < tol && ABS(dist[3] - fmindist) < tol)
	      mface->Rcorner = 3;
	    if(ABS(dist[1] - fmindist) < tol && ABS(dist[2] - fmindist) < tol)
	      mface->Rcorner = 1;
	    if(ABS(dist[2] - fmindist) < tol && ABS(dist[3] - fmindist) < tol)
	      mface->Rcorner = 2;
	  }
	  assert(mface->Rcorner != -1);
	  ndimensions = 2;
	  assert(f->macromesh.is2d);
	}
	break;
      case 4:
	ndimensions = 1;
	//assert(f->macromesh.is1d);
	// Not much to do here: there is no loop over the interace at
	// all, so it doesn't matter which corner we choose.
	mface->Rcorner = 0;
	break;
      default:
	assert(false);
    }

    //printf("mface->Rcorner: %d\n", mface->Rcorner);

    assert(mface->Rcorner != -1); 
    
    test_MacroFace_orientation(f, mface);
  }
}

void init_field_macrofaces(field *f)
{
  f->mface = calloc(f->macromesh.nbfaces, sizeof(MacroFace));
  for(int ifa = 0; ifa < f->macromesh.nbfaces; ++ifa) {
    MacroFace *mface = f->mface + ifa;
    mface->ifa = ifa;
    
    int *f2eifa = f->macromesh.face2elem + 4 * ifa;
    mface->ieL =    f2eifa[0];
    mface->locfaL = f2eifa[1];
    mface->ieR =    f2eifa[2];
    mface->locfaR = f2eifa[3];

    MacroCell *mcellL = f->mcell + mface->ieL;
    
    mface->npgf = NPGF(mcellL->raf, mcellL->deg, mface->locfaL);
  }
  
  init_field_macrointerfaces(f);
}

void setMacroCellmass(MacroCell *mcell, field *f)
{
  for(int ipg = 0; ipg < mcell->npg; ipg++) {
    real xpgref[3];
    real wpg;
    ref_pg_vol(mcell->raf, mcell->deg, ipg, xpgref, &wpg, NULL);

    real xphy[3];
    real dtau[3][3];
    real codtau[3][3];
    Ref2Phy(mcell->physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    
    real det = dot_product(dtau[0], codtau[0]);
    mcell->mass[ipg] = wpg * det;
  }
}

void init_field_macrocells(field *f)
{
  // Allocate and set MacroCells
  f->mcell = calloc(f->macromesh.nbelems, sizeof(MacroCell));
  int wcount = 0;
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    MacroCell *mcell = f->mcell + ie; 

    mcell->ie = ie;

    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      mcell->physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      mcell->physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      mcell->physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    mcell->raf[0] = f->interp_param[4];
    mcell->raf[1] = f->interp_param[5];
    mcell->raf[2] = f->interp_param[6];
    
    mcell->deg[0] = f->interp_param[1];
    mcell->deg[1] = f->interp_param[2];
    mcell->deg[2] = f->interp_param[3];

    mcell->nsubcell = mcell->raf[0] * mcell->raf[1]  * mcell->raf[2];

    mcell->npgsubcell
      = (mcell->deg[0] + 1) * (mcell->deg[1] + 1) * (mcell->deg[2] + 1);

    mcell->dnpg[0] = mcell->deg[0] + 1;
    mcell->dnpg[1] = mcell->deg[1] + 1;
    mcell->dnpg[2] = mcell->deg[2] + 1;
        
    mcell->npgdir[0] = mcell->raf[0] * (mcell->deg[0] + 1);
    mcell->npgdir[1] = mcell->raf[1] * (mcell->deg[1] + 1);
    mcell->npgdir[2] = mcell->raf[2] * (mcell->deg[2] + 1);
    
    mcell->npg = mcell->nsubcell * mcell->npgsubcell;

    assert(mcell->npgdir[0] * mcell->npgdir[1] *mcell->npgdir[2] == mcell->npg);

    mcell->nreal = f->model.m * mcell->npg;

    mcell->nrealsubcell = f->model.m * mcell->npgsubcell;

    mcell->woffset = wcount;

    mcell->mass = calloc(mcell->npg, sizeof(real));    
    setMacroCellmass(mcell, f);
      
    wcount += mcell->nreal;
  }
}

void Initfield(field *f)
{
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};

  // a copy for avoiding too much "->"
  for(int ip = 0; ip < 8; ip++)
    f->interp_param[ip] = f->interp.interp_param[ip];

  int raf[3] = {f->interp_param[4], f->interp_param[5], f->interp_param[6]};
  int deg[3] = {f->interp_param[1], f->interp_param[2], f->interp_param[3]};
  f->wsize = f->model.m * f->macromesh.nbelems * NPG(raf, deg);

  double g_memsize = f->wsize * sizeof(real) * 1e-9;
  if(sizeof(real) == sizeof(double))
    printf("Allocating %d doubles per array (%f GB).\n", f->wsize, g_memsize);
  else
    printf("Allocating %d floats per array (%f GB)\n", f->wsize, g_memsize);

  f->wn = calloc(f->wsize, sizeof(real));
  assert(f->wn);
  f->dtwn = calloc(f->wsize, sizeof(real));
  assert(f->dtwn);

  f->Diagnostics = NULL;
  f->pre_dtfield = NULL;
  f->post_dtfield = NULL;
  f->update_after_rk = NULL;
  f->model.Source = NULL;
  f->pic = NULL;

  // TODO: move this to the integrator code
  f->tnow = 0;
  f->itermax = 0;
  f->iter_time = 0;
  f->nb_diags = 0;

  //  printf("hmin=%f\n", f->hmin);

  init_field_macrocells(f);
  
  init_field_macrofaces(f);
  
  // Compute cfl parameter min_i vol_i/surf_i
  f->hmin = min_grid_spacing(f);

  init_data(f);
  
#ifdef _WITH_OPENCL
  if(!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable; OpenCL initialization disabled.\n");
  } else {
    init_field_cl(f);
  }
#endif
  
  printf("field init done\n");
}

// This is the destructor for a field
void free_field(field *f) 
{

#ifdef _WITH_OPENCL
  //  cl_int status;


#endif

  // FIXME: free mcells and mface contents.
  free(f->mcell);
  free(f->mface);
}

// Display the field on screen
void Displayfield(field *f)
{
  printf("Display field...\n");

  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    printf("elem %d\n", ie);
    MacroCell *mcell = f->mcell + ie;
    
    for(int ipg = 0; ipg < mcell->npg; ipg++) {
      real xref[3], wpg;
      int *raf = f->interp_param + 4;
      int *deg = f->interp_param + 1;
      ref_pg_vol(raf, deg, ipg, xref, &wpg, NULL);

      printf("Gauss point %d %f %f %f \n", ipg, xref[0], xref[1], xref[2]);
      printf("dtw= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ipg, iv) + mcell->woffset;
	printf("%f ", f->dtwn[imem]);
      }
      printf("\n");
      printf("w= ");
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ipg, iv) + mcell->woffset;
	printf("%f ", f->wn[imem]);
      }
      printf("\n");
    }
  }
}

// Save the results in a text file in order plot it with Gnuplot
void Gnuplot(field *f, int dir, real fixval, char *filename)
{
  FILE *gmshfile;
  gmshfile = fopen(filename, "w" );

  printf("Save for Gnuplot...\n");
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {

    MacroCell *mcell = f->mcell + ie;
    
    for(int ipg = 0; ipg < mcell->npg; ipg++) {

      real xref[3], xphy[3], wpg;
      real dtau[3][3];
      ref_pg_vol(mcell->raf, mcell->deg, ipg, xref, &wpg, NULL);

      Ref2Phy(mcell->physnode,
	      xref,
	      0, -1, // dphiref, ifa
	      xphy, dtau,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      if(xphy[dir] > -(fixval + 0.0001) && xphy[dir] < (fixval + 0.00001)){

	fprintf(gmshfile, "%f ",xphy[1-dir]);

	for(int iv = 0; iv < f->model.m; iv++) {
	  int imem = f->varindex(f->interp_param, ipg, iv) + mcell->woffset;
	  fprintf(gmshfile, "%f ",f->wn[imem]);
	}
	fprintf(gmshfile, "\n");

      }
    }
  }
  fclose(gmshfile);
}

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

  //int *elem2nodes = f->macromesh.elem2node;
  //real *node = f->macromesh.node;
  
  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int *raf = f->interp_param + 4;
  int *deg = f->interp_param + 1;

  // Refinement size in each direction
  real hh[3] = {1.0 / raf[0], 1.0 / raf[1], 1.0 / raf[2]};

  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(real));

  int nb_plotnodes = f->macromesh.nbelems * raf[0] * raf[1] * raf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);

  real *value = malloc(nb_plotnodes * sizeof(real));
  assert(value);
  int nodecount = 0;

  real tolerance;
  if(sizeof(real) == sizeof(double))
    tolerance = 1e-10;
  else
    tolerance = 1e-6;
  
  // Nodes
  int npgv = NPG(raf, deg);
  for(int i = 0; i < f->macromesh.nbelems; i++) {
    MacroCell *mcell = f->mcell + i;

    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < raf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < raf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < raf[2]; icL[2]++) {

	  for(int ino = 0; ino < 64; ino++) {
	    real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
			   hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
			   hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
	    
	    for(int ii = 0; ii < 3; ii++) {
	      assert(Xr[ii] < 1 +  1e-10);
	      assert(Xr[ii] > -1e-10);
	    }

	    real Xphy[3];
	    Ref2Phy(mcell->physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL,
		    NULL);

	    real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};

	    value[nodecount] = 0;
	    real testpsi = 0;
	    for(int ib = 0; ib < npgv; ib++) {
	      real psi;
	      psi_ref_subcell(f->interp_param + 1, icL, ib, Xr, &psi, NULL);
	      testpsi += psi;
	      int vi = f->varindex(f->interp_param, ib, typplot)
		+ mcell->woffset;
	      value[nodecount] += psi * f->wn[vi];
	    }
	    assert(fabs(testpsi-1) < tolerance);

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
	  f->macromesh.nbelems * raf[0] * raf[1] * raf[2]);

  int elm_type = 92;
  int num_tags = 0;

  // fwrite((char*) &elm_type, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_elm_follow, sizeof(int), 1, gmshfile);
  // fwrite((char*) &num_tags, sizeof(int), 1, gmshfile);

  for(int i = 0; i < f->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < raf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < raf[1]; icL[1]++) {
	for(icL[2] = 0; icL[2] < raf[2]; icL[2]++) {
	  // Get the subcell id
	  int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);

	  // Global subcell id
	  int numelem = ncL + i * raf[0] * raf[1] * raf[2] + 1;

	  //fwrite((char*) &numelem, sizeof(int), 1, gmshfile);
	  fprintf(gmshfile, "%d ", numelem);
	  fprintf(gmshfile, "%d ", elm_type);
	  fprintf(gmshfile, "%d ", num_tags);

	  for(int ii = 0; ii < 64; ii++) {
	    int numnoe = 64 * (i * raf[0] * raf[1] * raf[2] + ncL) + ii  + 1;
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
void DGSubCellInterface(MacroCell *mcell, field *f, real *wmc, real *dtwmc) 
{
  const int nraf[3] = {mcell->raf[0], mcell->raf[1], mcell->raf[2]};
  const int deg[3] = {mcell->deg[0], mcell->deg[1], mcell->deg[2]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
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

	    // non-coalescent access
	    //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;

	    // coalescent access (pairs of xyz are in the correct order).
	    const int altdim1[3] = {1, 0, 0};
	    const int altdim2[3] = {2, 2, 1};
	    int dim1 = altdim1[dim0];
	    int dim2 = altdim2[dim0];

	    int iL[3];
	    iL[dim0] = deg[dim0];
	    // Now loop on the left glops of the subface
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
		  int *raf = f->interp_param + 4;
		  int *deg = f->interp_param + 1;
		  ref_pg_vol(raf, deg, ipgL, xref, &wpg3, NULL);
		  // mapping from the ref glop to the physical glop
		  real dtau[3][3], codtau[3][3];
		  Ref2Phy(mcell->physnode,
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

		  real h1h2 = 1.0 / nraf[dim1] / nraf[dim2];
		  vnds[0] = codtau[0][dim0] * h1h2;
		  vnds[1] = codtau[1][dim0] * h1h2;
		  vnds[2] = codtau[2][dim0] * h1h2;
		}

		// numerical flux from the left and right state and
		// normal vector
		real wL[m], wR[m], flux[m];
		for(int iv = 0; iv < m; iv++) {
		  // TO DO change the varindex signature
		  int imemL = f->varindex(f->interp_param, ipgL, iv);
		  int imemR = f->varindex(f->interp_param, ipgR, iv);
		  // end TO DO
		  wL[iv] = wmc[imemL];
		  wR[iv] = wmc[imemR];
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
		  int imemL = f->varindex(f->interp_param, ipgL, iv);
		  int imemR = f->varindex(f->interp_param, ipgR, iv);
		  // end TO DO
		  dtwmc[imemL] -= flux[iv] * wpg;
		  dtwmc[imemR] += flux[iv] * wpg;
		}

	      }  // face yhat loop
	    } // face xhat loop
	  } // endif internal face
	} // dim loop
      } // subcell icl2 loop
    } // subcell icl1 loop
  } // subcell icl0 loop
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms
void DGMacroCellInterfaceSlow(MacroCell *mcell, field *f, real *w, real *dtw)
{
  // Local copy of the interpretation parameters
  int iparam[8];
  for(int ip = 0; ip < 8; ip++)
    iparam[ip] = f->interp_param[ip];

  int *raf = mcell->raf;
  int *deg = mcell->deg;
  
  int ie = mcell->ie;

  // loop on the 6 faces
  // or four faces for 2d computations
  int nbfa = 6;
  if (f->macromesh.is2d) nbfa = 4;
  if (f->macromesh.is1d) nbfa = 2;
  for(int nifa = 0; nifa < nbfa; nifa++) {
    // get the right elem or the boundary id
    int ifa = f->macromesh.is1d ?  2 * nifa + 1 : nifa;
    int ieR = f->macromesh.elem2elem[6*ie+ifa];

    // loop on the glops (numerical integration)
    // of the face ifa
    for(int ipgf = 0; ipgf < NPGF(raf, deg, ifa); ipgf++) {
      real xpgref[3], xpgref_in[3], wpg;
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipg = ref_pg_face(raf, deg, ifa, ipgf, xpgref, &wpg, xpgref_in);

      // get the left value of w at the gauss point
      real wL[f->model.m], wR[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(iparam, ipg, iv) + mcell->woffset;
	wL[iv] = f->wn[imem];
      }

      // The basis functions is also the gauss point index
      int ib = ipg;

      real dtau[3][3];
      real codtau[3][3];
      real xpg[3];
      real vnds[3]; // normal vector at gauss point ipg
      Ref2Phy(mcell->physnode,
	      xpgref,
	      NULL, ifa, // dpsiref, ifa
	      xpg, dtau,
	      codtau, NULL, vnds); // codtau, dpsi, vnds
      real flux[f->model.m];

      if (ieR >=0) {  // the right element exists
	MacroCell *mcellR = f->mcell + ieR;

	// find the corresponding point in the right elem
	real xpg_in[3];
	Ref2Phy(mcell->physnode,
		xpgref_in,
		NULL, ifa, // dpsiref, ifa
		xpg_in, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
	real xref[3];
	PeriodicCorrection(xpg_in, f->macromesh.period);
	Phy2Ref(mcellR->physnode, xpg_in, xref);
	int ipgR = ref_ipg(raf, deg, xref);
	real xpgR[3], xrefR[3], wpgR;
	int *raf = iparam + 4;
	int *deg = iparam + 1;
	ref_pg_vol(raf, deg, ipgR, xrefR, &wpgR, NULL);
	Ref2Phy(mcellR->physnode,
		xrefR,
		NULL, -1, // dphiref, ifa
		xpgR, NULL,
		NULL, NULL, NULL); // codtau, dphi, vnds
	  
	for(int iv = 0; iv < f->model.m; iv++) {
	  int imem = f->varindex(iparam, ipgR, iv) + mcellR->woffset;
	  wR[iv] = f->wn[imem];
	}
	// int_dL F(wL, wR, grad phi_ib )
	f->model.NumFlux(wL, wR, vnds, flux);

      } else { //the right element does not exist
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
      }

      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(iparam, ib, iv) + mcell->woffset;
	f->dtwn[imem] -= flux[iv] * wpg;
      }

    }

  }
}

// Compute the Discontinuous Galerkin macrocell boundary terms.
void DGMacroCellBoundary(MacroFace *mface, field *f, real *wmc, real *dtwmc) 
{
  MacroMesh *msh = &f->macromesh;
  const unsigned int m = f->model.m;

  // Assembly of the surface terms loop on the macrocells faces

  int ifa = mface->ifa;
  int ieL = mface->ieL;

  MacroCell *mcellL = f->mcell + ieL;
  int locfaL = mface->locfaL;
  
  // Loop over the points on a single macro cell interface.
  const int npgf = mface->npgf;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ipgfL = 0; ipgfL < npgf; ipgfL++) {

    real xpgref[3];
    real wpg;
    int ipgL = ref_pg_face(mcellL->raf, mcellL->deg,
			   locfaL, ipgfL, xpgref, &wpg, NULL);

    // Normal vector at gauss point ipgL
    real vnds[3];
    real xpg[3];
    {
      real dtau[3][3];
      real codtau[3][3];
      Ref2Phy(mcellL->physnode,
	      xpgref,
	      NULL, locfaL, // dpsiref, ifa
	      xpg, dtau,
	      codtau, NULL, vnds); // codtau, dpsi, vnds
    }

    real wL[m];
    for(int iv = 0; iv < m; iv++) {
      int imemL = f->varindex(f->interp_param, ipgL, iv);
      wL[iv] = wmc[imemL];
    }
    
    real flux[m];
    f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
    
    for(int iv = 0; iv < m; iv++) {
      int imemL = f->varindex(f->interp_param, ipgL, iv);
      dtwmc[imemL] -= flux[iv] * wpg;
    }
  }
}

// Compute the Discontinuous Galerkin inter-macrocells flux.
void DGMacroCellInterface(MacroFace *mface, field *f,
			  real *wmcL, real *wmcR, real *dtwmcL, real *dtwmcR) 
{
  const unsigned int m = f->model.m;

  // Assembly of the surface terms loop on the macrocells faces

  int ieL = mface->ieL;
  MacroCell *mcellL = f->mcell + ieL;
  int locfaL = mface->locfaL;

  int ieR = mface->ieR;
  MacroCell *mcellR = f->mcell + ieR;

  // Loop over the points on a single macro cell interface.
  const int npgf = mface->npgf;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int ipgfL = 0; ipgfL < npgf; ipgfL++) {

    real xpgref[3];
    real xpgref_in[3];
    real wpg;
    // Get the coordinates of the Gauss point and coordinates of a
    // point slightly inside the opposite element in xref_in
    int ipgL = ref_pg_face(mcellL->raf, mcellL->deg,
			   locfaL, ipgfL, xpgref, &wpg, xpgref_in);

    int ipgR;
    {
      real xpg_in[3];
      Ref2Phy(mcellL->physnode, xpgref_in, NULL, -1, xpg_in, NULL,
	      NULL, NULL, NULL);

      PeriodicCorrection(xpg_in, f->macromesh.period);

      real xrefR[3];
      Phy2Ref(mcellR->physnode, xpg_in, xrefR);

      ipgR = ref_ipg(mcellR->raf, mcellR->deg, xrefR);
    }
    
    real wL[m];
    real wR[m];
    for(int iv = 0; iv < m; iv++) {
      int imemL = f->varindex(f->interp_param, ipgL, iv);
      wL[iv] = wmcL[imemL];
      int imemR = f->varindex(f->interp_param, ipgR, iv);
      wR[iv] = wmcR[imemR];
    }

    // Normal vector at gauss point ipgL
    real vnds[3];
    {
      real dtau[3][3];
      real codtau[3][3];
      Ref2Phy(mcellL->physnode, xpgref, NULL, locfaL, NULL, dtau,
	      codtau, NULL, vnds);
    }
    
    // int_dL F(wL, wR, grad phi_ib)
    real flux[m];
    f->model.NumFlux(wL, wR, vnds, flux);
    
    for(int iv = 0; iv < m; iv++) {
      int imemL = f->varindex(f->interp_param, ipgL, iv);
      dtwmcL[imemL] -= flux[iv] * wpg;

      int imemR = f->varindex(f->interp_param, ipgR, iv);
      dtwmcR[imemR] += flux[iv] * wpg;
    }

  }
}

// Apply division by the mass matrix
void DGMass(MacroCell *mcell, field *f, real *dtwmc) 
{
  const int m = f->model.m;
  for(int ipg = 0; ipg < mcell->npg; ipg++) {
    const real overmass = 1.0 / mcell->mass[ipg];
    for(int iv = 0; iv < m; iv++) {
      int imem = f->varindex(f->interp_param, ipg, iv);
      dtwmc[imem] *= overmass;
    }
  }
}

// Apply the source term
void DGSource(MacroCell *mcell, field *f, real tnow, real *wmc, real *dtwmc) 
{
  if (f->model.Source == NULL) {
    return;
  }
  
  const int m = f->model.m;

  for(int ipg = 0; ipg < mcell->npg; ipg++) {
    const real mass = mcell->mass[ipg];

    real wpg;
    real xpgref[3], xphy[3];
    int *raf = f->interp_param + 4;
    int *deg = f->interp_param + 1;
    ref_pg_vol(raf, deg, ipg, xpgref, &wpg, NULL);
    Ref2Phy(mcell->physnode, // phys. nodes
	    xpgref, // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, NULL, // xphy, dtau
	    NULL, NULL, NULL); // codtau, dpsi, vnds

    real wL[m];
    for(int iv = 0; iv < m; ++iv){
      int imem = f->varindex(f->interp_param, ipg, iv);
      wL[iv] = wmc[imem];
    }

    real source[m];
    f->model.Source(xphy, tnow, wL, source, m);
      
    for(int iv = 0; iv < m; ++iv) {
      int imem = f->varindex(f->interp_param, ipg, iv);
      dtwmc[imem] += source[iv] * mass;
    }
  }
}

// Compute the Discontinuous Galerkin volume terms, fast version
void DGVolume(MacroCell *mcell, field *f, real *wmc, real *dtwmc) 
{
  const int m = f->model.m;
  const int deg[3] = {mcell->deg[0],mcell->deg[1],mcell->deg[2]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int raf[3] = {mcell->raf[0], mcell->raf[1], mcell->raf[2]};

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
  for(int icL0 = 0; icL0 < raf[0]; icL0++) {
    for(int icL1 = 0; icL1 < raf[1]; icL1++) {
      for(int icL2 = 0; icL2 < raf[2]; icL2++) {

	int icL[3] = {icL0, icL1, icL2};
	// get the L subcell id
	int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
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

	  int *raf = f->interp_param + 4;
	  int *deg = f->interp_param + 1;
	  ref_pg_vol(raf, deg, offsetL + p, xref, &tomega, NULL);
	  xref0[p] = xref[0];
	  xref1[p] = xref[1];
	  xref2[p] = xref[2];
	  omega[p] = tomega;

	  for(int im = 0; im < m; ++im) {
	    imems[pos++] = f->varindex(f_interp_param, offsetL + p, im);
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
		  wL[iv] = wmc[imems[m * (ipgL - offsetL) + iv]];
		}
		int q[3] = {p[0], p[1], p[2]};
		// loop on the direction dim0 on the "cross"
		for(int iq = 0; iq < npg[dim0]; iq++) {
		  q[dim0] = (p[dim0] + iq) % npg[dim0];
		  real dphiref[3] = {0, 0, 0};
		  // compute grad phi_q at glop p
		  dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) 
		    * raf[dim0];

		  real xrefL[3] = {xref0[ipgL - offsetL],
				   xref1[ipgL - offsetL],
				   xref2[ipgL - offsetL]};
		  real wpgL = omega[ipgL - offsetL];
		  /* real xrefL[3], wpgL; */
		  /* ref_pg_vol(f->interp_param+1,ipgL,xrefL, &wpgL, NULL); */

		  // mapping from the ref glop to the physical glop
		  real dtau[3][3], codtau[3][3], dphiL[3];
		  Ref2Phy(mcell->physnode,
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
		    //int imemR = f->varindex(f_interp_param, ipgR, iv);
		    int temp = m * (ipgR - offsetL) + iv;  
		    //assert(imemR == imems[temp]);
		    dtwmc[imems[temp]] += flux[iv] * wpgL;
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

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field
void dtfield(field *f, real tnow, real *w, real *dtw) {
  if(f->pre_dtfield != NULL) // FIXME: rename to before dtfield
    f->pre_dtfield(f, w);

#ifdef _OPENMP
#pragma omp parallel
#endif
  for(int iw = 0; iw < f->wsize; iw++)
    dtw[iw] = 0;

  bool facealgo = true;
  if(facealgo) {
    // Macrocell interfaces
    const int ninterfaces = f->macromesh.nmacrointerfaces;
    for(int i = 0; i < ninterfaces; ++i) {
      int ifa = f->macromesh.macrointerface[i];
      MacroFace *mface = f->mface + ifa;

      int ieL = mface->ieL;
      MacroCell *mcellL = f->mcell + ieL;
      real *wL = w + mcellL->woffset;
      real *dtwL = dtw + mcellL->woffset;

      int ieR = mface->ieR;
      MacroCell *mcellR = f->mcell + ieR;
      real *wR = w + mcellR->woffset;
      real *dtwR = dtw + mcellR->woffset;

      DGMacroCellInterface(f->mface + ifa, f, wL, wR, dtwL, dtwR);
    }

    // Macrocell boundaries
    const int nboundaryfaces = f->macromesh.nboundaryfaces;
    for(int i = 0; i < nboundaryfaces; ++i) {
      int ifa = f->macromesh.boundaryface[i];
      MacroFace *mface = f->mface + ifa;
      int ie = mface->ieL;
      MacroCell *mcell = f->mcell + ie;
      real *wmc = w + mcell->woffset;
      real *dtwmc = dtw + mcell->woffset;
      DGMacroCellBoundary(mface, f, wmc, dtwmc);
    }
  }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    MacroCell *mcell = f->mcell + ie;

    if(!facealgo)
      DGMacroCellInterfaceSlow(mcell, f, w, dtw);

    real *wmc = w + mcell->woffset;
    real *dtwmc = dtw + mcell->woffset;
    DGSubCellInterface(mcell, f, wmc, dtwmc);
    DGVolume(mcell, f, wmc, dtwmc);
    DGSource(mcell, f, tnow, wmc, dtwmc);
    DGMass(mcell, f, dtwmc);
  }

  if(f->post_dtfield != NULL) // FIXME: rename to after dtfield
    f->post_dtfield(f, w);
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
  int *raf = f->interp_param + 4;
  int *deg = f->interp_param + 1;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(raf, deg);
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

    dtfield(f, f->tnow, f->wn, f->dtwn);
    RK_out(wnp1, f->wn, f->dtwn, 0.5 * dt, sizew);

    f->tnow += 0.5 * dt;

    dtfield(f, f->tnow, wnp1, f->dtwn);
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
  int size_diags;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int *raf = f->interp_param + 4;
  int *deg = f->interp_param + 1;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(raf, deg);
  int iter = 0;

  // Allocate memory for RK time-stepping
  real *l1, *l2, *l3;
  l1 = calloc(sizew, sizeof(real));
  l2 = calloc(sizew, sizeof(real));
  l3 = calloc(sizew, sizeof(real));
  
  size_diags = f->nb_diags * f->itermax;
  f->iter_time = iter;
  
  if(f->nb_diags != 0)
    f->Diagnostics = malloc(size_diags * sizeof(real));
  
  while(f->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    // l_1 = w_n + 0.5dt * S(w_n, t_0)
    dtfield(f, f->tnow, f->wn, f->dtwn);
    RK_out(l1, f->wn, f->dtwn, 0.5 * dt, sizew);

    f->tnow += 0.5 * dt;

    // l_2 = w_n + 0.5dt * S(l_1, t_0 + 0.5 * dt)
    dtfield(f, f->tnow, l1, f->dtwn);
    RK_out(l2, f->wn, f->dtwn, 0.5 * dt, sizew);

    // l_3 = w_n + dt * S(l_2, t_0 + 0.5 * dt)
    dtfield(f, f->tnow, l2, f->dtwn);
    RK_out(l3, f->wn, f->dtwn, dt, sizew);

    f->tnow += 0.5 * dt;

    // Compute S(l_3, t_0 + dt)
    dtfield(f, f->tnow, l3, f->dtwn);
    RK4_final_inplace(f->wn, l1, l2, l3, f->dtwn, dt, sizew);

    
    if(f->update_after_rk != NULL)
      f->update_after_rk(f, f->wn);
    
    iter++;
    f->iter_time=iter;
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
    MacroCell *mcell = f->mcell + ie;

        // Loop on the glops (for numerical integration)
    const int npg = NPG(mcell->raf, mcell->deg);
    for(int ipg = 0; ipg < npg; ipg++) {
      real w[f->model.m];
      for(int iv = 0; iv < f->model.m; iv++) {
	int imem = f->varindex(f->interp_param, ipg, iv) + mcell->woffset;
	w[iv] = f->wn[imem];
      }

      real wex[f->model.m];
      real wpg, det;
      // Compute wpg, det, and the exact solution
      { 
	real xphy[3], xpgref[3];
	real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	int *raf = f->interp_param + 4;
	int *deg = f->interp_param + 1;
	ref_pg_vol(raf, deg, ipg, xpgref, &wpg, NULL);
	Ref2Phy(mcell->physnode, // phys. nodes
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


// FIXME: pass MacroCell*
void InterpField(field *f, int ie, real *xref, real *w){

  MacroCell *mcell = f->mcell + ie;

  const int nraf[3] = {mcell->raf[0], mcell->raf[1], mcell->raf[2]};
  //const int deg[3] = {mcell->deg[0], mcell->deg[1], mcell->deg[2]};

  for(int iv = 0; iv < f->model.m; iv++)
    w[iv] = 0;

  int is[3];

  for(int ii = 0; ii < 3; ii++){
    is[ii] = xref[ii] * nraf[ii];
    assert(is[ii] < nraf[ii] && is[ii]>= 0);
  }
  
  int npgv = NPG(mcell->raf, mcell->deg);
  // TODO: loop only on non zero basis function
  for(int ib = 0; ib < npgv; ib++) { 
    real psi;
    psi_ref_subcell(f->interp_param + 1, is, ib, xref, &psi, NULL);
    
    for(int iv=0;iv<f->model.m;iv++){
      int imem = f->varindex(f->interp_param, ib, iv) + mcell->woffset;
      w[iv] += psi * f->wn[imem];
    }
  }
}


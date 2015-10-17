#include "field.h"
#include "field_cl.h"
#include "clinfo.h"
#include "clutils.h"
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#ifdef _WITH_OPENCL
void CopyfieldtoCPU(field *f)
{
  cl_int status;

  // Ensure that all the buffers are mapped
  status = clFinish(f->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  void *chkptr;
  chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
			      f->dtwn_cl, // buffer to copy from
			      CL_TRUE, // block until the buffer is available
			      CL_MAP_READ, // we just want to see the results
			      0, // offset
			      f->wsize * sizeof(real), // buffersize
			      0, NULL, NULL, // events management
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == f->dtwn);

  chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
			      f->wn_cl, // buffer to copy from
			      CL_TRUE, // block until the buffer is available
			      CL_MAP_READ, // we just want to see the results
			      0, // offset
			      f->wsize * sizeof(real), // buffersize
			      0, NULL, NULL, // events management
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == f->wn);

  status = clFinish(f->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_source_CL(field *f, char *sourcename_cl) 
{
#ifdef _WITH_OPENCL
  f->use_source_cl = true;
  int namelength = strlen(sourcename_cl);
  if(namelength > 0) {
    f->sourcename_cl = calloc(sizeof(char), namelength + 1);
    strcat(f->sourcename_cl, sourcename_cl);
  } else {
    printf("The source name %s is empty!\n", sourcename_cl);
    assert(false);
  }
  printf("%s\n", f->sourcename_cl);
#endif
}

void init_DGBoundary_CL(field *f, 
			int ieL, int locfaL,
			cl_mem *wn_cl,
			size_t cachesize)
{
  cl_int status;
  cl_kernel kernel = f->dgboundary;

  MacroCell *mcell = f->mcell + ieL;
  
  unsigned int argnum = 0;
  
  // __constant int *param,        // 0: interp param
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // real tnow,                  // 1: current time
  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(real),
  			  &f->tnow);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(int),
  			  &mcell->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // int locfaL,                   // 3: left face index
  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(int),
  			  &locfaL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __constant real *physnodeL, // 4: left physnode
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &mcell->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *wn,          // 5: field
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *dtwn,        // 6: time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __local real *cache         // 7: local mem
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGBoundary_CL(MacroFace *mface, field *f, cl_mem *wn_cl,
		   cl_uint nwait, cl_event *wait, cl_event *done) 
{
  int *param = f->interp_param;
  cl_int status;

  int ifa = mface->ifa;
  int ieL =    f->macromesh.face2elem[4 * ifa];
  int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    
  size_t numworkitems = NPGF(f->interp_param + 1, locfaL);
  size_t cachesize = 1; // TODO make use of cache
  init_DGBoundary_CL(f, ieL, locfaL, wn_cl, cachesize);
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->dgboundary,
				  1, // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &numworkitems, // global_work_size, 
				  NULL, // size_t *local_work_size, 
				  nwait, 
				  wait, // *wait_list,
				  done); // *event
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Set the loop-dependant kernel arguments for DGMacroCellInterface_CL
void init_DGMacroCellInterface_CL(field *f, 
				  int ieL, int ieR, int locfaL, int locfaR,
				  cl_mem *wn_cl,
				  size_t cachesize)
{
  //printf("loop_init_DGMacroCellInterface_CL\n");
  cl_int status;
  cl_kernel kernel = f->dginterface;

  MacroCell *mcellL = f->mcell + ieL;
  MacroCell *mcellR = f->mcell + ieR;
  
  // Set kernel arguments
  unsigned int argnum = 0;

  status = clSetKernelArg(kernel,
			 argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &mcellL->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //int ieR,                      // 2: right macrocell
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &mcellR->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //int locfaL,                   // 3: left face index
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //int locfaR,                   // 4: right face index
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaR);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__constant real *physnodeL,
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
			  &mcellL->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__constant real *physnodeR,
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &mcellR->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real *wn,          // 7: field 
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real *dtwn,        // 8: time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__local real *cache         // 9: local mem
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGMacroCellInterface_CL(MacroFace *mface, field *f, cl_mem *wn_cl,
			     cl_uint nwait, cl_event *wait, cl_event *done) 
{
  int *param = f->interp_param;
  cl_int status;
  cl_kernel kernel = f->dginterface;
  status = clSetKernelArg(kernel,
                          7,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);


  int ifa = mface->ifa;

  int ieL =    f->macromesh.face2elem[4 * ifa];
  int locfaL = f->macromesh.face2elem[4 * ifa + 1];
  int ieR =    f->macromesh.face2elem[4 * ifa + 2];
  int locfaR = f->macromesh.face2elem[4 * ifa + 3];

  size_t numworkitems = NPGF(f->interp_param + 1, locfaL);
  if(ieR >= 0) {
    // Set the remaining loop-dependant kernel arguments
    size_t kernel_cachesize = 1;
    init_DGMacroCellInterface_CL(f, 
				 ieL, ieR, locfaL, locfaR, 
				 wn_cl, 
				 kernel_cachesize);

    status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				    kernel,
				    1, // cl_uint work_dim,
				    NULL, // global_work_offset,
				    &numworkitems, // global_work_size, 
				    NULL, // size_t *local_work_size, 
				    nwait,  // nwait, 
				    wait, // *wait_list,
				    done); // *event
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
  } else {
    // Set the event to completed status.
    if(done != NULL)
      clSetUserEventStatus(*done, CL_COMPLETE);
  }
}

// Set up kernel arguments, etc, for DGMass_CL.
void init_DGMass_CL(MacroCell *mcell, field *f)
{
  cl_int status;
  cl_kernel kernel = f->dgmass;
  int argnum = 0;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* int ie, // macrocel index */
  // Set in loop on call.

  int offset = mcell->woffset;
  status = clSetKernelArg(f->dgmass, 
			  argnum++, 
			  sizeof(int), 
			  &offset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  /* __constant real* physnode,  // macrocell nodes */
  status = clSetKernelArg(f->dgmass,           // kernel name
                          argnum++,              // arg num
                          sizeof(cl_mem),
			  &mcell->physnode_cl);     // opencl buffer
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* __global real* dtwn // time derivative */
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Apply division by the mass matrix OpenCL version
void DGMass_CL(MacroCell *mcell, field *f,
	       cl_uint nwait, cl_event *wait, cl_event *done) 
{
  //printf("DGMass_CL\n");
  int *param = f->interp_param;
  cl_int status;

  init_DGMass_CL(mcell, f);
  
  // The groupsize is the number of glops in a subcell
  size_t groupsize = (param[1] + 1) * (param[2] + 1) * (param[3] + 1);

  // The total work items number is (the number of glops in a
  // subcell) * (number of subcells)
  size_t numworkitems = param[4] * param[5] * param[6] * groupsize;

  unsigned int m = f->interp_param[0];
  unsigned int nreadsdgmass = m;
  unsigned int nmultsdgmass = 53 + 2 * m + 2601;
  f->flops_mass += numworkitems * nmultsdgmass;
  f->reads_mass += numworkitems * nreadsdgmass;
    
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->dgmass,
				  1, // cl_uint work_dim,
				  NULL, // size_t *global_work_offset, 
				  &numworkitems, // size_t *global_work_size,
				  &groupsize, // size_t *local_work_size, 
				  nwait, // cl_uint num_events_in_wait_list, 
				  wait, //*event_wait_list, 
				  done); // cl_event *event
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void init_DGFlux_CL(field *f, int ie, int dim0, cl_mem *wn_cl, 
		    size_t cachesize)
{
  //printf("DGFlux cachesize:%zu\n", cachesize);

  cl_kernel kernel = f->dgflux;
  cl_int status;
  int argnum = 0; 

  MacroCell *mcell = f->mcell + ie;
  
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel, 
			  argnum++, 
			  sizeof(int), 
			  &mcell->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // int dim0, // face direction
  status = clSetKernelArg(kernel, 
			  argnum++, 
			  sizeof(int), 
			  (void *)&dim0);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
    
  // __constant real *physnode, // macrocell nodes
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &mcell->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //  __global real *wn, // field values
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *dtwn // time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __local real* wnloc,       // 6: wn local memory
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGFlux_CL(field *f, int dim0, int ie, cl_mem *wn_cl,
	       cl_uint nwait, cl_event *wait, cl_event *done) 
{
  // Unpack the parameters
  int *param = f->interp_param;
  int m = param[0];
  int deg[3] = {param[1], param[2], param[3]};
  int nraf[3] = {param[4], param[5], param[6]};

  int dim[3];
  dim[0] = dim0;
  dim[1] = (dim[0] + 1) % 3;
  dim[2] = (dim[1] + 1) % 3;

  // Number of points per face
  int npgf = (deg[dim[1]] + 1) * (deg[dim[2]] + 1);

  // Number of faces
  int nf = (nraf[dim[0]] - 1) * nraf[dim[1]] * nraf[dim[2]];

  if(nf > 0) { // If there are faces to work on, launch the kernel
    // Set kernel args
    size_t numworkitems = nf * npgf;
    size_t groupsize = npgf;
    init_DGFlux_CL(f, ie, dim0, wn_cl, 4 * m * groupsize);

    unsigned int nreadsdgflux = 4 * m;
    unsigned int nmultsdgflux = 2601 + 92 + 28 * m;
    nmultsdgflux += 3 * m; // Using NUMFLUX = NumFlux
    f->flops_flux += numworkitems * nmultsdgflux;
    f->reads_flux += numworkitems * nreadsdgflux;
    
    // Launch the kernel
    cl_int status;
    status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				    f->dgflux,
				    1,
				    NULL,
				    &numworkitems,
				    &groupsize,
				    nwait,
				    wait,
				    done);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
  } else {
    if(done != NULL)
      clSetUserEventStatus(*done, CL_COMPLETE);
  }

}

// Set kernel argument for DGVolume_CL
void init_DGVolume_CL(MacroCell *mcell, field *f, cl_mem *wn_cl,
		      size_t cachesize)
{
  //printf("DGVolume cachesize:%zu\n", cachesize);

  cl_int status;
  int argnum = 0;
  cl_kernel kernel = f->dgvolume;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &mcell->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &mcell->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* status = clSetKernelArg(kernel, */
  /*                         argnum++, */
  /*                         sizeof(cl_real) * cachesize, */
  /*                         NULL); */
  /* if(status < CL_SUCCESS) printf("%s\n", clErrorString(status)); */
  /* assert(status >= CL_SUCCESS); */
}

// Apply division by the mass matrix OpenCL version
void DGVolume_CL(MacroCell *mcell, field *f, cl_mem *wn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done) 
{
  cl_kernel kernel = f->dgvolume;
  int *param = f->interp_param;

  cl_int status;
  int m = param[0];
  const int npgc[3] = {param[1] + 1, param[2] + 1, param[3] + 1};
  const int npg = npgc[0] * npgc[1] * npgc[2];
  size_t groupsize = npg;
  // The total work items number is the number of glops in a subcell
  // * number of subcells
  const int nraf[3] = {param[4], param[5], param[6]};
  size_t numworkitems = nraf[0] * nraf[1] * nraf[2] * groupsize;
  
  unsigned int nreadsdgvol = 2 * m; // read m from wn, write m to dtwn
  unsigned int nmultsdgvol = 2601 + (npgc[0] + npgc[1] + npgc[2]) * 6 * m;
  nmultsdgvol += (npgc[0] + npgc[1] + npgc[2]) * 54; 
  // Using NUMFLUX = NumFlux (3 * m multiplies):
  nmultsdgvol += (npgc[0] + npgc[1] + npgc[2]) * 6 * m; 
  
 

  //printf("DGVolume_CL loop: %d\n", end - start); // This is always 1!!!

  init_DGVolume_CL(mcell, f, wn_cl, 2 * groupsize * m);
  int ie = mcell->ie;
  
  f->flops_vol += numworkitems * nmultsdgvol;
  f->reads_vol += numworkitems * nreadsdgvol;

  // The groupsize is the number of glops in a subcell
  /* size_t groupsize = (param[1] + 1)* (param[2] + 1)*(param[3] + 1); */
  /* // The total work items number is the number of glops in a subcell */
  /* // * number of subcells */
  /* size_t numworkitems = param[4] * param[5] * param[6] * groupsize; */
  /* //printf("groupsize=%zd numworkitems=%zd\n", groupsize, numworkitems); */
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  kernel,
				  1,
				  NULL,
				  &numworkitems,
				  &groupsize,
				  nwait,
				  wait,
				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Set kernel argument for DGVolume_CL
void init_DGSource_CL(MacroCell *mcell, field *f,
		      real tnow, cl_mem *wn_cl, size_t cachesize)
{
  //printf("DGVolume cachesize:%zu\n", cachesize);

  cl_int status;
  int argnum = 0;
  cl_kernel kernel = f->dgsource;

  int ie = mcell->ie;
  
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &mcell->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &mcell->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(real),
  			  &tnow);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Apply division by the mass matrix OpenCL version
void DGSource_CL(MacroCell *mcell, field *f, real tnow, cl_mem *wn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done) 
{
  cl_kernel kernel = f->dgsource;
  int *param = f->interp_param;

  cl_int status;
  int m = param[0];
  size_t groupsize = (param[1] + 1) * (param[2] + 1) * (param[3] + 1);
  // The total work items number is the number of glops in a subcell
  // times the number of subcells
  size_t numworkitems = param[4] * param[5] * param[6] * groupsize;
  
  /* unsigned int nreadsdgvol = 2 * m; // read m from wn, write m to dtwn */
  /* unsigned int nmultsdgvol = 1296 + m; */
  /* nmultsdgvol += 3 * m; // Using NUMFLUX = NumFlux */
  
  /* f->nmults += numworkitems * nmultsdgvol; */
  /* f->nreads += numworkitems * nreadsdgvol; */
   
  init_DGSource_CL(mcell, f, tnow, wn_cl, 2 * groupsize * m);
  
  // The groupsize is the number of glops in a subcell
  /* size_t groupsize = (param[1] + 1)* (param[2] + 1)*(param[3] + 1); */
  /* // The total work items number is the number of glops in a subcell */
  /* // * number of subcells */
  /* size_t numworkitems = param[4] * param[5] * param[6] * groupsize; */
  /* //printf("groupsize=%zd numworkitems=%zd\n", groupsize, numworkitems); */
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  kernel,
				  1,
				  NULL,
				  &numworkitems,
				  &groupsize,
				  nwait,
				  wait,
				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_buf_to_zero_cl(cl_mem *buf, int size, field *f,
			cl_uint nwait, cl_event *wait,  cl_event *done)
{
  cl_int status;

  cl_kernel kernel = f->zero_buf;

  // associates the param buffer to the 0th kernel argument
  status = clSetKernelArg(kernel,
                          0, 
                          sizeof(cl_mem),
                          buf);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t numworkitems = size;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  kernel,
				  1, NULL,
				  &numworkitems,
				  NULL,
				  nwait, 
				  wait, 
				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field. OpenCL version.
void dtfield_CL(field *f, real tnow, cl_mem *wn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done)
{

  set_buf_to_zero_cl(&f->dtwn_cl, f->wsize, f,
  		     nwait, wait, &f->clv_zbuf);
  
  // Macrocell interfaces must be launched serially
  const int ninterfaces = f->macromesh.nmacrointerfaces;
  
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = f->macromesh.macrointerface[i];
    DGMacroCellInterface_CL((void*) (f->mface + ifa), f, wn_cl,
    			    1,
    			    i == 0 ? &f->clv_zbuf : f->clv_mci + i - 1,
    			    f->clv_mci + i);
  }

  cl_event *startboundary = ninterfaces > 0 ?
    f->clv_mci + ninterfaces - 1: &f->clv_zbuf;
  
  // Boundary terms may also need to be launched serially.
  const int nboundaryfaces = f->macromesh.nboundaryfaces;
  for(int i = 0; i < nboundaryfaces; ++i) {
    int ifa = f->macromesh.boundaryface[i];
    DGBoundary_CL((void*) (f->mface + ifa), f, wn_cl,
		  1,
		  i == 0 ?  startboundary : f->clv_boundary + i - 1,
		  f->clv_boundary + i);
  }

  // If there are no interfaces and no boundaries, just wait for zerobuf.
  cl_event *startmacroloop;
  if(ninterfaces == 0 && nboundaryfaces == 0) {
    startmacroloop = &f->clv_zbuf;
  } else {
    if(nboundaryfaces > 0) {
      startmacroloop = f->clv_boundary + nboundaryfaces - 1;
    } else {
      startmacroloop = f->clv_mci + ninterfaces - 1;
    }
  }

  cl_event *fluxdone = f->clv_flux2;
  unsigned int ndim = 3;
  if(f->macromesh.is2d) {
    ndim = 2;
    fluxdone = f->clv_flux1;
  }
  if(f->macromesh.is1d) {
    ndim = 1;
    fluxdone = f->clv_flux0;
  }
  
  cl_event *dtfielddone = f->use_source_cl ? f->clv_source : f->clv_mass;
  
  // The kernels for the intra-macrocell computations can be launched
  // in parallel between macrocells.
  const int nmacro = f->macromesh.nbelems;
  for(int ie = 0; ie < nmacro; ++ie) {
    //printf("ie: %d\n", ie);
    MacroCell *mcelli = f->mcell + ie;
    
    DGFlux_CL(f, 0, ie, wn_cl, 
	      1, startmacroloop,
	      f->clv_flux0 + ie);

    if(ndim > 1) {
      DGFlux_CL(f, 1, ie, wn_cl, 
		1, f->clv_flux0 + ie,
		f->clv_flux1 + ie);
    }
    
    if(ndim > 2) {
      DGFlux_CL(f, 2, ie, wn_cl, 
		1, f->clv_flux1 + ie, 
		f->clv_flux2 + ie);
    }
    
    DGVolume_CL(mcelli, f, wn_cl,
  		1,
  		fluxdone + ie,
  		f->clv_volume + ie);

    DGMass_CL(mcelli, f,
  	      1,
  	      f->clv_volume + ie,
  	      f->clv_mass + ie);

    if(f->use_source_cl) {
      DGSource_CL(mcelli, f, tnow, wn_cl, 
		  1,
		  f->clv_mass + ie,
		  f->clv_source + ie);
    }
  }
  clWaitForEvents(nmacro, dtfielddone);

  // Add times for sources after everything is finished
  f->zbuf_time += clv_duration(f->clv_zbuf);

  for(int i = 0; i < ninterfaces; ++i)
    f->minter_time += clv_duration(f->clv_mci[i]);

  for(int i = 0; i < nboundaryfaces; ++i)
    f->boundary_time += clv_duration(f->clv_boundary[i]);

  for(int ie = 0; ie < nmacro; ++ie) {
    if(f->use_source_cl)
      f->source_time += clv_duration(f->clv_source[ie]);
    f->flux_time += clv_duration(f->clv_flux0[ie]);
    if(ndim > 1) 
      f->flux_time += clv_duration(f->clv_flux1[ie]);
    if(ndim > 2)
      f->flux_time += clv_duration(f->clv_flux2[ie]);
    f->vol_time += clv_duration(f->clv_volume[ie]);
    f->mass_time += clv_duration(f->clv_mass[ie]);
  }

  if(done != NULL)
    clSetUserEventStatus(*done, CL_COMPLETE);
}

// Set kernel arguments for first stage of RK2
void init_RK2_CL_stage1(MacroCell *mcell, field *f,
			const real dt, cl_mem *wnp1_cl)
{
  cl_kernel kernel = f->RK_out_CL;
  cl_int status;
  int argnum = 0;

  //__global real *wnp1 // field values
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wnp1_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  // __global real *wn, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &f->wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real* dtwn // time derivative
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  real halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &halfdt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &mcell->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Launch first stage of RK2 integration
void RK2_CL_stage1(MacroCell *mcell, field *f,
		   size_t numworkitems, // FIXME: remove
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  size_t nwork = mcell->nreal;
  
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->RK_out_CL,
				  1, 
				  NULL,
				  &nwork,
				  NULL,
				  nwait, 
				  wait, 
				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Set kernel arguments for second stage of RK2
void init_RK2_CL_stage2(MacroCell *mcell, field *f, const real dt) 
{
  cl_kernel kernel = f->RK_in_CL;
  cl_int status;
  int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &f->wn_cl);
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &f->dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &mcell->woffset);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Launch second stage of RK2 integration
void RK2_CL_stage2(MacroCell *mcell, field *f,
		   size_t numworkitems, // FIXME: remove
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;

  size_t nwork = mcell->nreal;

  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->RK_in_CL,
				  1, NULL,
				  &nwork,
				  NULL,
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void RK4_CL_stageA(field *f, 
		   cl_mem *wnp1, cl_mem *wn, cl_mem *dtw, 
		   const real dt, const int sizew, size_t numworkitems, 
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  // l_1 = w_n + 0.5dt * S(w_n, t_0)
  
  cl_kernel kernel = f->RK4_first_stages;
  cl_int status;
  int argnum = 0;

  //__global real *wnp1 // field values
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wnp1);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  // __global real *wn, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wn);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real *dtwn // time derivative
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //real dt,
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
  				  kernel,
  				  1,
  				  NULL,
  				  &numworkitems,
  				  NULL,
  				  nwait,
  				  wait,
  				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void RK4_final_inplace_CL(field *f, 
			  cl_mem *w_cl, cl_mem *l1, cl_mem *l2, cl_mem *l3, 
			  cl_mem *dtw_cl, const real dt, 
			  const size_t numworkitems, 
			  cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = f->RK4_final_stage;
  cl_int status;
  int argnum = 0;

  // __global real *w,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *l1,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l1);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *l2,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l2);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *l3,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l3);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *dtw, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // const real dt
  real halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
  				  kernel,
  				  1,
  				  NULL,
  				  &numworkitems,
  				  NULL,
  				  nwait,
  				  wait,
  				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Time integration by a fourth-order Runge-Kutta algorithm, OpenCL
// version.
void RK4_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  if(dt <= 0)
    dt = set_dt(f);

  f->itermax = tmax / dt;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(f->interp_param + 1);
  int iter = 0;

  cl_int status;  
  cl_mem l1 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(real) * f->wsize,
			     NULL,
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  cl_mem l2 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(real) * f->wsize,
			     NULL,
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  cl_mem l3 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(real) * f->wsize,
			     NULL,
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t numworkitems = f->wsize;

  cl_mem *w = &f->wn_cl;
  cl_mem *dtw = &f->dtwn_cl;

  int nstages = 4;
  cl_event source[nstages];
  cl_event stage[nstages];
  for(int i = 0; i < nstages; ++i) {
    source[i] = clCreateUserEvent(f->cli.context, &status);
    stage[i] = clCreateUserEvent(f->cli.context, &status);
  }

  clWaitForEvents(nwait, wait);

  printf("Starting RK4_CL\n");
  
  status = clSetUserEventStatus(stage[3], CL_COMPLETE);

  struct timeval t_start;
  struct timeval t_end;
  gettimeofday(&t_start, NULL);
  while(f->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    // stage 0
    dtfield_CL(f, f->tnow, w, 
	       1, stage + 3, source);
    RK4_CL_stageA(f, &l1, w, dtw,
    		  0.5 * dt, sizew, numworkitems,
		  1, source, stage);

    f->tnow += 0.5 * dt;
    
    // stage 1
    dtfield_CL(f, f->tnow, &l1, 1, stage, source + 1);
    RK4_CL_stageA(f, &l2, w, dtw,
    		  0.5 * dt, sizew, numworkitems,
		  1, source + 1, stage + 1);
    
    // stage 2
    dtfield_CL(f, f->tnow, &l2, 
	       1, stage + 1, source + 2);
    RK4_CL_stageA(f, &l3, w, dtw,
    		  dt, sizew, numworkitems,
		  1, source + 2, stage + 2);

    f->tnow += 0.5 * dt;
    
    // stage 3
    dtfield_CL(f, f->tnow, &l3, 
	       1, stage + 2, source + 3);
    RK4_final_inplace_CL(f, w, &l1, &l2, &l3, 
			 dtw, dt, numworkitems,
			 1, source + 3, stage + 3);

    for(int i = 0; i < nstages; ++i)
      f->rk_time += clv_duration(stage[i]);
    
    iter++;
  }
  gettimeofday(&t_end, NULL); 

  if(done != NULL)
    status = clSetUserEventStatus(*done, CL_COMPLETE);

 double rkseconds = (t_end.tv_sec - t_start.tv_sec) * 1.0 // seconds
    + (t_end.tv_usec - t_start.tv_usec) * 1e-6; // microseconds
  printf("\nTotal RK time (s):\n%f\n", rkseconds);
  printf("\nTotal RK time per time-step (s):\n%f\n", rkseconds / iter );
 
  
  for(int i = 0; i < nstages; ++i) {
    clReleaseEvent(source[i]);
    clReleaseEvent(stage[i]);
  }
  
  printf("\nt=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
}


// Time integration by a second-order Runge-Kutta algorithm, OpenCL
// version.
void RK2_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  if(dt <= 0)
    dt = set_dt(f);

  f->itermax = tmax / dt;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int iter = 0;

  cl_int status;
  cl_mem wnp1_cl = clCreateBuffer(f->cli.context,
				  0, // no flags
				  sizeof(real) * f->wsize,
				  NULL, // do not use a host pointer
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Set up kernels



  const int nmacro = f->macromesh.nbelems;

  cl_event *stage1 =  calloc(nmacro, sizeof(cl_event));
  cl_event *stage2 =  calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    stage1[ie] = clCreateUserEvent(f->cli.context, &status);
    stage2[ie] = clCreateUserEvent(f->cli.context, &status);
    status = clSetUserEventStatus(stage2[ie], CL_COMPLETE);
  }
  
  cl_event source1 = clCreateUserEvent(f->cli.context, &status);
  //cl_event stage1 =  clCreateUserEvent(f->cli.context, &status);
  cl_event source2 = clCreateUserEvent(f->cli.context, &status);
  //cl_event stage2 =  clCreateUserEvent(f->cli.context, &status);

  if(nwait > 0)
    clWaitForEvents(nwait, wait);

  printf("Starting RK2_CL\n");
  
  struct timeval t_start;
  struct timeval t_end;
  gettimeofday(&t_start, NULL);
  while(f->tnow < tmax) {
    //printf("iter: %d\n", iter);
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
    
    dtfield_CL(f, f->tnow, &f->wn_cl, nmacro, stage2, &source1);
    
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      init_RK2_CL_stage1(mcell, f, dt, &wnp1_cl);
      RK2_CL_stage1(mcell, f, f->wsize, 1, &source1, stage1 + ie);
    }
    
    f->tnow += 0.5 * dt;

    dtfield_CL(f, f->tnow, &wnp1_cl, nmacro, stage1, &source2);
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      init_RK2_CL_stage2(mcell, f, dt);
      RK2_CL_stage2(mcell, f, f->wsize, 1, &source2, stage2 + ie);
    }
    
    f->tnow += 0.5 * dt;

    //f->rk_time += clv_duration(stage1); //FIXME: restore
    //f->rk_time += clv_duration(stage2);
    iter++;
  }
  gettimeofday(&t_end, NULL);
  if(done != NULL) 
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  //clReleaseEvent(stage2);
  //clReleaseEvent(stage1);
  
  double rkseconds = (t_end.tv_sec - t_start.tv_sec) * 1.0 // seconds
    + (t_end.tv_usec - t_start.tv_usec) * 1e-6; // microseconds
  printf("\nTotal RK time (s):\n%f\n", rkseconds);
  printf("\nTotal RK time per time-step (s):\n%f\n", rkseconds / iter );
 
  printf("\nt=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
}

void show_cl_timing(field *f)
{
  printf("\n");
  printf("Device characteristics:\n");
  printf("\tname:\t%s\n", f->cli.devicename);
  double dev_gflops = cl_dev_gflops(f->cli.devicename);
  double dev_bwidth = cl_dev_bwidth(f->cli.devicename);
  printf("\tgflops:   \t%f\n", dev_gflops);
  printf("\tbandwidth:\t%f\n", dev_bwidth);

  printf("\n");
  printf("Roofline counts:\n");
  unsigned long int flops_total = f->flops_vol + f->flops_flux + f->flops_mass;
  unsigned long int reads_total = f->reads_vol + f->reads_flux + f->reads_mass;
  printf("Number of real multiplies in kernels:               %lu\n", 
	 flops_total);
  printf("Number of reads/writes of reals from global memory: %lu\n", 
	 reads_total);
  //printf("NB: real multiplies assumes the flux involved 3m multiplies.\n");
  printf("Terms included in roofline: volume, flux, and mass.\n");

  cl_ulong roofline_time_ns = f->vol_time + f->flux_time + f->mass_time;
  double roofline_time_s = 1e-9 * roofline_time_ns;
  double roofline_flops = flops_total / roofline_time_s;
  double roofline_bw = sizeof(real) * reads_total / roofline_time_s;
  
  // Volume terms
  double vol_time_s = 1e-9 * f->vol_time;
  double vol_flops = f->flops_vol / vol_time_s;
  double vol_bw = sizeof(real) * f->reads_vol / vol_time_s;
  printf("DGVol:  GFLOP/s: %f\tbandwidth (GB/s): %f\n", 
	 1e-9 * vol_flops, 1e-9 * vol_bw);

  // Flux terms
  double flux_time_s = 1e-9 * f->flux_time;
  double flux_flops = f->flops_flux / flux_time_s;
  double flux_bw = sizeof(real) * f->reads_flux / flux_time_s;
  printf("DGFlux: GFLOP/s: %f\tbandwidth (GB/s): %f\n", 
	 1e-9 * flux_flops, 1e-9 * flux_bw);
  
  // Mass terms
  double mass_time_s = 1e-9 * f->mass_time;
  double mass_flops = f->flops_mass / mass_time_s;
  double mass_bw = sizeof(real) * f->reads_mass / mass_time_s;
  printf("DGMass: GFLOP/s: %f\tbandwidth (GB/s): %f\n", 
	 1e-9 * mass_flops, 1e-9 * mass_bw);

  printf("Total:  GFLOP/s: %f\tbandwidth (GB/s): %f\n", 
	 1e-9 * roofline_flops, 1e-9 * roofline_bw);
  
  printf("\n");

  printf("Kernel execution times:\n"); 
  cl_ulong total = 0;
  total += f->zbuf_time;
  total += f->minter_time;
  total += f->boundary_time;
  total += f->vol_time;
  total += f->mass_time;
  total += f->source_time;
  total += f->rk_time;
  total += f->flux_time;
  
  real N = 100.0 / total;

  cl_ulong ns;

  ns = f->zbuf_time;
  printf("set_buf_to_zero_cl time:      %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->minter_time;
  printf("DGMacroCellInterface_CL time: %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->boundary_time;
  printf("DGBoundary time:              %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->vol_time;
  printf("DGVolume_CL time:             %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->flux_time;
  printf("DGFlux_CL time:               %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->mass_time;
  printf("DGMass_CL time:               %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->source_time;
  printf("DGSource_CL time:             %f%% \t%luns \t%fs\n", 
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = f->rk_time;
  printf("rk time:                      %f%% \t%luns \t%fs\n", 
	 ns*N, (unsigned long) ns, 1e-9 * ns);

  printf("\n");
  
  ns = total;
  double total_time = 1e-9 * ns;
  printf("total time:                   %f%% \t%luns \t%fs\n", 
	 ns*N, (unsigned long) ns, total_time);

  printf("\n");
  print_kernel_perf(dev_gflops, dev_bwidth,
		    flops_total, reads_total, total);
  printf("\n");
}

#endif // _WITH_OPENCL

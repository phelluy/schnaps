#include "field.h"
#include "field_cl.h"
#include "clinfo.h"
#include "clutils.h"
#include <assert.h>
#include <string.h>

#ifdef _WITH_OPENCL
void CopyfieldtoCPU(field *f) {

  cl_int status;

  // ensures that all the buffers are mapped
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

cl_ulong clv_startime(cl_event clv) 
{
  cl_ulong time;
  clWaitForEvents(1, &clv);
  clGetEventProfilingInfo(clv,
			  CL_PROFILING_COMMAND_START,
			  sizeof(time),
			  &time, NULL);
  return time;
}

cl_ulong clv_endtime(cl_event clv) 
{
  cl_ulong time;
  clWaitForEvents(1, &clv);
  clGetEventProfilingInfo(clv,
			  CL_PROFILING_COMMAND_END,
			  sizeof(time),
			  &time, NULL);
  return time;
}

cl_ulong clv_duration(cl_event clv)
{
  return clv_endtime(clv) - clv_startime(clv);
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

void initDGBoundary_CL(field *f, 
		       int ieL, int locfaL,
		       cl_mem physnodeL_cl, 
		       size_t cachesize)
{

  cl_int status;
  cl_kernel kernel = f->dgboundary;

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

  // int ieL,                      // 2: left macrocell
  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(int),
  			  &ieL);
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
                          &physnodeL_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *wn,          // 5: field
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->wn_cl);
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


// Set the loop-dependant kernel arguments for DGMacroCellInterface_CL
void initDGMacroCellInterface_CL(field *f, 
				      int ieL, int ieR, int locfaL, int locfaR,
				      cl_mem physnodeL_cl, cl_mem physnodeR_cl,
				      size_t cachesize)
{
  //printf("loop_initDGMacroCellInterface_CL\n");
  cl_int status;
  cl_kernel kernel = f->dginterface;

  // Set kernel arguments
  unsigned int argnum = 0;

  //__constant int *param,        // 0: interp param
  status = clSetKernelArg(kernel,
			 argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //int ieL,                      // 1: left macrocell 
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &ieL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //int ieR,                      // 2: right macrocell
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &ieR);
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

  //__constant real *physnodeL, // 5: left physnode
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &physnodeL_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__constant real *physnodeR, // 6: right physnode
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &physnodeR_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real *wn,          // 7: field 
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->wn_cl);
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

  //printf("... loop_initDGMacroCellInterface_CL done.\n");
}

// Set up kernel arguments, etc, for DGMass_CL.
void init_DGMass_CL(field *f)
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
  argnum++;
      
  /* __constant real* physnode,  // macrocell nodes */
  status = clSetKernelArg(f->dgmass,           // kernel name
                          argnum++,              // arg num
                          sizeof(cl_mem),
                          &f->physnode_cl);     // opencl buffer
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

// Set kernel argument for DGVolume_CL
void init_DGVolume_CL(field *f, cl_mem *wn_cl, size_t cachesize)
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

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->physnode_cl);
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

// Set kernel argument for DGVolume_CL
void init_DGSource_CL(field *f, cl_mem *wn_cl, size_t cachesize)
{
  //printf("DGVolume cachesize:%zu\n", cachesize);

  cl_int status;
  int argnum = 0;
  cl_kernel kernel = f->dgsource;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(real),
  			  &f->tnow);
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

// Update the cl buffer with physnode data depending in the
// macroelement with index ie
void update_physnode_cl(field *f, int ie, cl_mem physnode_cl, real *physnode,
			cl_ulong *time,
			cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  
  void *chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
				    physnode_cl,
				    CL_TRUE,
				    CL_MAP_WRITE,
				    0, // offset
				    sizeof(real) * 60, // buffersize
				    nwait, 
				    wait, 
				    &f->clv_mapdone,
				    &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == physnode);

  if(time != NULL)
    *time += clv_duration(f->clv_mapdone);

  int ie20 = 20 * ie;
  for(int inoloc = 0; inoloc < 20; ++inoloc) {
    int ino = 3 * f->macromesh.elem2node[ie20 + inoloc];
    real *iphysnode = physnode + 3 * inoloc;
    real *nodeino = f->macromesh.node + ino;
    iphysnode[0] = nodeino[0];
    iphysnode[1] = nodeino[1];
    iphysnode[2] = nodeino[2];
  }

  status = clEnqueueUnmapMemObject(f->cli.commandqueue,
				   physnode_cl,
				   physnode,
				   1, &f->clv_mapdone, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  if(time != NULL)
    *time += clv_duration(f->clv_mapdone);
}

void DGMacroCellInterface_CL(void *mf, field *f, cl_mem *wn_cl,
			     cl_uint nwait, cl_event *wait, cl_event *done) 
{
  //printf("DGMacroCellInterface_CL\n");

  clWaitForEvents(nwait, wait);

  /* if(done != NULL) { */
  /*   status = clSetUserEventStatus(*done, CL_); */
  
  MacroFace *mface = (MacroFace*) mf;
  int *param = f->interp_param;
  cl_int status;
  cl_kernel kernel = f->dginterface;

  status = clSetKernelArg(kernel,
                          8,
                          sizeof(cl_mem),
                          wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Loop on the macro faces
  int start = mface->first;
  int end = mface->last_p1;
  for(int ifa = start; ifa < end; ++ifa) {
    //printf("ifa: %d\n", ifa);
        
    int ieL =    f->macromesh.face2elem[4 * ifa];
    int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    int ieR =    f->macromesh.face2elem[4 * ifa + 2];
    int locfaR = f->macromesh.face2elem[4 * ifa + 3];

    update_physnode_cl(f, ieL, f->physnode_cl, f->physnode, &f->minter_time,
		       0, 
		       NULL,
		       &f->clv_interupdate);

    if(ieR >= 0) {
      update_physnode_cl(f, ieR, f->physnodeR_cl, f->physnodeR, 
			 &f->minter_time,
			 1, 
			 &f->clv_interupdate,
			 &f->clv_interupdateR);
      f->minter_time += clv_duration(f->clv_interupdateR);
    } else {
      clWaitForEvents(1, &f->clv_interupdate);
      status = clSetUserEventStatus(f->clv_interupdateR, CL_COMPLETE);
    }

    size_t numworkitems = NPGF(f->interp_param + 1, locfaL);
    if(ieR >= 0) {
    
      // Set the remaining loop-dependant kernel arguments
      size_t kernel_cachesize = 1;
      initDGMacroCellInterface_CL(f, 
				  ieL, ieR, locfaL, locfaR, 
				  f->physnode_cl, f->physnodeR_cl, 
				  kernel_cachesize);

      status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				      kernel,
				      1, // cl_uint work_dim,
				      NULL, // global_work_offset,
				      &numworkitems, // global_work_size, 
				      NULL, // size_t *local_work_size, 
				      1,  // nwait, 
				      &f->clv_interupdateR, // *wait_list,
				      &f->clv_interkernel); // *event
      if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
      assert(status >= CL_SUCCESS);
      f->minter_time += clv_duration(f->clv_interkernel);
    } else {
      size_t cachesize = 1; // TODO make use of cache
      initDGBoundary_CL(f, ieL, locfaL, f->physnode_cl, cachesize);
      status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				      f->dgboundary,
				      1, // cl_uint work_dim,
				      NULL, // global_work_offset,
				      &numworkitems, // global_work_size, 
				      NULL, // size_t *local_work_size, 
				      1,  // nwait, 
				      &f->clv_interupdateR, // *wait_list,
				      &f->clv_interkernel); // *event
      if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
      assert(status >= CL_SUCCESS);
      f->boundary_time += clv_duration(f->clv_interkernel);
    }


  }
  
  clWaitForEvents(1, &f->clv_interkernel);

  if(done != NULL)
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  //printf("... DGMacroCellInterface_CL done.\n");
}

// Apply division by the mass matrix OpenCL version
void DGMass_CL(void *mc, field *f,
	       cl_uint nwait, cl_event *wait, cl_event *done) 
{
  //printf("DGMass_CL\n");

  MacroCell *mcell = (MacroCell*) mc;
  int *param = f->interp_param;
  cl_int status;

  init_DGMass_CL(f);
 
  clSetUserEventStatus(f->clv_mass, CL_COMPLETE);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Loop on the elements
  const int start = mcell->first;
  const int end = mcell->last_p1;
  for(int ie = start; ie < end; ie++) {
    
    status = clSetKernelArg(f->dgmass, 
			    1, 
			    sizeof(int), 
			    (void *)&ie);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

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
				    &f->clv_mass); // cl_event *event
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    f->mass_time += clv_duration(f->clv_mass);
  }
  
  clWaitForEvents(1, &f->clv_mass);
  if(done != NULL) {
    status = clSetUserEventStatus(*done, CL_COMPLETE);
  }
  //printf("... DGMass_CL done.\n");
}

void init_DGFlux_CL(field *f, int ie, int dim0, cl_mem *wn_cl, 
		    size_t cachesize)
{
  //printf("DGFlux cachesize:%zu\n", cachesize);

  cl_kernel kernel = f->dgflux;
  cl_int status;
  int argnum = 0; 
  
  // __constant int *param, // interp param
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // int ie, // macrocel index
  status = clSetKernelArg(kernel, 
			  argnum++, 
			  sizeof(int), 
			  (void *)&ie);
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
                          &f->physnode_cl);
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

  // debugging output
  /* printf("\ndim0: %d\n", dim0); */
  /* printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]); */
  /* printf("deg: %d %d %d\n", deg[0], deg[1], deg[2]); */
  /* printf("nraf: %d %d %d\n", nraf[0], nraf[1], nraf[2]); */
  /* printf("npgf: %d\n", npgf); */
  /* printf("nf: %d\n", nf); */

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

// Apply division by the mass matrix OpenCL version
void DGVolume_CL(void *mc, field *f, cl_mem *wn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done) 
{
  MacroCell *mcell = (MacroCell*) mc; // FIXME: just pass ie, it'll be easier.
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
  
  init_DGVolume_CL(f, wn_cl, 2 * groupsize * m);
  
  const int start = mcell->first;
  const int end = mcell->last_p1;
  
  //printf("DGVolume_CL loop: %d\n", end - start); // This is always 1!!!

  // Loop on the elements
  for(int ie = start; ie < end; ie++) {
    f->flops_vol += numworkitems * nmultsdgvol;
    f->reads_vol += numworkitems * nreadsdgvol;

    status = clSetKernelArg(kernel,
			    1,
			    sizeof(int),
			    &ie);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

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
    //f->vol_time += clv_duration(f->clv_volkernel);
  }
}

// Apply division by the mass matrix OpenCL version
void DGSource_CL(void *mc, field *f, cl_mem *wn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done) 
{
  MacroCell *mcell = (MacroCell*) mc; // FIXME: just pass ie, it'll be easier.
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
   
  init_DGSource_CL(f, wn_cl, 2 * groupsize * m);
  
  const int start = mcell->first;
  const int end = mcell->last_p1;
  
  // Loop on the elements
  for(int ie = start; ie < end; ie++) {
    status = clSetKernelArg(kernel,
			    1,
			    sizeof(int),
			    &ie);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

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
  if(done != NULL)
    f->zbuf_time += clv_duration(f->clv_zbuf);  
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field. OpenCL version.
void dtfield_CL(field *f, cl_mem *wn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done) {
  set_buf_to_zero_cl(&f->dtwn_cl, f->wsize, f,
  		     nwait, wait, &f->clv_zbuf);
  
  //printf("f->macromesh.nbfaces: %d\n", f->macromesh.nbfaces);
  for(int ifa = 0; ifa < f->macromesh.nbfaces; ++ifa) {
    //printf("ifa: %d\n", ifa);
    DGMacroCellInterface_CL((void*) (f->mface + ifa), f, wn_cl,
    			    1,
    			    &f->clv_zbuf,
    			    &f->clv_mci);
  }

  //printf("f->macromesh.nbelems: %d\n", f->macromesh.nbelems);
  
  // Start the loop
  if(f->use_source_cl) {
    clSetUserEventStatus(f->clv_source, CL_COMPLETE); 
  } else {
    clSetUserEventStatus(f->clv_mass, CL_COMPLETE);
  }
  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    //printf("ie: %d\n", ie);
    MacroCell *mcelli = f->mcell + ie;
    
    update_physnode_cl(f, ie, f->physnode_cl, f->physnode, &f->vol_time,
		       1, f->use_source_cl ? &f->clv_source : &f->clv_mass,
		       &f->clv_physnodeupdate);
    //f->???_time += clv_duration(f->clv_physnodeupdate);

    unsigned int ndim = f->macromesh.is2d ? 2 : 3;
    DGFlux_CL(f, 0, ie, wn_cl, 
	      1, &f->clv_physnodeupdate, 
	      f->clv_flux);
    f->flux_time += clv_duration(f->clv_flux[0]);
    for(unsigned int d = 1; d < ndim; ++d) {
      DGFlux_CL(f, d, ie, wn_cl, 
		1, f->clv_flux + (d - 1) % ndim,  f->clv_flux + d);
      f->flux_time += clv_duration(f->clv_flux[d]);
    }

    DGVolume_CL(mcelli, f, wn_cl,
  		1,
  		f->clv_flux + ndim - 1,
  		&f->clv_volume);
    f->vol_time += clv_duration(f->clv_volume);

    DGMass_CL(mcelli, f,
  	      1,
  	      &f->clv_volume,
  	      &f->clv_mass);
    f->mass_time += clv_duration(f->clv_mass);

    if(f->use_source_cl) {
      DGSource_CL(mcelli, f, wn_cl, 
		  1,
		  &f->clv_mass,
		  &f->clv_source);
      f->source_time += clv_duration(f->clv_source);
    }
  }
  clWaitForEvents(1, &f->clv_mass);
  if(done != NULL)
    clSetUserEventStatus(*done, CL_COMPLETE);
}

// Set kernel arguments for first stage of RK2
void init_RK2_CL_stage1(field *f, const real dt, cl_mem *wnp1_cl)
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
  
  //real dt, // time step for the stage
  real halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &halfdt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Launch first stage of RK2 integration
void RK2_CL_stage1(field *f, size_t numworkitems,
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->RK_out_CL,
				  1, 
				  NULL,
				  &numworkitems,
				  NULL,
				  nwait, 
				  wait, 
				  done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  if(done != NULL)
    f->rk_time += clv_duration(*done);
}

// Set kernel arguments for second stage of RK2
void init_RK2_CL_stage2(field *f, const real dt) 
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
}

// Launch second stage of RK2 integration
void RK2_CL_stage2(field *f, size_t numworkitems,
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->RK_in_CL,
				  1, NULL,
				  &numworkitems,
				  NULL,
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  if(done != NULL)
    f->rk_time += clv_duration(*done);
}

void RK4_CL_stageA(field *f, 
		   cl_mem *wnp1, cl_mem *wn, cl_mem *dtw, 
		   const real dt, const int sizew, size_t numworkitems, 
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  // l_1 = w_n + 0.5dt * S(w_n, t_0)
  
  cl_kernel kernel = f->RK_out_CL;
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

  if(done != NULL) 
    f->rk_time += clv_duration(*done);
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

  if(done != NULL)
    f->rk_time += clv_duration(*done);
}

// Time integration by a fourth-order Runge-Kutta algorithm, OpenCL
// version.
void RK4_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  if(dt <= 0)
    dt = set_dt(f);

  clWaitForEvents(nwait, wait);

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

  status = clSetUserEventStatus(stage[3], CL_COMPLETE);

  while(f->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    // stage 0
    dtfield_CL(f, w, 
	       1, stage + 3, source);
    RK4_CL_stageA(f, &l1, w, dtw,
    		  0.5 * dt, sizew, numworkitems,
		  1, source, stage);
    
    // stage 1
    dtfield_CL(f, &l1, 1, stage, source + 1);
    RK4_CL_stageA(f, &l2, w, dtw,
    		  0.5 * dt, sizew, numworkitems,
		  1, source + 1, stage + 1);
    
    // stage 2
    dtfield_CL(f, &l2, 
	       1, stage + 1, source + 2);
    RK4_CL_stageA(f, &l3, w, dtw,
    		  dt, sizew, numworkitems,
		  1, source + 2, stage + 2);
    
    // stage 3
    dtfield_CL(f, &l3, 
	       1, stage + 2, source + 3);
    RK4_final_inplace_CL(f, w, &l1, &l2, &l3, 
			 dtw, dt, numworkitems,
			 1, source + 3, stage + 3);

    f->tnow += dt;
    iter++;
  }
  if(done != NULL)
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
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
  init_RK2_CL_stage1(f, dt, &wnp1_cl);
  init_RK2_CL_stage2(f, dt);

  cl_event source1 = clCreateUserEvent(f->cli.context, &status);
  cl_event stage1 = clCreateUserEvent(f->cli.context, &status);
  cl_event source2 = clCreateUserEvent(f->cli.context, &status);
  cl_event stage2 = clCreateUserEvent(f->cli.context, &status);
  status = clSetUserEventStatus(stage2, CL_COMPLETE);

  if(nwait > 0)
    clWaitForEvents(nwait, wait);
  while(f->tnow < tmax) {
    //printf("iter: %d\n", iter);
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    dtfield_CL(f, &f->wn_cl, 1, &stage2, &source1);
    RK2_CL_stage1(f, f->wsize, 1, &source1, &stage1);

    f->tnow += 0.5 * dt;

    dtfield_CL(f, &wnp1_cl, 1, &stage1, &source2);
    RK2_CL_stage2(f, f->wsize, 1, &source2, &stage2);

    f->tnow += 0.5 * dt;
    iter++;
  }
  if(done != NULL) 
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
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
  
  // vol terms
  double vol_time_s = 1e-9 * f->vol_time;
  double vol_flops = f->flops_vol / vol_time_s;
  double vol_bw = sizeof(real) * f->reads_vol / vol_time_s;
  printf("DGVol:  GFLOP/s: %f\tbandwidth (GB/s): %f\n", 
	 1e-9 * vol_flops, 1e-9 * vol_bw);

  // flux terms
  double flux_time_s = 1e-9 * f->flux_time;
  double flux_flops = f->flops_flux / flux_time_s;
  double flux_bw = sizeof(real) * f->reads_flux / flux_time_s;
  printf("DGFlux: GFLOP/s: %f\tbandwidth (GB/s): %f\n", 
	 1e-9 * flux_flops, 1e-9 * flux_bw);
  
  // mass terms
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
  printf("Total kernel performance:\n");
  print_kernel_perf(dev_gflops, dev_bwidth,
		    flops_total, reads_total, total);
  
  printf("\n");
}

#endif // _WITH_OPENCL

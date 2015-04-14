#include "field.h"
#include "field_cl.h"
#include "clinfo.h"
#include "clutils.h"
#include <assert.h>

#ifdef _WITH_OPENCL
void CopyfieldtoCPU(field *f) {

  cl_int status;

  // ensures that all the buffers are mapped
  status = clFinish(f->cli.commandqueue);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  void *chkptr;
  chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
			      f->dtwn_cl, // buffer to copy from
			      CL_TRUE, // block until the buffer is available
			      CL_MAP_READ, // we just want to see the results
			      0, // offset
			      f->wsize * sizeof(double), // buffersize
			      0, NULL, NULL, // events management
			      &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == f->dtwn);

  chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
			      f->wn_cl, // buffer to copy from
			      CL_TRUE, // block until the buffer is available
			      CL_MAP_READ, // we just want to see the results
			      0, // offset
			      f->wsize * sizeof(double), // buffersize
			      0, NULL, NULL, // events management
			      &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == f->wn);

  status = clFinish(f->cli.commandqueue);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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

// Set OpenCL kernel arguments for DGMacroCellInterface
void initDGMacroCellInterface_CL(field *f, 
				 cl_mem physnodeL_cl, cl_mem physnodeR_cl)
{
  //printf("initDGMacroCellInterface_CL\n");

  cl_int status;
  cl_kernel kernel = f->dginterface;

  // Set kernel arguments
  unsigned int argnum = 0;

  // associates the param buffer to the 0th kernel argument
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->param_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // tnow
  argnum++;

  // ieL
  argnum++;

  // ieR
  argnum++;

  // locfaL
  argnum++;

  // locfaR
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &physnodeL_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &physnodeR_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // wn_cl
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  //printf("... initDGMacroCellInterface_CL done.\n");
}

// Set the loop-dependant kernel arguments for DGMacroCellInterface_CL
void loop_initDGMacroCellInterface_CL(field *f, 
				      int ieL, int ieR, int locfaL, int locfaR)
{
  //printf("loop_initDGMacroCellInterface_CL\n");
  cl_int status;
  cl_kernel kernel = f->dginterface;

  // Set kernel arguments
  unsigned int argnum = 1;

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &(f->tnow));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &ieL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &ieR);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaR);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
                          &(f->param_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* int ie, // macrocel index */
  // Set in loop on call.
  argnum++;
      
  /* __constant double* physnode,  // macrocell nodes */
  status = clSetKernelArg(f->dgmass,           // kernel name
                          argnum++,              // arg num
                          sizeof(cl_mem),
                          &f->physnode_cl);     // opencl buffer
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* __global double* dtwn // time derivative */
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Set kernel argument for DGVolume_CL
void init_DGVolume_CL(field *f, cl_mem *wn_cl)
{
  cl_int status;
  int argnum = 0;
  cl_kernel kernel = f->dgvolume;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->param_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->physnode_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Update the cl buffer with physnode data depending in the
// macroelement with index ie
void update_physnode_cl(field *f, int ie, cl_mem physnode_cl, double *physnode,
			cl_ulong *time,
			cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  
  void *chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
				    physnode_cl,
				    CL_TRUE,
				    CL_MAP_WRITE,
				    0, // offset
				    sizeof(cl_double) * 60, // buffersize
				    nwait, 
				    wait, 
				    &(f->clv_mapdone),
				    &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == physnode);

  if(time != NULL)
    *time += clv_duration(f->clv_mapdone);

  int ie20 = 20 * ie;
  for(int inoloc = 0; inoloc < 20; ++inoloc) {
    int ino = 3 * f->macromesh.elem2node[ie20 + inoloc];
    double *iphysnode = physnode + 3 * inoloc;
    double *nodeino = f->macromesh.node + ino;
    iphysnode[0] = nodeino[0];
    iphysnode[1] = nodeino[1];
    iphysnode[2] = nodeino[2];
  }

  status = clEnqueueUnmapMemObject(f->cli.commandqueue,
				   physnode_cl,
				   physnode,
				   1, &(f->clv_mapdone), done);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Loop on the macro faces
  int start = mface->first;
  int end = mface->last_p1;
  for(int ifa = start; ifa < end; ++ifa) {
    //printf("ifa: %d\n", ifa);
    
    initDGMacroCellInterface_CL(f, f->physnode_cl, f->physnodeR_cl);
    
    int ieL =    f->macromesh.face2elem[4 * ifa];
    int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    int ieR =    f->macromesh.face2elem[4 * ifa + 2];
    int locfaR = f->macromesh.face2elem[4 * ifa + 3];

    update_physnode_cl(f, ieL, f->physnode_cl, f->physnode, &(f->minter_time),
		       0, 
		       NULL,
		       &(f->clv_interupdate));

    if(ieR >= 0) {
      update_physnode_cl(f, ieR, f->physnodeR_cl, f->physnodeR, 
			 &(f->minter_time),
			 1, 
			 &(f->clv_interupdate),
			 &(f->clv_interupdateR));
      f->minter_time += clv_duration(f->clv_interupdateR);
    } else {
      clWaitForEvents(1, &(f->clv_interupdate));
      status = clSetUserEventStatus(f->clv_interupdateR, CL_COMPLETE);
    }

    // Set the remaining loop-dependant kernel arguments
    loop_initDGMacroCellInterface_CL(f, ieL, ieR, locfaL, locfaR);

    size_t numworkitems = NPGF(f->interp_param + 1, locfaL);
    status = clEnqueueNDRangeKernel(f->cli.commandqueue,
        			    kernel,
        			    1, // cl_uint work_dim,
				    NULL, // const size_t *global_work_offset,
        			    &numworkitems, // size_t *global_work_size, 
        			    NULL, // size_t *local_work_size, 
        			    1,  // cl_uint num_events_in_wait_list, 
				    &(f->clv_interupdateR), // *event_wait_list,
				    &(f->clv_interkernel)); // *event
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    f->minter_time += clv_duration(f->clv_interkernel);
  }
  
  clWaitForEvents(1, &(f->clv_interkernel));

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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Loop on the elements
  const int start = mcell->first;
  const int end = mcell->last_p1;
  for(int ie = start; ie < end; ie++) {
    
    status = clSetKernelArg(f->dgmass, 
			    1, 
			    sizeof(int), 
			    (void *)&ie);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    // The groupsize is the number of glops in a subcell
    size_t groupsize = (param[1] + 1) * (param[2] + 1) * (param[3] + 1);
    // The total work items number is (the number of glops in a
    // subcell) * (number of subcells)
    size_t numworkitems = param[4] * param[5] * param[6] * groupsize;

    status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				    f->dgmass,
				    1, // cl_uint work_dim,
				    NULL, // size_t *global_work_offset, 
				    &numworkitems, // size_t *global_work_size,
				    &groupsize, // size_t *local_work_size, 
				    nwait, // cl_uint num_events_in_wait_list, 
				    wait, //*event_wait_list, 
				    &(f->clv_mass)); // cl_event *event
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    f->mass_time += clv_duration(f->clv_mass);
  }
  
  clWaitForEvents(1, &(f->clv_mass));
  if(done != NULL) {
    status = clSetUserEventStatus(*done, CL_COMPLETE);
  }
  //printf("... DGMass_CL done.\n");
}

void init_DGFlux_CL(field *f, int ie, int dim0, cl_mem *wn_cl)
{
  cl_kernel kernel = f->dgflux;
  cl_int status;
  int argnum = 0; 
  
  // __constant int *param, // interp param
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->param_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // int ie, // macrocel index
  status = clSetKernelArg(kernel, 
			  argnum++, 
			  sizeof(int), 
			  (void *)&ie);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // int dim0, // face direction
  status = clSetKernelArg(kernel, 
			  argnum++, 
			  sizeof(int), 
			  (void *)&dim0);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
    
  // __constant double *physnode, // macrocell nodes
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->physnode_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //  __global double *wn, // field values
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global double *dtwn // time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGFlux_CL(field *f, int dim0, int ie, cl_mem *wn_cl,
	       cl_uint nwait, cl_event *wait, cl_event *done) 
{
  // Unpack the parameters
  int *param = f->interp_param;
  int deg[3] = {param[1], param[2], param[3]};
  int nraf[3] = {param[4], param[5], param[6]};


  int dim[3];
  dim[0] = dim0;
  dim[1] = (dim[0] + 1) % 3;
  dim[2] = (dim[1] + 1) % 3;


  // Number of points on per face
  int npgf = (deg[dim[1]] + 1) * (deg[dim[2]] + 1);

  // Number of faces
  int nf = (nraf[dim[0]] - 1) * nraf[dim[1]] * nraf[dim[2]];

  // FIXME: temp debugging output
  printf("\ndim0: %d\n", dim0);
  printf("dim: %d %d %d\n", dim[0], dim[1], dim[2]);
  printf("deg: %d %d %d\n", deg[0], deg[1], deg[2]);
  printf("nraf: %d %d %d\n", nraf[0], nraf[1], nraf[2]);
  printf("npgf: %d\n", npgf);
  printf("nf: %d\n", nf);

  if(nf > 0) { // If there are faces to work on, launch the kernel
    // Set kernel args
    init_DGFlux_CL(f, ie, dim0, wn_cl);
     
    // Launch the kernel
    size_t numworkitems = nf * npgf;
    size_t groupsize = npgf;
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
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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

  init_DGVolume_CL(f, wn_cl);
    
  const int start = mcell->first;
  const int end = mcell->last_p1;
  
  //printf("DGVolume_CL loop: %d\n", end - start); // This is always 1!!!

  // Loop on the elements
  for(int ie = start; ie < end; ie++) {
    status = clSetKernelArg(kernel,
			    1,
			    sizeof(int),
			    &ie);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    // The groupsize is the number of glops in a subcell
    size_t groupsize = (param[1] + 1)* (param[2] + 1)*(param[3] + 1);
    // The total work items number is the number of glops in a subcell
    // * number of subcells
    size_t numworkitems = param[4] * param[5] * param[6] * groupsize;
    //printf("groupsize=%zd numworkitems=%zd\n", groupsize, numworkitems);
    status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				    kernel,
				    1,
				    NULL,
				    &numworkitems,
				    &groupsize,
				    nwait,
				    wait,
				    done);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    //f->vol_time += clv_duration(f->clv_volkernel);
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  f->zbuf_time += clv_duration(f->clv_zbuf);  
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field. OpenCL version.
void dtfield_CL(field *f, cl_mem *wn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done) {
  set_buf_to_zero_cl(&(f->dtwn_cl), f->wsize, f,
  		     nwait, wait, &(f->clv_zbuf));
  
  //printf("f->macromesh.nbfaces: %d\n", f->macromesh.nbfaces);
  for(int ifa = 0; ifa < f->macromesh.nbfaces; ++ifa) {
    //printf("ifa: %d\n", ifa);
    DGMacroCellInterface_CL((void*) (f->mface + ifa), f, wn_cl,
    			    1,
    			    &(f->clv_zbuf),
    			    &(f->clv_mci));
  }

  //printf("f->macromesh.nbelems: %d\n", f->macromesh.nbelems);
  
  clSetUserEventStatus(f->clv_mass, CL_COMPLETE); // start the loop
  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    //printf("ie: %d\n", ie);
    MacroCell *mcelli = f->mcell + ie;
    
    update_physnode_cl(f, ie, f->physnode_cl, f->physnode, &(f->vol_time),
		       1, &(f->clv_mass),
		       &(f->clv_physnodeupdate));
    //f->???_time += clv_duration(f->clv_physnodeupdate);

    unsigned int ndim = f->macromesh.is2d ? 2 : 3;
    DGFlux_CL(f, 0, ie, wn_cl, 
	      1, &(f->clv_physnodeupdate), 
	      f->clv_flux);
    //f->flux_time += clv_duration(f->clv_flux[0]);
    for(unsigned int d = 1; d < ndim; ++d) {
      DGFlux_CL(f, d, ie, wn_cl, 
		1, f->clv_flux + (d - 1) % ndim,  f->clv_flux + d);
      //f->flux_time += clv_duration(f->clv_flux[d]);
    }

    DGVolume_CL(mcelli, f, wn_cl,
  		1,
  		f->clv_flux + ndim - 1,
  		&(f->clv_volume));
    f->vol_time += clv_duration(f->clv_volume);

    DGMass_CL(mcelli, f,
  	      1,
  	      &(f->clv_volume),
  	      &(f->clv_mass));
    f->mass_time += clv_duration(f->clv_mass);
  }
  clWaitForEvents(1, &(f->clv_mass));
  if(done != NULL)
    clSetUserEventStatus(*done, CL_COMPLETE);
}

// Set kernel arguments for first stage of RK2
void init_RK2_CL_stage1(field *f, const double dt, cl_mem *wnp1_cl)
{
  cl_kernel kernel = f->RK_out_CL;
  cl_int status;
  int argnum = 0;

  //__global double *wnp1 // field values
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wnp1_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  // __global double *wn, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &(f->wn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global double* dtwn // time derivative
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  //double dt, // time step for the stage
  double halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &halfdt);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  if(done != NULL)
    f->rk_time += clv_duration(*done);
}

// Set kernel arguments for second stage of RK2
void init_RK2_CL_stage2(field *f, const double dt) 
{
  cl_kernel kernel = f->RK_in_CL;
  cl_int status;
  int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &(f->wn_cl));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &(dt));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  if(done != NULL)
    f->rk_time += clv_duration(*done);
}

void RK4_CL_stageA(field *f, 
		   cl_mem *wnp1, cl_mem *wn, cl_mem *dtw, 
		   const double dt, const int sizew, size_t numworkitems, 
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  // l_1 = w_n + 0.5dt * S(w_n, t_0)
  
  cl_kernel kernel = f->RK_out_CL;
  cl_int status;
  int argnum = 0;

  //__global double *wnp1 // field values
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wnp1);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  // __global double *wn, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wn);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global double *dtwn // time derivative
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //double dt,
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &dt);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  if(done != NULL) 
    f->rk_time += clv_duration(*done);
}

void RK4_final_inplace_CL(field *f, 
			  cl_mem *w_cl, cl_mem *l1, cl_mem *l2, cl_mem *l3, 
			  cl_mem *dtw_cl, const double dt, 
			  const size_t numworkitems, 
			  cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = f->RK4_final_stage;
  cl_int status;
  int argnum = 0;

  // __global double *w,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          w_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global double *l1,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l1);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global double *l2,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l2);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global double *l3,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l3);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global double *dtw, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // const double dt
  double halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &dt);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  if(done != NULL)
    f->rk_time += clv_duration(*done);
}

// Time integration by a fourth-order Runge-Kutta algorithm, OpenCL
// version.
void RK4_CL(field *f, double tmax, 
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  clWaitForEvents(nwait, wait);

  f->itermax = tmax / f->dt;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int sizew = f->macromesh.nbelems * f->model.m * NPG(f->interp_param + 1);
  int iter = 0;

  cl_int status;  
  cl_mem l1 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(double) * f->wsize,
			     NULL,
			     &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  cl_mem l2 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(double) * f->wsize,
			     NULL,
			     &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  cl_mem l3 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(double) * f->wsize,
			     NULL,
			     &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t numworkitems = f->wsize;

  cl_mem *w = &(f->wn_cl);
  cl_mem *dtw = &(f->dtwn_cl);

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
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, f->dt);

    // stage 0
    dtfield_CL(f, w, 
	       1, stage + 3, source);
    RK4_CL_stageA(f, &l1, w, dtw,
    		  0.5 * f->dt, sizew, numworkitems,
		  1, source, stage);
    
    // stage 1
    dtfield_CL(f, &l1, 1, stage, source + 1);
    RK4_CL_stageA(f, &l2, w, dtw,
    		  0.5 * f->dt, sizew, numworkitems,
		  1, source + 1, stage + 1);
    
    // stage 2
    dtfield_CL(f, &l2, 
	       1, stage + 1, source + 2);
    RK4_CL_stageA(f, &l3, w, dtw,
    		  f->dt, sizew, numworkitems,
		  1, source + 2, stage + 2);
    
    // stage 3
    dtfield_CL(f, &l3, 
	       1, stage + 2, source + 3);
    RK4_final_inplace_CL(f, w, &l1, &l2, &l3, 
			 dtw, f->dt, numworkitems,
			 1, source + 3, stage + 3);

    f->tnow += f->dt;
    iter++;
  }
  if(done != NULL)
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, f->dt);
}


// Time integration by a second-order Runge-Kutta algorithm, OpenCL
// version.
void RK2_CL(field *f, double tmax,
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  f->itermax = tmax / f->dt;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int iter = 0;

  cl_int status;
  cl_mem wnp1_cl = clCreateBuffer(f->cli.context,
				  0, // no flags
				  sizeof(double) * f->wsize,
				  NULL, // do not use a host pointer
				  &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Set up kernels
  init_RK2_CL_stage1(f, f->dt, &wnp1_cl);
  init_RK2_CL_stage2(f, f->dt);

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
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, f->dt);

    dtfield_CL(f, &(f->wn_cl), 1, &stage2, &source1);
    RK2_CL_stage1(f, f->wsize, 1, &source1, &stage1);

    f->tnow += 0.5 * f->dt;

    dtfield_CL(f, &wnp1_cl, 1, &stage1, &source2);
    RK2_CL_stage2(f, f->wsize, 1, &source2, &stage2);

    f->tnow += 0.5 * f->dt;
    iter++;
  }
  if(done != NULL) 
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, f->dt);
}

void show_cl_timing(field *f)
{
  cl_ulong total = f->zbuf_time + f->minter_time + f->vol_time 
    + f->mass_time + f->rk_time;
  double N = 1.0 / total;

  cl_ulong ns;

  ns = f->zbuf_time;
  printf("set_buf_to_zero_cl time:      %f%% \t%luns \t%fs\n", 
	 ns*N, ns, 1e-9 * ns);

  ns = f->minter_time;
  printf("DGMacroCellInterface_CL time: %f%% \t%luns \t%fs\n", 
	 ns*N, ns, 1e-9 * ns);

  ns = f->vol_time;
  printf("DGVolume_CL time:             %f%% \t%luns \t%fs\n", 
	 ns*N, ns, 1e-9 * ns);

  ns = f->mass_time;
  printf("DGMass_CL time:               %f%% \t%luns \t%fs\n", 
	 ns*N, ns, 1e-9 * ns);

  ns = f->rk_time;
  printf("rk time:                      %f%% \t%luns \t%fs\n", 
	 ns*N, ns, 1e-9 * ns);

}

#endif // _WITH_OPENCL

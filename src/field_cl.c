#include "clinfo.h"
#include "clutils.h"
#include "field.h"
#include "field_cl.h"
#include "gyro.h"
#ifdef _WITH_PYTHONLIBS
#include "python_call.h"
#endif
#include <assert.h>
#include <string.h>
#include <sys/time.h>



#ifdef _WITH_OPENCL

void DGCharge_CL(Simulation *simu, cl_mem *w_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done);


void CopyfieldtoCPU(Simulation *simu)
{
  cl_int status;

  //DGCharge_CL(simu,&simu->w_cl,0,0,0);

  
  // Ensure that all the buffers are mapped
  status = clFinish(simu->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  void *chkptr;
  chkptr = clEnqueueMapBuffer(simu->cli.commandqueue,
			      simu->dtw_cl, // buffer to copy from
			      CL_TRUE, // block until the buffer is available
			      CL_MAP_READ, // we just want to see the results
			      0, // offset
			      simu->wsize * sizeof(schnaps_real), // buffersize
			      0, NULL, NULL, // events management
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == simu->dtw);

  chkptr = clEnqueueMapBuffer(simu->cli.commandqueue,
			      simu->w_cl, // buffer to copy from
			      CL_TRUE, // block until the buffer is available
			      CL_MAP_READ, // we just want to see the results
			      0, // offset
			      simu->wsize * sizeof(schnaps_real), // buffersize
			      0, NULL, NULL, // events management
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  assert(chkptr == simu->w);

  status = clFinish(simu->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
}


void CopyCPU2GPU(Simulation *simu)
{
  cl_int status;

  status = clFinish(simu->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  void *chkptr_dtwcl;
  void *chkptr_wcl;
  status = clEnqueueUnmapMemObject(simu->cli.commandqueue,
				simu->dtw_cl, // buffer to copy from
				chkptr_dtwcl,
				0, NULL, NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));

  status = clEnqueueUnmapMemObject(simu->cli.commandqueue,
				simu->w_cl, // buffer to copy from
				chkptr_wcl,
				0, NULL, NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));

  status = clFinish(simu->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
}

void set_source_CL(Simulation *simu, char *sourcename_cl)
{
#ifdef _WITH_OPENCL
  simu->use_source_cl = true;
  int namelength = strlen(sourcename_cl);
  if(namelength > 0) {
    simu->sourcename_cl = calloc(sizeof(char), namelength + 1);
    strcat(simu->sourcename_cl, sourcename_cl);
  } else {
    printf("The source name %s is empty!\n", sourcename_cl);
    assert(false);
  }
  printf("%s\n", simu->sourcename_cl);
#endif
}

void init_DGBoundary_CL(Simulation *simu,
			int ieL, int locfaL,
			cl_mem physnodeL_cl,
			cl_mem *w_cl,
			size_t cachesize)
{
  cl_int status;
  cl_kernel kernel = simu->dgboundary;

  unsigned int argnum = 0;

  // __constant int *param,        // 0: interp param
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // real tnow,                  // 1: current time
  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(schnaps_real),
  			  &simu->tnow);
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

  // __global real *w,          // 5: field
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *dtw,        // 6: time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __local real *cache         // 7: local mem
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(schnaps_real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGBoundary_CL(int ifa, Simulation *simu, cl_mem *w_cl,
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  int *param = simu->interp_param;
  cl_int status;

  int ieL =    simu->macromesh.face2elem[4 * ifa];
  int locfaL = simu->macromesh.face2elem[4 * ifa + 1];

  field *f = simu->fd;
  size_t numworkitems = NPGF(f->deg, f->raf, locfaL);
  size_t cachesize = 1; // TODO make use of cache
  init_DGBoundary_CL(simu, ieL, locfaL, simu->physnodes_cl, w_cl, cachesize);
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
				  simu->dgboundary,
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
void init_DGMacroCellInterface_CL(Simulation *simu,
				  int ieL, int ieR, int locfaL, int locfaR,
				  cl_mem physnodeL_cl,
				  cl_mem *w_cl,
				  size_t cachesize)
{
  //printf("loop_init_DGMacroCellInterface_CL\n");
  cl_int status;
  cl_kernel kernel = simu->dginterface;

  // Set kernel arguments
  unsigned int argnum = 0;

  //__constant int *param,        // 0: interp param
  status = clSetKernelArg(kernel,
			 argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
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
  /* status = clSetKernelArg(kernel, */
  /*                         argnum++, */
  /*                         sizeof(cl_mem), */
  /*                         &physnodeR_cl); */
  /* if(status < CL_SUCCESS) printf("%s\n", clErrorString(status)); */
  /* assert(status >= CL_SUCCESS); */

  //__global real *w,          // 7: field
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real *dtw,        // 8: time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__local real *cache         // 9: local mem
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(schnaps_real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGMacroCellInterface_CL(int ifa, Simulation *simu, cl_mem *w_cl,
			     cl_uint nwait, cl_event *wait, cl_event *done)
{
  int *param = simu->interp_param;
  cl_int status;
  cl_kernel kernel = simu->dginterface;
  status = clSetKernelArg(kernel,
                          7,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  int ieL =    simu->macromesh.face2elem[4 * ifa];
  int locfaL = simu->macromesh.face2elem[4 * ifa + 1];
  int ieR =    simu->macromesh.face2elem[4 * ifa + 2];
  int locfaR = simu->macromesh.face2elem[4 * ifa + 3];

  
  field *f = simu->fd;
  size_t numworkitems = NPGF(f->deg, f->raf, locfaL);
  if(ieR >= 0) {
    // Set the remaining loop-dependant kernel arguments
    size_t kernel_cachesize = 1;
    init_DGMacroCellInterface_CL(simu,
				 ieL, ieR, locfaL, locfaR,
				 simu->physnodes_cl,
				 w_cl,
				 kernel_cachesize);

    status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
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
void init_DGMass_CL(Simulation *simu)
{
  cl_int status;
  cl_kernel kernel = simu->dgmass;
  int argnum = 0;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* int ie, // macrocel index */
  // Set in loop on call.
  argnum++;

  /* __constant real* physnode,  // macrocell nodes */
  status = clSetKernelArg(simu->dgmass,           // kernel name
                          argnum++,              // arg num
                          sizeof(cl_mem),
                          &simu->physnodes_cl);     // opencl buffer
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  /* __global real* dtw // time derivative */
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Apply division by the mass matrix OpenCL version
void DGMass_CL(int ie, Simulation *simu,
	       cl_uint nwait, cl_event *wait, cl_event *done)
{

  int *param = simu->interp_param;
  cl_int status;

  init_DGMass_CL(simu);

  //clSetUserEventStatus(simu->clv_mass, CL_COMPLETE);
  //if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  //assert(status >= CL_SUCCESS);

  // Loop on the elements

  status = clSetKernelArg(simu->dgmass,
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

  unsigned int m = simu->interp_param[0];
  unsigned int nreadsdgmass = m;
  unsigned int nmultsdgmass = 53 + 2 * m + 2601;
  simu->flops_mass += numworkitems * nmultsdgmass;
  simu->reads_mass += numworkitems * nreadsdgmass;

  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
				  simu->dgmass,
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

void init_DGFlux_CL(Simulation *simu, int ie, int dim0, cl_mem *w_cl,
		    size_t cachesize)
{
  //printf("DGFlux cachesize:%zu\n", cachesize);

  cl_kernel kernel = simu->dgflux;
  cl_int status;
  int argnum = 0;

  // __constant int *param, // interp param
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //xxx
 
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
                          &simu->physnodes_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //  __global real *w, // field values
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *dtw // time derivative
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __local real* wloc,       // 6: w local memory
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(schnaps_real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void DGFlux_CL(Simulation *simu, int dim0, int ie, cl_mem *w_cl,
	       cl_uint nwait, cl_event *wait, cl_event *done)
{
  // Unpack the parameters
  int *param = simu->interp_param;
  int m = param[0];
  int deg[3] = {param[1], param[2], param[3]};
  int nraf[3] = {param[4], param[5], param[6]};
  size_t cachesize = 1; 

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
#if _CACHEMEMORY
  cachesize =  4 * m * groupsize;
#endif
    init_DGFlux_CL(simu, ie, dim0, w_cl, cachesize);

    unsigned int nreadsdgflux = 4 * m;
    unsigned int nmultsdgflux = 2601 + 92 + 28 * m;
    nmultsdgflux += 3 * m; // Using NUMFLUX = NumFlux
    simu->flops_flux += numworkitems * nmultsdgflux;
    simu->reads_flux += numworkitems * nreadsdgflux;

    // Launch the kernel
    cl_int status;
    status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
				    simu->dgflux,
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
void init_DGVolume_CL(Simulation *simu, cl_mem *w_cl, size_t cachesize)
{
  //printf("DGVolume cachesize:%zu\n", cachesize);

  cl_int status;
  int argnum = 0;
  cl_kernel kernel = simu->dgvolume;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->physnodes_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(schnaps_real) * cachesize,
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
void DGVolume_CL(int ie, Simulation *simu, cl_mem *w_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = simu->dgvolume;
  int *param = simu->interp_param;
  size_t cachesize = 1; 
  cl_int status;
  int m = param[0];
  const int npgc[3] = {param[1] + 1, param[2] + 1, param[3] + 1};
  const int npg = npgc[0] * npgc[1] * npgc[2];
  size_t groupsize = npg;
  // The total work items number is the number of glops in a subcell
  // * number of subcells
  const int nraf[3] = {param[4], param[5], param[6]};
  size_t numworkitems = nraf[0] * nraf[1] * nraf[2] * groupsize;

  unsigned int nreadsdgvol = 2 * m; // read m from w, write m to dtw
  unsigned int nmultsdgvol = 2601 + (npgc[0] + npgc[1] + npgc[2]) * 6 * m;
  nmultsdgvol += (npgc[0] + npgc[1] + npgc[2]) * 54;
  // Using NUMFLUX = NumFlux (3 * m multiplies):
  nmultsdgvol += (npgc[0] + npgc[1] + npgc[2]) * 6 * m;
 #if _CACHEMEMORY
  cachesize = 2 * groupsize * m;
#endif
  //printf("cache memory = :%zu \n",cachesize);
  init_DGVolume_CL(simu, w_cl, cachesize);

  simu->flops_vol += numworkitems * nmultsdgvol;
  simu->reads_vol += numworkitems * nreadsdgvol;

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
  //printf("groupsize=%zd numworkitems=%zd\n", groupsize, numworkitems);
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
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

// Set kernel argument for DGCharge_CL
void init_DGCharge_CL(Simulation *simu, cl_mem *w_cl)
{
  //printf("DGVolume cachesize:%zu\n", cachesize);

  cl_int status;
  int argnum = 0;
  cl_kernel kernel = simu->dgcharge;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Set kernel argument for DGVolume_CL
void init_DGSource_CL(Simulation *simu, cl_mem *w_cl, size_t cachesize)
{
  //printf("DGVolume cachesize:%zu\n", cachesize);

  cl_int status;
  int argnum = 0;
  cl_kernel kernel = simu->dgsource;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->physnodes_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
  			  argnum++,
  			  sizeof(schnaps_real),
  			  &simu->tnow);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(schnaps_real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// compute charge in opencl
void DGCharge_CL(Simulation *simu, cl_mem *w_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = simu->dgcharge;
  int *param = simu->interp_param;

  cl_int status;
  int m = param[0];
  size_t groupsize = (param[1] + 1) * (param[2] + 1) * (param[3] + 1);
  // The total work items number is the number of glops in a subcell
  // times the number of subcells
  size_t numworkitems = param[4] * param[5] * param[6] * groupsize;
  __constant KineticData *kd = &schnaps_kinetic_data; //xxx _constant
  init_DGCharge_CL(simu, w_cl);
  void* chkptr_rho[simu->macromesh.nbelems];
  void* chkptr_phi[simu->macromesh.nbelems];
  void* chkptr_elec[simu->macromesh.nbelems];
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    // Loop on the elements
    status = clSetKernelArg(kernel,
			    1,
			    sizeof(int),
			    &ie);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    // The groupsize is the number of glops in a subcell
    size_t groupsize = (param[1] + 1)* (param[2] + 1)*(param[3] + 1);
    // The total work items number is the number of glops in a subcell
    // * number of subcells
    size_t numworkitems = param[4] * param[5] * param[6] * groupsize;
    //printf("groupsize=%zd numworkitems=%zd\n", groupsize, numworkitems);
    status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
    				    kernel,
    				    1,
    				    NULL,
    				    &numworkitems,
    				    &groupsize,
    				    nwait,
    				    wait,
    				    //NULL);
                                    simu->clv_charge + ie); //bug here
    /* 	CopyfieldtoCPU(simu); */
    /* if(ie == 0) PlotFields(kd->index_rho,(1==0),simu,"rho","rho1.msh"); */
    /* assert(1==0); */
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    // accès à  la charge en Read Only -> cpu
    int woffset = ie * m * NPG(param + 1, param + 4);
    field * f = simu->fd + ie; 
    int start = f->varindex(param + 1, param + 4, m, 0, kd->index_rho) + woffset;
    int end = f->varindex(param + 1, param + 4, m, 0, kd->index_rho + 1) + woffset;
    size_t soffset = start * sizeof(schnaps_real);
    chkptr_rho[ie] = clEnqueueMapBuffer(simu->cli.commandqueue,
				*w_cl, // buffer to copy from
				CL_TRUE, // block until the buffer is available
				CL_MAP_READ, // we just want to see the results
				soffset, // offset
				(end - start) * sizeof(schnaps_real), // buffersize
				1, simu->clv_charge + ie, NULL, // events management
				&status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    assert(chkptr_rho[ie] == simu->w + start);
    // accès à  phi en write Only -> cpu
    start = f->varindex(param + 1, param + 4, m, 0, kd->index_phi) + woffset;
    end = f->varindex(param + 1, param + 4, m, 0, kd->index_phi + 1) + woffset;
    soffset = start * sizeof(schnaps_real);
    chkptr_phi[ie] = clEnqueueMapBuffer(simu->cli.commandqueue,
				*w_cl, // buffer to copy from
				CL_TRUE, // block until the buffer is available
				CL_MAP_WRITE, // we just want to see the results
				soffset, // offset
				(end - start) * sizeof(schnaps_real), // buffersize
				0, NULL, NULL, // events management
				&status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    assert(chkptr_phi[ie] == simu->w + start);
    // copier ex ey ez en Write Only
    start = f->varindex(param + 1, param + 4, m, 0, kd->index_ex) + woffset;
    end = f->varindex(param + 1, param + 4, m, 0, kd->index_ez + 1) + woffset;
    soffset = start * sizeof(schnaps_real);
    chkptr_elec[ie] = clEnqueueMapBuffer(simu->cli.commandqueue,
				*w_cl, // buffer to copy from
				CL_TRUE, // block until the buffer is available
				CL_MAP_WRITE, // we just want to see the results
				soffset, // offset
				(end - start) * sizeof(schnaps_real), // buffersize
				0, NULL, NULL, // events management
				&status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    assert(chkptr_elec[ie] == simu->w + start);
  }
  // résoudre poisson
  //status = clFinish(simu->cli.commandqueue);
  UpdateGyroPoisson(simu);
 
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++) {

    // libère accès à  la charge en Read Only -> cpu
 
    status = clEnqueueUnmapMemObject (simu->cli.commandqueue,
				*w_cl, // buffer to copy from
				chkptr_rho[ie],
				0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));

    status = clEnqueueUnmapMemObject (simu->cli.commandqueue,
				*w_cl, // buffer to copy from
				chkptr_phi[ie],
				0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));

    status = clEnqueueUnmapMemObject (simu->cli.commandqueue,
				*w_cl, // buffer to copy from
				chkptr_elec[ie],
				0, NULL, done);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    
  }

  
  status = clFinish(simu->cli.commandqueue);
}




// Apply division by the mass matrix OpenCL version
void DGSource_CL(int ie, Simulation *simu, cl_mem *w_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = simu->dgsource;
  int *param = simu->interp_param;

  cl_int status;
  int m = param[0];
  size_t groupsize = (param[1] + 1) * (param[2] + 1) * (param[3] + 1);
  // The total work items number is the number of glops in a subcell
  // times the number of subcells
  size_t numworkitems = param[4] * param[5] * param[6] * groupsize;

  /* unsigned int nreadsdgvol = 2 * m; // read m from w, write m to dtw */
  /* unsigned int nmultsdgvol = 1296 + m; */
  /* nmultsdgvol += 3 * m; // Using NUMFLUX = NumFlux */

  /* simu->nmults += numworkitems * nmultsdgvol; */
  /* simu->nreads += numworkitems * nreadsdgvol; */

  init_DGSource_CL(simu, w_cl, 2 * groupsize * m);


  // Loop on the elements
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
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
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

void set_buf_to_zero_cl(cl_mem *buf, int size, Simulation *simu,
			cl_uint nwait, cl_event *wait,  cl_event *done)
{
  cl_int status;

  cl_kernel kernel = simu->zero_buf;

  // associates the param buffer to the 0th kernel argument
  status = clSetKernelArg(kernel,
                          0,
                          sizeof(cl_mem),
                          buf);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t numworkitems = size;
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
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
void dtfield_CL(Simulation *simu, cl_mem *w_cl,
		cl_uint nwait, cl_event *wait, cl_event *done)
{

  set_buf_to_zero_cl(&simu->dtw_cl, simu->wsize, simu,
  		     nwait, wait, &simu->clv_zbuf);

   
  if (schnaps_ocl_getcharge) {
    DGCharge_CL(simu, &simu->w_cl, nwait, wait, simu->clv_poisson);
    __constant KineticData *kd = &schnaps_kinetic_data; //xxx _constant
    clWaitForEvents(1, simu->clv_poisson);
  }
  

  // Macrocell interfaces must be launched serially
  const int ninterfaces = simu->macromesh.nmacrointerfaces;
  
  for(int i = 0; i < ninterfaces; ++i) {
    int ifa = simu->macromesh.macrointerface[i];
    DGMacroCellInterface_CL(ifa, simu, w_cl,
    			    1,
    			    i == 0 ? &simu->clv_zbuf : simu->clv_mci + i - 1,
    			    simu->clv_mci + i);

  }
  cl_event *startboundary = ninterfaces > 0 ?
    simu->clv_mci + ninterfaces - 1: &simu->clv_zbuf;

  // Boundary terms may also need to be launched serially.
  const int nboundaryfaces = simu->macromesh.nboundaryfaces;
  for(int i = 0; i < nboundaryfaces; ++i) {
    int ifa = simu->macromesh.boundaryface[i];
    DGBoundary_CL(ifa, simu, w_cl,
		  1,
		  i == 0 ?  startboundary : simu->clv_boundary + i - 1,
		  simu->clv_boundary + i);
  }

  // If there are no interfaces and no boundaries, just wait for zerobuf.
  cl_event *startmacroloop;
  if(ninterfaces == 0 && nboundaryfaces == 0) {
    startmacroloop = &simu->clv_zbuf;
  } else {
    if(nboundaryfaces > 0) {
      startmacroloop = simu->clv_boundary + nboundaryfaces - 1;
    } else {
      startmacroloop = simu->clv_mci + ninterfaces - 1;
    }
  }

  cl_event *fluxdone = simu->clv_flux2;
  unsigned int ndim = 3;
  if(simu->macromesh.is2d) {
    ndim = 2;
    fluxdone = simu->clv_flux1;
  }
  if(simu->macromesh.is1d) {
    ndim = 1;
    fluxdone = simu->clv_flux0;
  }

  cl_event *dtfielddone = simu->use_source_cl ? simu->clv_source : simu->clv_mass;

  // The kernels for the intra-macrocell computations can be launched
  // in parallel between macrocells.
  const int nmacro = simu->macromesh.nbelems;
  for(int ie = 0; ie < nmacro; ++ie) {
    DGFlux_CL(simu, 0, ie, w_cl,
	      1, startmacroloop,
	      simu->clv_flux0 + ie);

    if(ndim > 1) {
      DGFlux_CL(simu, 1, ie, w_cl,
		1, simu->clv_flux0 + ie,
		simu->clv_flux1 + ie);
    }

    if(ndim > 2) {
      DGFlux_CL(simu, 2, ie, w_cl,
		1, simu->clv_flux1 + ie,
		simu->clv_flux2 + ie);
    }

    DGVolume_CL(ie, simu, w_cl,
  		1,
  		fluxdone + ie,
  		simu->clv_volume + ie);

    DGMass_CL(ie, simu,
  	      1,
  	      simu->clv_volume + ie,
  	      simu->clv_mass + ie);

    if(simu->use_source_cl) {
      DGSource_CL(ie, simu, w_cl,
		  1,
		  simu->clv_mass + ie,
		  simu->clv_source + ie);
     }
  }
  clWaitForEvents(nmacro, dtfielddone);

  
  // Add times for sources after everything is finished
  simu->zbuf_time += clv_duration(simu->clv_zbuf);

  for(int i = 0; i < ninterfaces; ++i)
    simu->minter_time += clv_duration(simu->clv_mci[i]);

  for(int i = 0; i < nboundaryfaces; ++i)
    simu->boundary_time += clv_duration(simu->clv_boundary[i]);

  for(int ie = 0; ie < nmacro; ++ie) {
    if(simu->use_source_cl)
      simu->source_time += clv_duration(simu->clv_source[ie]);
    simu->flux_time += clv_duration(simu->clv_flux0[ie]);
    if(ndim > 1)
      simu->flux_time += clv_duration(simu->clv_flux1[ie]);
    if(ndim > 2)
      simu->flux_time += clv_duration(simu->clv_flux2[ie]);
    simu->vol_time += clv_duration(simu->clv_volume[ie]);
    simu->mass_time += clv_duration(simu->clv_mass[ie]);
  }

  if(done != NULL)
    clSetUserEventStatus(*done, CL_COMPLETE);

 
}

// Set kernel arguments for first stage of RK2
void init_RK2_CL_stage1(Simulation *simu, const schnaps_real dt, cl_mem *wp1_cl)
{
  cl_kernel kernel = simu->RK_out_CL;
  cl_int status;
  int argnum = 0;

  //__global real *wp1 // field values
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          wp1_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *w,
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &simu->w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real* dtw // time derivative
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //real dt, // time step for the stage
  schnaps_real halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(schnaps_real),
			  &halfdt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  
}

// Launch first stage of RK2 integration
void RK2_CL_stage1(Simulation *simu, size_t numworkitems,
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
				  simu->RK_out_CL,
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

// Set kernel arguments for second stage of RK2
void init_RK2_CL_stage2(Simulation *simu, const schnaps_real dt)
{
  cl_kernel kernel = simu->RK_in_CL;
  cl_int status;
  int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &simu->w_cl);
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &simu->dtw_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(schnaps_real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Launch second stage of RK2 integration
void RK2_CL_stage2(Simulation *simu, size_t numworkitems,
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_int status;
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
				  simu->RK_in_CL,
				  1, NULL,
				  &numworkitems,
				  NULL,
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void RK4_CL_stageA(Simulation *simu,
		   cl_mem *wp1, cl_mem *w, cl_mem *dtw,
		   const schnaps_real dt, const int sizew, size_t numworkitems,
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  // l_1 = w_n + 0.5dt * S(w_n, t_0)

  cl_kernel kernel = simu->RK_out_CL;
  cl_int status;
  int argnum = 0;

  //__global real *wp1 // field values
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          wp1);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // __global real *w,
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          w);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //__global real *dtw // time derivative
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          dtw);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //real dt,
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(schnaps_real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
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

void RK4_final_inplace_CL(Simulation *simu,
			  cl_mem *w_cl, cl_mem *l1, cl_mem *l2, cl_mem *l3,
			  cl_mem *dtw_cl, const schnaps_real dt,
			  const size_t numworkitems,
			  cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = simu->RK4_final_stage;
  cl_int status;
  int argnum = 0;
 __constant KineticData *kd = &schnaps_kinetic_data; //xxx _constant
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
   schnaps_real halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(schnaps_real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  status = clEnqueueNDRangeKernel(simu->cli.commandqueue,
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
void RK4_CL(Simulation *simu, schnaps_real tmax, schnaps_real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done)
{
  simu->dt = Get_Dt_RK(simu);
  if(dt <= 0) dt = simu->dt;
  assert(simu->tnow == 0);
  //KineticData *kd = &schnaps_kinetic_data;
  simu->itermax_rk = tmax / dt;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  field *f = simu->fd;
  int sizew = simu->macromesh.nbelems * f->model.m * NPG(f->deg, f->raf);
  int iter = 0;
  cl_int status;
  cl_mem l1 = clCreateBuffer(simu->cli.context,
			     0,
			     sizeof(schnaps_real) * simu->wsize,
			     NULL,
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  cl_mem l2 = clCreateBuffer(simu->cli.context,
			     0,
			     sizeof(schnaps_real) * simu->wsize,
			     NULL,
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  cl_mem l3 = clCreateBuffer(simu->cli.context,
			     0,
			     sizeof(schnaps_real) * simu->wsize,
			     NULL,
			     &status);
  /* int size_diags;
  // size_diags = simu->nb_diags * simu->itermax_rk;
  //allouer pour Diagnostic
   if(simu->nb_diags != 0)
    simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
  */
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t numworkitems = simu->wsize;

  cl_mem *w = &simu->w_cl;
  cl_mem *dtw = &simu->dtw_cl;
  int nstages = 4;
  cl_event source[nstages];
  cl_event stage[nstages];
  for(int i = 0; i < nstages; ++i) {
    source[i] = clCreateUserEvent(simu->cli.context, &status);
    stage[i] = clCreateUserEvent(simu->cli.context, &status);
  }

  clWaitForEvents(nwait, wait);

  printf("Starting RK4_CL\n");

  status = clSetUserEventStatus(stage[3], CL_COMPLETE);



  struct timeval t_start;
  struct timeval t_end;
  gettimeofday(&t_start, NULL);
  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    // stage 0
    dtfield_CL(simu, w,
	       1, stage +3, source);
    RK4_CL_stageA(simu, &l1, w, dtw,
    		  0.5 * dt, sizew, numworkitems,
		  1, source, stage);

    //clFinish(simu->cli.commandqueue); //--
    //	assert(1==0);		 
    // stage 1
    dtfield_CL(simu, &l1, 1, stage, source + 1);
    RK4_CL_stageA(simu, &l2, w, dtw,
    		  0.5 * dt, sizew, numworkitems,
		  1, source + 1, stage + 1);
	
		  

    // stage 2
    dtfield_CL(simu, &l2,
	       1, stage + 1, source + 2);
    RK4_CL_stageA(simu, &l3, w, dtw,
    		  dt, sizew, numworkitems,
		  1, source + 2, stage + 2);

	/* CopyfieldtoCPU(simu);	 */
	/* PlotFields(schnaps_kinetic_data.index_rho,(1==0),simu,"sol","dens.msh"); */

			 
    // stage 3
    dtfield_CL(simu, &l3,
	       1, stage + 2, source + 3);
	       		  
    RK4_final_inplace_CL(simu, w, &l1, &l2, &l3,
			 dtw, dt, numworkitems,
			 1, source + 3, stage + 3);

	
	
    /*  //copy les valeur de GPU2CPU pour calculer les valeurs de Diagnostics
      CopyfieldtoCPU(simu);
    //copy les valeur dans Diagnostics
    if(simu->update_after_rk != NULL){
    simu->update_after_rk(simu, simu->w);
    }
    //unmap pour retourner les valeur en GPU
      CopyCPU2GPU(simu);
    */
      simu->tnow += dt;

    for(int i = 0; i < nstages; ++i)
      simu->rk_time += clv_duration(stage[i]);

    iter++;
  }
  gettimeofday(&t_end, NULL);

  if(done != NULL)
    status = clSetUserEventStatus(*done, CL_COMPLETE);

 schnaps_real rkseconds = (t_end.tv_sec - t_start.tv_sec) * 1.0 // seconds
    + (t_end.tv_usec - t_start.tv_usec) * 1e-6; // microseconds
  printf("\nTotal RK time (s):\n%f\n", rkseconds);
  printf("\nTotal RK time per time-step (s):\n%f\n", rkseconds / iter );


  for(int i = 0; i < nstages; ++i) {
    clReleaseEvent(source[i]);
    clReleaseEvent(stage[i]);
  }

  //printf("\nt=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
     
  
}


// Time integration by a second-order Runge-Kutta algorithm, OpenCL
// version.
void RK2_CL(Simulation *simu, schnaps_real tmax, schnaps_real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done)
{


  simu->dt = Get_Dt_RK(simu);

  if(dt <= 0)
    dt = simu->dt;

  assert(simu->tnow == 0);
  simu->itermax_rk = tmax / dt;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  cl_int status;
  cl_mem wp1_cl = clCreateBuffer(simu->cli.context,
				  0, // no flags
				  sizeof(schnaps_real) * simu->wsize,
				  NULL, // do not use a host pointer
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Set up kernels
  init_RK2_CL_stage1(simu, dt, &wp1_cl);
  init_RK2_CL_stage2(simu, dt);

  cl_event source1 = clCreateUserEvent(simu->cli.context, &status);
  cl_event stage1 =  clCreateUserEvent(simu->cli.context, &status);
  cl_event source2 = clCreateUserEvent(simu->cli.context, &status);
  cl_event stage2 =  clCreateUserEvent(simu->cli.context, &status);

  status = clSetUserEventStatus(stage2, CL_COMPLETE);

  if(nwait > 0)
    clWaitForEvents(nwait, wait);

  printf("Starting RK2_CL\n");

  struct timeval t_start;
  struct timeval t_end;
  gettimeofday(&t_start, NULL);
  while(simu->tnow < tmax) {
    //printf("iter: %d\n", iter);
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    dtfield_CL(simu, &simu->w_cl, 1, &stage2, &source1);
    RK2_CL_stage1(simu, simu->wsize, 1, &source1, &stage1);

    simu->tnow += 0.5 * dt;

    dtfield_CL(simu, &wp1_cl, 1, &stage1, &source2);
    RK2_CL_stage2(simu, simu->wsize, 1, &source2, &stage2);

    simu->tnow += 0.5 * dt;

    simu->rk_time += clv_duration(stage1);
    simu->rk_time += clv_duration(stage2);
    iter++;
  }
  status = clFinish(simu->cli.commandqueue);
  gettimeofday(&t_end, NULL);
  if(done != NULL)
    status = clSetUserEventStatus(*done, CL_COMPLETE);

  clReleaseEvent(stage2);
  clReleaseEvent(stage1);

  schnaps_real rkseconds = (t_end.tv_sec - t_start.tv_sec) * 1.0 // seconds
    + (t_end.tv_usec - t_start.tv_usec) * 1e-6; // microseconds
  printf("\nTotal RK time (s):\n%f\n", rkseconds);
  printf("\nTotal RK time per time-step (s):\n%f\n", rkseconds / iter );

  printf("\nt=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);


}

void show_cl_timing(Simulation *simu)
{
  printf("\n");
  printf("Device characteristics:\n");
  printf("\tname:\t%s\n", simu->cli.devicename);
  schnaps_real dev_gflops = cl_dev_gflops(simu->cli.devicename);
  schnaps_real dev_bwidth = cl_dev_bwidth(simu->cli.devicename);
  printf("\tgflops:   \t%f\n", dev_gflops);
  printf("\tbandwidth:\t%f\n", dev_bwidth);

  printf("\n");
  printf("Roofline counts:\n");
  unsigned long int flops_total = simu->flops_vol + simu->flops_flux + simu->flops_mass;
  unsigned long int reads_total = simu->reads_vol + simu->reads_flux + simu->reads_mass;
  printf("Number of real multiplies in kernels:               %lu\n",
	 flops_total);
  printf("Number of reads/writes of reals from global memory: %lu\n",
	 reads_total);
  //printf("NB: real multiplies assumes the flux involved 3m multiplies.\n");
  printf("Terms included in roofline: volume, flux, and mass.\n");

  cl_ulong roofline_time_ns = simu->vol_time + simu->flux_time + simu->mass_time;
  schnaps_real roofline_time_s = 1e-9 * roofline_time_ns;
  schnaps_real roofline_flops = flops_total / roofline_time_s;
  schnaps_real roofline_bw = sizeof(schnaps_real) * reads_total / roofline_time_s;

  // Volume terms
  schnaps_real vol_time_s = 1e-9 * simu->vol_time;
  schnaps_real vol_flops = simu->flops_vol / vol_time_s;
  schnaps_real vol_bw = sizeof(schnaps_real) * simu->reads_vol / vol_time_s;
  printf("DGVol:  GFLOP/s: %f\tbandwidth (GB/s): %f\n",
	 1e-9 * vol_flops, 1e-9 * vol_bw);

  // Flux terms
  schnaps_real flux_time_s = 1e-9 * simu->flux_time;
  schnaps_real flux_flops = simu->flops_flux / flux_time_s;
  schnaps_real flux_bw = sizeof(schnaps_real) * simu->reads_flux / flux_time_s;
  printf("DGFlux: GFLOP/s: %f\tbandwidth (GB/s): %f\n",
	 1e-9 * flux_flops, 1e-9 * flux_bw);

  // Mass terms
  schnaps_real mass_time_s = 1e-9 * simu->mass_time;
  schnaps_real mass_flops = simu->flops_mass / mass_time_s;
  schnaps_real mass_bw = sizeof(schnaps_real) * simu->reads_mass / mass_time_s;
  printf("DGMass: GFLOP/s: %f\tbandwidth (GB/s): %f\n",
	 1e-9 * mass_flops, 1e-9 * mass_bw);

  printf("Total:  GFLOP/s: %f\tbandwidth (GB/s): %f\n",
	 1e-9 * roofline_flops, 1e-9 * roofline_bw);

  printf("\n");

  printf("Kernel execution times:\n");
  cl_ulong total = 0;
  total += simu->zbuf_time;
  total += simu->minter_time;
  total += simu->boundary_time;
  total += simu->vol_time;
  total += simu->mass_time;
  total += simu->source_time;
  total += simu->rk_time;
  total += simu->flux_time;

  schnaps_real N = 100.0 / total;

  cl_ulong ns;

  ns = simu->zbuf_time;
  printf("set_buf_to_zero_cl time:      %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->minter_time;
  printf("DGMacroCellInterface_CL time: %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->boundary_time;
  printf("DGBoundary time:              %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->vol_time;
  printf("DGVolume_CL time:             %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->flux_time;
  printf("DGFlux_CL time:               %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->mass_time;
  printf("DGMass_CL time:               %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->source_time;
  printf("DGSource_CL time:             %f%% \t%luns \t%fs\n",
	 ns * N, (unsigned long) ns, 1e-9 * ns);

  ns = simu->rk_time;
  printf("rk time:                      %f%% \t%luns \t%fs\n",
	 ns*N, (unsigned long) ns, 1e-9 * ns);

  printf("\n");

  ns = total;
  schnaps_real total_time = 1e-9 * ns;
  printf("total time:                   %f%% \t%luns \t%fs\n",
	 ns*N, (unsigned long) ns, total_time);

  printf("\n");
  print_kernel_perf(dev_gflops, dev_bwidth,
		    flops_total, reads_total, total);
  printf("\n");
}


void init_field_cl(Simulation *simu)
{
  InitCLInfo(&simu->cli, nplatform_cl, ndevice_cl);
  cl_int status;

  simu->w_cl = clCreateBuffer(simu->cli.context,
			    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
			    sizeof(schnaps_real) * simu->wsize,
			    simu->w,
			    &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->dtw_cl = clCreateBuffer(simu->cli.context,
			      CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
			      sizeof(schnaps_real) * simu->wsize,
			      simu->dtw,
			      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  field *f = simu->fd;
  /* simu->interp_param[0]=f->model.m; */

  /* 		       f->deg[0], f->deg[1], f->deg[2], */
  /* 		       f->raf[0], f->raf[1], f->raf[2] }; */

  simu->param_cl = clCreateBuffer(simu->cli.context,
			       CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
			       sizeof(int) * 7,
			       simu->interp_param,
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  
  // Allocate one physnode buffer
  //f->physnode = calloc(60, sizeof(real));
  simu->physnode_cl = clCreateBuffer(simu->cli.context,
				  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				  sizeof(schnaps_real) * 60,
				  f->physnode,
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Allocate and fill buffer for all macrocell geometries.
  const int nmacro = simu->macromesh.nbelems;
  const size_t buf_size = sizeof(schnaps_real) * 60 * nmacro;
  simu->physnodes_cl = clCreateBuffer(simu->cli.context,
				   CL_MEM_READ_ONLY,
				   buf_size,
				   NULL,
				   &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  set_physnodes_cl(simu);

  // Allocate one physnode buffer for R macrocell
  simu->physnodeR = calloc(60, sizeof(schnaps_real));
  simu->physnodeR_cl = clCreateBuffer(simu->cli.context,
				   CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,
				   sizeof(schnaps_real) * 60,
				   f->physnode,
				   &status);

  // Program compilation
    char *strprog;
#ifdef _WITH_PYTHONLIBS
    Py_ssize_t ocl_buffer_size = 0;
    call_py_opencl_src_extract(&strprog,&ocl_buffer_size);
   
    // GetOpenCLCode();
    // ReadFile("schnaps.cl", &strprog);
    // printf("%s\n",strprog);
#else 
  GetOpenCLCode();
  ReadFile("schnaps.cl", &strprog);
#endif
  printf("\t%s\n", numflux_cl_name);
  //printf("\t%s\n", strprog);

  char buf[1000];

  // Set simulation parameters
  sprintf(buf, " -D _M=%d", simu->interp_param[0]);
  strcat(cl_buildoptions, buf);

  // If the source term is set (via set_source_CL) then add it to the
  // buildoptions and compile using the new buildoptions.
  if(simu->use_source_cl) {
    char *temp;
    int len0 = strlen(cl_buildoptions);
    char *D_SOURCE_FUNC = " -D_SOURCE_FUNC=";
    int len1 = strlen(D_SOURCE_FUNC);
    int len2 = strlen(simu->sourcename_cl);
    temp = calloc(sizeof(char), len0 + len1 + len2 + 3);
    printf("lens=%d %d %d \n",len0,len1,len2);
    strcat(temp, cl_buildoptions);
    strcat(temp, D_SOURCE_FUNC);
    strcat(temp, simu->sourcename_cl);
    strcat(temp, " ");
    BuildKernels(&simu->cli, strprog, temp);
  } else {
    printf("No source term\n");
    BuildKernels(&simu->cli, strprog, cl_buildoptions);
  }

  simu->dgmass = clCreateKernel(simu->cli.program,
			     "DGMass",
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->dgflux = clCreateKernel(simu->cli.program,
			     "DGFlux",
			     &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  
  simu->dgvolume = clCreateKernel(simu->cli.program,
			       "DGVolume",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->dgsource = clCreateKernel(simu->cli.program,
			       "DGSource",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

 if (schnaps_ocl_getcharge) {
  simu->dgcharge = clCreateKernel(simu->cli.program,
			       "DGCharge",
			       &status);
}
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->dginterface = clCreateKernel(simu->cli.program,
				  "DGMacroCellInterface",
				  &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->dgboundary = clCreateKernel(simu->cli.program,
				 "DGBoundary",
				 &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->RK_out_CL = clCreateKernel(simu->cli.program,
				"RK_out_CL",
				&status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->RK4_final_stage = clCreateKernel(simu->cli.program,
				      "RK4_final_stage",
				      &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->RK_in_CL = clCreateKernel(simu->cli.program,
			       "RK_in_CL",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  simu->zero_buf = clCreateKernel(simu->cli.program,
			       "set_buffer_to_zero",
			       &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // Initialize events. // FIXME: free on exit
  simu->clv_zbuf = clCreateUserEvent(simu->cli.context, &status);

  const int ninterfaces = simu->macromesh.nmacrointerfaces;
  if(ninterfaces > 0) {
    simu->clv_mci = calloc(ninterfaces, sizeof(cl_event));
    for(int ifa = 0; ifa < ninterfaces; ++ifa)
      simu->clv_mci[ifa] = clCreateUserEvent(simu->cli.context, &status);
  }

  const int nbound = simu->macromesh.nboundaryfaces;
  if(nbound > 0) {
    simu->clv_boundary = calloc(nbound, sizeof(cl_event));
    for(int ifa = 0; ifa < nbound; ++ifa)
      simu->clv_boundary[ifa] = clCreateUserEvent(simu->cli.context, &status);
  }

  simu->clv_mass = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie)
    simu->clv_mass[ie] = clCreateUserEvent(simu->cli.context, &status);

  simu->clv_flux0 = calloc(nmacro, sizeof(cl_event));
  simu->clv_flux1 = calloc(nmacro, sizeof(cl_event));
  simu->clv_flux2 = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    simu->clv_flux0[ie] = clCreateUserEvent(simu->cli.context, &status);
    simu->clv_flux1[ie] = clCreateUserEvent(simu->cli.context, &status);
    simu->clv_flux2[ie] = clCreateUserEvent(simu->cli.context, &status);
  }

  simu->clv_volume = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    simu->clv_volume[ie] = clCreateUserEvent(simu->cli.context, &status);
  }

  simu->clv_source = calloc(nmacro, sizeof(cl_event));
  for(int ie = 0; ie < nmacro; ++ie) {
    simu->clv_source[ie] = clCreateUserEvent(simu->cli.context, &status);
  }

 
  if (schnaps_ocl_getcharge) {
    simu->clv_charge = calloc(nmacro, sizeof(cl_event)); 
 	for(int ie = 0; ie < nmacro; ++ie) {
    	simu->clv_charge[ie] = clCreateUserEvent(simu->cli.context, &status);
  	}
  }

  // Set timers to zero
  simu->zbuf_time = 0;
  simu->mass_time = 0;
  simu->vol_time = 0;
  simu->flux_time = 0;
  simu->minter_time = 0;
  simu->boundary_time = 0;
  simu->source_time = 0;
  simu->rk_time = 0;

  // Set roofline counts to zero
  simu->flops_vol = 0;
  simu->flops_flux = 0;
  simu->flops_mass = 0;
  simu->reads_vol = 0;
  simu->reads_flux = 0;
  simu->reads_mass = 0;


  if (schnaps_ocl_getcharge) {
    status = clFinish(simu->cli.commandqueue);

    // just check that the good varindex is chosen
    field * f = simu->fd + 0; 
    int *param = simu->interp_param;
    int m = param[0];
    // ipg = 0
    int i1 = f->varindex(param + 1, param + 4, m, 0, 0);
    // ipg = 1
    int i2 = f->varindex(param + 1, param + 4, m, 1, 0);
    schnaps_real *addr1 = &simu->w[i1];
    schnaps_real *addr2 = &simu->w[i2];
    assert(addr1 + 1 == addr2);
    /////////////////////////////
    /* cl_event *wait = NULL; */
    /* DGCharge_CL(simu, &simu->w_cl, 0, wait, simu->clv_poisson); */
    /* status = clFinish(simu->cli.commandqueue); */
  }
  
}

#endif // _WITH_OPENCL

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
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);
  assert(chkptr == f->wn);
  status = clFinish(f->cli.commandqueue);
}

// Set OpenCL kernel arguments for DGMacroCellInterface
void initDGMacroCellInterface_CL(field *f, 
				 cl_mem physnodeL_cl, cl_mem physnodeR_cl)
{
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
  assert(status == CL_SUCCESS);

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
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &physnodeR_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // wn_cl
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
}

// Set the loop-dependant kernel arguments for DGMacroCellInterface_CL
void loop_initDGMacroCellInterface_CL(field *f, 
				      int ieL, int ieR, int locfaL, int locfaR)
{
  cl_int status;
  cl_kernel kernel = f->dginterface;

  // Set kernel arguments
  unsigned int argnum = 1;

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &(f->tnow));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &ieL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &ieR);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaR);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);

  /* int ie, // macrocel index */
  // Set in loop on call.
  argnum++;
      
  /* __constant double* physnode,  // macrocell nodes */
  status = clSetKernelArg(f->dgmass,           // kernel name
                          argnum++,              // arg num
                          sizeof(cl_mem),
                          &f->physnode_cl);     // opencl buffer
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  /* __global double* dtwn // time derivative */
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);

  // ie
  argnum++;

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->physnode_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
}

// Update the cl buffer with physnode data depending in the
// macroelement with index ie
void update_physnode_cl(field *f, int ie, cl_mem physnode_cl, double *physnode)
{
  cl_int status;
  /* status = clFinish(f->cli.commandqueue); */
  /* if(status != CL_SUCCESS) printf("%s\n", clErrorString(status)); */
  /* assert(status == CL_SUCCESS); */

  void *chkptr = clEnqueueMapBuffer(f->cli.commandqueue,
				    physnode_cl,
				    CL_TRUE,
				    CL_MAP_WRITE,
				    0, // offset
				    sizeof(cl_double) * 60, // buffersize
				    0, NULL, NULL, // events management
				    &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  assert(chkptr == physnode);

  for(int inoloc = 0; inoloc < 20; ++inoloc) {
    int ino = f->macromesh.elem2node[20 * ie + inoloc];
    for(unsigned int i = 0; i < 3; ++i)
      physnode[3 * inoloc + i] = f->macromesh.node[3 * ino + i];
  }

  status = clEnqueueUnmapMemObject(f->cli.commandqueue,
				   physnode_cl,
				   physnode,
				   0, NULL, NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
}

void DGMacroCellInterface_CL(void *mf, field *f, cl_mem *wn_cl,
			     cl_uint nwait, cl_event *wait, cl_event *done) 
{
  MacroFace *mface = (MacroFace*) mf;
  int *param = f->interp_param;

  cl_int status;
  cl_kernel kernel = f->dginterface;

  status = clSetKernelArg(kernel,
                          8,
                          sizeof(cl_mem),
                          wn_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // Set the kernel arguments
  initDGMacroCellInterface_CL(f, f->physnode_cl, f->physnodeR_cl);
  
  // Loop on the macro faces
  const unsigned int nbfaces = f->macromesh.nbfaces;
  for(int ifa = mface->first; ifa < mface->last_p1; ifa++) {
    int ieL =    f->macromesh.face2elem[4 * ifa + 0];
    int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    int ieR =    f->macromesh.face2elem[4 * ifa + 2];
    int locfaR = f->macromesh.face2elem[4 * ifa + 3];

    update_physnode_cl(f, ieL, f->physnode_cl, f->physnode);

    if(ieR >= 0)
      update_physnode_cl(f, ieR, f->physnodeR_cl, f->physnodeR);

    // Set the remaining, loop-dependant kernel arguments
    loop_initDGMacroCellInterface_CL(f, ieL, ieR, locfaL, locfaR);

    size_t numworkitems = NPGF(f->interp_param + 1, locfaL);
    status = clEnqueueNDRangeKernel(f->cli.commandqueue,
        			    kernel,
        			    1, // cl_uint work_dim,
				    NULL, // const size_t *global_work_offset,
        			    &numworkitems, // size_t *global_work_size, 
        			    NULL, // size_t *local_work_size, 
        			    nwait,  // cl_uint num_events_in_wait_list, 
				    wait, // cl_event *event_wait_list,
				    done); // cl_event *event
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);
  }
}

// Apply division by the mass matrix OpenCL version
void DGMass_CL(void *mc, field *f,
	       cl_uint nwait, cl_event *wait, cl_event *done) 
{
  MacroCell *mcell = (MacroCell*) mc;
  int *param = f->interp_param;
  cl_int status;

  init_DGMass_CL(f);

  // Loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {

    update_physnode_cl(f, ie, f->physnode_cl, f->physnode);

    status = clSetKernelArg(f->dgmass, 
			    1, 
			    sizeof(int), 
			    (void *)&ie);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);

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
				    wait, //  cl_event *event_wait_list, 
				    done); // cl_event *event
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);
  }
}

// Apply division by the mass matrix OpenCL version
void DGVolume_CL(void *mc, field *f, cl_mem *wn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done) {
  MacroCell *mcell = (MacroCell*) mc;
  cl_kernel kernel = f->dgvolume;
  int* param = f->interp_param;

  cl_int status;

  init_DGVolume_CL(f, wn_cl);

  // Loop on the elements
  for (int ie = mcell->first; ie < mcell->last_p1; ie++) {
    update_physnode_cl(f, ie, f->physnode_cl, f->physnode);

    status = clSetKernelArg(kernel,
			    1,
			    sizeof(int),
			    &ie);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);

    // Mass kernel launch
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
    //printf("%d\n", status);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);

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
  assert(status == CL_SUCCESS);
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field. OpenCL version.
void dtfield_CL(field *f, cl_mem *wn_cl,
		cl_uint nwait, cl_event *wait, cl_event *done) {
  
  set_buf_to_zero_cl(&(f->dtwn_cl), f->wsize, f, 
		     nwait, wait, &(f->clv_zbuf));

  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++) {
    DGMacroCellInterface_CL((void*) (f->mface + ifa), f, wn_cl,
			    1,
			    ifa == 0 ? &(f->clv_zbuf) : &(f->clv_mci[ifa - 1]),
			    &(f->clv_mci[ifa]));
  }

  int lastelem = f->macromesh.nbelems - 1;
  for(int ie = 0; ie <= lastelem; ++ie) {
    MacroCell *mcelli = f->mcell + ie;
    DGVolume_CL(mcelli, f, wn_cl,
		1,
		ie == 0 ? &(f->clv_mci[f->macromesh.nbelems]): 
		&(f->clv_mass[ie - 1]),
		&(f->clv_volume[ie]));
    DGMass_CL(mcelli, f,
	      1,
	      &(f->clv_volume[ie]),
	      ie < lastelem ? &(f->clv_mass[ie]) : done);
  }
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
  assert(status == CL_SUCCESS);
  
  // __global double *wn, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &(f->wn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  //__global double* dtwn // time derivative
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  
  //double dt, // time step for the stage
  double halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &halfdt);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &(f->dtwn_cl));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &(dt));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);
  
  // __global double *wn, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wn);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  //__global double *dtwn // time derivative
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  //double dt,
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &dt);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  clFinish(f->cli.commandqueue);

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
  assert(status == CL_SUCCESS);
  clFinish(f->cli.commandqueue);
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
  assert(status == CL_SUCCESS);

  // __global double *l1,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l1);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // __global double *l2,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l2);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // __global double *l3,
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l3);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // __global double *dtw, 
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw_cl);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // const double dt
  double halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(double),
			  &dt);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  clFinish(f->cli.commandqueue);

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
  assert(status == CL_SUCCESS);
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
  assert(status == CL_SUCCESS);

  cl_mem l2 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(double) * f->wsize,
			     NULL,
			     &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  cl_mem l3 = clCreateBuffer(f->cli.context,
			     0,
			     sizeof(double) * f->wsize,
			     NULL,
			     &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

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

  status = clSetUserEventStatus(stage[nstages - 1], CL_COMPLETE);

  while(f->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, f->dt);

    // stage 0
    dtfield_CL(f, w, 
	       1, &stage[3], &source[0]); // FIXME: use events
    RK4_CL_stageA(f, &l1, w, dtw,
    		  0.5 * f->dt, sizew, numworkitems,
		  1, &source[0], &stage[0]);
    
    // stage 1
    dtfield_CL(f, &l1, 1, &stage[0], &source[1]);
    RK4_CL_stageA(f, &l2, w, dtw,
    		  0.5 * f->dt, sizew, numworkitems,
		  1, &source[1], &stage[1]);
    
    // stage 2
    dtfield_CL(f, &l2, 
	       1, &stage[1], &source[2]);
    RK4_CL_stageA(f, &l3, w, dtw,
    		  f->dt, sizew, numworkitems,
		  1, &source[2], &stage[2]);
    
    // stage 3
    dtfield_CL(f, &l3, 
	       0, &stage[2], &source[3]);
    RK4_final_inplace_CL(f, w, &l1, &l2, &l3, 
			 dtw, f->dt, numworkitems,
			 1, &source[3], &stage[3]);

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
  clWaitForEvents(nwait, wait);

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
  assert(status == CL_SUCCESS);

  // Set up kernels
  init_RK2_CL_stage1(f, f->dt, &wnp1_cl);
  init_RK2_CL_stage2(f, f->dt);

  cl_event source1 = clCreateUserEvent(f->cli.context, &status);
  cl_event stage1 = clCreateUserEvent(f->cli.context, &status);
  cl_event source2 = clCreateUserEvent(f->cli.context, &status);
  cl_event stage2 = clCreateUserEvent(f->cli.context, &status);
  status = clSetUserEventStatus(stage2, CL_COMPLETE);

  while(f->tnow < tmax) {
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

#endif // _WITH_OPENCL

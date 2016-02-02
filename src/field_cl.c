#include "field.h"
#include "field_cl.h"
#include "clinfo.h"
#include "clutils.h"
#include <assert.h>
#include <string.h>
#include <sys/time.h>

void empty_kernel(field *f,
		  cl_uint nwait, cl_event *wait, cl_event *done) 
{
  cl_int status;
  status = clEnqueueTask(f->cli.commandqueue,
			 f->empty_kernel,
			 nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}
  
#ifdef _WITH_OPENCL


void CopyfieldtoCPU(field *f)
{
  // TODO: use map instead of read / writes?
  cl_int status;

  status = clFinish(f->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  const int nmacro = f->macromesh.nbelems;
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;

    status = clEnqueueReadBuffer(f->cli.commandqueue,
				 mcell->wn_cl,
				 CL_TRUE,
				 0,
				 mcell->nreal * sizeof(real),
				 f->wn + mcell->woffset,
				 0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    status = clEnqueueReadBuffer(f->cli.commandqueue,
				 mcell->dtwn_cl,
				 CL_TRUE,
				 0,
				 mcell->nreal * sizeof(real),
				 f->dtwn + mcell->woffset,
				 0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
  }

  status = clFinish(f->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void CopyfieldtoGPU(field *f)
{
  cl_int status;

  status = clFinish(f->cli.commandqueue);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  const int nmacro = f->macromesh.nbelems;

  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;

    status = clEnqueueWriteBuffer(f->cli.commandqueue,
				  f->wn_cl[ie],
				  CL_TRUE,
				  0,
				  mcell->nreal * sizeof(real),
				  f->wn + mcell->woffset,
				  0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    status = clEnqueueWriteBuffer(f->cli.commandqueue,
				  f->dtwn_cl[ie],
				  //mcell->dtwn_cl,
				  CL_TRUE,
				  0,
				  mcell->nreal * sizeof(real),
				  f->dtwn + mcell->woffset,
				  0, NULL, NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
  }

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

// wn_cl is an array of cl_mems, one per macrocell
void init_DGBoundary_CL(field *f, 
			int ieL, int locfaL,
			cl_mem *wn_cl,
			real tnow)
{
  cl_int status;
  cl_kernel kernel = f->dgboundary;

  MacroCell *mcell = f->mcell + ieL;
  
  unsigned int argnum = 0;
  
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
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
  			  sizeof(int),
  			  &locfaL);
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
                          wn_cl + ieL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          f->dtwn_cl + ieL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// wn_cl is an array of cl_mems, one per macrocell
void DGBoundary_CL(MacroFace *mface, field *f, cl_mem *wn_cl, real tnow,
		   cl_uint nwait, cl_event *wait, cl_event *done) 
{
  int *param = f->interp_param;
  cl_int status;

  int ifa = mface->ifa;
  int ieL = mface->ieL;
  int locfaL = mface->locfaL;
  int ieR = mface->ieR;

  assert(ieR == -1);
  
  size_t numworkitems = mface->npgf;
  init_DGBoundary_CL(f, ieL, locfaL, wn_cl, tnow);
  
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

void init_ExtractInterface_CL(MacroCell *mcell,
			      field *f,
			      int ifa,
			      cl_mem wn_cl)
{
  cl_int status;
  cl_kernel kernel = f->ExtractInterface;

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
                          &ifa);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &wn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
			  mcell->interface_cl + ifa);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void ExtractInterface_CL(MacroCell *mcell, field *f, int locfa, cl_mem wn_cl,
			 cl_uint nwait, cl_event *wait, cl_event *done)
{
  // Extract the faces with normals in direction i2 on both side of
  // the MacroCell.

  int *param = f->interp_param;
  int m = param[0];
  int *raf = mcell->raf;
  int *deg = mcell->deg;

  //const int faces[3][2] = { {2, 0}, {1, 3}, {4,5} };
  
  const int axis_permut[6][4] = { {0, 2, 1, 0},
				  {1, 2, 0, 1},
				  {2, 0, 1, 1},
				  {2, 1, 0, 0},
				  {0, 1, 2, 1},
				  {1, 0, 2, 0} };
  
  const int paxis[4] = {axis_permut[locfa][0],
			axis_permut[locfa][1],
			axis_permut[locfa][2],
			axis_permut[locfa][3] };

  const int d0 = paxis[0];
  const int d1 = paxis[1];
  
  init_ExtractInterface_CL(mcell, f, locfa, wn_cl);
  
  // Number of points on the face.
  int npgf = NPGF(raf, deg, locfa);
  size_t numworkitems[1] = {npgf};
  
  cl_int status;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->ExtractInterface,
				  1,            // cl_uint work_dim,
				  NULL,         // global_work_offset,
				  numworkitems, // global_work_size, 
				  NULL,          // size_t *local_work_size, 
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void init_InsertInterface_CL(MacroCell *mcell,
			     field *f,
			     int locfa,
			     cl_mem dtwn_cl)
{
  cl_int status;
  cl_kernel kernel = f->InsertInterface;

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
                          &locfa);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &dtwn_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
			  mcell->interface_cl + locfa);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void InsertInterface_CL(MacroCell *mcell, field *f, int locfa, cl_mem dtwn_cl,
			cl_uint nwait, cl_event *wait, cl_event *done)
{
  // Extract the faces with normals in direction i2 on both side of
  // the MacroCell.

  int *param = f->interp_param;
  int m = param[0];
  int *raf = mcell->raf;
  int *deg = mcell->deg;

  init_InsertInterface_CL(mcell, f, locfa, dtwn_cl);
  
  // Number of points on the face.
  int npgf = NPGF(raf, deg, locfa);
  size_t numworkitems[1] = {npgf};
  
  cl_int status;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->InsertInterface,
				  1,            // cl_uint work_dim,
				  NULL,         // global_work_offset,
				  numworkitems, // global_work_size, 
				  NULL,          // size_t *local_work_size, 
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void init_extracted_DGInterface_CL(MacroFace *mface, field *f)
{
  cl_int status;
  cl_kernel kernel = f->ExtractedDGInterfaceFlux;

  MacroCell *mcellL = f->mcell + mface->ieL;
  MacroCell *mcellR = f->mcell + mface->ieR;
    
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
			  &mface->locfaL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(int),
			  &mface->locfaR);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
			  &mcellL->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(int),
                          &mface->Rcorner);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
  			  argnum++,
                          sizeof(cl_mem),
			  &mface->wL_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
  			  argnum++,
                          sizeof(cl_mem),
			  &mface->wR_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void ExtractedDGInterface_CL(MacroFace *mface, field *f,
			     cl_uint nwait, cl_event *wait, cl_event *done)
{

  MacroCell *mcellL = f->mcell + mface->ieL;
  
  size_t numworkitems[1] = {mface->npgf};

  init_extracted_DGInterface_CL(mface, f);
  
  cl_int status;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->ExtractedDGInterfaceFlux,
				  1,            // cl_uint work_dim,
				  NULL,         // global_work_offset,
				  numworkitems, // global_work_size, 
				  NULL,         // size_t *local_work_size, 
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void init_extracted_DGBoundary_CL(MacroFace *mface, field *f, real tnow)
{
  cl_int status;
  cl_kernel kernel = f->ExtractedDGBoundaryFlux;

  MacroCell *mcellL = f->mcell + mface->ieL;
    
  unsigned int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
			  &mcellL->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
  			  argnum++,
                          sizeof(cl_mem),
			  mcellL->interface_cl + mface->locfaL);

  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
  			  argnum++,
                          sizeof(int),
			  &mface->locfaL);

  status = clSetKernelArg(kernel,
  			  argnum++,
                          sizeof(real),
			  &tnow);

  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void ExtractedDGBoundary_CL(MacroFace *mface, field *f, real tnow,
			    cl_uint nwait,
			    cl_event *wait,
			    cl_event *done)
{
  MacroCell *mcellL = f->mcell + mface->ieL;
  
  size_t numworkitems = mface->npgf;

  init_extracted_DGBoundary_CL(mface, f, tnow);
  
  cl_int status;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->ExtractedDGBoundaryFlux,
				  1,            // cl_uint work_dim,
				  NULL,         // global_work_offset,
				  &numworkitems, // global_work_size, 
				  NULL,         // size_t *local_work_size, 
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Set the loop-dependant kernel arguments for DGMacroCellInterface_CL
// wn_cl is an array of cl_mems, one per macrocell.
void init_DGMacroCellInterface_CL(field *f, 
				  int ieL, int ieR, int locfaL, int locfaR,
				  int Rcorner,
				  cl_mem *wn_cl)
{
  cl_int status;
  cl_kernel kernel = f->dginterface;

  MacroCell *mcellL = f->mcell + ieL;
  MacroCell *mcellR = f->mcell + ieR;
  
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
			  &locfaL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(int),
			  &locfaR);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
			  &mcellL->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &mcellR->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl + ieL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          f->dtwn_cl + ieL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));

  assert(status >= CL_SUCCESS);
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          wn_cl + ieR);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          f->dtwn_cl + ieR);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(int),
                          &Rcorner);
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

  int ieL = mface->ieL;
  int locfaL = mface->locfaL;
  int ieR = mface->ieR;
  int locfaR = mface->locfaR;

  size_t numworkitems = mface->npgf;
  assert(ieR >= 0);

  // Set the remaining loop-dependant kernel arguments
  init_DGMacroCellInterface_CL(f, 
			       ieL, ieR, locfaL, locfaR, mface->Rcorner,
			       wn_cl);

  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  kernel,
				  1,    // cl_uint work_dim,
				  NULL, // global_work_offset,
				  &numworkitems, // global_work_size, 
				  NULL, // size_t *local_work_size, 
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
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
  
  status = clSetKernelArg(f->dgmass,
                          argnum++,
                          sizeof(cl_mem),
			  &mcell->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(f->dgmass,
                          argnum++,
                          sizeof(cl_mem),
			  &mcell->mass_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          f->dtwn_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
			  sizeof(real) * mcell->nrealsubcell,
                          NULL);
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
  
  size_t global_work_size = mcell->npg;
  size_t local_work_size = mcell->npgsubcell;

  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
				  f->dgmass,
				  1,
				  NULL,
				  &global_work_size,
				  &local_work_size,
				  nwait, wait, done);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// wn_cl is a pointer to a macrocell's wn_cl
void init_DGFlux_CL(field *f, MacroCell *mcell, int dim0, cl_mem *wn_cl, 
		    size_t cachesize)
{
  //printf("DGFlux cachesize:%zu\n", cachesize);

  cl_kernel kernel = f->dgflux;
  cl_int status;
  int argnum = 0; 

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
                          &f->param_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel, 
			  argnum++, 
			  sizeof(int), 
			  (void *)&dim0);
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
                          f->dtwn_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * cachesize,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// wn_cl is a pointer to a macrocell's wn_cl
void DGFlux_CL(MacroCell *mcell, field *f, int dim0, cl_mem *wn_cl,
	       cl_uint nwait, cl_event *wait, cl_event *done)
{
  // FIXME: instead of ie, pass a MacroCell.
  
  // Unpack the parameters
  int *param = f->interp_param;
  int m = param[0];
  int deg[3] = {param[1], param[2], param[3]};
  int nraf[3] = {param[4], param[5], param[6]};

  int dim[3];
  dim[0] = dim0;
  dim[1] = (dim[0] + 1) % 3;
  dim[2] = (dim[1] + 1) % 3;

  // Number of points per subcell face
  int npgf = (deg[dim[1]] + 1) * (deg[dim[2]] + 1);

  // Number of subcell faces
  int nf = (nraf[dim[0]] - 1) * nraf[dim[1]] * nraf[dim[2]];
  assert(nf > 0);

  size_t numworkitems = nf * npgf;
  size_t groupsize = npgf;
  init_DGFlux_CL(f, mcell, dim0, wn_cl, 4 * m * groupsize);

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
}

// Set kernel argument for DGVolume_CL
void init_DGVolume_CL(MacroCell *mcell, field *f, cl_mem *wn_cl)
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
                          f->dtwn_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(real) * mcell->nrealsubcell,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
                          argnum++,
			  sizeof(real) * mcell->nrealsubcell,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

// Apply division by the mass matrix OpenCL version
void DGVolume_CL(MacroCell *mcell, field *f, cl_mem *wn_cl,
		 cl_uint nwait, cl_event *wait, cl_event *done) 
{
  cl_kernel kernel = f->dgvolume;

  int *param = f->interp_param;
  int m = param[0];

  size_t groupsize = mcell->npgsubcell;
  size_t numworkitems = mcell->npg;
  
  init_DGVolume_CL(mcell, f, wn_cl);
  int ie = mcell->ie;

  cl_int status;
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

void init_DGSource_CL(MacroCell *mcell, field *f,
		      real tnow, cl_mem *wn_cl)
{
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
                          sizeof(cl_mem),
                          &mcell->physnode_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
                          sizeof(cl_mem),
			  &mcell->mass_cl);
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
                          f->dtwn_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
                          argnum++,
			  sizeof(real) * mcell->nrealsubcell,
                          NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
                          argnum++,
			  sizeof(real) * mcell->nrealsubcell,
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

  size_t groupsize = mcell->npgsubcell;
  size_t numworkitems = mcell->npg;
     
  init_DGSource_CL(mcell, f, tnow, wn_cl);
  
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

void set_buf_to_zero_cl(cl_mem *buf, MacroCell *mcell, field *f,
			cl_uint nwait, cl_event *wait,  cl_event *done)
{
  cl_int status;

  cl_kernel kernel = f->zero_buf;

  int argnum = 0;
  
  status = clSetKernelArg(kernel,
                          argnum++, 
                          sizeof(cl_mem),
                          buf);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  size_t numworkitems = mcell->nreal;
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
// wn_cl is an array of pointers to the macrocell wn_cls.
void dtfield_noextract_CL(field *f, real tnow, cl_mem *wn_cl,
			  cl_uint nwait, cl_event *wait, cl_event *done)
{
  const int nmacro = f->macromesh.nbelems;

  cl_event macrofaces;
  cl_event macrobounds;
  cl_event zbuf;

  // We need to add an event for each wave so that the following
  // kernels have nwait=1, which reduces the overhead of launching
  // kernels.
  cl_event * cl_iwave = f->niwaves == 0 ? NULL :
    calloc(f->niwaves, sizeof(cl_event));
  cl_event * cl_bwave = f->nbwaves == 0 ? NULL :
    calloc(f->nbwaves, sizeof(cl_event));
  
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;
    set_buf_to_zero_cl(f->dtwn_cl + ie, mcell, f,
		       nwait, wait, f->clv_zbuf + ie);
  }

  empty_kernel(f, nmacro, f->clv_zbuf, &zbuf);
  
  if(f->niwaves > 0) {
    for(int i = 0; i < f->niwaves; ++i) {
      for(int j = 0; j < f->iwavecount[i]; ++j) {
	int ifa = f->iwave[i][j];
	MacroFace *mface = f->mface + ifa;
	DGMacroCellInterface_CL(mface, f, wn_cl,
				1,
				(i == 0) ?  &zbuf: cl_iwave + i - 1,
				f->clv_iwave[i] + j);
      }
      clWaitForEvents(f->iwavecount[i], f->clv_iwave[i]);
      empty_kernel(f, f->iwavecount[i], f->clv_iwave[i],
		   cl_iwave + i);
    }
    empty_kernel(f, 1, f->clv_iwave[f->niwaves - 1], &macrofaces);
  } else {
    empty_kernel(f, nmacro, f->clv_zbuf, &macrofaces);
  }
  
  if(f->nbwaves > 0) {
    for(int i = 0; i < f->nbwaves; ++i) {
      for(int j = 0; j < f->bwavecount[i]; ++j) {
	int ifa = f->bwave[i][j];
	MacroFace *mface = f->mface + ifa;
	DGBoundary_CL(mface, f, wn_cl, f->tnow,
		      1, (i == 0) ? &macrofaces : cl_bwave + i - 1,
		      f->clv_bwave[i] + j);
      }
      clWaitForEvents(f->bwavecount[i], f->clv_bwave[i]);
      empty_kernel(f, f->bwavecount[i], f->clv_bwave[i],
		   cl_bwave + i);
    }
    empty_kernel(f, 1, f->clv_bwave[f->nbwaves - 1], &macrobounds);
  } else {
    empty_kernel(f, 1, &macrofaces, &macrobounds);
  }

  unsigned int ndim = 3;
  if(f->macromesh.is2d)
    ndim = 2;
  if(f->macromesh.is1d)
    ndim = 1;
  
  // The kernels for the intra-macrocell computations can be launched
  // in parallel between macrocells.
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;

    int nflux = 0;
    int fluxlist[3];
    for(int d = 0; d < ndim; ++d) {
      if(mcell->raf[d] > 1) {
	fluxlist[nflux++] = d;
      }
    }

    DGFlux_CL(mcell, f, fluxlist[0], wn_cl + ie, 
	      1, &macrobounds,
	      f->clv_flux[fluxlist[0]] + ie);
    for(int d = 1; d < nflux; ++d) {
      DGFlux_CL(mcell, f, fluxlist[d], wn_cl + ie, 
		1, f->clv_flux[fluxlist[d - 1]] + ie,
		f->clv_flux[fluxlist[d]] + ie);
    }
    
    DGVolume_CL(mcell, f, wn_cl + ie,
  		1, f->clv_flux[nflux - 1] + ie, f->clv_volume + ie);

    if(f->use_source_cl) {
      DGSource_CL(mcell, f, tnow, wn_cl + ie, 
		  1, f->clv_volume + ie, f->clv_source + ie);
    }
    
    DGMass_CL(mcell, f,
  	      1,
	      f->use_source_cl ? f->clv_source + ie : f->clv_volume + ie,
  	      f->clv_mass + ie);
  }

  empty_kernel(f, nmacro, f->clv_mass, done);
  
  for(int i = 0; i < f->niwaves; ++i) {
    for(int j = 0; j < f->iwavecount[i]; ++j) {
      f->minter_calls++;
      f->minter_time += clv_duration(f->clv_iwave[i][j]);
    }
  }
  for(int i = 0; i < f->nbwaves; ++i) {
    for(int j = 0; j < f->bwavecount[i]; ++j) {
      f->boundary_calls++;
      f->boundary_time += clv_duration(f->clv_bwave[i][j]);
    }
  }
  for(int ie = 0; ie < nmacro; ++ie) {
    f->zbuf_calls++;
    f->zbuf_time += clv_duration(f->clv_zbuf[ie]);
    if(f->use_source_cl) {
      f->source_calls++;
      f->source_time += clv_duration(f->clv_source[ie]);
    }
    f->flux_calls++;
    f->flux_time += clv_duration(f->clv_flux[0][ie]);
    if(ndim > 1) {
      f->flux_time += clv_duration(f->clv_flux[1][ie]);
      f->flux_calls++;
    }
    if(ndim > 2) {
      f->flux_calls++;
      f->flux_time += clv_duration(f->clv_flux[2][ie]);
    }
    f->vol_calls++;
    f->vol_time += clv_duration(f->clv_volume[ie]);
    f->mass_calls++;
    f->mass_time += clv_duration(f->clv_mass[ie]);
  }

  for(int i = 0; i < f->niwaves; ++i) {
    clReleaseEvent(cl_iwave[i]);
  }
  if(f->niwaves > 0)
    free(cl_iwave);
  
  for(int i = 0; i < f->nbwaves; ++i) {
    clReleaseEvent(cl_bwave[i]);
  }
  if(f->nbwaves > 0)
    free(cl_bwave);
  
  clReleaseEvent(macrobounds);
  clReleaseEvent(macrofaces);
}

// Apply the Discontinuous Galerkin approximation for computing the
// time derivative of the field. OpenCL version with face extraction.
// wn_cl is an array of pointers to the macrocell wn_cls.
void dtfield_extract_CL(field *f, real tnow, cl_mem *wn_cl,
			cl_uint nwait, cl_event *wait, cl_event *done)
{
  const int nmacro = f->macromesh.nbelems;

  cl_event macrobounds;
  cl_event zbuf;

  // We need to add an event for each wave so that the following
  // kernels have nwait=1, which reduces the overhead of launching
  // kernels.

  // FIXME: deal with the cases where these can be zero.
  const int nboundary = f->macromesh.nboundaryfaces;
  const int ninterface = f->macromesh.nmacrointerfaces;
  
  unsigned int ndim = 3;
  if(f->macromesh.is2d)
    ndim = 2;
  if(f->macromesh.is1d)
    ndim = 1;
  
  const int ncellfaces = 2 * ndim;
  // 1d: face 0, 2; 2d: face 0, 1, 2, 3, etc.
  int ifalist[6] = {0, 2, 1, 3, 4, 5};

  // TODO: move to field struct.
  const int nextract = ncellfaces * nmacro;
  cl_event clv_extracts;
  cl_event clv_inserts;

  cl_event clv_bi[2];
  
  for(int ie = 0; ie < nmacro; ++ie) {
    set_buf_to_zero_cl(f->dtwn_cl + ie, f->mcell + ie, f,
		       nwait, wait, f->clv_zbuf + ie);
  }

  empty_kernel(f, nmacro, f->clv_zbuf, &zbuf);
  
  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    MacroCell *mcell = f->mcell + ie;
    int offset = ie * ncellfaces;
    for(int ifa = 0; ifa < ncellfaces; ++ifa) {
      ExtractInterface_CL(mcell, f, ifalist[ifa], wn_cl[ie],
			  1, &zbuf, f->clv_extract + offset + ifa);
    }
  }

  empty_kernel(f, nextract, f->clv_extract, &clv_extracts);
  
  if(nboundary > 0) {
    for(int i = 0; i < nboundary; ++i) {
      int ifa = f->macromesh.boundaryface[i];
      MacroFace *mface = f->mface + ifa;
      int ieL = mface->ieL;
      MacroCell *mcellL = f->mcell + ieL;
      ExtractedDGBoundary_CL(mface, f, tnow,
			     1, &clv_extracts, f->clv_boundary + i);
    }
    empty_kernel(f, nboundary, f->clv_boundary, &clv_bi[0]);
  } else {
    empty_kernel(f, 1, &clv_extracts, &clv_bi[0]);
  }

  if(ninterface > 0) {
    for(int i = 0; i < ninterface; ++i) {
      int ifa = f->macromesh.macrointerface[i];
      MacroFace *mface = f->mface + ifa;
      ExtractedDGInterface_CL(mface, f, 1, &clv_extracts, f->clv_interface + i);
    }
    empty_kernel(f, ninterface, f->clv_interface, &clv_bi[1]);
  } else {
    empty_kernel(f, 1, &clv_extracts, &clv_bi[1]);
  }

  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    MacroCell *mcell = f->mcell + ie;
    int offset = ie * ncellfaces;
    for(int ifa = 0; ifa < ncellfaces; ++ifa) {
      InsertInterface_CL(mcell, f, ifalist[ifa], f->dtwn_cl[ie],
			 ifa == 0 ? 2 : 1,
			 ifa == 0 ? clv_bi : f->clv_insert + offset + ifa - 1,
			 f->clv_insert + offset + ifa);
    }
  }

  empty_kernel(f, nextract, f->clv_insert, &macrobounds);
  
  // The kernels for the intra-macrocell computations can be launched
  // in parallel between macrocells.
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;

    int nflux = 0;
    int fluxlist[3];
    for(int d = 0; d < ndim; ++d) {
      if(mcell->raf[d] > 1) {
	fluxlist[nflux++] = d;
      }
    }

    DGFlux_CL(mcell, f, fluxlist[0], wn_cl + ie, 
	      1, &macrobounds,
	      f->clv_flux[fluxlist[0]] + ie);
    for(int d = 1; d < nflux; ++d) {
      DGFlux_CL(mcell, f, fluxlist[d], wn_cl + ie, 
		1, f->clv_flux[fluxlist[d - 1]] + ie,
		f->clv_flux[fluxlist[d]] + ie);
    }
    
    DGVolume_CL(mcell, f, wn_cl + ie,
  		1, f->clv_flux[nflux - 1] + ie, f->clv_volume + ie);

    if(f->use_source_cl) {
      DGSource_CL(mcell, f, tnow, wn_cl + ie, 
		  1, f->clv_volume + ie, f->clv_source + ie);
    }
    
    DGMass_CL(mcell, f,
  	      1,
	      f->use_source_cl ? f->clv_source + ie : f->clv_volume + ie,
  	      f->clv_mass + ie);
  }

  empty_kernel(f, nmacro, f->clv_mass, done);

  clWaitForEvents(nmacro, f->clv_mass);
  
  for(int ie = 0; ie < nmacro; ++ie) {
    f->zbuf_calls++;
    f->zbuf_time += clv_duration(f->clv_zbuf[ie]);
    if(f->use_source_cl) {
      f->source_calls++;
      f->source_time += clv_duration(f->clv_source[ie]);
    }
    f->flux_calls++;
    f->flux_time += clv_duration(f->clv_flux[0][ie]);
    if(ndim > 1) {
      f->flux_time += clv_duration(f->clv_flux[1][ie]);
      f->flux_calls++;
    }
    if(ndim > 2) {
      f->flux_calls++;
      f->flux_time += clv_duration(f->clv_flux[2][ie]);
    }
    f->vol_calls++;
    f->vol_time += clv_duration(f->clv_volume[ie]);
    f->mass_calls++;
    f->mass_time += clv_duration(f->clv_mass[ie]);
  }

  for(int ie = 0; ie < f->macromesh.nbelems; ++ie) {
    int offset = ie * ncellfaces;
    for(int ifa = 0; ifa < ncellfaces; ++ifa) {
      f->extractinsert_calls++;
      f->extractinsert_time += clv_duration(f->clv_extract[offset + ifa]);
      f->extractinsert_calls++;
      f->extractinsert_time += clv_duration(f->clv_insert[offset + ifa]);
    }
  }

  for(int i = 0; i < nboundary; ++i) {
    f->boundary_calls++;
    f->boundary_time += clv_duration(f->clv_boundary[i]);
  }

  for(int i = 0; i < ninterface; ++i) {
    f->minter_calls++;
    f->minter_time += clv_duration(f->clv_interface[i]);
  }

  clReleaseEvent(clv_extracts);
  clReleaseEvent(clv_inserts);
  
  clReleaseEvent(clv_bi[0]);
  clReleaseEvent(clv_bi[1]);
  
  clReleaseEvent(macrobounds);
}

void dtfield_CL(field *f, real tnow, cl_mem *wn_cl,
		  cl_uint nwait, cl_event *wait, cl_event *done)
{
  if(f->extract_cl)
    dtfield_extract_CL(f, tnow, wn_cl, nwait, wait, done);
  else
    dtfield_noextract_CL(f, tnow, wn_cl, nwait, wait, done);
}

// Set kernel arguments for first stage of RK2
// wnp1_cl is a pointer to an array of cl_mems, one per macrocell
void init_RK2_CL_stage1(MacroCell *mcell, field *f,
			const real dt, cl_mem *wnp1_cl)
{
  cl_kernel kernel = f->RK_out_CL;
  cl_int status;
  int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wnp1_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          f->wn_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          f->dtwn_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &dt);
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
                          f->wn_cl + mcell->ie);
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++,
                          sizeof(cl_mem),
                          f->dtwn_cl + mcell->ie);
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

// Time integration by a second-order Runge-Kutta algorithm, OpenCL version.
void RK2_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  if(dt <= 0)
    dt = set_dt(f);

  f->itermax = tmax / dt;
  int freq = (1 >= f->itermax / 10)? 1 : f->itermax / 10;
  int iter = 0;

  cl_int status;

  const int nmacro = f->macromesh.nbelems;
  
  cl_mem *wnp1_cl = calloc(nmacro, sizeof(cl_mem));
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;
    wnp1_cl[ie] = clCreateBuffer(f->cli.context,
				 0,
				 sizeof(real) * mcell->nreal,
				 NULL,
				 &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
  }

  cl_event start;
  cl_event *stage1 = calloc(nmacro, sizeof(cl_event));
  cl_event *stage2 = calloc(nmacro, sizeof(cl_event));
  cl_event source1;
  cl_event source2;

  printf("Starting RK2_CL\n");
  
  struct timeval t_start;
  gettimeofday(&t_start, NULL);

  empty_kernel(f, nwait, wait, &start);
  while(f->tnow < tmax) {
    
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    // stage 0
    dtfield_CL(f, f->tnow, f->wn_cl,
	       1, &start, &source1);
    
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      init_RK2_CL_stage1(mcell, f, 0.5 * dt, wnp1_cl);
      RK2_CL_stage1(mcell, f, f->wsize, 1, &source1, stage1 + ie);
    }
    
    f->tnow += 0.5 * dt;

    // stage 1
    dtfield_CL(f, f->tnow, wnp1_cl,
	       nmacro, stage1, &source2);
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      init_RK2_CL_stage2(mcell, f, dt);
      RK2_CL_stage2(mcell, f, f->wsize, 1, &source2, stage2 + ie);
    }
    
    f->tnow += 0.5 * dt;

    clWaitForEvents(nmacro, stage2);
    for(int ie = 0; ie < nmacro; ++ie) {
      f->rk_calls++;
      f->rk_time += clv_duration(stage1[ie]);
      f->rk_calls++;
      f->rk_time += clv_duration(stage2[ie]);
    }

    iter++;
  }

  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
  empty_kernel(f, nmacro, stage2, done);
  
  struct timeval t_end;
  gettimeofday(&t_end, NULL);

  // Release events
  clReleaseEvent(start);
  for(int ie = 0; ie < nmacro; ++ie) {
    clReleaseEvent(stage1[ie]);
    clReleaseEvent(stage2[ie]);
  }
  free(stage2);
  free(stage1);

  // Release memory
  for(int ie = 0; ie < nmacro; ++ie) {
    clReleaseMemObject(wnp1_cl[ie]);
  }
  free(wnp1_cl);
  
  double rkseconds = (t_end.tv_sec - t_start.tv_sec) * 1.0 // seconds
    + (t_end.tv_usec - t_start.tv_usec) * 1e-6; // microseconds
  printf("\nTotal RK time (s):\n%f\n", rkseconds);
  printf("\nTotal RK time per time-step (s):\n%f\n", rkseconds / iter );
}

void RK4_CL_stageA(MacroCell *mcell, field *f, 
		   cl_mem *wnp1, cl_mem *wn, cl_mem *dtw, 
		   const real dt, size_t numworkitems, 
		   cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = f->RK4_first_stages;
  cl_int status;
  int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wnp1);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          wn);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  //real dt,
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t nwork = mcell->nreal;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
  				  kernel,
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

void RK4_final_inplace_CL(MacroCell *mcell, field *f, 
			  cl_mem *w_cl, cl_mem *l1, cl_mem *l2, cl_mem *l3, 
			  cl_mem *dtw_cl, const real dt, 
			  const size_t numworkitems, 
			  cl_uint nwait, cl_event *wait, cl_event *done)
{
  cl_kernel kernel = f->RK4_final_stage;
  cl_int status;
  int argnum = 0;

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          w_cl);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l1);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l2);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          l3);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  status = clSetKernelArg(kernel,
			  argnum++, 
                          sizeof(cl_mem),
                          dtw_cl + mcell->ie);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  real halfdt = 0.5 * dt;
  status = clSetKernelArg(kernel,
			  argnum++,
			  sizeof(real),
			  &dt);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  size_t nwork = mcell->nreal;
  status = clEnqueueNDRangeKernel(f->cli.commandqueue,
  				  kernel,
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

// Time integration by a fourth-order Runge-Kutta algorithm, OpenCL version.
void RK4_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done) 
{
  if(dt <= 0)
    dt = set_dt(f);

  f->itermax = tmax / dt;
  int nouts = 10;

  int freq = (1 >= f->itermax / nouts)? 1 : f->itermax / nouts;
  int iter = 0;

  cl_int status;  

  const int nmacro = f->macromesh.nbelems;

  // Allocate memory for RK stages
  cl_mem *l1 = calloc(nmacro, sizeof(cl_mem));
  cl_mem *l2 = calloc(nmacro, sizeof(cl_mem));
  cl_mem *l3 = calloc(nmacro, sizeof(cl_mem));
  for(int ie = 0; ie < nmacro; ++ie) {
    MacroCell *mcell = f->mcell + ie;
    
    l1[ie] = clCreateBuffer(f->cli.context,
			    0,
			    sizeof(real) * mcell->nreal,
			    NULL,
			    &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    l2[ie] = clCreateBuffer(f->cli.context,
			    0,
			    sizeof(real) * mcell->nreal,
			    NULL,
			    &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);

    l3[ie] = clCreateBuffer(f->cli.context,
			    0,
			    sizeof(real) * mcell->nreal,
			    NULL,
			    &status);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
  }

  size_t numworkitems = f->wsize;

  cl_mem *w = f->wn_cl;
  cl_mem *dtw = f->dtwn_cl;

  // Set up events for source and stage.
  cl_event source0;
  cl_event source1;
  cl_event source2;
  cl_event source3;
  
  cl_event *stage0 = calloc(nmacro, sizeof(cl_event));
  cl_event *stage1 = calloc(nmacro, sizeof(cl_event));
  cl_event *stage2 = calloc(nmacro, sizeof(cl_event));
  cl_event *stage3 = calloc(nmacro, sizeof(cl_event));

  cl_event start;
  empty_kernel(f, nwait, wait, &start);
  
  printf("Starting RK4_CL\n");

  struct timeval t_start;
  gettimeofday(&t_start, NULL);

  while(f->tnow < tmax) {
    
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);

    
    // stage 0
    dtfield_CL(f, f->tnow, w, 
	       1, &start, &source0);
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      RK4_CL_stageA(mcell, f, l1 + ie, w + ie, dtw,
		    0.5 * dt, numworkitems,
		    1, &source0, stage0 + ie);
    }
    
    f->tnow += 0.5 * dt;
    
    // stage 1
    dtfield_CL(f, f->tnow, l1,
	       nmacro, stage0, &source1);
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      RK4_CL_stageA(mcell, f, l2 + ie, w + ie, dtw,
		    0.5 * dt, numworkitems,
		    1, &source1, stage1 + ie);
    }
    
    // stage 2
    dtfield_CL(f, f->tnow, l2, 
	       nmacro, stage1, &source2);
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      RK4_CL_stageA(mcell, f, l3 + ie, w + ie, dtw,
		    dt, numworkitems,
		    1, &source2, stage2 + ie);
    }
    
    f->tnow += 0.5 * dt;
        
    // stage 3
    dtfield_CL(f, f->tnow, l3, 
	       nmacro, stage2, &source3);
    for(int ie = 0; ie < nmacro; ++ie) { 
      MacroCell *mcell = f->mcell + ie;
      RK4_final_inplace_CL(mcell, f, w + ie, l1 + ie, l2 + ie, l3 + ie, 
			   dtw, dt, numworkitems,
			   1, &source3, stage3 + ie);
    }

    clWaitForEvents(nmacro, stage3);

    for(int ie = 0; ie < nmacro; ++ie) {
      f->rk_calls++;
      f->rk_time += clv_duration(stage0[ie]);
      f->rk_calls++;
      f->rk_time += clv_duration(stage1[ie]);
      f->rk_calls++;
      f->rk_time += clv_duration(stage2[ie]);
      f->rk_calls++;
      f->rk_time += clv_duration(stage3[ie]);
    }
    
    iter++;
  }

  empty_kernel(f, nmacro, stage3, done);
  printf("t=%f iter=%d/%d dt=%f\n", f->tnow, iter, f->itermax, dt);
  
  struct timeval t_end;
  gettimeofday(&t_end, NULL); 
  
  double rkseconds = (t_end.tv_sec - t_start.tv_sec) * 1.0 // seconds
    + (t_end.tv_usec - t_start.tv_usec) * 1e-6; // microseconds
  printf("\nTotal RK time (s):\n%f\n", rkseconds);
  printf("\nTotal RK time per time-step (s):\n%f\n", rkseconds / iter );
  
  // Release events
  clReleaseEvent(start);
  clReleaseEvent(source0);
  clReleaseEvent(source1);
  clReleaseEvent(source2);
  clReleaseEvent(source3);
  for(int ie = 0; ie < nmacro; ++ie) {
    clReleaseEvent(stage0[ie]);
    clReleaseEvent(stage1[ie]);
    clReleaseEvent(stage2[ie]);
    clReleaseEvent(stage3[ie]);
  }
  free(stage0);
  free(stage1);
  free(stage2);
  free(stage3);

  // Release buffers
  for(int ie = 0; ie < nmacro; ++ie) {
    clReleaseMemObject(l1[ie]);
    clReleaseMemObject(l2[ie]);
    clReleaseMemObject(l3[ie]);
  }
  free(l1);
  free(l2);
  free(l3);
}

void printtimes(char* name, cl_ulong ns, int ncall, double N)
{
  double s = 1e-9 * ns;
  double meantime = ncall > 0 ? s / ncall : 0.0;
  
  printf("%s%f%%\t%gs\t%d\t%g\n", 
	 name, ns * N, s, ncall, meantime);
}

void show_cl_timing(field *f)
{
  printf("\n");
 
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
  total += f->extractinsert_time;
  
  double N = 100.0 / total;

  cl_ulong ns;
  int ncall;
  
  printf("Kernel name:             percent\ttotal   \tncalls\tmean time\n");

  
  printtimes("set_buf_to_zero_cl:      ", f->zbuf_time, f->zbuf_calls, N);
  printtimes("DGMacroCellInterface_CL: ", f->minter_time, f->minter_calls, N); 
  printtimes("DGBoundary:              ", f->boundary_time, f->boundary_calls,
	     N); 
  printtimes("Insert / extract:        ", f->extractinsert_time,
	     f->extractinsert_calls, N);
  printtimes("DGVolume_CL:             ", f->vol_time, f->vol_calls, N);
  printtimes("DGFlux_CL:               ", f->flux_time, f->flux_calls, N);
  printtimes("DGMass_CL:               ", f->mass_time, f->mass_calls, N);
  printtimes("DGSource_CL:             ", f->source_time, f->source_calls, N);
  printtimes("rk:                      ", f->rk_time, f->rk_calls, N);

  printf("\n");
  
  ns = total;
  double total_time = 1e-9 * ns;
  printf("total:                   %f%%\t%fs\n", 
	 ns*N, total_time);

  printf("\n");
}

#endif // _WITH_OPENCL

#include "simulation_spu.h"
#include <stdlib.h>

//#define WORKER_ON_NODE

schnaps_real seconds() {
  struct timespec ts;
  schnaps_real res;
#ifdef __MACH__
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts.tv_sec = mts.tv_sec;
  ts.tv_nsec = mts.tv_nsec;
  res =  (schnaps_real)ts.tv_sec + 1e-9 * (schnaps_real)ts.tv_nsec;
#else
  clock_gettime(CLOCK_MONOTONIC, &ts);
  res = (schnaps_real)ts.tv_sec + 1e-9 * (schnaps_real)ts.tv_nsec;
#endif
  return res;
}



void DisplayHandle_SPU(starpu_data_handle_t handle,
                       const char* name) {
  starpu_task_wait_for_all();
  // Force update of memory on main node
  starpu_data_fetch_on_node(handle, 0, 0);
  // Warning: only schnaps_real arrays for the moment
  // Extention to every type could be made with starpu_vector_get_elemsize
  schnaps_real* ptr = (schnaps_real*) starpu_vector_get_local_ptr(handle);
  int size = starpu_vector_get_nx(handle);
  for (int i = 0; i < size; ++i)
    printf("%s[%d]: %f\n", name, i , ptr[i]);
}



void ZeroBuffer_SPU(starpu_data_handle_t handle) {
  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = ZeroBuffer_codelet();
  task->cl_arg = NULL;
  task->cl_arg_size = 0;
  task->handles[0] = handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}

void ZeroBuffer_SPU2(starpu_data_handle_t handle, int ie) {

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);
  
  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = ZeroBuffer_codelet();
  task->cl_arg = NULL;
  task->cl_arg_size = 0;
  task->handles[0] = handle;

  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[ie%nbgpu]);
  

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void ZeroBuffer_C(void *buffers[], void *cl_args);
void ZeroBuffer_OCL(void *buffers[], void *cl_args);

struct starpu_codelet* ZeroBuffer_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel zb_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "zerobuffer"
  };    

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "ZeroBuffer";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_RW;
    //codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &zb_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = ZeroBuffer_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = ZeroBuffer_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
        codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
      
    }
  }

  return &codelet;
}

void ZeroBuffer_C(void *buffers[], void *cl_args) {
  schnaps_real* buffer = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							       (struct starpu_vector_interface *)buffers[0]);
  size_t size = STARPU_VECTOR_GET_NX(buffers[0]);

  for (size_t i = 0; i < size; i++) buffer[i] = 0;
}

void ZeroBuffer_OCL(void *buffers[], void *cl_args) {
  cl_mem buffer = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  size_t size = STARPU_VECTOR_GET_NX(buffers[0]);

  int devid = starpu_worker_get_devid(starpu_worker_get_id());
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program,
                                            "set_buffer_to_zero",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &buffer);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  status = clEnqueueNDRangeKernel(queue, kernel,
                                  1, // work_dim
                                  NULL, // *global_work_offset
                                  &size, // *global_work_size
                                  NULL, // *local_work_size
                                  0, // num_events_in_wait_list
                                  NULL, // *wait_list
                                  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}



void AddBuffer_SPU(schnaps_real alpha,
                   starpu_data_handle_t handle_in,
                   starpu_data_handle_t handle_out) {
  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &alpha, sizeof(schnaps_real),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = AddBuffer_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = handle_in;
  task->handles[1] = handle_out;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void AddBuffer_SPU2(schnaps_real alpha,
		    starpu_data_handle_t handle_in,
                    starpu_data_handle_t handle_out,
                    int ie) {

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &alpha, sizeof(schnaps_real),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = AddBuffer_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = handle_in;
  task->handles[1] = handle_out;

  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[ie%nbgpu]);
  
  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}



void AddBuffer_C(void *buffers[], void *cl_args);
void AddBuffer_OCL(void *buffers[], void *cl_args);

struct starpu_codelet* AddBuffer_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel ab_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "addbuffer"
  }; 

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "AddBuffer";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;
    codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &ab_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = AddBuffer_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = AddBuffer_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
        codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
      
    }
  }

  return &codelet;
}

void AddBuffer_C(void *buffers[], void *cl_args) {
  schnaps_real alpha;

  starpu_codelet_unpack_args(cl_args, &alpha);
  free(cl_args);

  schnaps_real* buffer_in = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								  (struct starpu_vector_interface *)buffers[0]);
  schnaps_real* buffer_out = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								   (struct starpu_vector_interface *)buffers[1]);

  size_t size_in = STARPU_VECTOR_GET_NX(buffers[0]);
  size_t size_out = STARPU_VECTOR_GET_NX(buffers[1]);

  assert(size_in == size_out);

  for (size_t i = 0; i < size_in; i++) buffer_out[i] += alpha * buffer_in[i];
}

void AddBuffer_OCL(void *buffers[], void *cl_args) {
  schnaps_real alpha;

  starpu_codelet_unpack_args(cl_args, &alpha);
  free(cl_args);

  cl_mem buffer_in = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem buffer_out = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);

  size_t size_in = STARPU_VECTOR_GET_NX(buffers[0]);
  size_t size_out = STARPU_VECTOR_GET_NX(buffers[1]);

  assert(size_in == size_out);

  int devid = starpu_worker_get_devid(starpu_worker_get_id());
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program,
                                            "AddBuffer",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(schnaps_real), &alpha);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &buffer_in);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &buffer_out);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  status = clEnqueueNDRangeKernel(queue, kernel,
                                  1, // work_dim
                                  NULL, // *global_work_offset
                                  &size_in, // *global_work_size
                                  NULL, // *local_work_size
                                  0, // num_events_in_wait_list
                                  NULL, // *wait_list
                                  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}


//xxx begin init stuff




void DGInitField_C(void *buffers[], void *cl_args);
//void DGInitField_OCL(void *buffers[], void *cl_args);

struct starpu_codelet* DGInitField_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgif_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dginitfield"
  };

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGInitField";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_W;
    //codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &dgif_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGInitField_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    /* if (starpu_ocl_use) { */
    /*   LOAD_OPENCL_PROGRAM_SPU(); */

    /*   size_t ncodelet = 0; */
    /*   codelet.opencl_funcs[ncodelet++] = DGInitField_OCL; */
    /*   codelet.opencl_funcs[ncodelet++] = NULL; */

    /*   for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++) */
    /* 	codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC; */
    /* } // end if use openCL */
  } // end if not is_init

  return &codelet;
}

void DGInitField_C(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  initdataptr InitData;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex, &InitData);
  free(cl_args);

  schnaps_real* w = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							  (struct starpu_vector_interface *)buffers[0]);

  

  for (int ipg = 0; ipg < NPG(deg, raf); ++ipg) {
    schnaps_real xref[3], xphy[3];
    ref_pg_vol(deg, raf, ipg, xref, NULL, NULL);
    schnaps_ref2phy(physnode,
		    xref,
		    NULL,  // dpsiref
		    -1,  // ifa
		    xphy,
		    NULL,
		    NULL,
		    NULL,  // dpsi
		    NULL);  // vnds
    
    schnaps_real wL[m];

    InitData(xphy, wL);
    
    for (int iv = 0; iv < m; ++iv) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      w[imem] = wL[iv];
    }

  }
}

void DGInitField_SPU(field* f) {
  assert (f->model.InitData != NULL);

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.InitData, sizeof(initdataptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGInitField_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}



//xxx end init stuff




void DGSubCellInterface_SPU(field *f) {
  if (f->raf[0] * f->raf[1] * f->raf[2] == 1) return;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.NumFlux, sizeof(fluxptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGSubCellInterface_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void DGSubCellInterface_SPU2(field *f, int ie) {

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);
  
  
  if (f->raf[0] * f->raf[1] * f->raf[2] == 1) return;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.NumFlux, sizeof(fluxptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGSubCellInterface_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;

  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[ie%nbgpu]);  

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}



void DGSubCellInterface_C(void *buffers[], void *cl_args);
void DGSubCellInterface_OCL(void *buffers[], void *cl_args);
void DGSubCellInterface3D_OCL(void *buffers[], void *cl_args);

struct starpu_codelet* DGSubCellInterface_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgsci_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgsubcellinterface"
  }; 
  

  if (!is_init){
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGSubCellInterface";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;
    codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &dgsci_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGSubCellInterface_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGSubCellInterface_OCL;
      codelet.opencl_funcs[ncodelet++] = DGSubCellInterface3D_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
        codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
      
    }
  }

  return &codelet;
}

void DGSubCellInterface_C(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
                             &Varindex, &NumFlux);
  free(cl_args);

  schnaps_real* w = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							  (struct starpu_vector_interface *)buffers[0]);
  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[1]);


  int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};

  // Loop over the subcells
  for(int icL0 = 0; icL0 < raf[0]; icL0++) {
    for(int icL1 = 0; icL1 < raf[1]; icL1++) {
      for(int icL2 = 0; icL2 < raf[2]; icL2++) {

	int icL[3] = {icL0, icL1, icL2};

	// Get the left subcell id
	int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);

	// First glop index in the subcell
	int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	// Loop over subcell faces in the three directions
	for (int dim0 = 0; dim0 < 3; ++dim0) {
	  // Compute the subface flux only if we do not touch the macrocell boundary
	  if (icL[dim0] != raf[dim0] - 1) {
	    // The right cell index corresponds to an increment in the dim0 direction
	    int icR[3] = {icL[0], icL[1], icL[2]};
	    icR[dim0]++;
	    int ncR = icR[0] + raf[0] * (icR[1] + raf[1] * icR[2]);
	    int offsetR = npg[0] * npg[1] * npg[2] * ncR;

	    // FIXME: only write to L-values (and do both faces) to parallelize better

	    const int altdim1[3] = {1, 0, 0};
	    const int altdim2[3] = {2, 2, 1};

	    // Loop on the left glops of the subface
	    int dim1 = altdim1[dim0];
	    int dim2 = altdim2[dim0];
	    int iL[3];
	    iL[dim0] = deg[dim0];
	    for (iL[dim2] = 0; iL[dim2] < npg[dim2]; ++iL[dim2]) {
	      for (iL[dim1] = 0; iL[dim1] < npg[dim1]; ++iL[dim1]) {
		// Find the right and left glops volume indices
		int iR[3] = {iL[0], iL[1], iL[2]};
		iR[dim0] = 0;

		int ipgL = offsetL + iL[0] + (deg[0] + 1) * (iL[1] + (deg[1] + 1) * iL[2]);
		int ipgR = offsetR + iR[0] + (deg[0] + 1) * (iR[1] + (deg[1] + 1) * iR[2]);
		//printf("ipgL=%d ipgR=%d\n", ipgL, ipgR);

		// Compute the normal vector for integrating on the facee
		schnaps_real vnds[3];
		{
		  schnaps_real xref[3], wpg_vol;
		  ref_pg_vol(deg, raf, ipgL, xref, &wpg_vol, NULL);
		  schnaps_real dtau[3][3], codtau[3][3];
		  schnaps_ref2phy(physnode,
                                  xref,
                                  NULL, // dphiref
                                  -1,  // ifa
                                  NULL, // xphy
                                  dtau,
                                  codtau,
                                  NULL, // dphi
                                  NULL);  // vnds

		  // Compute the normal vector to take into account the subcell surface
		  schnaps_real h1h2 = 1. / raf[dim1] / raf[dim2];
		  vnds[0] = codtau[0][dim0] * h1h2;
		  vnds[1] = codtau[1][dim0] * h1h2;
		  vnds[2] = codtau[2][dim0] * h1h2;
		}

		// Numerical flux from the left and right state
		schnaps_real wL[m], wR[m], flux[m];
		for (int iv = 0; iv < m; ++iv) {
		  int imemL = Varindex(deg, raf, m, ipgL, iv);
		  int imemR = Varindex(deg, raf, m, ipgR, iv);
		  wL[iv] = w[imemL];
		  wR[iv] = w[imemR];
		}
		NumFlux(wL, wR, vnds, flux);

		// Subcell surface glop weight
		schnaps_real wpg = wglop(deg[dim1], iL[dim1]) *
		  wglop(deg[dim2], iL[dim2]);

		/* printf("wL %f wR %f vnds %f %f %f flux %f wpg %f\n", */
                /*        wL[0], wR[0], vnds[0], vnds[1], vnds[2], flux[0], wpg); */

		// Distribute the flux on the two sides
		for (int iv = 0; iv < m; ++iv) {
		  int imemL = Varindex(deg, raf, m, ipgL, iv);
		  int imemR = Varindex(deg, raf, m, ipgR, iv);
		  res[imemL] -= flux[iv] * wpg;
		  res[imemR] += flux[iv] * wpg;
		}
	      }  // face glop dim1
	    }  // face glop dim2
	  }  // internal face
	}  // dim
      }  // icL2
    }  // icL1
  }  // icL0
}


void DGSubCellInterface_OCL(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_args,
                             &m, deg, raf, physnode,
                             &Varindex, &NumFlux);
  free(cl_args);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);


  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program,
                                            "DGFlux",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
                                      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      60 * sizeof(schnaps_real),
                                      physnode,
                                      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
                                   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   7 * sizeof(int),
                                   param,
                                   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_event events[3];
  cl_event *event = events;
  cl_event *wait_event = NULL;
  int nb_wait_events = 0;

  // Loop over the dimensions
  for (int dim0 = 0; dim0 < 3; ++dim0) {
    // More than one subcell according to that dimension
    if (raf[dim0] > 1) {
      // Other dimensions
      int dim1 = (dim0 + 1) % 3;
      int dim2 = (dim0 + 2) % 3;

      int ie = 0;

      // Number of glops on a subface orthogonal to that direction
      size_t lsize = (deg[dim1] + 1) * (deg[dim2] + 1);
      // Number glops on the subfaces of the macrocell
      size_t wsize = lsize * (raf[dim0] - 1) * raf[dim1] * raf[dim2];

      int narg = 0;
      status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
      status |= clSetKernelArg(kernel, narg++, sizeof(int), &dim0);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
      status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real) * 4 * lsize * m, NULL);
      if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

      status = clEnqueueNDRangeKernel(queue, kernel,
                                      1, // work_dim
                                      NULL, // *global_work_offset
                                      &wsize, // *global_work_size
                                      &lsize, // *local_work_size
                                      nb_wait_events, // num_events_in_wait_list
                                      wait_event, // *wait_list
                                      event); // *event
      if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

      nb_wait_events = 1;
      wait_event = event;
      event++;
    }
  }
}

void DGSubCellInterface3D_OCL(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_args,
                             &m, deg, raf, physnode,
                             &Varindex, &NumFlux);
  free(cl_args);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);


  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
                                            &opencl_program,
                                            "DGFlux3D",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
                                      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                      60 * sizeof(schnaps_real),
                                      physnode,
                                      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
                                   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   7 * sizeof(int),
                                   param,
                                   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_event events[3];
  cl_event *event = events;
  cl_event *wait_event = NULL;
  int nb_wait_events = 0;

  // Loop over the dimensions
  for (int dim0 = 0; dim0 < 3; ++dim0) {
    // More than one subcell according to that dimension
    if (raf[dim0] > 1) {
      // Other dimensions
      int dim1 = (dim0 + 1) % 3;
      int dim2 = (dim0 + 2) % 3;

      int ie = 0;

      // Number of glops on a subface orthogonal to that direction
      size_t lsize[3];
      lsize[dim0] = 2;
      lsize[dim1] = deg[dim1] + 1;
      lsize[dim2] = deg[dim2] + 1;
      // Number glops on the subfaces of the macrocell
      size_t wsize[3];
      wsize[dim0] = (raf[dim0] - 1) * lsize[dim0];
      wsize[dim1] = raf[dim1] * lsize[dim1];
      wsize[dim2] = raf[dim2] * lsize[dim2];

      int narg = 0;
      status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
      status |= clSetKernelArg(kernel, narg++, sizeof(int), &dim0);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
      status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
      status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real) *
                               lsize[0] * lsize[1] * lsize[2] * m, NULL);
      if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

      status = clEnqueueNDRangeKernel(queue, kernel,
                                      3, // work_dim
                                      NULL, // *global_work_offset
                                      wsize, // *global_work_size
                                      lsize, // *local_work_size
                                      nb_wait_events, // num_events_in_wait_list
                                      wait_event, // *wait_list
                                      event); // *event
      if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

      nb_wait_events = 1;
      wait_event = event;
      event++;
    }
  }
}



void DGVolume_SPU(field* f) {

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.NumFlux, sizeof(fluxptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGVolume_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}

void DGVolume_SPU2(field* f, int ie) {

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);  
  
  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.NumFlux, sizeof(fluxptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGVolume_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;

  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[ie%nbgpu]);


  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void DGVolume_C(void *buffers[], void *cl_args);
void DGVolume_OCL(void *buffers[], void *cl_args);
void DGVolume3D_OCL(void *buffers[], void *cl_args);

void DGVolume_CUDA(void *buffers[], void *cl_args);


struct starpu_codelet* DGVolume_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgv_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgvolume"
  };

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGVolume";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;
    codelet.modes[codelet.nbuffers++] = STARPU_RW;
    //codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &dgv_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGVolume_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGVolume_OCL;
      codelet.opencl_funcs[ncodelet++] = DGVolume3D_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
        codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
      
    }

    if (starpu_cuda_use) {
      
      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGVolume_CUDA;;
      codelet.opencl_funcs[ncodelet++] = NULL;
      
      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
        codelet.opencl_flags[i] = STARPU_CUDA_ASYNC;
      
    }
    
  }

  return &codelet;
}

void DGAverage_C(void *buffers[], void *cl_args);
struct starpu_codelet* DGAverage_codelet();

void DGAverage_SPU(field* f) {

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGAverage_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


struct starpu_codelet* DGAverage_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dga_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgaverage"
  };

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGAverage";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_RW;
    //codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &dga_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGAverage_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    /* if (starpu_ocl_use) { */
    /*   LOAD_OPENCL_PROGRAM_SPU(); */

    /*   size_t ncodelet = 0; */
    /*   codelet.opencl_funcs[ncodelet++] = DGAverage_OCL; */
    /*   codelet.opencl_funcs[ncodelet++] = NULL; */

    /*   for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++) */
    /*     codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC; */
      
    /* } */

    /* if (starpu_cuda_use) { */
      
    /*   size_t ncodelet = 0; */
    /*   codelet.opencl_funcs[ncodelet++] = DGAverage_CUDA;; */
    /*   codelet.opencl_funcs[ncodelet++] = NULL; */
      
    /*   for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++) */
    /*     codelet.opencl_flags[i] = STARPU_CUDA_ASYNC; */
      
    /* } */
    
  }

  return &codelet;
}



void DGAverage_C(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
                             &Varindex); // ! fluxnum has been removed !!!
  free(cl_args);

  schnaps_real* w = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							  (struct starpu_vector_interface *)buffers[0]);


  int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const unsigned int npg_subcell = npg[0] * npg[1] * npg[2];

  // Loop over the subcells
  for (int icL0 = 0; icL0 < raf[0]; ++icL0) {
    for (int icL1 = 0; icL1 < raf[1]; ++icL1) {
      for (int icL2 = 0; icL2 < raf[2]; ++icL2) {
	int icL[3] = {icL0, icL1, icL2};

	// Get the left subcell id
	int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);

	// First glop index in the subcell
	int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	// compute all of the xref for the subcell
	schnaps_real *xref0 = malloc(npg_subcell * sizeof(schnaps_real));
	schnaps_real *xref1 = malloc(npg_subcell * sizeof(schnaps_real));
	schnaps_real *xref2 = malloc(npg_subcell * sizeof(schnaps_real));
	schnaps_real *omega = malloc(npg_subcell * sizeof(schnaps_real));
	int *imems = malloc(m * npg_subcell * sizeof(int));
	int pos = 0;
	for (int p = 0; p < npg_subcell; ++p) {
	  schnaps_real xref[3];
	  schnaps_real wpg;

	  ref_pg_vol(deg, raf, offsetL + p, xref, &wpg, NULL);
	  xref0[p] = xref[0];
	  xref1[p] = xref[1];
	  xref2[p] = xref[2];
	  omega[p] = wpg;

	  for (int iv = 0; iv < m; ++iv)
            imems[pos++] = Varindex(deg, raf, m, offsetL + p, iv);
	}

	schnaps_real vol_subcell = 0;
	schnaps_real wmean[m];
	for (int iv = 0; iv < m; ++iv) wmean[iv] = 0;
	
	// Loop over the subcell gauss points
	for (int p0 = 0; p0 < npg[0]; ++p0) {
	  for (int p1 = 0; p1 < npg[1]; ++p1) {
	    for (int p2 = 0; p2 < npg[2]; ++p2) {
	      schnaps_real wL[m];
	      int p[3] = {p0, p1, p2};
	      int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
	      
	      for (int iv = 0; iv < m; ++iv)
		wL[iv] = w[imems[m * (ipgL - offsetL) + iv]];
	      

	      schnaps_real xrefL[3] = {xref0[ipgL - offsetL],
				       xref1[ipgL - offsetL],
				       xref2[ipgL - offsetL]};
	      schnaps_real wpgL = omega[ipgL - offsetL];

	      schnaps_real dtau[3][3], codtau[3][3], dphiL[3];
	      schnaps_ref2phy(physnode,
			      xrefL,
			      NULL,
			      -1,  // ifa
			      NULL,  // xphy
			      dtau,
			      codtau,
			      NULL,
			      NULL);  // vnds

	      schnaps_real det = dot_product(dtau[0], codtau[0]);
	      
	      for (int iv = 0; iv < m; ++iv) {
		wmean[iv] += wL[iv] * wpgL * det;
	      }
	      vol_subcell += wpgL * det;
	    }  // p2
	  }  // p1
	}  // p0

	// Loop over the subcell gauss points
	for (int p0 = 0; p0 < npg[0]; ++p0) {
	  for (int p1 = 0; p1 < npg[1]; ++p1) {
	    for (int p2 = 0; p2 < npg[2]; ++p2) {
	      schnaps_real wL[m];
	      int p[3] = {p0, p1, p2};
	      int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
	      
	      for (int iv = 0; iv < m; ++iv)
		w[imems[m * (ipgL - offsetL) + iv]] = wmean[iv] / vol_subcell;
	    }  // p2
	  }  // p1
	}  // p0


	
	free(omega);
	free(xref0);
	free(xref1);
	free(xref2);
	free(imems);

      }  // icL2
    }  // icL1
  }  // icL0
}


void DGVolume_C(void *buffers[], void *cl_args) {

  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex, &NumFlux);
  free(cl_args);

  schnaps_real* w = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							  (struct starpu_vector_interface *)buffers[0]);
  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[1]);


  int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const unsigned int npg_subcell = npg[0] * npg[1] * npg[2];

  // Loop over the subcells
  for (int icL0 = 0; icL0 < raf[0]; ++icL0) {
    for (int icL1 = 0; icL1 < raf[1]; ++icL1) {
      for (int icL2 = 0; icL2 < raf[2]; ++icL2) {
	int icL[3] = {icL0, icL1, icL2};

	// Get the left subcell id
	int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);

	// First glop index in the subcell
	int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	// compute all of the xref for the subcell
	schnaps_real *xref0 = malloc(npg_subcell * sizeof(schnaps_real));
	schnaps_real *xref1 = malloc(npg_subcell * sizeof(schnaps_real));
	schnaps_real *xref2 = malloc(npg_subcell * sizeof(schnaps_real));
	schnaps_real *omega = malloc(npg_subcell * sizeof(schnaps_real));
	int *imems = malloc(m * npg_subcell * sizeof(int));
	int pos = 0;
	for (int p = 0; p < npg_subcell; ++p) {
	  schnaps_real xref[3];
	  schnaps_real wpg;

	  ref_pg_vol(deg, raf, offsetL + p, xref, &wpg, NULL);
	  xref0[p] = xref[0];
	  xref1[p] = xref[1];
	  xref2[p] = xref[2];
	  omega[p] = wpg;

	  for (int iv = 0; iv < m; ++iv)
	    imems[pos++] = Varindex(deg, raf, m, offsetL + p, iv);
	}

	// Loop over the "cross" in the three directions
	for (int dim0 = 0; dim0 < 3; ++dim0) {
	  for (int p0 = 0; p0 < npg[0]; ++p0) {
	    for (int p1 = 0; p1 < npg[1]; ++p1) {
	      for (int p2 = 0; p2 < npg[2]; ++p2) {
		schnaps_real wL[m], flux[m];
		int p[3] = {p0, p1, p2};
		int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);

		for (int iv = 0; iv < m; ++iv)
		  wL[iv] = w[imems[m * (ipgL - offsetL) + iv]];

		int q[3] = {p[0], p[1], p[2]};

		// Loop over the direction dim0 on the "cross"
		for (int iq = 0; iq < npg[dim0]; ++iq) {
		  q[dim0] = (p[dim0] + iq) % npg[dim0];

		  // Compute grad phi_q at glop p
		  schnaps_real dphiref[3] = {0, 0, 0};
		  dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) * raf[dim0];

		  schnaps_real xrefL[3] = {xref0[ipgL - offsetL],
					   xref1[ipgL - offsetL],
					   xref2[ipgL - offsetL]};
		  schnaps_real wpgL = omega[ipgL - offsetL];

		  schnaps_real dtau[3][3], codtau[3][3], dphiL[3];
		  schnaps_ref2phy(physnode,
				  xrefL,
				  dphiref,
				  -1,  // ifa
				  NULL,  // xphy
				  dtau,
				  codtau,
				  dphiL,
				  NULL);  // vnds

		  NumFlux(wL, wL, dphiL, flux);

		  int ipgR = offsetL + q[0] + npg[0] * (q[1] + npg[1] * q[2]);

		  for (int iv = 0; iv < m; ++iv) {
		    int imemR = Varindex(deg, raf, m, ipgR, iv);
		    int temp = m * (ipgR - offsetL) + iv;
		    assert(imemR == imems[temp]);
		    res[imems[temp]] += flux[iv] * wpgL;
		  }
		}  // iq
	      }  // p2
	    }  // p1
	  }  // p0
	}  // dim

	free(omega);
	free(xref0);
	free(xref1);
	free(xref2);
	free(imems);

      }  // icL2
    }  // icL1
  }  // icL0
}

void DGVolume_OCL(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex, &NumFlux);
  free(cl_args);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);

  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
					    &opencl_program,
					    "DGVolume",
					    devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
				      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				      60 * sizeof(schnaps_real),
				      physnode,
				      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
				   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   7 * sizeof(int),
				   param,
				   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;

  size_t lsize = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1);
  size_t wsize = lsize * raf[0] * raf[1] * raf[2];

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real) * 2 * lsize * m, NULL);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  status = clEnqueueNDRangeKernel(queue, kernel,
				  1, // work_dim
				  NULL, // *global_work_offset
				  &wsize, // *global_work_size
				  &lsize, // *local_work_size
				  0, // num_events_in_wait_list
				  NULL, // *wait_list
				  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}

void DGVolume_CUDA(void *buffers[], void *cl_args) {
  int m;
  //  int deg[3];
  //  int raf[3];
  //  schnaps_real physnode[20][3];
  //  varindexptr Varindex;
  //  fluxptr NumFlux;
  //
  //  starpu_codelet_unpack_args(cl_args,
  //			     &m, deg, raf, physnode,
  //                             &Varindex, &NumFlux);
  //  free(cl_args);
  //
  //  schnaps_real* w = (schnaps_real *)STARPU_VECTOR_GET_PTR(
  //      (struct starpu_vector_interface *)buffers[0]);
  //  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
  //      (struct starpu_vector_interface *)buffers[1]);
  //
  //
  //  CUfunction kernel;
  //  
  //  cudaError_t status = cudaGetDevice(wid);
  //  if (status != cudaSuccess) STARPU_CUDA_REPORT_ERROR(status);
  //
  //
  //  CUresult ret;
  //  
  //  CUdeviceptr physnode_cl;
  //  ret = cuMemAlloc(&physnode_cl,
  //                            60 * sizeof(schnaps_real));
  //  ret  = cuMemcpyHtoD (physnode_cl,
  //                                 physenode,
  //                                 60 * sizeof(schnaps_real));
  //  
  //  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  //  CUdeviceptr param_cl;
  //  ret = cuMemAlloc(&param_cl,
  //                             7 * sizeof(int));
  //  ret  = cuMemcpyHtoD(param_cl,
  //                                param,
  //                                7 * sizeof(int));
  //
  //  
  //  int ie = 0;
  //
  //  size_t lsize = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1);
  //  size_t wsize = lsize * raf[0] * raf[1] * raf[2];
  //
  //  int narg = 0;
  //
  //  ret = cuParamSetv(kernel, narg++, param_cl, sizeof(CUdeviceptr)); 	
  //  ret |= cuParamSeti(kernel, narg++ , ie);
  //  ret |= cuParamSetv(kernel, narg++, physnode_cl, sizeof(CUdeviceptr));
  //  ret |= cuParamSetv(kernel, narg++, w, sizeof(CUdeviceptr));
  //  ret |= cuParamSetv(kernel, narg++, res, sizeof(CUdeviceptr));
  //  ret |= cuParamSetv(kernel, narg++, NULL, sizeof(schnaps_real) * 2 * lsize * m);
  //
  //  
  //  ret = cuFuncSetBlockShape(kernel, lsize, 1, 1);
  //  ret = cuLaunchGrid(kernel,
  //                     wsize,
  //                     1); 	
}


void DGVolume3D_OCL(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  fluxptr NumFlux;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex, &NumFlux);
  free(cl_args);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);

  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
					    &opencl_program,
					    "DGVolume3D",
					    devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
				      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				      60 * sizeof(schnaps_real),
				      physnode,
				      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
				   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   7 * sizeof(int),
				   param,
				   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;

  size_t lsize[3] = {
    deg[0] + 1,
    deg[1] + 1,
    deg[2] + 1,
  };
  size_t wsize[3] = {
    raf[0] * lsize[0],
    raf[1] * lsize[1],
    raf[2] * lsize[2],
  };

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real) *
			   lsize[0] * lsize[1] * lsize[2] * m, NULL);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  status = clEnqueueNDRangeKernel(queue, kernel,
				  3, // work_dim
				  NULL, // *global_work_offset
				  wsize, // *global_work_size
				  lsize, // *local_work_size
				  0, // num_events_in_wait_list
				  NULL, // *wait_list
				  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}



void DGSource_SPU(field* f) {
  if (f->model.Source == NULL) return;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.Source, sizeof(sourceptr),
			   STARPU_VALUE, &f->tnow, sizeof(schnaps_real),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGSource_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}

void DGSource_SPU2(field* f, int ie) {

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER,
					    gpuids, STARPU_NMAXWORKERS);

  if (f->model.Source == NULL) return;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   STARPU_VALUE, &f->model.Source, sizeof(sourceptr),
			   STARPU_VALUE, &f->tnow, sizeof(schnaps_real),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGSource_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->wn_handle;
  task->handles[1] = f->res_handle;

  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[ie%nbgpu]);  

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void DGSource_C(void *buffers[], void *cl_args);
void DGSource_OCL(void *buffers[], void *cl_args);

struct starpu_codelet* DGSource_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgs_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgsource"
  };

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGSource";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;
    codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &dgs_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGSource_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGSource_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
	codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
    } // end if use openCL
  } // end if not is_init

  return &codelet;
}

void DGSource_C(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  sourceptr Source;
  schnaps_real tnow;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex, &Source, &tnow);
  free(cl_args);

  schnaps_real* w = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							  (struct starpu_vector_interface *)buffers[0]);
  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[1]);


  for (int ipg = 0; ipg < NPG(deg, raf); ++ipg) {
    schnaps_real dtau[3][3], codtau[3][3], xref[3], xphy[3], wpg;
    ref_pg_vol(deg, raf, ipg, xref, &wpg, NULL);
    schnaps_ref2phy(physnode,
		    xref,
		    NULL,  // dpsiref
		    -1,  // ifa
		    xphy,
		    dtau,
		    codtau,
		    NULL,  // dpsi
		    NULL);  // vnds
    schnaps_real det = dot_product(dtau[0], codtau[0]);  //// temp !!!!!!!!!!!!!!!
    schnaps_real wL[m], source[m];
    for (int iv = 0; iv < m; ++iv) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      wL[iv] = w[imem];
    }

    Source(xphy, tnow, wL, source);
    // printf("tnow=%f\n",tnow);

    for (int iv = 0; iv < m; ++iv) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      res[imem] += source[iv] * det * wpg;
    }
  }
}

void DGSource_OCL(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;
  sourceptr Source;
  schnaps_real tnow;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex, &Source, &tnow);
  free(cl_args);

  cl_mem w = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);

  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
					    &opencl_program,
					    "DGSource",
					    devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
				      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				      60 * sizeof(schnaps_real),
				      physnode,
				      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
				   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   7 * sizeof(int),
				   param,
				   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;

  size_t lsize = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1);
  size_t wsize = lsize * raf[0] * raf[1] * raf[2];

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real), &tnow);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &w);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real) * 2 * lsize * m, NULL);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  status = clEnqueueNDRangeKernel(queue, kernel,
				  1, // work_dim
				  NULL, // *global_work_offset
				  &wsize, // *global_work_size
				  &lsize, // *local_work_size
				  0, // num_events_in_wait_list
				  NULL, // *wait_list
				  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}



void DGMass_SPU(field* f) {
  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGMass_codelet(); 
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->res_handle;
  task->handles[1] = f->dtwn_handle;

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void DGMass_SPU2(field* f, int ie) {

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);
  
  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, &f->model.m, sizeof(int),
			   STARPU_VALUE, f->deg, 3 * sizeof(int),
			   STARPU_VALUE, f->raf, 3 * sizeof(int),
			   STARPU_VALUE, f->physnode, 60 * sizeof(schnaps_real),
			   STARPU_VALUE, &f->varindex, sizeof(varindexptr),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGMass_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  task->handles[0] = f->res_handle;
  task->handles[1] = f->dtwn_handle;

  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[ie%nbgpu]);  

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void DGMass_C(void *buffers[], void *cl_args);
void DGMass_OCL(void *buffers[], void *cl_args);

struct starpu_codelet* DGMass_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgm_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgmass"
  };

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGMass";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;
    codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;
    codelet.model = &dgm_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGMass_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGMass_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
	codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;

    }
  }

  return &codelet;
}

void DGMass_C(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex);
  free(cl_args);

  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[0]);
  schnaps_real* dtw = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[1]);


  for (int ipg = 0; ipg < NPG(deg, raf); ++ipg) {
    schnaps_real dtau[3][3], codtau[3][3], xref[3], xphy[3], wpg;
    ref_pg_vol(deg, raf, ipg, xref, &wpg, NULL);
    schnaps_ref2phy(physnode,
		    xref,
		    NULL,  // dpsiref
		    -1,  // ifa
		    xphy,
		    dtau,
		    codtau,
		    NULL,  // dpsi
		    NULL);  // vnds
    schnaps_real det = dot_product(dtau[0], codtau[0]);
    for (int iv = 0; iv < m; ++iv) {
      int imem = Varindex(deg, raf, m, ipg, iv);
      dtw[imem] = res[imem] / (wpg * det) - dtw[imem] / 2 ;
    }
  }
}

void DGMass_OCL(void *buffers[], void *cl_args) {
  int m;
  int deg[3];
  int raf[3];
  schnaps_real physnode[20][3];
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, physnode,
			     &Varindex);
  free(cl_args);

  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[0]);
  cl_mem dtw = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[1]);

  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
					    &opencl_program,
					    "DGMassRes",
					    devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  cl_mem physnode_cl = clCreateBuffer(context,
				      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				      60 * sizeof(schnaps_real),
				      physnode,
				      &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
				   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   7 * sizeof(int),
				   param,
				   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;


  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &physnode_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &dtw);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  size_t wsize = NPG(deg, raf);

  status = clEnqueueNDRangeKernel(queue, kernel,
				  1, // work_dim
				  NULL, // *global_work_offset
				  &wsize, // *global_work_size
				  NULL, // *local_work_size
				  0, // num_events_in_wait_list
				  NULL, // *wait_list
				  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}



void DGMacroCellInterface_SPU(Interface* inter, int side) {
  if (inter->fL == NULL || inter->fR == NULL) return;

  field *f;
  int locfa;
  starpu_data_handle_t index;
  starpu_data_handle_t wn_in;
  starpu_data_handle_t wn_ext;

  if (side == 0) {
    f = inter->fL;
    locfa = inter->locfaL;
    index = inter->vol_indexL_handle;
    wn_in = inter->wL_handle;
    wn_ext = inter->wR_handle;

  } else if (side == 1) {
    f = inter->fR;
    locfa = inter->locfaR;
    index = inter->vol_indexR_handle;
    wn_in = inter->wR_handle;
    wn_ext = inter->wL_handle;

  } else {
    assert(1==2);
  }

  const int sign = 1 - 2 * side;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, f, sizeof(field),
			   STARPU_VALUE, &locfa, sizeof(int),
			   STARPU_VALUE, &sign, sizeof(int),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGMacroCellInterface_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  int nhandle = 0;
  task->handles[nhandle++] = index;
  task->handles[nhandle++] = inter->vnds_handle;
  task->handles[nhandle++] = inter->wpg_handle;
  task->handles[nhandle++] = wn_in;
  task->handles[nhandle++] = wn_ext;
  if (f->solver != NULL) {
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->rhs_handle;
  } else {
    task->handles[nhandle++] = f->res_handle;
  }

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}


void DGMacroCellInterface_SPU2(Interface* inter, int side) {
  if (inter->fL == NULL || inter->fR == NULL) return;

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);
  
  field *f;
  int locfa;
  starpu_data_handle_t index;
  starpu_data_handle_t wn_in;
  starpu_data_handle_t wn_ext;

  if (side == 0) {
    f = inter->fL;
    locfa = inter->locfaL;
    index = inter->vol_indexL_handle;
    wn_in = inter->wL_handle;
    wn_ext = inter->wR_handle;

  } else if (side == 1) {
    f = inter->fR;
    locfa = inter->locfaR;
    index = inter->vol_indexR_handle;
    wn_in = inter->wR_handle;
    wn_ext = inter->wL_handle;

  } else {
    assert(1==2);
  }

  const int sign = 1 - 2 * side;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, f, sizeof(field),
			   STARPU_VALUE, &locfa, sizeof(int),
			   STARPU_VALUE, &sign, sizeof(int),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGMacroCellInterface_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  int nhandle = 0;
  task->handles[nhandle++] = index;
  task->handles[nhandle++] = inter->vnds_handle;
  task->handles[nhandle++] = inter->wpg_handle;
  task->handles[nhandle++] = wn_in;
  task->handles[nhandle++] = wn_ext;
  if (f->solver != NULL) {
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->rhs_handle;
  } else {
    task->handles[nhandle++] = f->res_handle;
  }

  int sendto;
  sendto = side*(inter->ieR) + (1-side)*(inter->ieL);
  
  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[sendto%nbgpu]);  

  
  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}



void DGMacroCellInterface_C(void* buffers[], void* cl_args);
void DGMacroCellInterface_OCL(void* buffers[], void* cl_args);

struct starpu_codelet* DGMacroCellInterface_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgmci_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgmacrocellinterface"
  }; 
  

  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGMacroCellInterface";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // index
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // vnds
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // wpg
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // wn_in
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // wn_ext
    codelet.modes[codelet.nbuffers++] = STARPU_RW;  // res
    //codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;  // res
    codelet.model = &dgmci_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGMacroCellInterface_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGMacroCellInterface_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
	codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
      
    }
  }

  return &codelet;
}

void DGMacroCellInterface_C(void* buffers[], void* cl_args) {
  field f0;
  field *f = &f0;
  int locfa;
  int sign;

  starpu_codelet_unpack_args(cl_args,
			     f, &locfa, &sign);
  free(cl_args);

  int nbuffer = 0;
  int* index = (int *)STARPU_VECTOR_GET_PTR(
					    (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* vnds_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								 (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* wpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								(struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* wn_in = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							      (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* wn_ext = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							       (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[nbuffer++]);


  // Warning: refinement has to be the same for the moment
  const int npgf = NPGF(f->deg, f->raf, locfa);
  const int m = f->model.m;

  for (int ipgf = 0; ipgf < npgf; ++ipgf) {
    schnaps_real wL[m];
    schnaps_real wR[m];
    for (int iv = 0; iv < m; ++iv) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      wL[iv] = wn_in[imemf];
      wR[iv] = wn_ext[imemf];
    }

    schnaps_real vnds[3];
    vnds[0] = sign * vnds_buf[3 * ipgf + 0];
    vnds[1] = sign * vnds_buf[3 * ipgf + 1];
    vnds[2] = sign * vnds_buf[3 * ipgf + 2];

    schnaps_real flux[m];
    f->model.NumFlux(wL, wR, vnds, flux);

    // Add flux to the selected side
    int ipgL = index[ipgf];
    const schnaps_real wpg = wpg_buf[ipgf];
    for (int iv = 0; iv < m; ++iv) {
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      // Warning: no wpg because already applied in vnds buffer
      res[imemL] -= flux[iv] * wpg;
    }
    /* printf("ipgf=%d ipgL=%d wL=%f wR=%f vnds=%f %f %f flux=%f\n", */
    /*        ipgf, ipgL, wL[0], wR[0], vnds[0], vnds[1], vnds[2], flux[0]); */
  }
}

void DGMacroCellInterface_OCL(void* buffers[], void* cl_args) {
  field f0;
  field *f = &f0;
  int locfa;
  int sign;

  starpu_codelet_unpack_args(cl_args,
			     f, &locfa, &sign);
  free(cl_args);

  int nbuffer = 0;
  cl_mem index = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem vnds_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem wpg_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem wn_in = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem wn_ext = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);

  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
					    &opencl_program,
					    "DGMacroCellInterfaceRes",
					    devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {f->model.m, f->deg[0], f->deg[1], f->deg[2], f->raf[0], f->raf[1], f->raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
				   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   7 * sizeof(int),
				   param,
				   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &locfa);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &sign);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &index);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn_in);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn_ext);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &vnds_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wpg_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  // Number of glops on the macrocell interface
  size_t wsize = NPGF(f->deg, f->raf, locfa);

  status = clEnqueueNDRangeKernel(queue, kernel,
				  1, // work_dim
				  NULL, // *global_work_offset
				  &wsize, // *global_work_size
				  NULL, // *local_work_size
				  0, // num_events_in_wait_list
				  NULL, // *wait_list
				  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}



void DGMacroCellBoundaryFlux_SPU(Interface* inter) {
  if (inter->fR != NULL) return;

  field *f = inter->fL;
  int locfa = inter->locfaL;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, f, sizeof(field),
			   STARPU_VALUE, &locfa, sizeof(int),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGMacroCellBoundaryFlux_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  int nhandle = 0;
  task->handles[nhandle++] = inter->vol_indexL_handle;
  task->handles[nhandle++] = inter->vnds_handle;
  task->handles[nhandle++] = inter->xpg_handle;
  task->handles[nhandle++] = inter->wpg_handle;
  task->handles[nhandle++] = inter->wL_handle;
  if (f->solver != NULL) {
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->rhs_handle;
  }
  else {
    task->handles[nhandle++] = f->res_handle;
  }

  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}

void DGMacroCellBoundaryFlux_SPU2(Interface* inter) {
  if (inter->fR != NULL) return;

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);
  
  field *f = inter->fL;
  int locfa = inter->locfaL;

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
			   STARPU_VALUE, f, sizeof(field),
			   STARPU_VALUE, &locfa, sizeof(int),
			   NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = DGMacroCellBoundaryFlux_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  int nhandle = 0;
  task->handles[nhandle++] = inter->vol_indexL_handle;
  task->handles[nhandle++] = inter->vnds_handle;
  task->handles[nhandle++] = inter->xpg_handle;
  task->handles[nhandle++] = inter->wpg_handle;
  task->handles[nhandle++] = inter->wL_handle;
  if (f->solver != NULL) {
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->rhs_handle;
  }
  else {
    task->handles[nhandle++] = f->res_handle;
  }

  int sendto;
  sendto = inter->ieL;
  
  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[sendto%nbgpu]);  

  
  STARPU_CHECK_RETURN_VALUE(
			    starpu_task_submit(task),
			    "starpu_task_submit");
}

void DGMacroCellBoundaryFlux_C(void* buffers[], void* cl_args);
void DGMacroCellBoundaryFlux_OCL(void* buffers[], void* cl_args);

struct starpu_codelet* DGMacroCellBoundaryFlux_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel dgmcbf_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dgmacrocellboundaryflux"
  }; 


  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "DGMacroCellBoundaryFlux";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // index
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // vnds
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // xpg
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // wpg
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // wn
    codelet.modes[codelet.nbuffers++] = STARPU_RW;  // res
    //codelet.modes[codelet.nbuffers++] = STARPU_RW|STARPU_COMMUTE;  // res
    codelet.model = &dgmcbf_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = DGMacroCellBoundaryFlux_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = DGMacroCellBoundaryFlux_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;

      for(int i=0; i< STARPU_MAXIMPLEMENTATIONS; i++)
	codelet.opencl_flags[i] = STARPU_OPENCL_ASYNC;
      
    }
  }

  return &codelet;
}

void DGMacroCellBoundaryFlux_C(void* buffers[], void* cl_args) {
  field f0;
  field *f = &f0;
  int locfa;

  starpu_codelet_unpack_args(cl_args,
			     f, &locfa);
  free(cl_args);

  int nbuffer = 0;
  int* index = (int *)STARPU_VECTOR_GET_PTR(
					    (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* vnds_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								 (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* xpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								(struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* wpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(
								(struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* wn = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							   (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* res = (schnaps_real *)STARPU_VECTOR_GET_PTR(
							    (struct starpu_vector_interface *)buffers[nbuffer++]);


  const int npgf = NPGF(f->deg, f->raf, locfa);
  const int m = f->model.m;

  for (int ipgf = 0; ipgf < npgf; ++ipgf) {
    schnaps_real wL[m];
    for (int iv = 0; iv < m; ++iv) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      wL[iv] = wn[imemf];
    }

    schnaps_real* xpg = xpg_buf + 3 * ipgf;
    schnaps_real* vnds = vnds_buf + 3 * ipgf;

    schnaps_real flux[m];
    f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);

    // Add flux
    int ipgL = index[ipgf];
    const schnaps_real wpg = wpg_buf[ipgf];
    for (int iv = 0; iv < m; ++iv) {
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv);
      res[imemL] -= flux[iv] * wpg;
    }
  }
}

void DGMacroCellBoundaryFlux_OCL(void* buffers[], void* cl_args) {
  field f0;
  field *f = &f0;
  int locfa;

  starpu_codelet_unpack_args(cl_args,
			     f, &locfa);
  free(cl_args);

  int nbuffer = 0;
  cl_mem index = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem vnds_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem xpg_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem wpg_buf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem wn = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem res = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);

  int wid = starpu_worker_get_id();
  int devid = starpu_worker_get_devid(wid);
  cl_context context;
  starpu_opencl_get_context(devid, &context);

  // Usefull to determine device type on runtime
  //char buffer[100];
  //starpu_worker_get_name(wid, buffer, 100);
  //switch (GetDeviceTypeFromName(buffer)) {
  //  case CL_DEVICE_TYPE_CPU:
  //
  //    break;
  //  default:
  //
  //    break;
  //}

  cl_kernel kernel;
  cl_command_queue queue;
  cl_int status = starpu_opencl_load_kernel(&kernel, &queue,
					    &opencl_program,
					    "DGBoundaryRes",
					    devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {f->model.m, f->deg[0], f->deg[1], f->deg[2], f->raf[0], f->raf[1], f->raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
				   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				   7 * sizeof(int),
				   param,
				   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int ie = 0;

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(schnaps_real), &(f->tnow));
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &ie);
  status |= clSetKernelArg(kernel, narg++, sizeof(int), &locfa);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &index);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &vnds_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &xpg_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wpg_buf);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &res);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  // Number of glops on the macrocell interface
  size_t wsize = NPGF(f->deg, f->raf, locfa);

  status = clEnqueueNDRangeKernel(queue, kernel,
				  1, // work_dim
				  NULL, // *global_work_offset
				  &wsize, // *global_work_size
				  NULL, // *local_work_size
				  0, // num_events_in_wait_list
				  NULL, // *wait_list
				  NULL); // *event
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);
}







void InterfaceExplicitFlux_bis(Interface* inter, int side){


  field *f;
  field *fext;
  int locfa;
  int *index_ext;
  int *index;
  schnaps_real* wn;
  schnaps_real* wn_ext;

  const int sign = 1 - 2 * side;

  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    index_ext = inter->vol_indexR;
    index = inter->vol_indexL;
    wn = inter->wL;
    wn_ext = inter->wR;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    index_ext = inter->vol_indexL;
    index = inter->vol_indexR;
    wn = inter->wR;
    wn_ext = inter->wL;
  }
  else {
    assert(1==2);
  }

  int npgf = NPGF(f->deg, f->raf, locfa);

  if (f != NULL){

    const unsigned int m = f->model.m;

    //printf("locfa=%d \n",locfa);

    for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {


      schnaps_real flux[m];
      schnaps_real wL[m];
      int ipgL = index[ipgf];
      for(int iv = 0; iv < m; iv++) {
	int imemf = VarindexFace(npgf, m, ipgf, iv);
	wL[iv] = wn[imemf];
	//wL[iv] = 0;
      }

      const schnaps_real wpg = inter->wpg[ipgf];

      if (fext != NULL) {  // the right element exists
	schnaps_real wR[m];
	int ipgR = index_ext[ipgf];
	//int ipgL = index[ipgf];
	for(int iv = 0; iv < m; iv++) {
	  int imemf = VarindexFace(npgf, m, ipgf, iv);
	  wR[iv] = wn_ext[imemf];
	  //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}

	// int_dL F(wL, wR, grad phi_ib)
	schnaps_real vndsloc[3];

	vndsloc[0] = sign * inter->vnds[3 * ipgf + 0];
	vndsloc[1] = sign * inter->vnds[3 * ipgf + 1];
	vndsloc[2] = sign * inter->vnds[3 * ipgf + 2];

	f->model.NumFlux(wL, wR, vndsloc, flux);
	//printf("flux=%f %f\n",flux[0],flux[0]);

	// Add flux  to the selected side
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  f->dtwn[imemL] -= flux[iv] * wpg;
	  //printf("imem=%d res=%f\n",imemL,f->dtwn[imemL]);
	}
      } else { // The point is on the boundary.


	assert(sign == 1);
	schnaps_real* xpg = inter->xpg + 3 * ipgf;
	schnaps_real* vnds = inter->vnds + 3 * ipgf;

	//printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
	//printf("wL=%f, ipgf=%d\n",wL[0], ipgf);
	//printf("flux=%f, ipgf=%d\n",flux[0], ipgf);
	int ipgL = index[ipgf];
	/* printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d \n", */
	/*        xpg[0], xpg[1], xpg[2], */
	/*        vnds[0], vnds[1],vnds[2], ipgL); */
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  f->dtwn[imemL] -= flux[iv] * sign * wpg;
	  /* printf("boundary flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      }

    }


  }
}


void DtFields_bis(Simulation *simu,
		  schnaps_real* w,
		  schnaps_real* dtw){

  if(simu->pre_dtfields != NULL) {
    simu->pre_dtfields(simu);
  }


#ifdef _OPENMP
#pragma omp parallel
#endif

  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int iw = 0; iw < simu->wsize; iw++)
    dtw[iw] = 0;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }

  // the field pointers must be updated
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].wn = w + ie * fsize;
    simu->fd[ie].dtwn = dtw + ie * fsize;
  }


  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    // left = 0  right = 1
    ExtractInterface(inter, 0);
    ExtractInterface(inter, 1);
    if (inter->fR != NULL) {
      InterfaceExplicitFlux_bis(inter, 0);
      InterfaceExplicitFlux_bis(inter, 1);
    }
    else{
      InterfaceExplicitFlux_bis(inter, 0);
      //InterfaceBoundaryFlux(inter);
    }
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

  if(simu->post_dtfields != NULL) {
    simu->post_dtfields(simu);
  }

}


void DtFields_SPU(Simulation *simu,
		  starpu_data_handle_t* w_handle,
		  starpu_data_handle_t* dtw_handle)
{

  /* if(simu->pre_dtfields != NULL) { */
  /*   simu->pre_dtfields(simu, w); */
  /* } */


  int fsize =  simu->wsize / simu->macromesh.nbelems;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
#ifndef WORKER_ON_NODE
    ZeroBuffer_SPU(simu->fd[ie].res_handle);
#else
    ZeroBuffer_SPU2(simu->fd[ie].res_handle, ie);
#endif
  }

 

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    simu->fd[ie].tnow = simu->tnow;
  }

  // the field pointers must be updated
  if (w_handle != NULL){
    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      simu->fd[ie].wn_handle = w_handle[ie];
      simu->fd[ie].dtwn_handle = dtw_handle[ie];
      simu->fd[ie].res_handle = simu->res_handle[ie];
    }
  }
  else {
    assert(dtw_handle == NULL);
  }

  
#ifndef WORKER_ON_NODE
  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    // left = 0  right = 1
    ExtractInterface_SPU(inter, 0);
    ExtractInterface_SPU(inter, 1);

    if (inter->fR != NULL) {
      DGMacroCellInterface_SPU(inter, 0);
      DGMacroCellInterface_SPU(inter, 1);
    }
    else{
      DGMacroCellBoundaryFlux_SPU(inter);
    }
  }
#else
  assert(1==2);
  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++){
    Interface* inter = simu->interface + ifa;
    // left = 0  right = 1
    ExtractInterface_SPU2(inter, 0);
    ExtractInterface_SPU2(inter, 1);

    if (inter->fR != NULL) {
      DGMacroCellInterface_SPU2(inter, 0);
      DGMacroCellInterface_SPU2(inter, 1);
    }
    else{
      DGMacroCellBoundaryFlux_SPU2(inter);
    }
  }
#endif

  
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
#ifndef WORKER_ON_NODE
    DGSubCellInterface_SPU(simu->fd + ie);
    DGVolume_SPU(simu->fd + ie);
    DGSource_SPU(simu->fd + ie);
    DGMass_SPU(simu->fd + ie);
#else
    assert(1==2);
    DGSubCellInterface_SPU2(simu->fd + ie, ie);
    DGVolume_SPU2(simu->fd + ie, ie);
    DGSource_SPU2(simu->fd + ie, ie);
    DGMass_SPU2(simu->fd + ie, ie);
#endif
  }

}




void RK2_SPU(Simulation *simu, schnaps_real tmax){

  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;

  simu->tmax = tmax;

  simu->itermax_rk = tmax / simu->dt + 1;
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  int iter = 0;

  // FIXME: remove
  //size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

  /* if(simu->nb_diags != 0) { */
  /*   simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real)); */
  /* } */

  simu->tnow = 0;

  assert(starpu_use);
  RegisterSimulation_SPU(simu);
  
  /* printf("Init macro elems data"); */
  /* for(int ie = 0; ie < simu->macromesh.nbelems; ++ie){ */
  /*   printf(" %d",ie); */
  /*   DGInitField_SPU(simu->fd + ie); */
  /* } */
  /* printf("\n"); */

  /* starpu_task_wait_for_all() ; */
  printf("Start time steps...\n");
  
  // top chrono
  schnaps_real tps_debut = seconds();
  starpu_profiling_init();
  while(simu->tnow < tmax) {
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
#ifndef WORKER_ON_NODE
      ZeroBuffer_SPU(simu->fd[ie].dtwn_handle);
#else
      ZeroBuffer_SPU2(simu->fd[ie].dtwn_handle, ie);
#endif
    }

    DtFields_SPU(simu, NULL, NULL);
    //UnregisterSimulation_SPU(simu);
    /* for(int i=0 ; i<simu->wsize; i++) */
    /*   printf("i=%d dtw=%f\n",i,simu->dtw[i]); */
    /* assert(1==2); */

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      schnaps_real alpha = dt/2;
#ifndef WORKER_ON_NODE
      AddBuffer_SPU(alpha, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle);
#else
      AddBuffer_SPU2(alpha, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle, ie);
#endif
    }

    //UnregisterSimulation_SPU(simu);
    /* for(int i=0 ; i<simu->wsize; i++) */
    /*   //printf("i=%d w=%f\n",i,simu->w[i]+dt/2*simu->dtw[i]); */
    /*   printf("i=%d w=%f\n",i,simu->w[i]); */
    /* assert(1==2); */
    simu->tnow += 0.5 * dt;

    DtFields_SPU(simu, NULL, NULL);

    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
#ifndef WORKER_ON_NODE
      AddBuffer_SPU(simu->dt, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle);
#else
      AddBuffer_SPU2(simu->dt, simu->fd[ie].dtwn_handle, simu->fd[ie].wn_handle, ie);
      assert(1==2);
#endif
    }

    simu->tnow += 0.5 * dt;

    /* if(simu->update_after_rk != NULL){  */
    /*   simu->update_after_rk(simu, simu->w);  */
    /* } */

    iter++;
    simu->iter_time_rk = iter;

    //starpu_task_wait_for_all() ;

  }
starpu_task_wait_for_all() ;
  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);



  // top chrono
  schnaps_real tps_fin = seconds();
  printf("Temps total (no memory transfer) =%f\n", tps_fin - tps_debut);

  // The memory transfers are not counted
  UnregisterSimulation_SPU(simu);

}



void SmartPrefetch_SPU(Simulation *simu){

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);

  int nbdev = starpu_worker_get_count();

  for(int i=0; i<nbgpu; i++){
    printf("nbgpu=%d nbdev=%d - i=%d gpuids=%d \n", nbgpu, nbdev, i, gpuids[i]);
  }
  
  if (!starpu_is_init && starpu_use){
    int ret;
    ret = starpu_init(NULL);
    assert(ret != -ENODEV) ;
    starpu_is_init = true;
    init_global_arbiter();
  }

  for(int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    if (starpu_use && !simu->fd[ie].starpu_registered){
      
      starpu_vector_data_register(&(simu->fd[ie].wn_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->fd[ie].wn), // vector location
				  simu->fd[ie].wsize,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->fd[ie].wn_handle);
      
      starpu_vector_data_register(&(simu->fd[ie].dtwn_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->fd[ie].dtwn), // vector location
				  simu->fd[ie].wsize,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->fd[ie].dtwn_handle);

      starpu_vector_data_register(&(simu->fd[ie].res_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->fd[ie].res), // vector location
				  simu->fd[ie].wsize,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->fd[ie].res_handle);
            
    
      simu->fd[ie].starpu_registered = true;
    
    }

    int sendto = ie%nbgpu;
    starpu_data_prefetch_on_node(simu->fd[ie].wn_handle,
				 starpu_worker_get_memory_node(gpuids[sendto]),
				 0); // blocking
    starpu_data_prefetch_on_node(simu->fd[ie].dtwn_handle,
				 starpu_worker_get_memory_node(gpuids[sendto]),
				 0); // blocking

    starpu_data_prefetch_on_node(simu->fd[ie].res_handle,
				 starpu_worker_get_memory_node(gpuids[sendto]),
				 0); // blocking
    
  }



  for(int ifa = 0; ifa < simu->macromesh.nbfaces; ifa++) {

    if (starpu_use && !simu->interface[ifa].starpu_registered){

      starpu_vector_data_register(&(simu->interface[ifa].vol_indexL_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->interface[ifa].vol_indexL), // vector location
				  simu->interface[ifa].npgL,  // size
				  sizeof(int));  // type
	register_data_arbiter(simu->interface[ifa].vol_indexL_handle);
      
      if (simu->interface[ifa].vol_indexR != NULL) {
	starpu_vector_data_register(&(simu->interface[ifa].vol_indexR_handle), // mem handle
				    0, // location: CPU
				    (uintptr_t)(simu->interface[ifa].vol_indexR), // vector location
				    simu->interface[ifa].npgR,  // size
				    sizeof(int));  // type
	register_data_arbiter(simu->interface[ifa].vol_indexR_handle);
      
      }
      
      starpu_vector_data_register(&(simu->interface[ifa].wL_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->interface[ifa].wL), // vector location
				  simu->interface[ifa].wsizeL,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->interface[ifa].wL_handle);
      
      if (simu->interface[ifa].wR != NULL) {
	starpu_vector_data_register(&(simu->interface[ifa].wR_handle), // mem handle
				    0, // location: CPU
				    (uintptr_t)(simu->interface[ifa].wR), // vector location
				    simu->interface[ifa].wsizeR,  // size
				    sizeof(schnaps_real));  // type
	register_data_arbiter(simu->interface[ifa].wR_handle);
      
      }
      
      starpu_vector_data_register(&(simu->interface[ifa].vnds_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->interface[ifa].vnds), // vector location
				  simu->interface[ifa].npgL * 3,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->interface[ifa].vnds_handle);
      
      starpu_vector_data_register(&(simu->interface[ifa].xpg_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->interface[ifa].xpg), // vector location
				  simu->interface[ifa].npgL * 3,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->interface[ifa].xpg_handle);
      
      starpu_vector_data_register(&(simu->interface[ifa].wpg_handle), // mem handle
				  0, // location: CPU
				  (uintptr_t)(simu->interface[ifa].wpg), // vector location
				  simu->interface[ifa].npgL,  // size
				  sizeof(schnaps_real));  // type
	register_data_arbiter(simu->interface[ifa].wpg_handle);
      
      
      
      simu->interface[ifa].starpu_registered = true;
    }

    for(int gpu=0; gpu<nbgpu; gpu++){

      //starpu_data_prefetch_on_node(simu->interface[ifa].vol_indexL_handle,
      //                             starpu_worker_get_memory_node(gpu),
      //                             0); // blocking
      //
      //starpu_data_prefetch_on_node(simu->interface[ifa].vol_indexR_handle,
      //                             starpu_worker_get_memory_node(gpu),
      //                             0); // blocking
      //
      //starpu_data_prefetch_on_node(simu->interface[ifa].wL_handle,
      //                             starpu_worker_get_memory_node(gpu),
      //                             0); // blocking
      //
      //starpu_data_prefetch_on_node(simu->interface[ifa].wR_handle,
      //                             starpu_worker_get_memory_node(gpu),
      //                             0); // blocking

      
      starpu_data_prefetch_on_node(simu->interface[ifa].vnds_handle,
				   starpu_worker_get_memory_node(gpu),
				   0); // blocking

      starpu_data_prefetch_on_node(simu->interface[ifa].xpg_handle,
				   starpu_worker_get_memory_node(gpu),
				   0); // blocking

      starpu_data_prefetch_on_node(simu->interface[ifa].wpg_handle,
				   starpu_worker_get_memory_node(gpu),
				   0); // blocking
    }
    
  }

}

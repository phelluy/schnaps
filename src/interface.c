#include "interface.h"
#include <assert.h>
#include <stdio.h>

#pragma start_opencl
int VarindexFace(int npg, int m, int ipgf, int iv){

  return ipgf * m + iv;

}
#pragma end_opencl

void RegisterInterface_SPU(Interface* inter){

  if (!starpu_is_init && starpu_use){
    assert(starpu_init(NULL) != -ENODEV);
    starpu_is_init = true;
    init_global_arbiter();
  }

  if (starpu_use && !inter->starpu_registered) {

    starpu_vector_data_register(&(inter->vol_indexL_handle), // mem handle
        			0, // location: CPU
        			(uintptr_t)(inter->vol_indexL), // vector location
        			inter->npgL,  // size
        			sizeof(int));  // type
	register_data_arbiter(inter->vol_indexL_handle);
    
    if (inter->vol_indexR != NULL) {
      starpu_vector_data_register(&(inter->vol_indexR_handle), // mem handle
        			  0, // location: CPU
        			  (uintptr_t)(inter->vol_indexR), // vector location
        			  inter->npgR,  // size
        			  sizeof(int));  // type
	register_data_arbiter(inter->vol_indexR_handle);
    }

    starpu_vector_data_register(&(inter->wL_handle), // mem handle
        			0, // location: CPU
        			(uintptr_t)(inter->wL), // vector location
        			inter->wsizeL,  // size
    				sizeof(schnaps_real));  // type
	register_data_arbiter(inter->wL_handle);

    if (inter->wR != NULL) {
      starpu_vector_data_register(&(inter->wR_handle), // mem handle
        			  0, // location: CPU
        			  (uintptr_t)(inter->wR), // vector location
        			  inter->wsizeR,  // size
        			  sizeof(schnaps_real));  // type
	register_data_arbiter(inter->wR_handle);
    }

    starpu_vector_data_register(&(inter->vnds_handle), // mem handle
        			0, // location: CPU
        			(uintptr_t)(inter->vnds), // vector location
        			inter->npgL * 3,  // size  !!!!!!!!!!! same for left and right ????
        			sizeof(schnaps_real));  // type
	register_data_arbiter(inter->vnds_handle);

    starpu_vector_data_register(&(inter->xpg_handle), // mem handle
        			0, // location: CPU
        			(uintptr_t)(inter->xpg), // vector location
        			inter->npgL * 3,  // size  !!!!!!!!!!! same for left and right ????
        			sizeof(schnaps_real));  // type
	register_data_arbiter(inter->xpg_handle);

    starpu_vector_data_register(&(inter->wpg_handle), // mem handle
        			0, // location: CPU
        			(uintptr_t)(inter->wpg), // vector location
        			inter->npgL,  // size  !!!!!!!!!!! same for left and right ????
        			sizeof(schnaps_real));  // type
	register_data_arbiter(inter->wpg_handle);

    inter->starpu_registered = true;
  }
  
}


void UnregisterInterface_SPU(Interface* inter) {
  if (starpu_use && inter->starpu_registered) {
    starpu_data_unregister(inter->vol_indexL_handle);
    /* if (inter->vol_indexR_handle != NULL) */
    /*   starpu_data_unregister(inter->vol_indexR_handle); */
    starpu_data_unregister(inter->wL_handle);
    /* if (inter->wR_handle != NULL) */
    /*   starpu_data_unregister(inter->wR_handle); */
    starpu_data_unregister(inter->vnds_handle);
    starpu_data_unregister(inter->xpg_handle);
    starpu_data_unregister(inter->wpg_handle);

    inter->starpu_registered = false;
  }
}


void ExtractInterface_SPU_with_prio(Interface* inter, int side, int prio) {
    assert(side == 1 || side == 0);

    field *f;
    int locfa, npgf;
    starpu_data_handle_t w_handle;
    starpu_data_handle_t vol_index_handle;

    if (side == 0) {
      if (inter->fL == NULL) return;
      f = inter->fL;
      locfa = inter->locfaL;
      npgf = inter->npgL;
      w_handle = inter->wL_handle;
      vol_index_handle = inter->vol_indexL_handle;
    } else  {
      if (inter->fR == NULL) return;
      f = inter->fR;
      locfa = inter->locfaR;
      npgf = inter->npgR;
      w_handle = inter->wR_handle;
      vol_index_handle = inter->vol_indexR_handle;
    }

    void* args;
    size_t arg_size;

    starpu_codelet_pack_args(&args, &arg_size,
                             STARPU_VALUE, &f->model.m, sizeof(int),
                             STARPU_VALUE, f->deg, 3 * sizeof(int),
                             STARPU_VALUE, f->raf, 3 * sizeof(int),
                             STARPU_VALUE, &npgf, sizeof(int),
                             STARPU_VALUE, &f->varindex, sizeof(varindexptr),
                             NULL);

    struct starpu_task *task;

    task = starpu_task_create();
    task->cl = ExtractInterface_codelet();
    task->cl_arg = args;
    task->cl_arg_size = arg_size;
    int nhandle = 0;
    task->handles[nhandle++] = w_handle;
    if (f->solver != NULL){
      Skyline_SPU* sky_spu = f->solver->matrix;
      task->handles[nhandle++] = sky_spu->sol_handle;
    } else {
      task->handles[nhandle++] = f->wn_handle;
    }
    task->handles[nhandle++] = vol_index_handle;

    task->priority = prio;

    STARPU_CHECK_RETURN_VALUE(
        starpu_task_submit(task),
        "starpu_task_submit");
}

void ExtractInterface_SPU(Interface* inter, int side) {
  ExtractInterface_SPU_with_prio(inter, side, 0);
}

void ExtractInterface_SPU2(Interface* inter, int side) {
  assert(side == 1 || side == 0);

  int gpuids[STARPU_NMAXWORKERS];
  int nbgpu = starpu_worker_get_ids_by_type(STARPU_OPENCL_WORKER, gpuids, STARPU_NMAXWORKERS);
  
  field *f;
  int locfa, npgf;
  starpu_data_handle_t w_handle;
  starpu_data_handle_t vol_index_handle;

  if (side == 0) {
    if (inter->fL == NULL) return;
    f = inter->fL;
    locfa = inter->locfaL;
    npgf = inter->npgL;
    w_handle = inter->wL_handle;
    vol_index_handle = inter->vol_indexL_handle;
  } else  {
    if (inter->fR == NULL) return;
    f = inter->fR;
    locfa = inter->locfaR;
    npgf = inter->npgR;
    w_handle = inter->wR_handle;
    vol_index_handle = inter->vol_indexR_handle;
  }

  void* args;
  size_t arg_size;

  starpu_codelet_pack_args(&args, &arg_size,
                           STARPU_VALUE, &f->model.m, sizeof(int),
                           STARPU_VALUE, f->deg, 3 * sizeof(int),
                           STARPU_VALUE, f->raf, 3 * sizeof(int),
                           STARPU_VALUE, &npgf, sizeof(int),
                           STARPU_VALUE, &f->varindex, sizeof(varindexptr),
                           NULL);

  struct starpu_task *task;

  task = starpu_task_create();
  task->cl = ExtractInterface_codelet();
  task->cl_arg = args;
  task->cl_arg_size = arg_size;
  int nhandle = 0;
  task->handles[nhandle++] = w_handle;
  if (f->solver != NULL){
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->sol_handle;
  } else {
    task->handles[nhandle++] = f->wn_handle;
  }
  task->handles[nhandle++] = vol_index_handle;


  int sendto;
  sendto = side*(inter->ieR) + (1-side)*(inter->ieL);
  
  task->execute_on_a_specific_worker = 1;
  task->workerid = starpu_worker_get_by_type(STARPU_OPENCL_WORKER, gpuids[sendto%nbgpu]);  

  
  STARPU_CHECK_RETURN_VALUE(
      starpu_task_submit(task),
      "starpu_task_submit");
}



void ExtractInterface_C(void* buffers[], void* cl_args);
void ExtractInterface_OCL(void* buffers[], void* cl_args);

struct starpu_codelet* ExtractInterface_codelet() {
  static bool is_init = false;
  static struct starpu_codelet codelet;
  static struct starpu_perfmodel ei_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "extractinterface"
  }; 


  if (!is_init) {
    is_init = true;
    starpu_codelet_init(&codelet);

    codelet.name = "ExtractInterface";
    printf("init codelet %s...\n", codelet.name);

    codelet.nbuffers = 0;
    codelet.modes[codelet.nbuffers++] = STARPU_W;  // wface
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // wvol
    codelet.modes[codelet.nbuffers++] = STARPU_R;  // vol_index
    codelet.model = &ei_model;

    if (starpu_c_use) {
      size_t ncodelet = 0;
      codelet.cpu_funcs[ncodelet++] = ExtractInterface_C;
      codelet.cpu_funcs[ncodelet++] = NULL;
    }

    if (starpu_ocl_use) {
      LOAD_OPENCL_PROGRAM_SPU();

      size_t ncodelet = 0;
      codelet.opencl_funcs[ncodelet++] = ExtractInterface_OCL;
      codelet.opencl_funcs[ncodelet++] = NULL;
    }
  }

  return &codelet;
}

void ExtractInterface_C(void* buffers[], void* cl_args){
  int m;
  int deg[3];
  int raf[3];
  int npgf;
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, &npgf,
                             &Varindex);
  free(cl_args);

  int nbuffer = 0;
  schnaps_real* wf = (schnaps_real *)STARPU_VECTOR_GET_PTR(
      (struct starpu_vector_interface *)buffers[nbuffer++]);
  schnaps_real* wn = (schnaps_real *)STARPU_VECTOR_GET_PTR(
      (struct starpu_vector_interface *)buffers[nbuffer++]);
  int* index = (int *)STARPU_VECTOR_GET_PTR(
      (struct starpu_vector_interface *)buffers[nbuffer++]);


  for (int ipgf = 0; ipgf < npgf; ++ipgf) {
    int ipgv = index[ipgf];
    for (int iv = 0; iv < m; ++iv) {
      int imem = Varindex(deg, raf, m, ipgv, iv);
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      wf[imemf] = wn[imem];
    }
  }
}


void ExtractInterface_OCL(void* buffers[], void* cl_args){
  int m;
  int deg[3];
  int raf[3];
  int npgf;
  varindexptr Varindex;

  starpu_codelet_unpack_args(cl_args,
			     &m, deg, raf, &npgf,
                             &Varindex);
  free(cl_args);

  int nbuffer = 0;
  cl_mem wf = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem wn = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);
  cl_mem index = (cl_mem)STARPU_VECTOR_GET_DEV_HANDLE(buffers[nbuffer++]);

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
                                            "ExtractInterface",
                                            devid);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int param[7] = {m, deg[0], deg[1], deg[2], raf[0], raf[1], raf[2]};
  cl_mem param_cl = clCreateBuffer(context,
                                   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                   7 * sizeof(int),
                                   param,
                                   &status);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  int narg = 0;
  status = clSetKernelArg(kernel, narg++, sizeof(cl_mem), &param_cl);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &index);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wn);
  status |= clSetKernelArg(kernel, narg++, sizeof(cl_mem), &wf);
  if (status != CL_SUCCESS) STARPU_OPENCL_REPORT_ERROR(status);

  // Number of glops on the macrocell interface
  size_t wsize = npgf;

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


void ExtractInterface(Interface* inter, int side){

  field *fd;
  int locfa,npgf;
  schnaps_real *wf;
  int mod_ie;
  int* vol_index;

  if (side == 0) {
    fd = inter->fL;
    locfa = inter->locfaL;
    mod_ie = inter->ieL;
    npgf = inter->npgL;
    wf = inter->wL;
    vol_index = inter->vol_indexL;
  }

  if (side == 1) {
    fd = inter->fR;
    locfa = inter->locfaR;
    mod_ie = inter->ieR;
    npgf = inter->npgR;
    wf = inter->wR;
    vol_index = inter->vol_indexR;
  }

  assert(side ==1 || side == 0);

  if (fd !=NULL){
    for(int ipgf = 0; ipgf < npgf; ipgf++){
      int ipgv = vol_index[ipgf];
      for(int iv=0; iv < fd->model.m; iv++){
	int imem = fd->varindex(fd->deg, fd->raf, fd->model.m,
				ipgv, iv);
	int imemf = VarindexFace(npgf, fd->model.m, ipgf, iv);
	//printf("extract to=%d \n",mod_ie);
	wf[imemf] = fd->wn[imem];
      }
    }
  }


}

void InterfaceExplicitFlux_C(void* buffers[], void* cl_args);

void InterfaceExplicitFlux_SPU(Interface* inter, int side)
{


  field *f;
  field *fext;
  int locfa;
  starpu_data_handle_t index_ext;
  starpu_data_handle_t index;
  starpu_data_handle_t wn_ext;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet InterfaceExplicitFlux...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = InterfaceExplicitFlux_C;
    codelet.nbuffers = 7;
    codelet.modes[0] = STARPU_R;  // index_ext
    codelet.modes[1] = STARPU_R; // index
    codelet.modes[2] = STARPU_R; // vnds
    codelet.modes[3] = STARPU_R; // xpg
    codelet.modes[4] = STARPU_R; // wpg
    codelet.modes[5] = STARPU_R; // wn_ext
    codelet.modes[6] = STARPU_RW; // rhs
    codelet.name="InterfaceExplicitFlux";
  }

  const int sign = 1 - 2 * side;

  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    index_ext = inter->vol_indexR_handle;
    wn_ext = inter->wR_handle;
    index = inter->vol_indexL_handle;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    index_ext = inter->vol_indexL_handle;
    wn_ext = inter->wL_handle;
    index = inter->vol_indexR_handle;
  }
  else {
    assert(1==2);
  }


  if (f != NULL && fext != NULL){

    const unsigned int m = f->model.m;

    void* arg_buffer;
    size_t arg_buffer_size;

    starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			     STARPU_VALUE, f, sizeof(field),
			     STARPU_VALUE, fext, sizeof(field),
			     STARPU_VALUE, &locfa, sizeof(int),
			     STARPU_VALUE, &sign, sizeof(int),
			     0);
    //printf("sizeof_field=%d\n", (int) sizeof(field));

    task = starpu_task_create();
    task->cl = &codelet;
    task->cl_arg = arg_buffer;
    task->cl_arg_size = arg_buffer_size;
    int nhandle = 0;
    task->handles[nhandle++] = index_ext;
    task->handles[nhandle++] = index;
    task->handles[nhandle++] = inter->vnds_handle;
    task->handles[nhandle++] = inter->xpg_handle;
    task->handles[nhandle++] = inter->wpg_handle;
    task->handles[nhandle++] = wn_ext;
    if (f->solver != NULL){
      Skyline_SPU* sky_spu = f->solver->matrix;
      task->handles[nhandle++] = sky_spu->rhs_handle;
    }
    else {
      task->handles[nhandle++] = f->res_handle;
    }
    int ret = starpu_task_submit(task);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");


  }
}


void InterfaceExplicitFlux_C(void* buffers[], void* cl_args){


  field f0;
  field fext0;
  field *f = &f0;
  field *fext = &fext0;
  int locfa;
  int sign;

  starpu_codelet_unpack_args(cl_args,f,fext,&locfa,&sign);
  free(cl_args);


  int buf_num=0;

  struct starpu_vector_interface *index_ext_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  int* index_ext = (int *)STARPU_VECTOR_GET_PTR(index_ext_v);

  struct starpu_vector_interface *index_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  int* index = (int *)STARPU_VECTOR_GET_PTR(index_v);

  struct starpu_vector_interface *vnds_buf_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* vnds_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* xpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);

  struct starpu_vector_interface *wpg_buf_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* wpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(wpg_buf_v);

  struct starpu_vector_interface *wn_ext_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* wn_ext = (schnaps_real *)STARPU_VECTOR_GET_PTR(wn_ext_v);

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* rhs = (schnaps_real *)STARPU_VECTOR_GET_PTR(rhs_v);


  int npgf = NPGF(f->deg, f->raf, locfa);
  int m = fext->model.m;

  for(int ipgf = 0; ipgf < npgf; ipgf++) {


    schnaps_real flux[m];
    schnaps_real wL[m];
    for(int iv = 0; iv < m; iv++) {
      wL[iv] = 0;
    }

    schnaps_real wR[m];
    int ipgR = index_ext[ipgf];
    int ipgL = index[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int imemf = VarindexFace(npgf, m, ipgf, iv);
      //int imemR = fext->varindex(fext->deg, fext->raf,m , ipgR, iv);
      wR[iv] = wn_ext[imemf];
      //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    }

    // int_dL F(wL, wR, grad phi_ib)
    schnaps_real vndsloc[3];

    vndsloc[0] = sign * vnds_buf[3 * ipgf + 0];
    vndsloc[1] = sign * vnds_buf[3 * ipgf + 1];
    vndsloc[2] = sign * vnds_buf[3 * ipgf + 2];

    //printf("sign=%d ipgL=%d ipgR=%d vndsloc=%f %f\n",sign,ipgL,ipgR,vndsloc[0],vndsloc[1]);

    f->model.NumFlux(wL, wR, vndsloc, flux);
    //printf("flux=%f %f\n",flux[0],flux[0]);

    // Add flux  to the selected side
    const schnaps_real wpg = wpg_buf[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int ipgL = index[ipgf];
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv] * f->dt * wpg;
      //printf("imem=%d res=%f\n",imemL,rhs[imemL]);
      /* real* xpg = inter->xpg + 3 * ipgf; */
      /* real* vnds = inter->vnds + 3 * ipgf; */
      /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
      /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
    }

  }





}


void InterfaceBoundaryFlux_C(void* buffers[], void* cl_args);

void InterfaceBoundaryFlux_SPU(Interface* inter)
{


  field *f;
  int locfa;
  starpu_data_handle_t index;

  static bool is_init = false;
  static struct starpu_codelet codelet;
  struct starpu_task *task;

  if (!is_init){
    printf("init codelet InterfaceBoundaryFlux...\n");
    is_init = true;
    starpu_codelet_init(&codelet);
    codelet.cpu_funcs[0] = InterfaceBoundaryFlux_C;
    codelet.nbuffers = 5;
    codelet.modes[0] = STARPU_R;  // wface
    codelet.modes[1] = STARPU_R; // wvol
    codelet.modes[2] = STARPU_R; // vol_index
    codelet.modes[3] = STARPU_R; // wpg
    codelet.modes[4] = STARPU_RW; // vol_index
    codelet.name="InterfaceBoundaryFlux";
  }


  f = inter->fL;
  locfa = inter->locfaL;
  index = inter->vol_indexL_handle;



  void* arg_buffer;
  size_t arg_buffer_size;

  starpu_codelet_pack_args(&arg_buffer, &arg_buffer_size,
			   STARPU_VALUE, f, sizeof(field),
			   STARPU_VALUE, &locfa, sizeof(int),
			   0);

  task = starpu_task_create();
  task->cl = &codelet;
  task->cl_arg = arg_buffer;
  task->cl_arg_size = arg_buffer_size;
  int nhandle = 0;
  task->handles[nhandle++] = index;
  task->handles[nhandle++] = inter->vnds_handle;
  task->handles[nhandle++] = inter->xpg_handle;
  task->handles[nhandle++] = inter->wpg_handle;
  if (f->solver != NULL){
    Skyline_SPU* sky_spu = f->solver->matrix;
    task->handles[nhandle++] = sky_spu->rhs_handle;
  }
  else {
    task->handles[nhandle++] = f->res_handle;
  }
  int ret = starpu_task_submit(task);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");


}

void InterfaceBoundaryFlux_C(void* buffers[], void* cl_args){


  field f0;
  field *f = &f0;
  int locfa;

  starpu_codelet_unpack_args(cl_args,f,&locfa);
  free(cl_args);

  int m = f->model.m;

  int buf_num=0;

  struct starpu_vector_interface *index_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  int* index = (int *)STARPU_VECTOR_GET_PTR(index_v);

  struct starpu_vector_interface *vnds_buf_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* vnds_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(vnds_buf_v);

  struct starpu_vector_interface *xpg_buf_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* xpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(xpg_buf_v);

  struct starpu_vector_interface *wpg_buf_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* wpg_buf = (schnaps_real *)STARPU_VECTOR_GET_PTR(wpg_buf_v);

  struct starpu_vector_interface *rhs_v =
    (struct starpu_vector_interface *) buffers[buf_num++];
  schnaps_real* rhs = (schnaps_real *)STARPU_VECTOR_GET_PTR(rhs_v);

  for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {


    schnaps_real flux[m];
    schnaps_real wL[m];
    for(int iv = 0; iv < m; iv++) {
      wL[iv] = 0;
    }

    schnaps_real* xpg = xpg_buf + 3 * ipgf;
    schnaps_real* vnds = vnds_buf + 3 * ipgf;

    //printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
    f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
    int ipgL = index[ipgf];
    /* printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d \n", */
    /*        xpg[0], xpg[1], xpg[2], */
    /*        vnds[0], vnds[1],vnds[2], ipgL); */
    const schnaps_real wpg = wpg_buf[ipgf];
    for(int iv = 0; iv < m; iv++) {
      int ipgL = index[ipgf];
      // The basis functions is also the gauss point index
      int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
      rhs[imemL] -= flux[iv] * f->dt * wpg;
    }


  }








}


void InterfaceExplicitFlux(Interface* inter, int side){


  field *f;
  field *fext;
  int locfa;
  int mod_ie;
  int *index_ext;
  int *index;

  const int sign = 1 - 2 * side;

  if (side == 0) {
    f = inter->fL;
    fext = inter->fR;
    locfa = inter->locfaL;
    mod_ie = inter->ieL;
    index_ext = inter->vol_indexR;
    index = inter->vol_indexL;
  }
  else if (side == 1) {
    f = inter->fR;
    fext = inter->fL;
    locfa = inter->locfaR;
    mod_ie = inter->ieR;
    index_ext = inter->vol_indexL;
    index = inter->vol_indexR;
  }
  else {
    assert(1==2);
  }


  if (f != NULL){
    schnaps_real* res;
    if (f->solver != NULL){
      res = f->solver->rhs;
    } else {
      res = f->res;
    }

    const unsigned int m = f->model.m;

    //printf("locfa=%d \n",locfa);

    for(int ipgf = 0; ipgf < NPGF(f->deg, f->raf, locfa); ipgf++) {


      schnaps_real flux[m];
      schnaps_real wL[m];
      for(int iv = 0; iv < m; iv++) {
	wL[iv] = 0;
      }

      const schnaps_real wpg = inter->wpg[ipgf];

      if (fext != NULL) {  // the right element exists
	schnaps_real wR[m];
	int ipgR = index_ext[ipgf];
	int ipgL = index[ipgf];
	for(int iv = 0; iv < m; iv++) {
	  int imemR = fext->varindex(fext->deg, fext->raf,fext->model.m, ipgR, iv);
	  wR[iv] = fext->wn[imemR];
	  //wR[iv] = 1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}

	// int_dL F(wL, wR, grad phi_ib)
	schnaps_real vndsloc[3];

	vndsloc[0] = sign * inter->vnds[3 * ipgf + 0];
	vndsloc[1] = sign * inter->vnds[3 * ipgf + 1];
	vndsloc[2] = sign * inter->vnds[3 * ipgf + 2];

	//printf("sign=%d ipgL=%d ipgR=%d vndsloc=%f %f\n",sign,ipgL,ipgR,vndsloc[0],vndsloc[1]);

	//assert(f != NULL);

	f->model.NumFlux(wL, wR, vndsloc, flux);
	//printf("flux=%f %f\n",flux[0],flux[0]);

	// Add flux  to the selected side
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  res[imemL] -= flux[iv] * f->dt * wpg;
	  printf("imem=%d res=%f\n",imemL,res[imemL]);

	  /* real* xpg = inter->xpg + 3 * ipgf; */
	  /* real* vnds = inter->vnds + 3 * ipgf; */
	  /* printf("interface flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      } else { // The point is on the boundary.


	assert(sign == 1);
	schnaps_real* xpg = inter->xpg + 3 * ipgf;
	schnaps_real* vnds = inter->vnds + 3 * ipgf;

	//printf("tnow=%f wL=%f\n",f->tnow,wL[0]);
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);
	int ipgL = index[ipgf];
	printf("xpg=%f %f %f vnds=%f %f %f ipgL=%d &rhs=%p\n",
	       xpg[0], xpg[1], xpg[2],
	       vnds[0], vnds[1],vnds[2], ipgL,res);
	for(int iv = 0; iv < m; iv++) {
	  int ipgL = index[ipgf];
	  // The basis functions is also the gauss point index
	  int imemL = f->varindex(f->deg, f->raf,f->model.m, ipgL, iv);
	  res[imemL] -= flux[iv] * sign * f->dt * wpg;
	  /* printf("boundary flux=%f xpg=%f %f vnds=%f %f\n", */
	  /* 	 flux[0],xpg[0],xpg[1],vnds[0],vnds[1]); */
	}
      }

    }


  }
}

#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

bool submit_task() {
  bool test = true;

  // Data buffer
  const int size = 10000;
  schnaps_real* buffer = calloc(size, sizeof(schnaps_real));

  // Create data handle (init and register)
  for (int i = 0; i < size; ++i) buffer[i] = i;
  starpu_data_handle_t handle;
  starpu_vector_data_register(&handle, 0, (uintptr_t) buffer,
			      size, sizeof(schnaps_real));

  // Task
  ZeroBuffer_SPU(handle);

  // Check output
  starpu_task_wait_for_all();
  starpu_data_unregister(handle);
  for (int i = 0; i < size; ++i) test &= (fabs(buffer[i]) < _VERY_SMALL);

  if (test) printf(" OK\n");
  else printf(" KO !\n");

  return test;
}


int TestCodelet_ZeroBuffer_SPU() {
  bool test = true;


  // Kernel compilation options
  sprintf(cl_buildoptions, "%s", "");
  char buf[1000];
#ifdef _DOUBLE_PRECISION
  sprintf(buf, "-D schnaps_real=double -D _M=1");
#else
  sprintf(buf, "-D schnaps_real=float -D _M=1 -cl-single-precision-constant");
#endif
  strcat(cl_buildoptions, buf);


  // Init StarPU
  struct starpu_conf conf;
  test &= (starpu_conf_init(&conf) == 0);
  assert(test);
  test &= (starpu_init(&conf) != -ENODEV);
  assert(test);
  init_global_arbiter();

  // Interesting environment variables
  printf("Environment variables...\n");
  printf("STARPU_NCPU                : %s\n", getenv("STARPU_NCPU"));
  printf("STARPU_NOPENCL             : %s\n", getenv("STARPU_NOPENCL"));
  printf("STARPU_NCUDA               : %s\n", getenv("STARPU_NCUDA"));
  printf("STARPU_NMIC                : %s\n", getenv("STARPU_NMIC"));
  printf("STARPU_OPENCL_ON_CPUS      : %s\n", getenv("STARPU_OPENCL_ON_CPUS"));
  printf("STARPU_OPENCL_ONLY_ON_CPUS : %s\n", getenv("STARPU_OPENCL_ONLY_ON_CPUS"));

  // Workers
  const int nb_workers = starpu_worker_get_count();
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("Available workers...\n");
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  printf("OPENCL workers             : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());

  if (nb_workers == 0) {
    printf("No available worker.\n");
    return test;
  }

  // Regenerate opencl code
  if (nb_ocl > 0) GetOpenCLCode();

  // Codelet
  starpu_c_use = true;
  starpu_ocl_use = true;
  struct starpu_codelet* codelet = AddBuffer_codelet();

  // Empty codelet for function selection
  struct starpu_codelet codelet_backup = *codelet;
  for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
    codelet->cpu_funcs[i] = NULL;
    codelet->opencl_funcs[i] = NULL;
    codelet->cuda_funcs[i] = NULL;
    codelet->mic_funcs[i] = NULL;
  }

  // TRICK Need at least one C function
  codelet->cpu_funcs[0] = codelet_backup.cpu_funcs[0];

  // Loop over workers to submit tasks
  for (int wid = 0; wid < nb_workers; ++wid) {
    // Create context with a single worker
    unsigned int ctxid = starpu_sched_ctx_create(
        &wid, 1, "ctx",
        // TRICK Default scheduler:
        // avoids seg fault in _starpu_push_task_on_specific_worker
        STARPU_SCHED_CTX_POLICY_NAME, conf.sched_policy_name,
        NULL);
    starpu_sched_ctx_set_context(&ctxid);

    char name[100];
    starpu_worker_get_name(wid, name, 100);
    printf("StarPU worker name: %s\n", name);

    // Loop over the codelet implementations
    switch (starpu_worker_get_type(wid)) {
      case STARPU_CPU_WORKER:
        if (codelet_backup.cpu_funcs[0] == NULL) {
          printf("No C codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.cpu_funcs[i] != NULL) {
              codelet->cpu_funcs[0] = codelet_backup.cpu_funcs[i];
              printf("Submit C codelet %d...", i);
              test &= submit_task();
            }
          }
        }
        break;

      case STARPU_OPENCL_WORKER:
        if (codelet_backup.opencl_funcs[0] == NULL) {
          printf("No OpenCL codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.opencl_funcs[i] != NULL) {
              codelet->opencl_funcs[0] = codelet_backup.opencl_funcs[i];
              printf("Submit OpenCL codelet %d...", i);
              test &= submit_task();
            }
          }
        }
        break;

      case STARPU_CUDA_WORKER:
        if (codelet_backup.cuda_funcs[0] == NULL) {
          printf("No CUDA codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.cuda_funcs[i] != NULL) {
              codelet->cuda_funcs[0] = codelet_backup.cuda_funcs[i];
              printf("Submit CUDA codelet %d...", i);
              test &= submit_task();
            }
          }
        }
        break;

      case STARPU_MIC_WORKER:
        if (codelet_backup.mic_funcs[0] == NULL) {
          printf("No MIC codelet implementation.\n");
        } else {
          for (int i = 0; i < STARPU_MAXIMPLEMENTATIONS; ++i) {
            if (codelet_backup.mic_funcs[i] != NULL) {
              codelet->mic_funcs[0] = codelet_backup.mic_funcs[i];
              printf("Submit MIC codelet %d...", i);
              test &= submit_task();
            }
          }
        }
        break;

      default:
        printf("Untreated worker type.\n");
    }

    // Delete context
    starpu_sched_ctx_delete(ctxid);
  }

  // Delete opencl program if it has been created
  if (nb_ocl > 0) {
    int ret = starpu_opencl_unload_opencl(&opencl_program);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_opencl_unload_opencl");
  }

  destroy_global_arbiter();
  starpu_shutdown();


  return test;
}

int main(void) {
  // Unit tests
  int resu = TestCodelet_ZeroBuffer_SPU();
  if (resu) printf("StarPU ZeroBuffer Codelet test OK !\n");
  else printf("StarPU ZeroBuffer Codelet test failed !\n");
  return !resu;
}

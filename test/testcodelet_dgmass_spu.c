#include "test.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

bool submit_task(Simulation* simu, schnaps_real* buffer) {
  bool test = true;

  const int fsize =  simu->wsize / simu->macromesh.nbelems;
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field* f = simu->fd + ie;

    // Create data handle (init and register)
    for (int i = 0; i < fsize; ++i) {
      f->res[i] = sin(i + 1.);
      f->dtwn[i] = 2;
    }
    starpu_vector_data_register(&f->res_handle, 0, (uintptr_t)f->res,
                                fsize, sizeof(schnaps_real));
    starpu_vector_data_register(&f->dtwn_handle, 0, (uintptr_t)f->dtwn,
                                fsize, sizeof(schnaps_real));

    // Submit task
    DGMass_SPU(f);
  }

  starpu_task_wait_for_all();

  // Check output
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field* f = simu->fd + ie;

    starpu_data_unregister(f->res_handle);
    starpu_data_unregister(f->dtwn_handle);

    for (int i = 0; i < fsize; ++i)
      assert(fabs(f->res[i] - sin(i + 1.)) < _SMALL);
    for (int i = 0; i < fsize; ++i) {
      if (!(fabs(buffer[ie * fsize + i] - f->dtwn[i] - 1) < _SMALL))
        printf("field: %d  reference[%d]: %f  result[%d]: %f err=%e\n",
               ie, i, buffer[ie * fsize + i], i,
	       f->dtwn[i],fabs(buffer[ie * fsize + i] - f->dtwn[i] -1));
      test &= (fabs(buffer[ie * fsize + i] - f->dtwn[i] -1) < 20 * _SMALL);
    }
  }

  if (test) printf(" OK\n");
  else printf(" KO !\n");

  return test;
}


int TestCodelet_DGMass_SPU() {
  bool test = true;

  int deg[]={2, 2, 2};
  int raf[]={1, 1, 1};

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testdisque.msh");
  //ReadMacroMesh(&mesh,"../test/testcubetordu.msh");
  BuildConnectivity(&mesh);
  CheckMacroMesh(&mesh, deg, raf);

  Model model;
  model.m = 6;
  model.NumFlux = Maxwell3DNumFlux_upwind;
  model.BoundaryFlux = Maxwell3DBoundaryFlux_upwind;
  model.InitData = Maxwell3DInitData;
  model.ImposedData = Maxwell3DImposedData;
  model.Source = NULL;

  Simulation simu;
  EmptySimulation(&simu);


  // Kernel compilation options
  char buf[1000];
#ifdef _DOUBLE_PRECISION
  sprintf(buf, " -D schnaps_real=double -D _M=%d",model.m);
#else
  sprintf(buf, " -cl-single-precision-constant -D schnaps_real=float -D _M=%d",model.m);
#endif
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D NUMFLUX=%s", "Maxwell3DNumFlux_upwind");
  strcat(cl_buildoptions, buf);
  InitSimulation(&simu, &mesh, deg, raf, &model);

  printf("StarPU compilation options: %s\n", cl_buildoptions);


  // Init StarPU
  struct starpu_conf conf;
  test &= (starpu_conf_init(&conf) == 0);
  assert(test);
  test &= (starpu_init(&conf) != -ENODEV);
  assert(test);

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


  // Compute the comparision buffer
  schnaps_real* buffer = calloc(simu.wsize, sizeof(schnaps_real));
  const int fsize =  simu.wsize / simu.macromesh.nbelems;
  for (int ie = 0; ie < simu.macromesh.nbelems; ++ie) {
    field* f = simu.fd + ie;

    // Init data
    for (int i = 0; i < fsize; ++i)
      f->dtwn[i] = sin(i + 1.);

    // Compute mass term
    DGMass(f, NULL, f->dtwn);

    // Store data for comparision
    for (int i = 0; i < fsize; ++i)
      buffer[ie * fsize + i] = f->dtwn[i];
  }


  // Codelet
  starpu_c_use = true;
  starpu_ocl_use = true;
  struct starpu_codelet* codelet = DGMass_codelet();

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
              test &= submit_task(&simu, buffer);
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
              test &= submit_task(&simu, buffer);
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
              test &= submit_task(&simu, buffer);
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
              test &= submit_task(&simu, buffer);
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

  starpu_shutdown();

  free(buffer);

  FreeMacroMesh(&mesh);


  return test;
}

int main(void) {
  // Unit tests
  int resu = TestCodelet_DGMass_SPU();
  if (resu) printf("StarPU DGMass Codelet test OK !\n");
  else printf("StarPU DGMass Codelet test failed !\n");
  return !resu;
}

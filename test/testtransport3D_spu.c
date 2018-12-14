#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "acoustic.h"
#include "getopt.h"
#include "global.h"
#include "io.h"
#include "schnaps.h"
#include "../test/test.h"

/* DEBUG xdmf filename and format*/
const int  xdmf_filename_size = sizeof("test_transport3D_000000.xdmf");

const char xdmf_filename_fmt[sizeof("stest_transport3D__%06d.xdmf")] 
                            = "test_transport3D_%06d.xdmf"; 


int test_sn_transport(){

#ifdef _WITH_OPENCL
    if (!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
        printf("OpenCL device not acceptable.\n");
        return true;
    }
#endif 

    putenv("STARPU_PROFILING=1");
    putenv("STARPU_SCHED=dmda");
    
    
    // putenv("STARPU_NOPENCL=10");
    // putenv("STARPU_NCPU=0");
    // putenv("STARPU_OPENCL_ON_CPUS=1");
    // putenv("STARPU_OPENCL_ONLY_ON_CPUS=1"); 
    // putenv("STARPU_GENERATE_TRACE=1"); 
    
        
    // putenv("STARPU_WORKER_STATS=1");

    /* schnaps STARPU boolean*/
    starpu_use = true;
    starpu_c_use = true;
    starpu_c_use = true;
    starpu_ocl_use = true;
    starpu_cuda_use = false;

    /*Read and init mesh struct*/
    MacroMesh mesh;
    ReadMacroMesh(&mesh, "../test/testcube.msh");
    BuildConnectivity(&mesh);

    /*DG parameters*/
    int thedeg  = 2;
    int therafx = 8;
    int therafy = 8;
    int therafz = 8;
    int deg[] = {thedeg, thedeg, thedeg};
    int raf[] = {therafx, therafy, therafz};
    CheckMacroMesh(&mesh, deg, raf);

    Model model;
    Simulation simu;
    EmptySimulation(&simu);
   
    strcpy(model.name,"testtransport3D_spu");
    
    model.m = 1;
    model.InitData = ocl_sin3D_init_data;
    model.ImposedData = ocl_sin3D_imposed_data;
    model.NumFlux = ocl_sin3d_upwind_numflux;
    model.BoundaryFlux = ocl_sin3d_boundary_numflux;
    model.Source = NULL;
    
    simu.cfl = 0.95;
    schnaps_real tmax = 0.2;
    
    /* OpenCL compilation options */
    char buf[1000];
#ifndef _DOUBLE_PRECISION
    sprintf(buf, " -D schnaps_real=float -D _M=%d -cl-single-precision-constant",
    model.m);
#else
    sprintf(buf, " -D schnaps_real=double -D _M=%d", model.m);
#endif
    strcat(cl_buildoptions, buf);
    sprintf(numflux_cl_name, "%s", "ocl_sin3d_upwind_numflux");
    sprintf(buf," -D NUMFLUX=");
    strcat(buf, numflux_cl_name);
    strcat(cl_buildoptions, buf);

    sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math",
    "ocl_sin3d_boundary_numflux");
    strcat(cl_buildoptions, buf);

    InitSimulation(&simu, &mesh, deg, raf, &model);
    
    /* Solve */
    RK2_SPU(&simu, tmax);
    
    /* Diagnostic */
    schnaps_real error = 0;
    error = L2error_onefield(&simu,0);
    
    // printf("Error L2 = %f\n",error);

    // char xdmf_data_filename[xdmf_filename_size];
    // sprintf(xdmf_data_filename,xdmf_filename_fmt,simu.iter_time_rk);
    // schnaps_simulation_xdmf_plot_fields_xml_structured_glops(&simu, xdmf_data_filename);
    destroy_global_arbiter();
    starpu_shutdown(); 
    

    if(error < 0.03) return true;
    else return false;
}

int main(void) {
  bool resu = test_sn_transport();
  return !resu;
}




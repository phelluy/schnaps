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
const int  xdmf_filename_size = sizeof("sn_spec_gaussian_cube_000000.xdmf");

const char xdmf_filename_fmt[sizeof("sn_spec_gaussian_cube_%06d.xdmf")] 
                            = "sn_spec_gaussian_cube_%06d.xdmf"; 
/**/


int test_sn_transport(){
    int test = true;

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
    
        
    putenv("STARPU_WORKER_STATS=1");

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


    /*Init velocity space*/
    schnaps_acoustic_data = lebedev_6;

    /* Model parameters */
    Model model;
    Simulation simu;
    EmptySimulation(&simu);

    simu.cfl = 0.95;

    
    strcpy(model.name,"test_transport_sn_lebedev");
    
    model.m = 6;
    model.InitData     = ocl_s6_init_gaussian_density;
    model.ImposedData  = ocl_s6_imposed_gaussian_density;
    model.NumFlux      = ocl_s6_upwind_numflux;
    model.BoundaryFlux = ocl_s6_free_boundary;
    
   
    model.Source       = NULL;

    schnaps_real tmax = 0.1;
    
    /* OpenCL compilation options */
    char buf[1000];
#ifndef _DOUBLE_PRECISION
    sprintf(buf, " -D schnaps_real=float -D _M=%d -cl-single-precision-constant",
    model.m);
#else
    sprintf(buf, " -D schnaps_real=double -D _M=%d", model.m);
#endif
    strcat(cl_buildoptions, buf);
    sprintf(numflux_cl_name, "%s", "ocl_s6_upwind_numflux");
    sprintf(buf," -D NUMFLUX=");
    strcat(buf, numflux_cl_name);
    strcat(cl_buildoptions, buf);

    sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math",
    "ocl_s6_free_boundary");
    strcat(cl_buildoptions, buf);

    InitSimulation(&simu, &mesh, deg, raf, &model);
    
    /* Debug */
    /*

    
    // schnaps_real error = 0;
    error = sn_l2_error_all_field(&simu);
    printf("Error L2 = %f\n",error);
    return test;
    */
    
    /* Solve */
    RK2_SPU(&simu, tmax);
    /* Diagnostic */
    schnaps_real error = 0;
    error = L2error(&simu);
    printf("Error L2 = %f\n",error);

    // CopyfieldtoCPU(&simu);
    // schnaps_real rec_rho=0;
    // schnaps_real rec_rho_ar[1];
    // sn_compute_macro(&simu, &rec_rho_ar[0]);
    // char xdmf_data_filename[xdmf_filename_size];
    // sprintf(xdmf_data_filename,xdmf_filename_fmt,simu.iter_time_rk);
    // schnaps_simulation_xdmf_plot_fields_xml_structured_glops(&simu, xdmf_data_filename);
    // destroy_global_arbiter();
    // starpu_shutdown(); 
    

    if(error < 1e-2) return true;
    else return false;
}

int main(void) {
  int resu = test_sn_transport();
  return !resu;
}




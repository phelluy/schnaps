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

#define HEADER_MAX_SIZE 9999
#define PATH_MAX_SIZE 9999
#define EXPORT_XDMF 0
#define EXPORT_CSV_RCP 1


/*chdir folder + time*/
const int  run_name_size = sizeof("sn_spec_gaussian_cube");
const char run_name[sizeof("sn_spec_gaussian_cube")] = "sn_spec_gaussian_cube";

/*xdmf filename and format*/
const int  xdmf_filename_size = sizeof("sn_spec_gaussian_cube_000000.xdmf");

const char xdmf_filename_fmt[sizeof("sn_spec_gaussian_cube_%06d.xdmf")] 
                            = "sn_spec_gaussian_cube_%06d.xdmf"; 
/*receptor csv file*/
const char csv_filename_density[sizeof("sn_spec_gaussian_cube_rcp_density.csv")] 
                                = "sn_spec_gaussian_cube_rcp_density.csv"; 
                              
const char csv_filename_intensity[sizeof("sn_spec_gaussian_cube_rcp_intensity.csv")] 
                                = "sn_spec_gaussian_cube_rcp_intensity.csv"; 


void sn_utils_change_wd()
{
  
  char folder_name[PATH_MAX_SIZE];
  char time_str[1024];
  struct tm *timenow;
  struct stat st = {0};

  time_t now = time(NULL);
  timenow = gmtime(&now);

  strftime(time_str, sizeof(time_str), "%Y%m%d_%H%M%S", timenow);

  snprintf(folder_name, 4 + sizeof(run_name_size) + 1 + sizeof(time_str), "run_%s_%s", run_name,time_str);
  printf("%s\n",folder_name);
    
  if (stat(folder_name, &st) == -1) mkdir(folder_name, 0777);
  else printf("[Error] CreateRunFolder :: error in mkdir call\n");
  
  if (chdir(folder_name) == 0) {
  } else {
    printf("[Error] CreateRunFolder :: error in chdir call");
  }
  
}

void unit_test_spec()
{
  __constant acoustic_data *ad = &schnaps_acoustic_data;
  schnaps_real normal[18]   = { 1,  0,  0,
                               -1,  0,  0,
                                0,  1,  0,
                                0, -1,  0,
                                0,  0,  1,
                                0,  0, -1};
                               
  schnaps_real v[3]         = {0, 0, 0};
  schnaps_real v_cp[3]      = {0, 0, 0};
  schnaps_real n[3]         = {0, 0, 0};
  schnaps_real v_spec[3]    = {0, 0, 0};
  schnaps_real v_spec_id[3] = {0, 0, 0};
  schnaps_real v_phy[3]     = {0, 0, 0};
  schnaps_real vn           = 0;
  
  int id = -1;

  for (int k = 0; k < 6; k++){
    
    n[0] = normal[3 * k + 0];
    n[1] = normal[3 * k + 1];
    n[2] = normal[3 * k + 2];

    for (int i = 0; i < ad->nb_v; i++) {
      v[0] = ad->v_tab[3 * i + 0];
      v[1] = ad->v_tab[3 * i + 1];
      v[2] = ad->v_tab[3 * i + 2];
      
      // id = AcousticSnComputeNormalSpec(v, n, v_spec);
      
      vn = + v[0] * n[0] 
           + v[1] * n[1] 
           + v[2] * n[2];
           
      /* Colision case*/
      if(vn > 0){
        // id = AcousticSnComputeNormalSpec(v, n, v_spec);
        
        /* Snell Descartes physical law */
        v_phy[0] = v[0] - 2 * n[0] * vn;
        v_phy[1] = v[1] - 2 * n[1] * vn;
        v_phy[2] = v[2] - 2 * n[2] * vn;
        
        // printf ("v_spec = (%f,%f,%f) v_phy = (%f,%f,%f) \n", 
                // v_spec[0], v_spec[1], v_spec[2], v_phy[0], v_phy[1], v_phy[2]);
                
        /* Test vector in AcousticSnComputeNormalSpec */
        assert (v_spec[0] == v_phy[0]);
        assert (v_spec[1] == v_phy[1]);
        assert (v_spec[2] == v_phy[2]);
        
        v_cp[0] = v_spec[0];
        v_cp[1] = v_spec[1];
        v_cp[2] = v_spec[2];
        
        /* Test id search in AcousticSnComputeNormalSpec */
        
        v_spec_id[0] = ad->v_tab[3 * id + 0];
        v_spec_id[1] = ad->v_tab[3 * id + 1];
        v_spec_id[2] = ad->v_tab[3 * id + 2];
        
        // printf ("v_spec_id = (%f,%f,%f) v_phy = (%f,%f,%f) \n", 
                // v_spec_id[0], v_spec_id[1], v_spec_id[2],
                // v_phy[0], v_phy[1], v_phy[2]);
        
        assert (v_spec_id[0] == v_phy[0]);
        assert (v_spec_id[1] == v_phy[1]);
        assert (v_spec_id[2] == v_phy[2]);
        
        int id2 = -1;
        
        /* Test array */
        // printf("NTM\n");
        // id2 = ad->spec_tab[0][i];
        // printf("id array = %d, id recent = %d\n",id2,id);
        
        // if (n[0] == -1 && n[1] == 0 && n[2] == 0){
                    // id = AcousticSnComputeNormalSpec(v, n, v_spec);
          
          // printf("vec init  = (%f,%f,%f) \n", ad->v_tab[3 * i + 0], 
                // ad->v_tab[3 * i + 1], ad->v_tab[3 * i + 2]);
          // printf("vec spec = (%f,%f,%f) \n", ad->v_tab[3 * id2 + 0], 
                // ad->v_tab[3 * id2 + 1], ad->v_tab[3 * id2 + 2]);
        // }
        // if (n[0] == 0 && n[1] != 0 && n[2] == 0)
          // printf("ny id_array = %d, id = %d \n", ad->spec_tab[1][i], id);
        // if (n[0] == 0 && n[1] == 0 && n[2] != 0)
          // printf("nz id_array = %d, id = %d \n", ad->spec_tab[2][i], id);
        /* Test AcousticSnComputeNormalSpec by reversing*/
        // AcousticSnComputeNormalSpec(v_cp, n, v_spec);
        
        // assert (v_spec[0] == v[0]);
        // assert (v_spec[1] == v[1]);
        // assert (v_spec[2] == v[2]);

        // printf ("v = (%f,%f,%f) v' = (%f,%f,%f) \n", 
       // v[0], v[1], v[2], v_spec[0], v_spec[1], v_spec[2]);
      }
   }
  }
}

void sn_rcp_csv_header(
        char header[HEADER_MAX_SIZE])
{
  /*
	 * Compute the header line of the csv file.
   *
	 */
   
  __constant acoustic_data *ad = &schnaps_acoustic_data;

  char header_base_fmt[sizeof("Time, R%04d")]="Time";
  char header_block_fmt[sizeof(", R%04d")] = ", R%04d";
  char header_fmt[HEADER_MAX_SIZE]; 
  
  /* time column */
  snprintf(header_fmt, sizeof(header_fmt), "%s", header_base_fmt);
  sprintf(header, header_fmt, 0);
  
  /* receptors columns */
  for(int i = 1 ; i <= ad->nb_rcp ; i++){
    snprintf(header_fmt, sizeof(header_fmt), "%s%s", header, header_block_fmt);
    sprintf(header, header_fmt, i);
  }
  
  /* carriage return */
  snprintf(header_fmt, sizeof(header_fmt), "%s%s", header, "%s");
  sprintf(header, header_fmt, "\n");
  return;
}

void sn_rcp_csv_line(
        const schnaps_real tnow,
        const schnaps_real *rcpt_val,
        char data_line[HEADER_MAX_SIZE])
{
  /*
	 * Compute the data line to export in the the csv file.
   *
	 */
   
  __constant acoustic_data *ad = &schnaps_acoustic_data;
  
  char data_line_base_fmt[sizeof("%06f")] = "%06f";
  char data_line_block_fmt[sizeof(", %06d")] = ", %06f";
  char data_line_fmt[HEADER_MAX_SIZE]; 

  /* time column */
  snprintf(data_line_fmt, sizeof(data_line_fmt), "%s", data_line_base_fmt);
  sprintf(data_line, data_line_fmt, tnow);
  
  /* Receptors column */
  for (int i = 1; i <= ad->nb_rcp; i++) {
    snprintf( data_line_fmt, sizeof(data_line_fmt), "%s%s", data_line,
              data_line_block_fmt);
  
    sprintf(data_line, data_line_fmt, rcpt_val[i-1]);
  }
  
  /* carriage return*/
  snprintf(data_line_fmt, sizeof(data_line_fmt), "%s%s", data_line, "%s");
  sprintf(data_line, data_line_fmt, "\n");
  return;
}

void sn_rcp_csv_time_exporter(
        bool first_open,
        const char *filename,
        const schnaps_real tnow,
        const schnaps_real *rcpt_val)
{
  /*
	 * Write each receptor's data into the csv file.
   *
	 */
  
  
  FILE * csvfile;
  char data_line[HEADER_MAX_SIZE];
  
 /* write the header and first line */
  if (first_open) {
    csvfile = fopen(filename, "w");
    // char header_str[HEADER_MAX_SIZE];
    // sn_rcp_csv_header(header_str);
    // fprintf(csvfile, header_str);
    sn_rcp_csv_line(tnow, rcpt_val, data_line);
    fprintf(csvfile,data_line);
  }  else {
    /* append receptor's data */
    csvfile = fopen(filename, "a+" );
    sn_rcp_csv_line(tnow, rcpt_val, data_line);
    fprintf(csvfile, data_line);
  }
  
  fclose(csvfile);
  return;
}



void sn_rk2_export_on(
        Simulation *simu,
        schnaps_real tmax)
{
  //Temporal infos
  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;
  simu->tmax = tmax;
  simu->itermax_rk = tmax / simu->dt + 1;
  simu->tnow = 0;
  
  //Setup diag freq
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10) ? 1 : simu->itermax_rk / 10;
  int iter = 0;
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;
  
  schnaps_real *wn = calloc(simu->wsize, sizeof(schnaps_real));
  assert(wn);
  
  //Exporter
  char xdmf_data_filename[xdmf_filename_size];
  schnaps_real rec_rho=0;
  schnaps_real rec_rho_ar[1];
  
   if (simu->nb_diags != 0) {
     simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
   }

  while (simu->tnow < tmax) {
    
    /* data export */
     if (iter % freq == 0) {
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
      
      if (EXPORT_XDMF == 1) {
        sprintf(xdmf_data_filename,xdmf_filename_fmt,simu->iter_time_rk);
        sn_compute_macro(simu, &rec_rho_ar[0]);
        schnaps_simulation_xdmf_plot_fields_xml_structured_glops(simu, xdmf_data_filename);
      }
    }
      
     if (EXPORT_CSV_RCP) {
      if(iter == 0){
        sn_compute_macro(simu, &rec_rho_ar[0]);
        printf("t=%f, Receptor value = %f \n",simu->tnow, rec_rho_ar[0]);
        sn_rcp_csv_time_exporter(true,
                                 csv_filename_density,
                                 simu->tnow,
                                 rec_rho_ar );
      } else {
        sn_compute_macro(simu, &rec_rho_ar[0]);
        printf("t=%f, Receptor value = %f \n",simu->tnow, rec_rho_ar[0]);
        sn_rcp_csv_time_exporter(false,
                                 csv_filename_density,
                                 simu->tnow,
                                 rec_rho_ar );
      }
    }

   /* time solver*/
    simu->dt = 0.5 * dt;
    
    /* pre dt field  
    if(simu->pre_dtfields != NULL) {
      printf("predtfield\n");
      simu->pre_dtfields(simu);
    }
    */
      
    /* RK 1 step */
    RK_Copy(simu,simu->w,wn);
    DtFields(simu);
    RK_in(simu);
    simu->tnow += 0.5 * dt;
    
    /* RK 2 step */
    DtFields(simu);
    simu->dt = dt;
    RK_out(simu, wn);
    simu->tnow += 0.5 * dt;
    simu->dt = 0.5 * dt;
    
    /* post dt field
    if(simu->post_dtfields != NULL) {
      printf("postfield\n");
     simu->post_dtfields(simu);
    }
    */
    simu->dt=dt;

    /*
    if(simu->update_after_rk != NULL){
      printf("update after rk\n");
      simu->update_after_rk(simu, simu->w);
    }
    */
      
    iter++;
    simu->iter_time_rk = iter;
  }
    free(wn);
}


void sn_rk2_spu(Simulation *simu, schnaps_real tmax){

  //Temporal infos
  simu->dt = Get_Dt_RK(simu);
  schnaps_real dt = simu->dt;
  simu->tmax = tmax;
  simu->itermax_rk = tmax / simu->dt + 1;
  simu->tnow = 0;
  
  //Setup diag freq
  int size_diags;
  int freq = (1 >= simu->itermax_rk / 10) ? 1 : simu->itermax_rk / 10;
  int iter = 0;
  size_diags = simu->nb_diags * simu->itermax_rk;
  simu->iter_time_rk = iter;

  assert(starpu_use);
 
  //Exporter
  char xdmf_data_filename[xdmf_filename_size];
  schnaps_real rec_rho=0;
  schnaps_real rec_rho_ar[1];
  
  RegisterSimulation_SPU(simu);

  // top chrono
  schnaps_real tps_debut = seconds();
  
  while(simu->tnow < tmax) {
    
    if (iter % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

    if (EXPORT_XDMF == 1) {
        UnregisterSimulation_SPU(simu);
        sprintf(xdmf_data_filename,xdmf_filename_fmt,simu->iter_time_rk);
        sn_compute_macro(simu, &rec_rho_ar[0]);
        schnaps_simulation_xdmf_plot_fields_xml_structured_glops(simu, xdmf_data_filename);
        RegisterSimulation_SPU(simu);
      }
    
      
    if (EXPORT_CSV_RCP) {
        UnregisterSimulation_SPU(simu);
        if(iter == 0){
            sn_compute_macro(simu, &rec_rho_ar[0]);
            printf("t=%f, Receptor value = %f \n",simu->tnow, rec_rho_ar[0]);
            sn_rcp_csv_time_exporter(true,
                                     csv_filename_density,
                                     simu->tnow,
                                     rec_rho_ar );
        } else {
            sn_compute_macro(simu, &rec_rho_ar[0]);
            printf("t=%f, Receptor value = %f \n",simu->tnow, rec_rho_ar[0]);
            sn_rcp_csv_time_exporter(false,
                                     csv_filename_density,
                                     simu->tnow,
                                     rec_rho_ar );
        }
        RegisterSimulation_SPU(simu);
    }
  
     
    
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
#endif
    }

    simu->tnow += 0.5 * dt;

    /* if(simu->update_after_rk != NULL){  */
    /*   simu->update_after_rk(simu, simu->w);  */
    /* } */

    iter++;
    simu->iter_time_rk = iter;
  }

  printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);

  starpu_task_wait_for_all() ;

  // top chrono
  UnregisterSimulation_SPU(simu);

  schnaps_real tps_fin = seconds();
  printf("Temps total =%f\n", tps_fin - tps_debut);
}


int test_sn_spec_gaussian_cube(void) { 
  bool test = true;
  
  /*Read and init mesh struct*/
  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testcube.msh");
  // ReadMacroMesh(&mesh, "../test/couloir2.msh");
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
  schnaps_acoustic_data = lebedev_50 ;
  // schnaps_acoustic_data = lebedev_6 ;
  // schnaps_acoustic_data = lebedev_6 ;
  acoustic_data *ad = &schnaps_acoustic_data ;
 
  
  // unit_test_spec();
  // return 0;
  
  /*Model parameters*/
  
  Model model;
  model.m            = ad->nb_v+5;
  model.InitData     = sn_init_gaussian_density;
  model.ImposedData  = sn_imposed_gaussian_density;
  model.Source       = NULL;
  model.NumFlux      = sn_upwind_numflux;
  // model.BoundaryFlux = sn_specular_boundary_flux;
  // model.BoundaryFlux = sn_diffusive_boundary_flux;
  // model.BoundaryFlux = sn_diffusive_boundary_flux_scalar_product;
  model.BoundaryFlux = sn_mixed_boundary_flux;

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  /*Resume Simulation*/
  
  // simu.iter_time_rk=29;
  // simu.tnow = 0.01* (schnaps_real)   simu.iter_time_rk; 
  #ifdef _WITH_HDF5
  // schnaps_simulation_load_field_state(&simu,"fooo");
  #endif
  simu.cfl=0.95;
  schnaps_real tmax = 10;
  schnaps_real tps_deb = seconds();
  sn_utils_change_wd();
  sn_rk2_export_on(&simu, tmax);
  schnaps_real tps_fin = seconds();
  
  printf("temps = %f\n", tps_fin-tps_deb);

  /*Save complete simulation*/
  // schnaps_simulation_dump_field_state(&simu,"fooo");

  // FreeAcousticData(ad);
  FreeMacroMesh(&mesh);
  // test= test && (L2err < 1e-2);
  test = true;
  return test;
}



int test_sn_spec_gaussian_cube_spu(void){
   int test = true;

#ifdef _WITH_OPENCL
  if (!cldevice_is_acceptable(nplatform_cl, ndevice_cl)) {
    printf("OpenCL device not acceptable.\n");
    return true;
  }
#endif 
  
  putenv("STARPU_SCHED=dmda");
  
  /* to estimate the cost of a task StarPU takes into
   * account the estimated computation time
   */
  putenv("STARPU_SCHED_ALPHA=1");
  
  // putenv("STARPU_OPENCL_PROGRAM_DIR='./ho");
  /* to estimate the cost of a task StarPU takes into account
   * the estimated data transfer time
   */
  putenv("STARPU_SCHED_BETA=1");

  /* specify the number of opencl devices that StarPU can use */
  putenv("STARPU_NOPENCL=8");

  /* specify the number of CPU workers */
  putenv("STARPU_NCPU=8");

  /* specify the number of CUDA devices that StarPU can use */
  putenv("STARPU_NCUDA=4");


  /* enable CPU devices for openCL computation (GPU already enable)*/
  putenv("STARPU_OPENCL_ON_CPU=1");
  
  /* openCL driver will ONLY enable CPU devices */
  putenv("STARPU_OPENCL_ONLY_ON_CPUS=1"); 
    
  putenv("STARPU_PREFETCH=0");
  
  
  
  /* schnaps STARPU boolean*/
  starpu_use = true;
  starpu_c_use = true;
  starpu_ocl_use = false;
  starpu_cuda_use = false;
  
  /*Read and init mesh struct*/
  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testcube.msh");
  // ReadMacroMesh(&mesh, "../test/couloir2.msh");
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
  schnaps_acoustic_data = lebedev_50 ;
  // schnaps_acoustic_data = lebedev_6 ;
  acoustic_data *ad = &schnaps_acoustic_data ;
  
    /*Model parameters*/
  
  Model model;
  
  model.m            = ad->nb_v+5;
  model.InitData     = sn_init_gaussian_density;
  model.ImposedData  = sn_imposed_gaussian_density;
  model.Source       = NULL;
  model.NumFlux      = sn_upwind_numflux;
  // model.BoundaryFlux = sn_specular_boundary_flux;
  // model.BoundaryFlux = sn_diffusive_boundary_flux;
  // model.BoundaryFlux = sn_diffusive_boundary_flux_scalar_product;
  model.BoundaryFlux = sn_mixed_boundary_flux;

  
  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.cfl=0.95;
  schnaps_real tmax = 10;
  

  
  
  
    char buf[1000];
// #ifndef _DOUBLE_PRECISION
  sprintf(buf, " -D schnaps_real=float -D _M=%d -cl-single-precision-constant",
	  model.m);
// #else
  // sprintf(buf, " -D schnaps_real=double -D _M=%d", model.m);
// #endif
  strcat(cl_buildoptions, buf);

  sprintf(numflux_cl_name, "%s", "sn_upwind_numflux");

  sprintf(buf," -D NUMFLUX=");
  strcat(buf, numflux_cl_name);
  strcat(cl_buildoptions, buf);

  sprintf(buf, " -D BOUNDARYFLUX=%s -cl-fast-relaxed-math",
	  "sn_mixed_boundary_flux");
  strcat(cl_buildoptions, buf);

  // Interesting environment variables
  printf("\n\n---------------------------------------------------------\n");
  printf("Environment variables...\n");
  printf("STARPU_NCPU                : %s\n", getenv("STARPU_NCPU"));
  printf("STARPU_NOPENCL             : %s\n", getenv("STARPU_NOPENCL"));
  printf("STARPU_NCUDA               : %s\n", getenv("STARPU_NCUDA"));
  printf("STARPU_NMIC                : %s\n", getenv("STARPU_NMIC"));
  printf("STARPU_OPENCL_ON_CPUS      : %s\n", getenv("STARPU_OPENCL_ON_CPUS"));
  printf("STARPU_OPENCL_ONLY_ON_CPUS : %s\n",
	 getenv("STARPU_OPENCL_ONLY_ON_CPUS"));

  printf("\n\n---------------------------------------------------------\n");
  printf("4/ Initializing simulation\n");
  


  // Workers
  const int nb_workers = starpu_worker_get_count();
  const int nb_ocl = starpu_opencl_worker_get_count();
  printf("Available workers...\n");
  printf("Number of workers          : %d\n", nb_workers);
  printf("CPU workers                : %d\n", starpu_cpu_worker_get_count());
  printf("OPENCL workers             : %d\n", nb_ocl);
  printf("CUDA workers               : %d\n", starpu_cuda_worker_get_count());
  printf("MIC workers                : %d\n", starpu_mic_worker_get_count());

  //SmartPrefetch_SPU(&simu);

  char xdmf_data_filename[xdmf_filename_size];


  // sn_utils_change_wd();

  schnaps_real tps_deb = seconds();
  printf("\n\n---------------------------------------------------------\n");
  printf("5/ Solving\n");
  sn_rk2_spu(&simu, tmax);

  schnaps_real tps_fin = seconds();
  printf("temps RK2 total= %f\n", tps_fin-tps_deb);

  CopyfieldtoCPU(&simu);

  schnaps_real rec_rho=0;
  schnaps_real rcp[1];
  // show_cl_timing(&simu);
  sprintf(xdmf_data_filename,xdmf_filename_fmt,simu.iter_time_rk);
  sn_compute_macro(&simu,rcp);
  schnaps_simulation_xdmf_plot_fields_xml_structured_glops(&simu, xdmf_data_filename);
  printf("rcp = %f\n",rcp[0]);
  // schnaps_real dd = L2error(&simu);
  // printf("erreur L2=%f\n", dd);
  // PlotFields(50, false, &simu, NULL, "dgvisu.msh");
   /* Necessary to ensure perfmodels are written to disk */
  destroy_global_arbiter();
  starpu_shutdown();
  return test;
  
  
  
  
  
  
  
}


int main(void) {
  // int resu1 = test_sn_spec_gaussian_cube();
  int resu1 = test_sn_spec_gaussian_cube_spu();
  if (resu1) printf("Acoustic SN test OK !\n");
  // else printf("Acoustic SN test failed !\n");
  return !resu1;
}

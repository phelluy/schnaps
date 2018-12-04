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
#include "skyline.h"
#include "getopt.h"
#include "global.h"
#include "io.h"
#include "schnaps.h"
#include "../test/test.h"

#define HEADER_MAX_SIZE 9999
#define PATH_MAX_SIZE 9999
#define EXPORT_XDMF 0
#define EXPORT_CSV_RCP 1

#define TICK(X) clock_t X = clock()
#define TOCK(X) printf("time %s: %g sec.\n", (#X), (double)(clock() - (X)) / CLOCKS_PER_SEC)



/*chdir folder + time*/
const int  run_name_size = sizeof("M1_spec_gaussian_cube");
const char run_name[sizeof("m1_spec_gaussian_cube")] = "m1_spec_gaussian_cube";

/*xdmf filename and format*/
const int  xdmf_filename_size = sizeof("m1_spec_gaussian_cube_000000.xdmf");

const char xdmf_filename_fmt[sizeof("m1_spec_gaussian_cube_%06d.xdmf")] 
                            = "m1_spec_gaussian_cube_%06d.xdmf"; 
/*receptor csv file*/
const char csv_filename_density[sizeof("m1_spec_gaussian_cube_rcp_density.csv")] 
                                = "m1_spec_gaussian_cube_rcp_density.csv"; 
                              
const char csv_filename_intensity[sizeof("m1_spec_gaussian_cube_rcp_intensity.csv")] 
                                = "m1_spec_gaussian_cube_rcp_intensity.csv"; 


void m1_utils_change_wd()
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

void m1_rcp_csv_header(
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

void m1_rcp_csv_line(
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

void m1_rcp_csv_time_exporter(
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
    m1_rcp_csv_line(tnow, rcpt_val, data_line);
    fprintf(csvfile,data_line);
  }  else {
    /* append receptor's data */
    csvfile = fopen(filename, "a+" );
    m1_rcp_csv_line(tnow, rcpt_val, data_line);
    fprintf(csvfile, data_line);
  }
  
  fclose(csvfile);
  return;
}



void m1_rk2_export_on(
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
  schnaps_real rec_rho_ar[1]={0};
  
   if (simu->nb_diags != 0) {
     simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
   }

  int fsize =  simu->wsize / simu->macromesh.nbelems;
  while (simu->tnow < tmax ) {
    // sprintf(xdmf_data_filename,xdmf_filename_fmt,simu->iter_time_rk);
    // schnaps_simulation_xdmf_plot_fields_xml_structured_glops(simu, xdmf_data_filename);
    
    
    for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
      DGAverage(simu->fd + ie, simu->w + ie * fsize);
    }

    
    if (EXPORT_CSV_RCP) {
      if(iter == 0){
        m1_compute_rcp(simu, &rec_rho_ar[0]);
        printf("t = %f Receptor value = %f \n", simu->tnow, rec_rho_ar[0]);
        m1_rcp_csv_time_exporter(true,
                                 csv_filename_density,
                                 simu->tnow,
                                 rec_rho_ar );
      } else {
        m1_compute_rcp(simu, &rec_rho_ar[0]);
        printf("t = %f Receptor value = %f \n", simu->tnow, rec_rho_ar[0]);
        m1_rcp_csv_time_exporter(false,
                                 csv_filename_density,
                                 simu->tnow,
                                 rec_rho_ar );
      }
    }
    
    /* data export */
     if (iter % freq == 0) {
      printf("t=%f iter=%d/%d dt=%f\n", simu->tnow, iter, simu->itermax_rk, dt);
      
      if (EXPORT_XDMF == 1) {
        if (EXPORT_CSV_RCP == 0)  m1_compute_rcp(simu, &rec_rho_ar[0]);
        sprintf(xdmf_data_filename,xdmf_filename_fmt,simu->iter_time_rk);
        schnaps_simulation_xdmf_plot_fields_xml_structured_glops(simu, xdmf_data_filename);
      }
    }

    /* time solver*/
    simu->dt = 0.5 * dt;
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
    simu->dt=dt;
    iter++;
    simu->iter_time_rk = iter;
  }
    free(wn);
}

int test_Newton(int N, schnaps_real xmin, schnaps_real xmax)
{
  schnaps_real dx = fabs(xmax - xmin) / N;
  schnaps_real R1 = xmin;
  schnaps_real b = 0;
  schnaps_real b_2 = 0;
  schnaps_real b_3 = 0;
  schnaps_real b_4 = 0;
  int iter = 0;
  int iter_2 = 0;
  
  for(int i = 0; i < N; i++){
    TICK(TIME_A);
    b = m1_get_b_householder(R1,&iter);
    TOCK(TIME_A);
    b_2 = m1_get_b_newton(R1,&iter_2);
#ifdef _WITH_GSL
    TICK(TIME_B);
    b_3 =  m1_get_b_brendt_gsl(R1);
    TOCK(TIME_B);
    TICK(TIME_C);
    b_4 = m1_get_b_newton_gsl(R1);
    TOCK(TIME_C);
#endif
    printf("R1 = %f, b = (%f, %f, %f, %f)\n",R1, b, b_2, b_3, b_4);
    // if( iter > iter_2) printf(" fail for R1 = %f\n", R1);
    R1 = R1 + dx;
  }
  return 1;
}



int test_lebedev_integration_pdf()
{
 
  schnaps_real R1 = 0e0;
  schnaps_real w[4] = {1, 0, 0.90, 0};
  schnaps_real vnorm[3] = {1 , 0, 0};
  schnaps_real int_leb_L[4] = {0,0,0,0};
  schnaps_real int_leb_R[4] = {0,0,0,0};
  schnaps_real error[4] = {0, 0, 0, 0};
  schnaps_real b_norm = 0;
  schnaps_real I_norm = 0;
  schnaps_real bx = 0;
  schnaps_real by = 0;
  schnaps_real bz = 0;
  schnaps_real a;
  schnaps_real fi;
  
  schnaps_real vin;
  schnaps_real vix;
  schnaps_real viy;
  schnaps_real viz;
  schnaps_real wi;
  schnaps_real wifivin = 0;
  int iter;
  

  struct acoustic_data lebedev_set[5] = { lebedev_6, lebedev_14,
                                          lebedev_26, lebedev_38,
                                          lebedev_50 };

                                        
  /* Lebedev set */ 
  for (int k = 0; k < 5; k++) { 
    /* Maxwelian reconstruction */
    I_norm = sqrt(w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
    
#ifdef _WITH_GSL
    b_norm = m1_get_b_mixed_gsl(I_norm / w[0]);
#else
    b_norm = m1_get_b_newton(I_norm / w[0], &iter);
#endif

    if(b_norm <= 10e-12) {
      
      bx = b_norm;
      by = b_norm;
      bz = b_norm;
      schnaps_real t1 = b_norm * b_norm;
      
      /* Use 5-th Order taylor for a*/
      a = w[0] * (0.1e1 - (0.1e1 / 0.6e1) * t1
         + (0.7e1 / 360e0) * t1 * t1);
         
    } else {
      bx = b_norm * w[1] / I_norm;
      by = b_norm * w[2] / I_norm;
      bz = b_norm * w[3] / I_norm;
      a = w[0] * b_norm / sinh(b_norm);
    }
        

    
    /* Left and right integrals*/
    for (int p = 0; p < 4; p++) {
      int_leb_R[p] = 0;
      int_leb_L[p] = 0;
    }
    
    
    
    /* Velocity integration */
    for (int i = 0; i < lebedev_set[k].nb_v; i++) {
      
      vix = lebedev_set[k].v_tab[3 * i];
      viy = lebedev_set[k].v_tab[3 * i + 1];
      viz = lebedev_set[k].v_tab[3 * i + 2];
      wi = lebedev_set[k].w_tab[i];
      
      // vin = vix * vnorm[0] + viy * vnorm[1] + viz * vnorm[2];

      /* Right half sphere */
      // if (vin >= 0) {
      fi = a * exp(bx * vix + by * viy + bz * viz);
      wifivin = wi * fi; 

      int_leb_R[0] += wifivin; 
      int_leb_R[1] += wifivin * vix; 
      int_leb_R[2] += wifivin * viy; 
      int_leb_R[3] += wifivin * viz; 

      
    /* Left half sphere */
    // } else {
      // fi = a * exp(bx * vix + by * viy + bz * viz);
      // wifivin = wi * fi; 

      // int_leb_L[0] += wifivin; 
      // int_leb_L[1] += wifivin * vix; 
      // int_leb_L[2] += wifivin * viy; 
      // int_leb_L[3] += wifivin * viz; 
      // }
    } /* Velocity */

  /* Compute error */

    error[0] = fabs( w[0] - ( int_leb_R[0] + int_leb_L[0] ));
    error[1] = fabs( w[1] - ( int_leb_R[1] + int_leb_L[1] ));
    error[2] = fabs( w[2] - ( int_leb_R[2] + int_leb_L[2] ));
    error[3] = fabs( w[3] - ( int_leb_R[3] + int_leb_L[3] ));

    printf("#Error computation for R1 = %.15f lebedev_%d\n",I_norm / w[0], lebedev_set[k].nb_v);
    printf(" err(w[0]) = %.15f \n", error[0]);
    printf(" err(w[1]) = %.15f \n", error[1]);
    printf(" err(w[2]) = %.15f \n", error[2]);
    printf(" err(w[3]) = %.15f \n", error[3]);
  } /* Lebedev set */
    return 1;
}

int TestDGAverage(void)
{

  int test = true;

  Model model;
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = NULL;
  model.BoundaryFlux = NULL;
  model.InitData = m1_init_test_average;
  model.ImposedData = m1_imposed_test_average;
  model.Source = NULL ;


  int deg[]={1, 1, 1};
  int raf[]={16, 16, 16};
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  BuildConnectivity(&mesh);

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
        // exit(-1);
  simu.tnow = 0;
  schnaps_simulation_xdmf_plot_fields_xml_structured_glops(&simu, "test_dgaverage_pre.xdmf");
  const int fsize =  simu.wsize / simu.macromesh.nbelems;
  schnaps_real w_res[fsize];
  // for (int i = 0; i < fsize; i++) {
    // w_res[i] = 1;
  // }
  
  
  for (int ie = 0; ie < simu.macromesh.nbelems; ie++) {
    field *f = simu.fd + ie;
    DGAverage(simu.fd + ie, simu.w + ie * fsize);
  }
  
  schnaps_simulation_xdmf_plot_fields_xml_structured_glops(&simu, "test_dgaverage_post.xdmf");

  FreeMacroMesh(&mesh);
  
  return test;
};


int test_m1_spec_gaussian_cube(void) { 
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
  // schnaps_acoustic_data = lebedev_50 ;
  schnaps_acoustic_data = lebedev_50;
  // schnaps_acoustic_data = lebedev_6 ;
  acoustic_data *ad = &schnaps_acoustic_data ;
 
  /*Model parameters*/
  
  Model model;
  model.m = 5;
  model.InitData = m1_init_gaussian_density;
  model.ImposedData =  m1_imposed_gaussian_density;
  model.Source      = NULL;
  // model.NumFlux = m1_num_flux_kinetic_discrete;
  model.NumFlux = m1_num_flux_kinetic;
  model.BoundaryFlux = m1_boundary_flux; 
  // model.BoundaryFlux = m1_boundary_flux_discrete;

  Simulation simu;
  EmptySimulation(&simu);
  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  /*Resume Simulation*/
  
  // simu.iter_time_rk=29;
  // simu.tnow = 0.01* (schnaps_real)   simu.iter_time_rk; 
  #ifdef _WITH_HDF5
  // schnaps_simulation_load_field_state(&simu,"fooo");
  #endif

  simu.cfl = 0.95;
  schnaps_real tmax = 30;
  schnaps_real tps_deb = seconds();
  m1_utils_change_wd();
  // schnaps_simulation_xdmf_plot_fields_xml_structured_glops(&simu, "t=0.xdmf");
  // return 1;
  m1_rk2_export_on(&simu, tmax);
  schnaps_real tps_fin = seconds();
  printf("temps = %f\n", tps_fin-tps_deb);
  FreeMacroMesh(&mesh);
  test = true;
  return test;
}


int main(void) {
  // TestPLU();
  // TestDGAverage();
  // schnaps_acoustic_data = lebedev_38;
  // schnaps_real w[4] = {1.0, 0.57, 0.57, 0.57};
  // schnaps_real b[3] = {0.4, 0.4, 0.4};
  // m1_discrete_newton(w,b);
  // printf("b = '%f, %f, %f\n",b[0],b[1],b[2]);
  // schnaps_real bnorm = sqrt( b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  // printf("b = '%f, %f, %f\n",b[0]/bnorm,b[1]/bnorm,b[2]/bnorm);
  // schnaps_real R1 = 1e-12;
  // schnaps_real r = m1_get_b_brendt_gsl(R1);
  // printf("r = %.15f\n",r);
  // r = m1_get_b_newton_gsl(R1); 
  // printf("r = %.15f\n", r);
  
  // int resu = test_lebedev_integration_pdf();
  
  
  int resu = test_m1_spec_gaussian_cube();
  // int resu = test_Newton(2*10000, 1e-14, 0.1);
  // if (resu) printf("Acoustic M1 test OK !\n");
  // else printf("Acoustic M1 test failed !\n");
  // return !resu;
}



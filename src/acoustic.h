/* M1 and SN routines for acoustic simulation.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 and
 * only version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not see <http://www.gnu.org/licenses/>.
 */

#ifndef _ACOUSTIC_
#define _ACOUSTIC_

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#include "schnaps.h"
#include "skyline.h"

#ifdef _WITH_GSL
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#endif


#pragma start_opencl

#define _DERIVATIVE_STEP 		1E-6
#define _NewtonTol 					1E-7
#define _NewtonIter 					30
#define _NewtonInitGuess			0.40
#pragma end_opencl

#define HEADER_MAX_SIZE 9999
#define PATH_MAX_SIZE 9999

#pragma start_opencl
typedef struct acoustic_data{
  int nb_v; //number of discrete SN points
  bool is_cst_w;
  schnaps_real v_tab[3*100]; //!< nodes coordinates array
  schnaps_real w_tab[3*100];
  int spec_tab[3][100]; //Normal array 
  schnaps_real v_sound;
  int nb_rcp;
  schnaps_real rcp_data[4*10];
} acoustic_data;
#pragma end_opencl

// #pragma start_opencl
extern __constant acoustic_data schnaps_acoustic_data;
extern __constant acoustic_data lebedev_6;
extern __constant acoustic_data lebedev_14;
extern __constant acoustic_data lebedev_26;
extern __constant acoustic_data lebedev_38;
extern __constant acoustic_data lebedev_50;
// #pragma end_opencl

/*
================================================================================
                    OCL TEST ROUTINES 
================================================================================
*/

#pragma start_opencl
void ocl_sin3D_imposed_data(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w);
    
void ocl_sin3D_init_data(
    schnaps_real *x,
    schnaps_real *w);
    
void ocl_sin3d_upwind_numflux(
    schnaps_real *wL, 
    schnaps_real *wR, 
    schnaps_real* vnorm, 
    schnaps_real* flux);
        
void ocl_sin3d_boundary_numflux(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL, 
    schnaps_real *vnorm, 
    schnaps_real *flux);

schnaps_real ocl_test_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t);

schnaps_real gaussian_density_v(
    const schnaps_real *x,
    const schnaps_real t,
    const schnaps_real vx,
    const schnaps_real vy,
    const schnaps_real vz);
        
        
void ocl_s6_imposed_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w);
    
void ocl_s6_init_gaussian_density(
    schnaps_real *x, 
    schnaps_real *w);
    
void ocl_s6_upwind_numflux_unroll(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux);

void ocl_s6_free_boundary_unroll(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux);
    
void ocl_s6_upwind_numflux(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux);
    
void ocl_s6_free_boundary(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux);   
 
void ocl_s50_imposed_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w);
    
void ocl_s50_init_gaussian_density(
    schnaps_real *x, 
    schnaps_real *w);
 
void ocl_s50_imposed_gaussian_density_unroll(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w);
    

void ocl_s50_upwind_numflux(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux);
    
void ocl_s50_free_boundary(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux);     
    
void ocl_s50_upwind_numflux_unroll(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux);
    
void ocl_s50_free_boundary_unroll(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux);  
        
void ocl_test_sn_mixed_boundary_flux(
        schnaps_real *x, 
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm, 
        schnaps_real *flux);
#pragma end_opencl



#pragma start_opencl
schnaps_real gaussian_density(
        const schnaps_real *x,
        const schnaps_real t);   
      
void sn_imposed_gaussian_density(
        const schnaps_real *x,
        const schnaps_real t,
        schnaps_real *w);

void sn_init_gaussian_density(
        schnaps_real *x, 
        schnaps_real *w);
        
        
void sn_imposed_constant_density(
        const schnaps_real *x, 
        const schnaps_real t,
        schnaps_real *w);

void sn_init_constant_density(
        schnaps_real *x, 
        schnaps_real *w);        
      
void sn_upwind_numflux(
        schnaps_real *wL,
        schnaps_real *wR,
        schnaps_real *vnorm,
        schnaps_real *flux);
     
void sn_specular_boundary_flux(
        schnaps_real *x, 
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm,
        schnaps_real *flux);
     
void sn_diffusive_boundary_flux(
        schnaps_real *x,
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm, 
        schnaps_real *flux);     
      
void sn_diffusive_boundary_flux_scalar_product(
        schnaps_real *x,
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm, 
        schnaps_real *flux);
      
void sn_diffusive_boundary_flux_equal(
        schnaps_real *x,
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm, 
        schnaps_real *flux);   

void sn_diffusive_boundary_flux_lebedev(
        schnaps_real *x,
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm, 
        schnaps_real *flux);
     
void sn_mixed_boundary_flux(
        schnaps_real *x, 
        schnaps_real  t, 
        schnaps_real *wL,
        schnaps_real *vnorm,
        schnaps_real *flux);
#pragma end_opencl


schnaps_real sn_rcp_sphere_mask( 
        const schnaps_real *x, 
        const schnaps_real rcp_data[4]);
                                 

schnaps_real sn_integrate_field(
        Simulation *simu, 
        int id_field);

schnaps_real sn_l2_error_all_field(
        Simulation *simu);

schnaps_real sn_l1_error_field(
        Simulation *simu, 
        int id_field);

schnaps_real sn_linfinity_error_field(
        Simulation *simu, 
        int id_field);
      
void sn_compute_macro(
        Simulation *simu,
        schnaps_real *rcp_density); 

#pragma start_foo

schnaps_real m1_func_to_min(
        const schnaps_real x, 
        schnaps_real r1);

schnaps_real m1_func_to_min_der(
        const schnaps_real x);

schnaps_real m1_get_b_newton(
        const schnaps_real r1, 
        int *iter);

schnaps_real m1_get_b_householder(
        const schnaps_real r1, 
        int *iter);        



schnaps_real m1_get_r2(
        const schnaps_real r1);
#pragma end_foo

schnaps_real m1_get_r2_der(
        const schnaps_real r1, 
        const schnaps_real fmid, 
        const int order, 
        schnaps_real h);

        
#ifdef _WITH_GSL
struct root_finding_1arg{
  schnaps_real R1;
};  
  
 
// Used for passing 2 (int/schnaps_real) arguments to a gsl_function
struct func_args{
  int argI1;
  schnaps_real argD1;
};  

//Used for passing 4 schnaps_reals arguments to a gsl_function
struct func_args2{
  schnaps_real argD1;
  schnaps_real argD2;
  schnaps_real argD3;
  schnaps_real argD4;
};
  

schnaps_real m1_get_b_mixed_gsl(
  const schnaps_real R1);
  
  
schnaps_real m1_get_b__func_f(
      const schnaps_real x, 
      void *args);
      
      
schnaps_real m1_get_b_func_df(
      const schnaps_real x,
      void *args);
      
      
void m1_get_b_func_fdf(
      schnaps_real x, 
      void *args,
      schnaps_real *f,
      schnaps_real *df);
     
schnaps_real m1_get_b_brendt_gsl(
      const schnaps_real R1);

schnaps_real m1_get_b_newton_gsl(
      const schnaps_real R1);
      
schnaps_real m1_get_dr2_gsl(
      const schnaps_real R1, 
      const int type, 
      const schnaps_real h);
      
      
schnaps_real m1_get_r2_gsl(
      const schnaps_real R1,
      void * params);
      
int m1_get_eigen(
      const schnaps_real R1,
      schnaps_real *R2, 
      schnaps_real *VP1,
      schnaps_real *VP2);
      
      
schnaps_real m1_get_eigen_i(
      const int i12,
      const schnaps_real R1);
      
int m1_exact_riemann_solver_gsl(
      const schnaps_real xi,
      const schnaps_real uL[2],
      const schnaps_real uR[2],
      schnaps_real u[2]);
      
schnaps_real m1_riemann_mid_state_to_min_gsl(
      const schnaps_real x,
      void *args);
      
schnaps_real m1_riemann_get_mid_state_gsl(
      const schnaps_real uL[2],
      const schnaps_real uR[2],
      const schnaps_real tol);
      
schnaps_real m1_riemann_li(
      const int i12,
      const schnaps_real R1m,
      const schnaps_real R1p);
      
schnaps_real m1_riemann_shock(
      const int i12,
      const schnaps_real R1m, 
      const schnaps_real R1p);
      
schnaps_real m1_riemann_rarefaction(
      const int i12,
      const schnaps_real R1m,
      const schnaps_real R1p);
      
schnaps_real m1_riemann_get_r1_gsl(
      const int i12,
      const schnaps_real xi,
      const schnaps_real tol);
      
schnaps_real m1_riemann_to_min_gsl(
      const schnaps_real x,
      void *args);
      
schnaps_real m1_riemann_invariant_to_int_gsl(
      schnaps_real R1,
      void * args);
      
schnaps_real m1_riemann_get_invariant_gsl(
      int i12,
      const schnaps_real R1m,
      const schnaps_real R1p);
#endif


#pragma start_foo
        
void m1_imposed_test_average(
        const schnaps_real *x,
        const schnaps_real t,
        schnaps_real *w);
        
        
void m1_init_test_average(
        schnaps_real *x,
        schnaps_real *w);

void m1_imposed_gaussian_density(
        const schnaps_real *x,
        const schnaps_real t,
        schnaps_real *w);

void m1_init_gaussian_density(
        schnaps_real *x, 
        schnaps_real *w);

void m1_discrete_compute_fk_jac(
        schnaps_real xk[3], 
        schnaps_real w[4],
        schnaps_real fk[3],
        schnaps_real jac[9]);
        
void m1_discrete_newton(
        schnaps_real w[4],
        schnaps_real b[3]);

void m1_num_flux_kinetic(
        schnaps_real * wL,
        schnaps_real * wR,
        schnaps_real * vnorm,
        schnaps_real * flux);
        
void m1_num_flux_rusanov(
        schnaps_real * wL,
        schnaps_real * wR,
        schnaps_real * vnorm,
        schnaps_real * flux);
        
void m1_num_flux_kinetic_discrete(
        schnaps_real * wL,
        schnaps_real * wR,
        schnaps_real * vnorm,
        schnaps_real * flux);
        
        
void m1_boundary_flux(
        schnaps_real *x,
        schnaps_real t,
        schnaps_real *wL,
        schnaps_real *vnorm,
        schnaps_real *flux);
        
void m1_boundary_flux_discrete(
        schnaps_real *x,
        schnaps_real t,
        schnaps_real *wL,
        schnaps_real *vnorm,
        schnaps_real *flux);
        
#pragma end_foo

void m1_compute_rcp(
        Simulation *simu,
        schnaps_real *rcp_density);

#endif

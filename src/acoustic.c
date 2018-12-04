/* M1 and SN routines for acoustic numerical simulation.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3 and
 * only version 3 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not see <http://www.gnu.org/licenses/>.
 */
 
#include "acoustic.h"

#pragma start_opencl
__constant acoustic_data schnaps_acoustic_data = {0};
#pragma end_opencl

#pragma start_opencl
__constant acoustic_data lebedev_6 = {
 .nb_v = 6, 
 .is_cst_w = false,
 .v_tab = { 
      1, 0, 0,
      -1, 0, 0,
      0, 1, 0,
      0, -1, 0,
      0, 0, 1,
      0, 0, -1,
      },
 .w_tab = {
      0.166666666666667,
      0.166666666666667,
      0.166666666666667,
      0.166666666666667,
      0.166666666666667,
      0.166666666666667,
      },
  .spec_tab = { 
         {1,0,2,3,4,5},
         {0,1,3,2,4,5},
         {0,1,2,3,5,4},
        },
 .v_sound = 1.0,
 .nb_rcp = 1,
 .rcp_data = { 
        0.8,
        0.5,
        0.5,
        0.05,
       }
};
#pragma end_opencl

#pragma start_opencl
__constant acoustic_data lebedev_14 = {
 .nb_v = 14, 
 .is_cst_w = false,
 .v_tab = { 
      1, 0, 0,
      -1, 0, 0,
      0, 1, 0,
      0, -1, 0,
      0, 0, 1,
      0, 0, -1,
      0.57735026919, 0.57735026919, 0.57735026919,
      0.57735026919, 0.57735026919, -0.57735026919,
      0.57735026919, -0.57735026919, 0.57735026919,
      0.57735026919, -0.57735026919, -0.57735026919,
      -0.57735026919, 0.57735026919, 0.57735026919,
      -0.57735026919, 0.57735026919, -0.57735026919,
      -0.57735026919, -0.57735026919, 0.57735026919,
      -0.57735026919, -0.57735026919, -0.57735026919,
      },
 .w_tab = {
      0.066666666666667,
      0.066666666666667,
      0.066666666666667,
      0.066666666666667,
      0.066666666666667,
      0.066666666666667,
      0.075000000000000,
      0.075000000000000,
      0.075000000000000,
      0.075000000000000,
      0.075000000000000,
      0.075000000000000,
      0.075000000000000,
      0.075000000000000,
      },
  .spec_tab = {
         {1,0,2,3,4,5,10,11,12,13,6,7,8,9},
         {0,1,3,2,4,5,8,9,6,7,12,13,10,11},
         {0,1,2,3,5,4,7,6,9,8,11,10,13,12},
        },
 .v_sound = 1.0,
 .nb_rcp = 1,
 .rcp_data = { 
        0.8,
        0.5,
        0.5,
        0.05,
       }
};
#pragma end_opencl

#pragma start_opencl
__constant acoustic_data lebedev_26 = {
 .nb_v = 26, 
 .is_cst_w = false,
 .v_tab = { 
      1, 0, 0,
      -1, 0, 0,
      0, 1, 0,
      0, -1, 0,
      0, 0, 1,
      0, 0, -1,
      0,        0.707106781187,  0.707106781187,
      0,        0.707106781187, -0.707106781187,
      0,        -0.707106781187,  0.707106781187,
      0,        -0.707106781187, -0.707106781187,
      0.707106781187,  0,        0.707106781187,
      0.707106781187,  0,        -0.707106781187,
      -0.707106781187,  0,        0.707106781187,
      -0.707106781187,  0,        -0.707106781187,
      0.707106781187,  0.707106781187,  0,
      0.707106781187, -0.707106781187,  0,
      -0.707106781187,  0.707106781187,  0,
      -0.707106781187, -0.707106781187,  0,
      0.57735026919,  0.57735026919,  0.57735026919,
      0.57735026919,  0.57735026919,  -0.57735026919,
      0.57735026919,  -0.57735026919,  0.57735026919,
      0.57735026919,  -0.57735026919,  -0.57735026919,
      -0.57735026919,  0.57735026919,  0.57735026919,
      -0.57735026919,  0.57735026919,  -0.57735026919,
      -0.57735026919,  -0.57735026919,  0.57735026919,
      -0.57735026919,  -0.57735026919,  -0.57735026919,
      },
 .w_tab = {
      0.047619047619048,
      0.047619047619048,
      0.047619047619048,
      0.047619047619048,
      0.047619047619048,
      0.047619047619048,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.038095238095238,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      },
  .spec_tab = {
         {1,0,2,3,4,5,6,7,8,9,12,13,10,11,16,
         17,14,15,22,23,24,25,18,19,20,21},
         {0,1,3,2,4,5,8,9,6,7,10,11,12,13,15,
         14,17,16,20,21,18,19,24,25,22,23},
         {0,1,2,3,5,4,7,6,9,8,11,10,13,12,14,
         15,16,17,19,18,21,20,23,22,25,24},
        },
 .v_sound = 1.0,
 .nb_rcp  = 1,
 .rcp_data = { 
        0.8,
        0.5,
        0.5,
        0.05,
       }
};
#pragma end_opencl

#pragma start_opencl
__constant acoustic_data lebedev_38 = {
 .nb_v = 38, 
 .is_cst_w = false,
 .v_tab = { 
      1, 0, 0,
      -1, 0, 0,
      0, 1, 0,
      0, -1, 0,
      0, 0, 1,
      0, 0, -1,
      0.57735026919,  0.57735026919,  0.57735026919,
      0.57735026919,  0.57735026919,  -0.57735026919,
      0.57735026919,  -0.57735026919,  0.57735026919,
      0.57735026919,  -0.57735026919,  -0.57735026919,
      -0.57735026919,  0.57735026919,  0.57735026919,
      -0.57735026919,  0.57735026919,  -0.57735026919,
      -0.57735026919,  -0.57735026919,  0.57735026919,
      -0.57735026919,  -0.57735026919,  -0.57735026919,
      0.459700843381,  0.888073833977,  0,
      0.459700843381, -0.888073833977,  0,
      -0.459700843381,  0.888073833977,  0,
      -0.459700843381, -0.888073833977,  0,
      0.888073833977,  0.459700843381,  0,
      0.888073833977, -0.459700843381,  0,
      -0.888073833977,  0.459700843381,  0,
      -0.888073833977, -0.459700843381,  0,
      0.459700843381,  0,        0.888073833977,
      0.459700843381,  0,        -0.888073833977,
      -0.459700843381,  0,        0.888073833977,
      -0.459700843381,  0,        -0.888073833977,
      0.888073833977,  0,        0.459700843381,
      0.888073833977,  0,        -0.459700843381,
      -0.888073833977,  0,        0.459700843381,
      -0.888073833977,  0,        -0.459700843381,
      0,        0.459700843381,  0.888073833977,
      0,        0.459700843381, -0.888073833977,
      0,        -0.459700843381,  0.888073833977,
      0,        -0.459700843381, -0.888073833977,
      0,        0.888073833977,  0.459700843381,
      0,        0.888073833977, -0.459700843381,
      0,        -0.888073833977,  0.459700843381,
      0,        -0.888073833977, -0.459700843381,
      },
 .w_tab = {
      0.009523809523810,
      0.009523809523810,
      0.009523809523810,
      0.009523809523810,
      0.009523809523810,
      0.009523809523810,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.032142857142857,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      0.028571428571429,
      },
  .spec_tab = {
         {1,0,2,3,4,5,10,11,12,13,6,7,8,9,16,17,14,15,20,21,
         18,19,24,25,22,23,28,29,26,27,30,31,32,33,34,35,36,37},
         {0,1,3,2,4,5,8,9,6,7,12,13,10,11,15,14,17,16,19,18,
         21,20,22,23,24,25,26,27,28,29,32,33,30,31,36,37,34,35},
         {0,1,2,3,5,4,7,6,9,8,11,10,13,12,14,15,16,17,18,19,
         20,21,23,22,25,24,27,26,29,28,31,30,33,32,35,34,37,36},
        },
 .v_sound = 1.0,
 .nb_rcp  = 1,
 .rcp_data = { 
        0.8,
        0.5,
        0.5,
        0.05,
       }
};
#pragma end_opencl

#pragma start_opencl
__constant acoustic_data lebedev_50 = {
 .nb_v = 50, 
 .is_cst_w = false,
 .v_tab = { 
      1, 0, 0,
      -1, 0, 0,
      0, 1, 0,
      0, -1, 0,
      0, 0, 1,
      0, 0, -1,
      0,        0.707106781187, 0.707106781187,
      0,        0.707106781187, -0.707106781187,
      0,       -0.707106781187, 0.707106781187,
      0,       -0.707106781187, -0.707106781187,
      0.707106781187, 0,        0.707106781187,
      0.707106781187, 0,       -0.707106781187,
      -0.707106781187, 0,        0.707106781187,
      -0.707106781187, 0,       -0.707106781187,
      0.707106781187, 0.707106781187, 0,       
      0.707106781187, -0.707106781187, 0,       
      -0.707106781187, 0.707106781187, 0,       
      -0.707106781187, -0.707106781187, 0,       
      0.57735026919,  0.57735026919,  0.57735026919, 
      0.57735026919,  0.57735026919, -0.57735026919, 
      0.57735026919, -0.57735026919,  0.57735026919, 
      0.57735026919, -0.57735026919, -0.57735026919, 
      -0.57735026919,  0.57735026919,  0.57735026919, 
      -0.57735026919,  0.57735026919, -0.57735026919, 
      -0.57735026919, -0.57735026919,  0.57735026919, 
      -0.57735026919, -0.57735026919, -0.57735026919, 
      0.301511344578, 0.301511344578, 0.904534033733,
      0.301511344578, 0.301511344578, -0.904534033733,
      0.301511344578, -0.301511344578, 0.904534033733,
      0.301511344578, -0.301511344578, -0.904534033733,
      -0.301511344578, 0.301511344578, 0.904534033733,
      -0.301511344578, 0.301511344578, -0.904534033733,
      -0.301511344578, -0.301511344578, 0.904534033733,
      -0.301511344578, -0.301511344578, -0.904534033733,
      0.301511344578, 0.904534033733, 0.301511344578,
      0.301511344578, -0.904534033733, 0.301511344578,
      0.301511344578, 0.904534033733, -0.301511344578,
      0.301511344578, -0.904534033733, -0.301511344578,
      -0.301511344578, 0.904534033733, 0.301511344578,
      -0.301511344578, -0.904534033733, 0.301511344578,
      -0.301511344578, 0.904534033733, -0.301511344578,
      -0.301511344578, -0.904534033733, -0.301511344578,
      0.904534033733, 0.301511344578, 0.301511344578,
      -0.904534033733, 0.301511344578, 0.301511344578,
      0.904534033733, 0.301511344578, -0.301511344578,
      -0.904534033733, 0.301511344578, -0.301511344578,
      0.904534033733, -0.301511344578, 0.301511344578,
      -0.904534033733, -0.301511344578, 0.301511344578,
      0.904534033733, -0.301511344578, -0.301511344578,
      -0.904534033733, -0.301511344578, -0.301511344578,
      },
 .w_tab = {
      0.012698412698413,
      0.012698412698413,
      0.012698412698413,
      0.012698412698413,
      0.012698412698413,
      0.012698412698413,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.022574955908289,
      0.021093750000000,
      0.021093750000000,
      0.021093750000000,
      0.021093750000000,
      0.021093750000000,
      0.021093750000000,
      0.021093750000000,
      0.021093750000000,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      0.020173335537919,
      },
  .spec_tab = { {1,0,2,3,4,5,6,7,8,9,12,13,10,11,
         16,17,14,15,22,23,24,25,18,19,20,
         21,30,31,32,33,26,27,28,29,38,39,
         40,41,34,35,36,37,43,42,45,44,47,
         46,49,48},
         {0,1,3,2,4,5,8,9,6,7,10,11,12,13,
         15,14,17,16,20,21,18,19,24,25,22,
         23,28,29,26,27,32,33,30,31,35,34,
         37,36,39,38,41,40,46,47,48,49,42,
         43,44,45},
         {0,1,2,3,5,4,7,6,9,8,11,10,13,12,
         14,15,16,17,19,18,21,20,23,22,25,
         24,27,26,29,28,31,30,33,32,36,37,
         34,35,40,41,38,39,44,45,42,43,48,
         49,46,47},
        },
 .v_sound = 1.0,
 .nb_rcp  = 1,
 .rcp_data = { 
        0.8,
        0.5,
        0.5,
        0.05,
       }
};
#pragma end_opencl

/*
================================================================================
          OCL TEST ROUTINES 
================================================================================
*/



#pragma start_opencl
void ocl_sin3D_imposed_data(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w)
{
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  // const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
  const schnaps_real transport_v[] = {1,0,0};
  schnaps_real vx
    = transport_v[0] * x[0]
    + transport_v[1] * x[1]
    + transport_v[2] * x[2];
  schnaps_real xx = vx - t;
  w[0] = sin(10*xx);
  // w[0] = xx * xx;

}
#pragma end_opencl


#pragma start_opencl
void ocl_sin3D_init_data(
    schnaps_real *x,
    schnaps_real *w)
{
    ocl_sin3D_imposed_data(x, 0, w);
}
#pragma end_opencl


#pragma start_opencl
void ocl_sin3d_upwind_numflux(
        schnaps_real *wL, 
        schnaps_real *wR, 
        schnaps_real* vnorm, 
        schnaps_real* flux)
{
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  // const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
  const schnaps_real transport_v[] = {1, 0, 0};
  schnaps_real vn
    = transport_v[0] * vnorm[0]
    + transport_v[1] * vnorm[1]
    + transport_v[2] * vnorm[2];
  schnaps_real vnp = vn > 0 ? vn : 0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
}
#pragma end_opencl


#pragma start_opencl
void ocl_sin3d_boundary_numflux(
        schnaps_real *x, 
        schnaps_real t, 
        schnaps_real *wL, 
        schnaps_real *vnorm, 
        schnaps_real *flux)
{
  schnaps_real wR[1];
  ocl_sin3D_imposed_data(x, t, wR);
  ocl_sin3d_upwind_numflux(wL, wR, vnorm, flux);
}
#pragma end_opencl


#pragma start_opencl
schnaps_real ocl_test_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t)
{
 schnaps_real sigx = 0.10;
 schnaps_real sigy = 0.10;
 schnaps_real sigz = 0.10;
 schnaps_real A = 1.0 / (15.7496099458 * sigx * sigy * sigz);
 schnaps_real center[3] = {0.5, 0.5, 0.5};
 
 schnaps_real X = (x[0] - center[0]) * (x[0] - center[0]) / (sigx * sigx);
 schnaps_real Y = (x[1] - center[1]) * (x[1] - center[1]) / (sigy * sigy);
 schnaps_real Z = (x[2] - center[2]) * (x[2] - center[2]) / (sigz * sigz);
 
 return A * exp(- 0.5*(X + Y + Z));
}
#pragma end_opencl


#pragma start_opencl
schnaps_real gaussian_density_v(
    const schnaps_real *x,
    const schnaps_real t,
    const schnaps_real vx,
    const schnaps_real vy,
    const schnaps_real vz)
{
/*
* 3D spatial gaussian .
*
*/
 schnaps_real sigx = 0.2;
 schnaps_real sigy = 0.2;
 schnaps_real sigz = 0.2;
 schnaps_real A = 1.0 / (15.7496099458 * sigx * sigy * sigz);
 schnaps_real center[3] = {0.5, 0.5, 0.5};
 
 schnaps_real x0 = (x[0] - vx * t - center[0]) / sigx;
 schnaps_real x1 = (x[1] - vy * t - center[1]) / sigy;
 schnaps_real x2 = (x[2] - vz * t - center[2]) / sigz;
 
 return A * exp(- 0.5*(x0*x0 + x1*x1 + x2*x2)); 

}
#pragma end_opencl


#pragma start_opencl
void ocl_s6_imposed_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w)
{ 
    __constant acoustic_data *ad = &lebedev_6;
    
    #pragma unroll 6
    for (unsigned int i = 0; i < ad->nb_v ; i++) {
        w[i] = gaussian_density_v(x, 
                                 t,
                                 ad->v_tab[3 * i + 0],
                                 ad->v_tab[3 * i + 1],
                                 ad->v_tab[3 * i + 2]);
    }
}
#pragma end_opencl


#pragma start_opencl
void ocl_s6_init_gaussian_density(
    schnaps_real *x, 
    schnaps_real *w)
{
 ocl_s6_imposed_gaussian_density(x, 0, w);
}
#pragma end_opencl


#pragma start_opencl
void ocl_s6_upwind_numflux_unroll(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux)
{

    schnaps_real vx = vnorm[0];
    schnaps_real vy = vnorm[1];
    schnaps_real vz = vnorm[2];
    schnaps_real vn;
    schnaps_real vnp;
    schnaps_real vnm;
 
    vn = vx;
    vnp = vn > 0 ? vn : 0;
    vnm = vn - vnp;
    flux[0] = vnp * wL[0] + vnm * wR[0];
    vn = -vx;
    vnp = vn > 0 ? vn : 0;
    vnm = vn - vnp;
    flux[1] = vnp * wL[1] + vnm * wR[1];
    vn = vy;
    vnp = vn > 0 ? vn : 0;
    vnm = vn - vnp;
    flux[2] = vnp * wL[2] + vnm * wR[2];
    vn = -vy;
    vnp = vn > 0 ? vn : 0;
    vnm = vn - vnp;
    flux[3] = vnp * wL[3] + vnm * wR[3];
    vn = vz;
    vnp = vn > 0 ? vn : 0;
    vnm = vn - vnp;
    flux[4] = vnp * wL[4] + vnm * wR[4];
    vn = -vz;
    vnp = vn > 0 ? vn : 0;
    vnm = vn - vnp;
    flux[5] = vnp * wL[5] + vnm * wR[5];
}
#pragma end_opencl


#pragma start_opencl
void ocl_s6_free_boundary_unroll(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
    schnaps_real wR[6]; 
    ocl_s6_imposed_gaussian_density(x, t, wR);
    ocl_s6_upwind_numflux_unroll(wL, wR, vnorm, flux);
}
#pragma end_opencl


#pragma start_opencl
void ocl_s6_upwind_numflux(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
 
    __constant acoustic_data *ad = &lebedev_6;

    for (unsigned int i = 0; i < ad->nb_v; i++) {
        schnaps_real vn = + vnorm[0] * ad->v_tab[3 * i + 0]
          + vnorm[1] * ad->v_tab[3 * i + 1]
          + vnorm[2] * ad->v_tab[3 * i + 2];
        schnaps_real vnp = vn > 0 ? vn : 0;
        schnaps_real vnm = vn - vnp;
        flux[i] = vnp * wL[i] + vnm * wR[i];
    }
}
#pragma end_opencl


#pragma start_opencl
void ocl_s6_free_boundary(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
    schnaps_real wR[6]; 
    ocl_s6_imposed_gaussian_density(x, t, wR);
    ocl_s6_upwind_numflux(wL, wR, vnorm, flux);
}
#pragma end_opencl


#pragma start_opencl
void ocl_s50_imposed_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w)
{ 
    
    __constant acoustic_data *ad = &lebedev_50;
    schnaps_real rho ;
    // = ocl_test_gaussian_density(x, t) / (4.0 * M_PI);
#pragma unroll 10
    for (unsigned int i = 0; i < 50 ; i++) {
        w[i] = gaussian_density_v(x, 
                                 t,
                                 ad->v_tab[3 * i + 0],
                                 ad->v_tab[3 * i + 1],
                                 ad->v_tab[3 * i + 2]);
    }
}
#pragma end_opencl




#pragma start_opencl
void ocl_s50_init_gaussian_density(
    schnaps_real *x, 
    schnaps_real *w)
{
 ocl_s50_imposed_gaussian_density(x, 0, w);
}
#pragma end_opencl

#pragma start_opencl
void ocl_s50_upwind_numflux_unroll(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
  schnaps_real vx = vnorm[0];
  schnaps_real vy = vnorm[1];
  schnaps_real vz = vnorm[2];
  schnaps_real vn;
  schnaps_real vnp;
  schnaps_real vnm;
  
  schnaps_real t0 = 0.707106781187 * vx;
  schnaps_real t1 = 0.707106781187 * vy;
  schnaps_real t2 = 0.707106781187 * vz;
  schnaps_real t3 = 0.57735026919 * vx;
  schnaps_real t4 = 0.57735026919 * vy;
  schnaps_real t5 = 0.57735026919 * vz;
  schnaps_real t6 = 0.301511344578 * vx;
  schnaps_real t7 = 0.301511344578 * vy;
  schnaps_real t8 = 0.301511344578 * vz;
  schnaps_real t9 = 0.904534033733 * vx;
  schnaps_real t10 = 0.904534033733 * vy;
  schnaps_real t11 = 0.904534033733 * vz;

  vn = vx;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  vn = -vx;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[1] = vnp * wL[1] + vnm * wR[1];
  vn = vy ;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[2] = vnp * wL[2] + vnm * wR[2];
  vn = -vy ;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[3] = vnp * wL[3] + vnm * wR[3];
  vn = 0 * vx + vz;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[4] = vnp * wL[4] + vnm * wR[4];
  vn = 0 * vx + -vz;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[5] = vnp * wL[5] + vnm * wR[5];
  vn = t1 + t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[6] = vnp * wL[6] + vnm * wR[6];
  vn = t1 + -t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[7] = vnp * wL[7] + vnm * wR[7];
  vn = -t1 + t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[8] = vnp * wL[8] + vnm * wR[8];
  vn = -t1 + -t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[9] = vnp * wL[9] + vnm * wR[9];
  vn = t0 + t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[10] = vnp * wL[10] + vnm * wR[10];
  vn = t0 + -t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[11] = vnp * wL[11] + vnm * wR[11];
  vn = -t0 + t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[12] = vnp * wL[12] + vnm * wR[12];
  vn = -t0 + -t2;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[13] = vnp * wL[13] + vnm * wR[13];
  vn = t0 + t1 ;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[14] = vnp * wL[14] + vnm * wR[14];
  vn = t0 + -t1 ;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[15] = vnp * wL[15] + vnm * wR[15];
  vn = -t0 + t1 ;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[16] = vnp * wL[16] + vnm * wR[16];
  vn = -t0 + -t1 ;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[17] = vnp * wL[17] + vnm * wR[17];
  vn = t3 + t4 + t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[18] = vnp * wL[18] + vnm * wR[18];
  vn = t3 + t4 + -t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[19] = vnp * wL[19] + vnm * wR[19];
  vn = t3 + -t4 + t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[20] = vnp * wL[20] + vnm * wR[20];
  vn = t3 + -t4 + -t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[21] = vnp * wL[21] + vnm * wR[21];
  vn = -t3 + t4 + t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[22] = vnp * wL[22] + vnm * wR[22];
  vn = -t3 + t4 + -t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[23] = vnp * wL[23] + vnm * wR[23];
  vn = -t3 + -t4 + t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[24] = vnp * wL[24] + vnm * wR[24];
  vn = -t3 + -t4 + -t5;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[25] = vnp * wL[25] + vnm * wR[25];
  vn = t6 + t7 + t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[26] = vnp * wL[26] + vnm * wR[26];
  vn = t6 + t7 + -t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[27] = vnp * wL[27] + vnm * wR[27];
  vn = t6 + -t7 + t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[28] = vnp * wL[28] + vnm * wR[28];
  vn = t6 + -t7 + -t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[29] = vnp * wL[29] + vnm * wR[29];
  vn = -t6 + t7 + t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[30] = vnp * wL[30] + vnm * wR[30];
  vn = -t6 + t7 + -t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[31] = vnp * wL[31] + vnm * wR[31];
  vn = -t6 + -t7 + t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[32] = vnp * wL[32] + vnm * wR[32];
  vn = -t6 + -t7 + -t11;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[33] = vnp * wL[33] + vnm * wR[33];
  vn = t6 + t10 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[34] = vnp * wL[34] + vnm * wR[34];
  vn = t6 + -t10 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[35] = vnp * wL[35] + vnm * wR[35];
  vn = t6 + t10 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[36] = vnp * wL[36] + vnm * wR[36];
  vn = t6 + -t10 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[37] = vnp * wL[37] + vnm * wR[37];
  vn = -t6 + t10 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[38] = vnp * wL[38] + vnm * wR[38];
  vn = -t6 + -t10 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[39] = vnp * wL[39] + vnm * wR[39];
  vn = -t6 + t10 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[40] = vnp * wL[40] + vnm * wR[40];
  vn = -t6 + -t10 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[41] = vnp * wL[41] + vnm * wR[41];
  vn = t9 + t7 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[42] = vnp * wL[42] + vnm * wR[42];
  vn = -t9 + t7 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[43] = vnp * wL[43] + vnm * wR[43];
  vn = t9 + t7 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[44] = vnp * wL[44] + vnm * wR[44];
  vn = -t9 + t7 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[45] = vnp * wL[45] + vnm * wR[45];
  vn = t9 + -t7 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[46] = vnp * wL[46] + vnm * wR[46];
  vn = -t9 + -t7 + t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[47] = vnp * wL[47] + vnm * wR[47];
  vn = t9 + -t7 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[48] = vnp * wL[48] + vnm * wR[48];
  vn = -t9 + -t7 + -t8;
  vnp = vn > 0 ? vn : 0;
  vnm = vn - vnp;
  flux[49] = vnp * wL[49] + vnm * wR[49];
}
#pragma end_opencl


#pragma start_opencl
void ocl_s50_free_boundary_unroll(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
    schnaps_real wR[50]; 
    ocl_s50_imposed_gaussian_density(x, t, wR);
    ocl_s50_upwind_numflux_unroll(wL, wR, vnorm, flux);
  
}
#pragma end_opencl

#pragma start_opencl
void ocl_s50_upwind_numflux(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
 
    __constant acoustic_data *ad = &lebedev_50;
    schnaps_real vx = vnorm[0];
    schnaps_real vy = vnorm[1];
    schnaps_real vz = vnorm[2];
    #pragma unroll 5
    for (unsigned int i = 0; i < ad->nb_v; i++) {
        schnaps_real vn = vx * ad->v_tab[3 * i + 0] 
                        + vy * ad->v_tab[3 * i + 1]
                        + vz * ad->v_tab[3 * i + 2];
                        
        schnaps_real vnp = vn > 0 ? vn : 0;
        schnaps_real vnm = vn - vnp;
        flux[i] = vnp * wL[i] + vnm * wR[i];
    }
}
#pragma end_opencl


#pragma start_opencl
void ocl_s50_free_boundary(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
    schnaps_real wR[50]; 
    ocl_s50_imposed_gaussian_density(x, t, wR);
    ocl_s50_upwind_numflux(wL, wR, vnorm, flux);
}
#pragma end_opencl






#pragma start_opencl
void ocl_test_sn_mixed_boundary_flux(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm, 
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;


 schnaps_real flux_out = 0;
 schnaps_real flux_diff = 0;
 schnaps_real flux_spec = 0;
 schnaps_real w_sum   = 0; 
 schnaps_real vn    = 0;
 schnaps_real tol = 10e-3;
 schnaps_real absorb = 0.1;
 // schnaps_real alpha = (0.40)*(1-absorb);
 schnaps_real alpha = 0.10;
 // schnaps_real beta = (1 - alpha)*(1-absorb);
 schnaps_real beta = 1-alpha;
 int count = 0;
 int count2 = 0;
 
 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
 }

 
 for (int i = 0; i < ad->nb_v; i++) {
  vn = + vnorm[0] * ad->v_tab[3 * i + 0]
     + vnorm[1] * ad->v_tab[3 * i + 1]
     + vnorm[2] * ad->v_tab[3 * i + 2];
    
  /*Colision Case */
  if (vn > 0 ) {
   /*Compute outflux flux*/
   flux_out += vn * wL[i];
   flux[i] += vn * wL[i];

   /*Compute incoming flux*/
   /*X faces*/
   if (vnorm[0] !=0 && fabs(vnorm[1]) <= tol && fabs(vnorm[2]) <= tol) {
    flux[ad->spec_tab[0][i]] -= (1-absorb) * alpha * vn * wL[i];
    flux_spec += alpha * vn * wL[i];
   }
   /*Y faces*/
   if (fabs(vnorm[0]) <= tol && vnorm[1] !=0 && fabs(vnorm[2]) <= tol) {
    flux[ad->spec_tab[1][i]] -= (1-absorb) *alpha * vn * wL[i];
    flux_spec += alpha * vn * wL[i];
   }
   /*Z faces*/
   if (fabs(vnorm[0]) <= tol && fabs(vnorm[1]) <= tol && vnorm[2] !=0) {
    flux[ad->spec_tab[2][i]] -= (1-absorb) * alpha * vn * wL[i];
    flux_spec += alpha * vn * wL[i];
   }
  } else {
   w_sum += vn;
   count2 += 1;// Compute weights
  }   
 }

 for (int k = 0; k < ad->nb_v; k++) {
   vn = + vnorm[0] * ad->v_tab[3 * k + 0]
      + vnorm[1] * ad->v_tab[3 * k + 1]
      + vnorm[2] * ad->v_tab[3 * k + 2];

   if (vn < 0) {
    flux[k]  +=  - beta * (1.0 - absorb) * flux_out * vn / w_sum;
   }
 }
}
#pragma end_opencl


/*
================================================================================
            END OCL TEST ROUTINES
================================================================================
*/


#pragma start_opencl
schnaps_real gaussian_density(
    const schnaps_real *x,
    const schnaps_real t)
{
 /*
	 * 3D spatial gaussian .
  *
	 */
 schnaps_real sigx = 0.1;
 schnaps_real sigy = 0.1;
 schnaps_real sigz = 0.1;
 schnaps_real A = 1.0 / (15.7496099458 * sigx * sigy * sigz);
 schnaps_real center[3] = {0.5, 0.5, 0.5};
 
 schnaps_real X = (x[0] - center[0]) * (x[0] - center[0]) / (sigx * sigx);
 schnaps_real Y = (x[1] - center[1]) * (x[1] - center[1]) / (sigy * sigy);
 schnaps_real Z = (x[2] - center[2]) * (x[2] - center[2]) / (sigz * sigz);
 
 // schnaps_real T=exp(-5*t);
 return A * exp(- 0.5*(X + Y + Z));
}
#pragma end_opencl






/* Sn model routines */
#pragma start_opencl
void sn_imposed_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 schnaps_real w_src_watt = 0.01;
 schnaps_real dt_mc = 0.1;
 /*Case of isotrop intensity*/
 schnaps_real rho = gaussian_density(x, t)*dt_mc*w_src_watt ;
 for (int i = 0; i < ad->nb_v; i++) {
  w[i] = rho / (4.0 * M_PI);
 }

 /*We set macro quantities to 0*/
 w[ad->nb_v] = rho;
 w[ad->nb_v + 1] = 0;
 w[ad->nb_v + 2] = 0;
 w[ad->nb_v + 3] = 0;
 w[ad->nb_v + 4] = 0;
 return;
}
#pragma end_opencl

#pragma start_opencl
void sn_init_gaussian_density(
    schnaps_real *x, 
    schnaps_real *w)
{
 sn_imposed_gaussian_density(x, 0, w);
}
#pragma end_opencl

#pragma start_opencl
void sn_imposed_constant_density(
    const schnaps_real *x, 
    const schnaps_real t,
    schnaps_real *w)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 
 for (int i = 0; i < ad->nb_v; i++) {
   w[i] = 1.0/(4.0*M_PI);
  }
 
 /*We set macro quantities to 0*/
 w[ad->nb_v] = 1.0;
 w[ad->nb_v + 1] = 0;
 w[ad->nb_v + 2] = 0;
 w[ad->nb_v + 3] = 0;
 w[ad->nb_v + 4] = 0;
 }
#pragma end_opencl
 
#pragma start_opencl
void sn_init_constant_density(
    schnaps_real *x, 
    schnaps_real *w)
{
 sn_imposed_constant_density(x, 0, w);
}
#pragma end_opencl

#pragma start_opencl
void sn_upwind_numflux(
    schnaps_real *wL,
    schnaps_real *wR,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 for (int i = 0; i < ad->nb_v; i++) {
  schnaps_real vn = + vnorm[0] * ad->v_tab[3 * i + 0]
            + vnorm[1] * ad->v_tab[3 * i + 1]
            + vnorm[2] * ad->v_tab[3 * i + 2];

  schnaps_real vnp = vn > 0 ? vn : 0;
  schnaps_real vnm = vn - vnp;
  flux[i] = vnp * wL[i] + vnm * wR[i];
 }
}
#pragma end_opencl

#pragma start_opencl
void sn_specular_boundary_flux(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 schnaps_real tol = 1e-4; 

 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
 }

 for (int i = 0; i < ad->nb_v; i++) {
  schnaps_real vn = + vnorm[0] * ad->v_tab[3 * i + 0]
            + vnorm[1] * ad->v_tab[3 * i + 1]
            + vnorm[2] * ad->v_tab[3 * i + 2];
  /*Colision Case */
  if (vn > tol) {
   /*Compute incoming flux*/
   flux[i] += vn * wL[i];
   /*X faces*/
   if (vnorm[0] !=0 && fabs(vnorm[1]) <= tol && fabs(vnorm[2]) <= tol)
    flux[ad->spec_tab[0][i]] = - flux[i];
   /*Y faces*/
   if (fabs(vnorm[0]) <= tol && vnorm[1] !=0 && fabs(vnorm[2]) <= tol)
    flux[ad->spec_tab[1][i]] = - flux[i];
   /*Z faces*/
   if (fabs(vnorm[0]) <= tol && fabs(vnorm[1]) <= tol && vnorm[2] !=0)
    flux[ad->spec_tab[2][i]] = - flux[i];
  }
 //Find NN for flux redistribution
 }
}
#pragma end_opencl

#pragma start_opencl
void sn_diffusive_boundary_flux(
    schnaps_real *x,
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm, 
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 schnaps_real flux_out = 0;
 schnaps_real flux_in  = 0;
 schnaps_real w_sum   = 0; 
 schnaps_real vn    = 0;

 int weight_type    = 2;
 
 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
 }

 for (int i = 0; i < ad->nb_v; i++) {
  vn = + vnorm[0] * ad->v_tab[3 * i + 0]
     + vnorm[1] * ad->v_tab[3 * i + 1]
     + vnorm[2] * ad->v_tab[3 * i + 2];
    
  /*Colision Case */
  if (vn > 0 ) {
   /*Compute outflux flux*/
   flux_out += vn * wL[i];
   flux[i] += vn * wL[i];
  } 
 }
 
   
 /*Redistribution on the normal sphere*/
 for (int k = 0; k < ad->nb_v; k++) {
  vn = + vnorm[0] * ad->v_tab[3 * k + 0]
     + vnorm[1] * ad->v_tab[3 * k + 1]
     + vnorm[2] * ad->v_tab[3 * k + 2];
    
   /* Type of weights */ 
   if (vn < 0) {
   switch (weight_type) {
    case 0 : // Lebedev
     w_sum += ad->w_tab[k];
     break;

    case 1 : // Equal
     w_sum += 1;
     break;
    
    case 2 : // Scalar product
     w_sum -= vn; 
     break;
   }
  }
 }

 for (int k = 0; k < ad->nb_v; k++) {
   vn = + vnorm[0] * ad->v_tab[3 * k + 0]
      + vnorm[1] * ad->v_tab[3 * k + 1]
      + vnorm[2] * ad->v_tab[3 * k + 2];

   if (vn < 0) {
    switch (weight_type) {
     case 0 : // Lebedev
      flux_in += flux_out / w_sum;
      flux[k] = - flux_out /w_sum;
      break;

     case 1 : // Equal
      flux_in += flux_out / w_sum;
      flux[k] = - flux_out /w_sum;
      break;
     
     case 2 : // Scalar product
      flux_in -= flux_out * vn / w_sum;
      flux[k] = flux_out *vn /w_sum;
      break;
    }
   }
 }
 /*assert for debug*/
 // assert(abs(flux_in - flux_out) == 0);
 // assert(vnorm[0]*vnorm[0] + vnorm[1]*vnorm[1] + vnorm[2]*vnorm[2] == 1);
}
#pragma end_opencl

#pragma start_opencl
void sn_diffusive_boundary_flux_scalar_product(
    schnaps_real *x,
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm, 
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 schnaps_real flux_out = 0;
 schnaps_real flux_in  = 0;
 schnaps_real w_sum   = 0; 
 schnaps_real vn    = 0;
 schnaps_real absorb  = 0.1;
 
 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
  // printf("v = (%f,%f,%f)\n",ad->v_tab[3 * i + 0], ad->v_tab[3 * i + 1], ad->v_tab[3 * i + 2]);
 }
 for (int i = 0; i < ad->nb_v; i++) {
  vn = + vnorm[0] * ad->v_tab[3 * i + 0]
     + vnorm[1] * ad->v_tab[3 * i + 1]
     + vnorm[2] * ad->v_tab[3 * i + 2];
    
  /*Colision Case */
  if (vn > 0 ) {
   /*Compute outflux flux*/
   flux_out += vn * wL[i];
   flux[i] += vn * wL[i];
  } else w_sum -= vn; // Compute weights
 }

 for (int k = 0; k < ad->nb_v; k++) {
   vn = + vnorm[0] * ad->v_tab[3 * k + 0]
      + vnorm[1] * ad->v_tab[3 * k + 1]
      + vnorm[2] * ad->v_tab[3 * k + 2];

   if (vn < 0) {
    flux_in -= flux_out * vn / w_sum;
    flux[k] = absorb*flux_out * vn /w_sum;
   }
 }
}
#pragma end_opencl

#pragma start_opencl
void sn_diffusive_boundary_flux_equal(
    schnaps_real *x,
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm, 
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 schnaps_real flux_out = 0;
 schnaps_real flux_in  = 0;
 schnaps_real w_sum   = 0; 
 schnaps_real vn    = 0;
 
 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
 }

 for (int i = 0; i < ad->nb_v; i++) {
  vn = + vnorm[0] * ad->v_tab[3 * i + 0]
     + vnorm[1] * ad->v_tab[3 * i + 1]
     + vnorm[2] * ad->v_tab[3 * i + 2];
    
  /*Colision Case */
  if (vn > 0 ) {
   /*Compute outflux flux*/
   flux_out += vn * wL[i];
   flux[i] += vn * wL[i];
  } else w_sum += 1; // Compute weights
 }

 for (int k = 0; k < ad->nb_v; k++) {
   vn = + vnorm[0] * ad->v_tab[3 * k + 0]
      + vnorm[1] * ad->v_tab[3 * k + 1]
      + vnorm[2] * ad->v_tab[3 * k + 2];

   if (vn < 0) {
    flux_in += flux_out / w_sum;
    flux[k] = - flux_out / w_sum;
   }
 }
}
#pragma end_opencl

#pragma start_opencl
void sn_diffusive_boundary_flux_lebedev(
    schnaps_real *x,
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm, 
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 schnaps_real flux_out = 0;
 schnaps_real flux_in  = 0;
 schnaps_real w_sum   = 0; 
 schnaps_real vn    = 0;
 
 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
 }

 for (int i = 0; i < ad->nb_v; i++) {
  vn = + vnorm[0] * ad->v_tab[3 * i + 0]
     + vnorm[1] * ad->v_tab[3 * i + 1]
     + vnorm[2] * ad->v_tab[3 * i + 2];
    
  /*Colision Case */
  if (vn > 0 ) {
   /*Compute outflux flux*/
   flux_out += vn * wL[i];
   flux[i] += vn * wL[i];
  } else w_sum += ad->w_tab[i]; // Compute weights
 }

 for (int k = 0; k < ad->nb_v; k++) {
   vn = + vnorm[0] * ad->v_tab[3 * k + 0]
      + vnorm[1] * ad->v_tab[3 * k + 1]
      + vnorm[2] * ad->v_tab[3 * k + 2];

   if (vn < 0) {
    flux_in += flux_out / w_sum;
    flux[k] = - flux_out / w_sum;
   }
 }
}
#pragma end_opencl

#pragma start_opencl
void sn_mixed_boundary_flux(
    schnaps_real *x, 
    schnaps_real t, 
    schnaps_real *wL,
    schnaps_real *vnorm, 
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;

 schnaps_real flux_out = 0;
 schnaps_real flux_diff = 0;
 schnaps_real flux_spec = 0;
 schnaps_real w_sum   = 0; 
 schnaps_real vn    = 0;
 schnaps_real tol = 10e-3;
 schnaps_real absorb = 0.1;
 // schnaps_real alpha = (0.40)*(1-absorb);
 schnaps_real alpha = 0.10;
 // schnaps_real beta = (1 - alpha)*(1-absorb);
 schnaps_real beta = 1-alpha;
 int count = 0;
 int count2 = 0;
 
 for (int i = 0; i < ad->nb_v; i++) {
  flux[i] = 0;
 }

 
 for (int i = 0; i < ad->nb_v; i++) {
  vn = + vnorm[0] * ad->v_tab[3 * i + 0]
     + vnorm[1] * ad->v_tab[3 * i + 1]
     + vnorm[2] * ad->v_tab[3 * i + 2];
    
  /*Colision Case */
  if (vn > 0 ) {
   /*Compute outflux flux*/
   flux_out += vn * wL[i];
   flux[i] += vn * wL[i];

   /*Compute incoming flux*/
   /*X faces*/
   if (vnorm[0] !=0 && fabs(vnorm[1]) <= tol && fabs(vnorm[2]) <= tol) {
    flux[ad->spec_tab[0][i]] -= (1-absorb) * alpha * vn * wL[i];
    flux_spec += alpha * vn * wL[i];
   }
   /*Y faces*/
   if (fabs(vnorm[0]) <= tol && vnorm[1] !=0 && fabs(vnorm[2]) <= tol) {
    flux[ad->spec_tab[1][i]] -= (1-absorb) *alpha * vn * wL[i];
    flux_spec += alpha * vn * wL[i];
   }
   /*Z faces*/
   if (fabs(vnorm[0]) <= tol && fabs(vnorm[1]) <= tol && vnorm[2] !=0) {
    flux[ad->spec_tab[2][i]] -= (1-absorb) * alpha * vn * wL[i];
    flux_spec += alpha * vn * wL[i];
   }
  } else {
   w_sum += vn;
   // count2 += 1;// Compute weights
  }   
 }

 // printf("ratio = %f, count=%d, count2=%d, wsum = %f\n",(flux_out/flux_spec),count,count2,w_sum);
 // assert((flux_out/flux_spec) == (1.0 / alpha));
 
 // schnaps_real flux_tmp = alpha*flux_out - flux_spec;
 // printf("flux_test = %f\n",flux_tmp);
 for (int k = 0; k < ad->nb_v; k++) {
   vn = + vnorm[0] * ad->v_tab[3 * k + 0]
      + vnorm[1] * ad->v_tab[3 * k + 1]
      + vnorm[2] * ad->v_tab[3 * k + 2];

   if (vn < 0) {
    // flux_diff += - (beta*(1-absorb)*flux_out) * vn / w_sum;
    flux[k]  +=  - beta * (1.0 - absorb) * flux_out * vn / w_sum;
   }
 }
 // assert(flux_tmp/flux_diff == 1.0);
 // if (fabs(flux_tmp-flux_diff) != 0)
 // printf("diff = %.15f\n",fabs(flux_tmp - flux_diff));

 /*assert for debug*/
 // assert(fabs(flux_out - flux_diff - flux_spec) == 0);
 // if (fabs(flux_out + flux_diff + flux_spec) != 0)
 // printf("diff = %.15f\n",fabs(flux_out - flux_diff - flux_spec));

 // printf("flux_out = %f, flux_diff = %f flux_spec=%f\n",flux_out, flux_diff, flux_spec);
 
 
 
 
 
 
 
}
#pragma end_opencl


/* Diagnosis routines for SN model */


void sn_compute_macro(
    Simulation *simu,
    schnaps_real *rcp_density) 
{
 /*
	 * Compute macro quantities at fixed t=tnow on each node.
  * Macro quantities are stored as conservatives variables in the field struct.
	 */
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 
 schnaps_real density; 
 schnaps_real intensity_x; 
 schnaps_real intensity_y;
 schnaps_real intensity_z;
 schnaps_real v_g = 0;
 schnaps_real v_r = 0;
 schnaps_real wpg;
 schnaps_real det;
 schnaps_real xphy[3];
 schnaps_real xpgref[3];
 schnaps_real dtau[3][3];
 schnaps_real codtau[3][3];
  *rcp_density = 0;
 
 int npg = 0;
 int imem_density = 0;
 int imem_intensity_x = 0;
 int imem_intensity_y = 0;
 int imem_intensity_z = 0;
 int imem_intensity_norm = 0;
 int imem = 0;
 
 /* macromesh */
 for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
  field *f = simu->fd + ie;
  
  /* gauss point */
  npg = NPG(f->deg, f->raf);
  for (int ipg = 0; ipg < npg; ipg++) {
   
   /* get memory index */
   imem_density = f->varindex(f->deg, 
                 f->raf,
                 f->model.m,
                 ipg,
                 ad->nb_v);
                
   imem_intensity_x = f->varindex(f->deg, 
                   f->raf, 
                   f->model.m, 
                   ipg, 
                   ad->nb_v+1);
                  
   imem_intensity_y = f->varindex(f->deg,
                   f->raf,
                   f->model.m,
                   ipg,
                   ad->nb_v+2);
                  
   imem_intensity_z = f->varindex(f->deg,
                   f->raf,
                   f->model.m,
                   ipg,
                   ad->nb_v+3);
                  
   imem_intensity_norm = f->varindex(f->deg,
                     f->raf, 
                     f->model.m,
                     ipg, 
                     ad->nb_v+4);
   
   density = 0;
   intensity_x = 0;
   intensity_y = 0;
   intensity_z = 0; 
   
   /* compute fluid quantities through velocity integration */
   for (int i = 0; i < ad->nb_v; i++) {
    
    imem = f->varindex(f->deg,
              f->raf,
              f->model.m,
              ipg,
              i);
    
    density += 4 * M_PI * f->wn[imem] * ad->w_tab[i];
    
    intensity_x += (4 * M_PI * f->wn[imem] * ad->w_tab[i] 
            * ad->v_tab[3 * i + 0]);
           
    intensity_y += (4 * M_PI * f->wn[imem] * ad->w_tab[i] 
            * ad->v_tab[3 * i + 1]);
           
    intensity_z += (4 * M_PI * f->wn[imem] * ad->w_tab[i]
            * ad->v_tab[3 * i + 2]);
   }

   /* get GP coordinates in ref FE */
   ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
  
   /* get GE physical coordinates */
   schnaps_ref2phy(f->physnode,
           xpgref,
           0, 
           -1,
           xphy,
           dtau,
           codtau,
           NULL,
           NULL); 
          
   det = dot_product(dtau[0], codtau[0]);
      
   /* compute receptors values */
   v_r += wpg * det * sn_rcp_sphere_mask(xphy, ad->rcp_data);
   // v_g += wpg * det * gaussian_density(xphy,0);
   *rcp_density = (*rcp_density + density * wpg * det 
           * sn_rcp_sphere_mask(xphy, ad->rcp_data));
   
   /* store macro values in some field entries*/
    
   f->wn[imem_density] = density;
   f->wn[imem_intensity_x] = intensity_x;
   f->wn[imem_intensity_y] = intensity_y;
   f->wn[imem_intensity_z] = intensity_z;
   f->wn[imem_intensity_norm] = sqrt(+ intensity_x*intensity_x 
                     + intensity_y*intensity_y
                     + intensity_z*intensity_z);
  }
 }
 *rcp_density = *rcp_density/v_r;
}

#pragma start_opencl
schnaps_real sn_rcp_sphere_mask(
    const schnaps_real *x, 
    const schnaps_real rcp_data[4]) {
 schnaps_real X = x[0] - rcp_data[0];
 schnaps_real Y = x[1] - rcp_data[1];
 schnaps_real Z = x[2] - rcp_data[2];

 if (X*X + Y*Y + Z*Z <= rcp_data[3]*rcp_data[3] ) {
  return 1;
 }
 else {
  return 0;
 }
}
#pragma end_opencl


schnaps_real sn_integrate_field(
    Simulation *simu, 
    int id_field)
{
 /*
	 * integral over whole x, of a conservative variable which id is id_field.
  * 
	 */
 
 schnaps_real wpg;
 schnaps_real det;
 schnaps_real xphy[3];
 schnaps_real xpgref[3];
 schnaps_real dtau[3][3];
 schnaps_real codtau[3][3];
 schnaps_real value = 0;
 int imem;
 
 /* macromesh */
 for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
  field *f = simu->fd + ie;
  
  /* gauss points */
  for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
   schnaps_real w[f->model.m];
   
   /* conservatives variabe*/
   for (int iv = 0; iv < f->model.m; iv++) {
    imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
    w[iv] = f->wn[imem];
   }

   /* get GP coordinates in ref FE */
   ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
  
   /* get GE physical coordinates */
   schnaps_ref2phy(f->physnode,
           xpgref,
           0,
           -1,
           xphy,
           dtau,
           codtau,
           NULL,
           NULL); 
          
   det = dot_product(dtau[0], codtau[0]);
   value += w[id_field] * wpg * det;
  }
 }
 return value;
}



schnaps_real sn_l2_error_all_field(
Simulation *simu)
{
  /*
  * L2 norm of the error between the approximate and exact 
* solution provided by ImposedData function.
* 
  */
  schnaps_real wpg;
  schnaps_real det;
  schnaps_real xphy[3];
  schnaps_real xpgref[3];
  schnaps_real dtau[3][3];
  schnaps_real codtau[3][3];
  schnaps_real value = 0;

  schnaps_real error_cumul = 0;
  int imem;

  

  
  
  /* macromesh */
  for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    field *f = simu->fd + ie;
    schnaps_real error[f->model.m];
    
    for (int iv = 0; iv < f->model.m; iv++) {
      error[iv] = 0;
    }
    
    
    /* gauss points */
    for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      schnaps_real w[f->model.m];
      schnaps_real wex[f->model.m];
      
      /* conservatives variables*/
      for (int iv = 0; iv < f->model.m; iv++) {
        imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
        w[iv] = f->wn[imem];
      }

      /* get GP coordinates in ref FE */
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      
      /* get GE physical coordinates */
      schnaps_ref2phy(f->physnode,
      xpgref,
      NULL,
      -1,
      xphy,
      dtau,
      codtau,
      NULL,
      NULL); 
      
      det = dot_product(dtau[0], codtau[0]);
      f->model.ImposedData(xphy, f->tnow, wex);
      // printf("tnow = %f\n", f->tnow);
      
      /* L2 error */
      for (int iv = 0; iv < f->model.m; iv++) {
        schnaps_real diff = w[iv] - wex[iv]; 
        // printf("iv = %d f[approx] = %f , %f = f[ex]\n",iv, w[iv], wex[iv]);
        error[iv] += diff * diff * wpg * det;
      }
      // value += diff * diff * wpg * det;
      // schnaps_real diff = w[id_field] - wex[id_field];

    }
    for (int iv = 0; iv < f->model.m; iv++) {
        error_cumul += sqrt(error[iv]);
    }
  }


  return error_cumul;
}



schnaps_real sn_l1_error_field(
    Simulation *simu,
    int id_field)
{
 /*
	 * L1 norm of the error between the approximate and exact 
  * solution provided by ImposedData function.
  * 
	 */
  
 schnaps_real wpg;
 schnaps_real det;
 schnaps_real xphy[3];
 schnaps_real xpgref[3];
 schnaps_real dtau[3][3];
 schnaps_real codtau[3][3];
 schnaps_real value = 0;
 int imem;
 
 /* macromesh */
 for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
  field *f = simu->fd + ie;
  
  /* gauss points */
  for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
   schnaps_real w[f->model.m];
   schnaps_real wex[f->model.m];
   
   /* conservatives variabe*/
   for (int iv = 0; iv < f->model.m; iv++) {
    imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
    w[iv] = f->wn[imem];
   }

   /* get GP coordinates in ref FE */
   ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
  
   /* get GE physical coordinates */
   schnaps_ref2phy(f->physnode,
           xpgref,
           0,
           -1,
           xphy,
           dtau,
           codtau,
           NULL,
           NULL); 
          
   det = dot_product(dtau[0], codtau[0]);
   f->model.ImposedData(xphy, f->tnow, wex);
   
   /* L1 error */
   value += fabs(w[id_field] - wex[id_field]) * wpg * det;
  }
 }
 return value;
}



schnaps_real sn_linfinity_error_field(
    Simulation *simu, 
    int id_field)
{
 /*
	 * Loo norm of the error between the approximate and exact 
  * solution provided by ImposedData function.
  * 
	 */
  
 schnaps_real wpg;
 schnaps_real det;
 schnaps_real xphy[3];
 schnaps_real xpgref[3];
 schnaps_real dtau[3][3];
 schnaps_real codtau[3][3];
 schnaps_real value = 10e8;
 schnaps_real current_val = 0;
 int imem;
 
 /* macromesh */
 for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
  field *f = simu->fd + ie;
  
  /* gauss points */
  for (int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
   schnaps_real w[f->model.m];
   schnaps_real wex[f->model.m];
   
   /* conservatives variabe*/
   for (int iv = 0; iv < f->model.m; iv++) {
    imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
    w[iv] = f->wn[imem];
   }

   /* get GP coordinates in ref FE */
   ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
  
   /* get GE physical coordinates */
   schnaps_ref2phy(f->physnode,
           xpgref,
           0,
           -1,
           xphy,
           dtau,
           codtau,
           NULL,
           NULL); 
          
   det = dot_product(dtau[0], codtau[0]);
   f->model.ImposedData(xphy, f->tnow, wex);
   
   /* L_infinity error */
   current_val = fabs(w[id_field] - wex[id_field]) * wpg * det ;
   if (current_val > value) {
    value = current_val;
   }
  }
 return value;
 }
}


/* M1 model routines */

#pragma start_opencl
schnaps_real m1_func_to_min(
    const schnaps_real x, 
    schnaps_real r1)
{
 return 1. / tanh(x) -1. / x - r1;
};
#pragma end_opencl

#pragma start_opencl
schnaps_real m1_func_to_min_der(
    const schnaps_real x)
{
 return -1. / (sinh(x) * sinh(x)) + 1. / (x * x);
};
#pragma end_opencl

#pragma start_opencl
schnaps_real m1_get_b_newton(
    const schnaps_real r1, 
    int *iter)
{	
 /* Newton solver for vectorial Lagrange multiplier
   Research only on the positive side */
   
 schnaps_real abs_r1 = fabs(r1);
 schnaps_real xk = _NewtonInitGuess;
 schnaps_real xkp1;
 schnaps_real error = 0;
 int i = 1;
 
 /* Deal I = 0 case */
 if (abs_r1 < _NewtonTol) return 0;
 
 /* Initial guess optimization */
 if (abs_r1 > _NewtonInitGuess) {
  xk = 1. / (1. - abs_r1);
 }
 xkp1 = xk;
 while (i < _NewtonIter) {
  // xkp1 = xk - m1_func_to_min(xk, abs_r1) / m1_func_to_min_der(xk);
  xkp1 = xk - (+ 1. / tanh(xk) - 1. / xk - abs_r1) / 
        (- 1. / (sinh(xk) * sinh(xk)) + 1. / (xk * xk));
  error = fabs(xkp1 - xk);
  
  /*Minimum found*/
  if (error < _NewtonTol) {
   if (r1 > 0) {
    *iter = i;
    return xkp1;
   } else {	
    *iter = i;
    return - xkp1;
   }
  }
  xk = xkp1;
  i++;
 }
 *iter = i;
 printf("[ERROR} Newton solver failed to converge after %d iteration\n",_NewtonIter);
 printf("R r1 = %.15f, \n",r1);
 return - 1000;
};
#pragma end_opencl

#pragma start_opencl
schnaps_real m1_get_b_householder(
    const schnaps_real r1, 
    int *iter)
{	
 /* Householdersolver for vectorial Lagrange multiplier
   Research only on the positive side */
 schnaps_real f_0;
 schnaps_real f_1;
 schnaps_real f_2;
 schnaps_real f_3;
 schnaps_real t1;
 schnaps_real t2;
 schnaps_real t3;
 schnaps_real t4;
 schnaps_real t5;
 schnaps_real t6;
 schnaps_real hk;
 schnaps_real xk = _NewtonInitGuess;
 schnaps_real xkp1;
 schnaps_real error = 0;
 int i = 1;
 
 /* Deal I = 0 case */
 if (r1 <= 10e-8) {
  return 0;
 }
 
 /* Initial guess optimization */
 if (r1 > _NewtonInitGuess) {
  xk = 1. / (1. - r1);
 }
 
 xkp1 = xk;
 while (i < _NewtonIter) {
  // printf("NEwton\n");
  t1 = 0.1e1 / tanh(xk);
  t2 = t1 * t1;
  t3 = t2 * t2;
  t4 = xk * xk;
  t5 = t4 * t4;

  f_0 = t1 - 0.1e1 / xk - r1;
  f_1 = 0.1e1 - t2 + 0.1e1 / t4;
  f_2 = 0.2e1 * t2 * t1 - 0.2e1 * t1 - 0.2e1 / xk / t4;
  f_3 = 0.1e1 / t5 * (0.8e1 * t5 * t2 - 0.6e1 * t5 * t3 - 0.2e1 * t5 + 0.6e1);
  
  hk = - f_0 / f_1;
  t6 = (f_2 / f_1) * hk;
  
  /* Householder Order 1 : Newton-Raphson */
  //xkp1 = xk + hk;
  
  /* Householder Order 2 */
  xkp1 = xk + hk / (1 + 0.5 * t6);

  /* Householder Order 3 */
  // xkp1 = xk + hk * (0.1e1 + 0.5e0 * t6 ) / (0.1e1 + t6 + (0.1e1 / 0.6e1) * (f_3 / f_1) * hk *hk);
  error = fabs(xkp1 - xk);
  
  /*Minimum found*/
  if (error < _NewtonTol) {
    *iter = i;
    return xkp1;
  }
  xk = xkp1;
  i++;
 }
 printf("[ERROR] HouseHolder solver failed to converge after %d iteration\n",_NewtonIter);
 printf("R r1 = %.15f, \n",r1);
 };
#pragma end_opencl

#pragma start_opencl
schnaps_real m1_get_r2(
    const schnaps_real r1)
{
 schnaps_real b1,b2,b, cothb;
 int iter;
 volatile schnaps_real h1 = 7E-5;
 volatile schnaps_real h2 = 8E-5;
 if (fabs(r1) <= h1) {
  b = m1_get_b_newton(h2,&iter);
 } else {
  b = m1_get_b_newton(r1,&iter);
 }
 cothb = 1./tanh(b);
 return 1. - 2.*(cothb -1./b)/b;
};
#pragma end_opencl

#pragma start_opencl
schnaps_real m1_get_r2_der(
    const schnaps_real r1, 
    const schnaps_real fmid, 
    const int order, 
    schnaps_real h)
{
 volatile schnaps_real r1ph,r1mh;
 schnaps_real fph, fmh,resu;
 /*Backward Derivative*/
 
 if (1-fabs(r1)<=h) h=1E-9;
 if (fabs(r1)<=h) h=h*2;
 if (order == -1) {
  r1mh=r1-h;
  h=r1-r1mh;
  fmh = m1_get_r2(r1mh);
  resu= (fmid-fmh)/h;
 }
 
 /*Centred Derivative*/
 else if (order == 0) {
  r1mh=r1-h;
  h=r1-r1mh;
  r1ph=r1+h;
  fmh = m1_get_r2(r1mh);
  fph = m1_get_r2(r1ph);
  resu= (fph - fmh)/(2*h);
 }
 
 /*Forward Derivative*/
 else if (order ==1) {
  r1ph = r1 +h;
  h=r1ph-r1;
  fph = m1_get_r2(r1+h);
  resu= (fph-fmid)/h;
 }

 // 
 if (1-fabs(r1)<=1E-10) resu=2;
 
 if (fabs(r1)<1e-6) {
  if (r1<0) {
   if (resu <0) {
    return resu;
   }
   else{
    return -resu;
   }
  }
  if (r1>0) {
   if (resu>0) {
    return resu;
   }
   else return -resu;
  }
 }
 return resu;
};
#pragma end_opencl

#ifdef _WITH_GSL

schnaps_real m1_get_b_func_f(
   const schnaps_real x, 
   void *args)
{
 struct root_finding_1arg *p  = (struct root_finding_1arg *) args;
 schnaps_real R1 = p->R1;
 return 0.1e1 / tanh(x) - 0.1e1 / x - R1;
}

schnaps_real m1_get_b_func_df(
   const schnaps_real x,
   void *args)
{
 return 0.1e1 - (0.1e1 / tanh(x)) * (0.1e1 / tanh(x)) + 0.1e1 / (x * x);
}


void m1_get_b_func_fdf(
   schnaps_real x, 
   void *args,
   schnaps_real *f,
   schnaps_real *df)
{
 struct root_finding_1arg *p  = (struct root_finding_1arg *) args;
 schnaps_real R1 = p->R1;
 *f = 0.1e1 / tanh(x) - 0.1e1 / x - R1;
 *df = 0.1e1 - (0.1e1 / tanh(x)) * (0.1e1 / tanh(x)) + 0.1e1 / (x * x);
}

schnaps_real m1_get_b_mixed_gsl(
 const schnaps_real R1)
{
/* HouseHolder Order 2 */
 bool status = false;
 schnaps_real f_0;
 schnaps_real f_1;
 schnaps_real f_2;
 schnaps_real f_3;
 schnaps_real t1;
 schnaps_real t2;
 schnaps_real t3;
 schnaps_real t4;
 schnaps_real t5;
 schnaps_real t6;
 schnaps_real hk;
 schnaps_real step;
 schnaps_real xk = _NewtonInitGuess;
 schnaps_real xkp1;
 schnaps_real error = 0;
 int iter = 1;
 if (R1 == 0) { 
  // printf("Null R1\n");
  return 0;
 }
 /* Initial guess optimization */
 if (R1 > _NewtonInitGuess) xk = 0.1e1 / (0.1e1 - R1);

 while (iter < _NewtonIter) {
  
  /* Compute f(xk), f'(xk) and f''(xk) */
  t1 = 0.1e1 / tanh(xk);
  t2 = t1 * t1;
  t3 = t2 * t2;
  t4 = xk * xk;
  t5 = t4 * t4;

  f_0 = t1 - 0.1e1 / xk - R1;
  f_1 = 0.1e1 - t2 + 0.1e1 / t4;
  f_2 = 0.2e1 * t2 * t1 - 0.2e1 * t1 - 0.2e1 / xk / t4;
  
  hk = - f_0 / f_1;

  /* Householder Order 1 : Newton-Raphson */
  // xkp1 = xk + hk;
  
  /* Householder Order 2 */
  t6 = (f_2 / f_1) * hk;
  xkp1 = xk + hk / (0.1e1 + 0.5e0 * t6);

  /* Householder Order 3 */
  // f_3 = 0.1e1 / t5 * (0.8e1 * t5 * t2 - 0.6e1 * t5 * t3 - 0.2e1 * t5 + 0.6e1);
  // xkp1 = xk + hk * (0.1e1 + 0.5e0 * t6 ) / (0.1e1 + t6 + (0.1e1 / 0.6e1) * (f_3 / f_1) * hk *hk);
  
  error = fabs(xkp1 - xk);
  
  /*Minimum found*/
  if (error < _NewtonTol) {
   // printf("[Info] Newton converged \n");
   return xkp1;
  }
  xk = xkp1;
  iter++;
 }
 // printf("[Warning] HouseHolder solver failed to converge after %d iteration for R1 = %.15f\n", _NewtonIter, R1);
 
 /* Call GSL Brendt Algorithm */
 xkp1 = m1_get_b_brendt_gsl(R1);
 return xkp1;
}



schnaps_real m1_get_b_brendt_gsl(
   const schnaps_real R1)
{
 int status = GSL_CONTINUE;
 int iter = 0, max_iter = 100;
 schnaps_real r = 0;
 schnaps_real x_lo = 1e-10;
 schnaps_real x_hi = 1000;
 schnaps_real tol = 1e-12;
 /* GSL Solver */
 const gsl_root_fsolver_type *T;
 struct root_finding_1arg args = {R1};

 gsl_root_fsolver *s;
 gsl_function F;
 F.function = &m1_get_b_func_f;
 F.params = &args;
 T = gsl_root_fsolver_brent;
 s = gsl_root_fsolver_alloc (T);

 gsl_root_fsolver_set (s, &F, x_lo, x_hi);

 while (status == GSL_CONTINUE && iter < max_iter) {
  iter++;
  status = gsl_root_fsolver_iterate (s);
  r = gsl_root_fsolver_root (s);
  x_lo = gsl_root_fsolver_x_lower (s);
  x_hi = gsl_root_fsolver_x_upper (s);
  status = gsl_root_test_interval (x_lo, x_hi,0,tol);
 }
 // if (iter > 15) printf("Brent GSL iter = %d\n",iter);
 gsl_root_fsolver_free (s);
 return r;	
};



schnaps_real m1_get_b_newton_gsl(
   const schnaps_real R1)
{
 int status = GSL_CONTINUE;
 int iter = 0, max_iter = 100;
 schnaps_real xk = 0.4, xkp1;
 schnaps_real tol = 1e-12;
 struct root_finding_1arg args = {R1};
 
 /*GSL Solver Newton*/
 const gsl_root_fdfsolver_type *T;
 gsl_root_fdfsolver *s;
 gsl_function_fdf FDF;
 FDF.f = &m1_get_b_func_f;
 FDF.df = &m1_get_b_func_df;
 FDF.fdf = &m1_get_b_func_fdf;
 FDF.params = &args;

 // T = gsl_root_fdfsolver_newton;
 T = gsl_root_fdfsolver_steffenson;
 s = gsl_root_fdfsolver_alloc (T);
 gsl_root_fdfsolver_set (s, &FDF, xk);

 // printf ("using %s method\n",
     // gsl_root_fdfsolver_name (s));

 // printf ("%-5s %10s %10s %10s\n",
     // "iter", "root", "err", "err(est)");
     
 while (status == GSL_CONTINUE && iter < max_iter) {
   iter++;
   status = gsl_root_fdfsolver_iterate(s);
   xk = xkp1;
   xkp1 = gsl_root_fdfsolver_root(s);
   status = gsl_root_test_delta(xkp1, xk, 0, tol);

   // if (status == GSL_SUCCESS)
   // printf ("%5d %10.7f %+10.7f %10.7f\n",
       // iter, x, x - r_expected, x - x0);
  }
 printf("Newton-Raphson GSL iter = %d\n", iter);
 gsl_root_fdfsolver_free (s);
 return xkp1;
}



schnaps_real m1_get_dr2_gsl(
   const schnaps_real R1,
   const int type,
   const schnaps_real h)
{
 gsl_function F;
 schnaps_real result, abserr;
 F.function = &m1_get_r2_gsl;
 F.params = 0;
 //Backward
 if (type == -1) {
  gsl_deriv_backward(&F, R1, h, &result, &abserr);
  return result;
 }
 //Centered
 if (type == 0) {
  gsl_deriv_central(&F, R1, h, &result, &abserr);
  return result;
 }
 // Forward
 if (type == 1) {
  gsl_deriv_forward(&F, R1,h, &result, &abserr);
  return result;
 }
};



schnaps_real m1_get_r2_gsl(
   const schnaps_real R1,
   void * params)
{
 (void)(params);
 return (schnaps_real) m1_get_r2(R1);
};



int m1_get_eigen(
   const schnaps_real R1,
   schnaps_real *R2,
   schnaps_real *VP1,
   schnaps_real *VP2)
{
 schnaps_real DR2, sqdelta_2, R2temp,err;
 R2temp = m1_get_r2(R1);
 DR2 = m1_get_r2_der(R1,R2temp,0,_DERIVATIVE_STEP);
 // DR2 = m1_get_r2_derGSL(R1,0, 1e-4);
 sqdelta_2 = sqrt(DR2*DR2 - 4*R1*DR2 + 4*R2temp)/2;
 *R2 = R2temp;
 DR2=DR2/2;
 *VP1=DR2 - sqdelta_2;
 *VP2=DR2 + sqdelta_2;
 return 0;
};



schnaps_real m1_get_eigen_i(
   const int i12, 
   const schnaps_real R1)
{
 schnaps_real DR2,sqdelta_2, R2,err;
 R2 = m1_get_r2(R1);
 DR2=m1_get_r2_der(R1,R2,0,_DERIVATIVE_STEP);
 // DR2 = m1_get_r2_der_gsl(R1,0, 1e-4);
 sqdelta_2 = 0.5*sqrt(DR2*DR2 - 4*R1*DR2 + 4*R2);
 DR2=0.5*DR2;
 if (i12==1) return DR2-sqdelta_2;
 if (i12==2) return DR2+sqdelta_2;
 // In case there was an error and i12 is not 1 nor 2:
 perror("In m1_get_eigen_i : i12 should always be 1 or 2.");
 return -1;
}

int m1_exact_riemann_solver_gsl(
   const schnaps_real xi,
   const schnaps_real uL[2],
   const schnaps_real uR[2],
   schnaps_real u[2])
{

 schnaps_real R1, R1L, R1R, R2, R1M, uM[2];
 schnaps_real VP1L, VP2L, VP1R, VP2R, VP1M, VP2M, vp1m, vp2m, vp1p, vp2p;
 schnaps_real X1,X2,Rf1,Rf2;
 schnaps_real tol = 1e-5;
  
 R1L = uL[1] / uL[0];
 R1R = uR[1] / uR[0];
 
 /* Get middle state uM */
 R1M = m1_riemann_get_mid_state_gsl(uL,uR,tol);
 uM[0] = uL[0] * m1_riemann_li(1, R1L, R1M);

 /*Compute caracteristics speed*/
 if (m1_get_eigen(R1L,&R2,&VP1L,&VP2L) != 0) return -1;
 if (m1_get_eigen(R1R,&R2,&VP1R,&VP2R) != 0) return -1;
 if (m1_get_eigen(R1M,&R2,&VP1M,&VP2M) != 0) return -1;
  

 if (R1M < R1L) {
  /* 1-shock entropic*/
  X1 = m1_riemann_shock(1, R1L, R1M);
  /* Compute shock's speed*/
  vp1m = (X1 * R1M - R1L) / (X1 - 1);
  vp1p = vp1m;
 } else {
 /* 1-Rarefaction case */
  vp1m = VP1L;
  vp1p = VP1M;
 }
 
 if (R1M > R1R) {
  /* 2-shock entropic*/
  X2 = m1_riemann_shock(2, R1M, R1R);
  vp2m = (X2 * R1R - R1M) / (X2 - 1);
  /* Compute shock's speed*/
  vp2p = vp2m;
 } else { 
 /* 2-Rarefaction case */
  vp2m = VP2M;
  vp2p = VP2R;
 }
 
 
 /* Solution computation */

 /* Left state */
 if (xi < vp1m) {
  u[0] = uL[0];
  u[1] = uL[1];
  return 0;
 }
 
 /* 1-Rarefaction connecting left and middle states */
 if ((xi >= vp1m) && (xi <= vp1p) ) {
  /* Get intersection point with the caracteristic in the rarefaction wave*/
  R1 = m1_riemann_get_r1_gsl(1, xi, tol);
  u[0] = uL[0] * m1_riemann_rarefaction(1, R1L, R1);
  u[1] = R1 * u[0];
  return 0;
 }
 
 /* Middle state */
 if ((xi >= vp1p) && (xi < vp2m) ) {
  u[0] = uM[0];
  u[1] = u[0] * R1M;
  return 0;
 }
 
 /* 2-Rarefaction connecting middle and right states */
 if ((xi >= vp2m) && (xi <= vp2p)) {
  /* Get intersection point with the caracteristic in the rarefaction wave*/
  R1 = m1_riemann_get_r1_gsl(2, xi, tol);
  u[0] = uR[0] * m1_riemann_rarefaction(2, R1R, R1);
  u[1] = R1 * u[0];
  return 0;
 }
 
 /* Right state */
 if (xi > vp2p) {
  u[0] = uR[0];
  u[1] = uR[1];
  return 0;
 }
 
}



schnaps_real m1_riemann_get_mid_state_gsl(
   const schnaps_real uL[2],
   const schnaps_real uR[2],
   const schnaps_real tol) 
{
 int status = GSL_CONTINUE;
 int iter = 0;
 int max_iter = 100;
 schnaps_real x_lo;
 schnaps_real x_hi;
 schnaps_real r = 0;
 schnaps_real dx = 1e-5;
 schnaps_real R1L;
 schnaps_real R1R;
 R1L = uL[1] / uL[0];
 R1R = uR[1] / uR[0];
 x_lo = -1 + dx;
 x_hi =1 - dx;

 const gsl_root_fsolver_type *T;
 struct func_args2 args = {uL[0], uR[0], R1L, R1R};

 gsl_root_fsolver *s;
 gsl_function F;
 F.function = &m1_riemann_mid_state_to_min_gsl;
 F.params = &args;
 T = gsl_root_fsolver_brent;
 s = gsl_root_fsolver_alloc(T);

 gsl_root_fsolver_set(s, &F, x_lo, x_hi);

 while (status == GSL_CONTINUE && iter < max_iter) {
  iter++;
  status = gsl_root_fsolver_iterate(s);
  r = gsl_root_fsolver_root(s);
  x_lo = gsl_root_fsolver_x_lower(s);
  x_hi = gsl_root_fsolver_x_upper(s);
  status = gsl_root_test_interval(x_lo, x_hi, 0, tol);
 }
 gsl_root_fsolver_free(s);
 return r;	
};



schnaps_real m1_riemann_mid_state_to_min_gsl(
   const schnaps_real x,
   void *args)
{
 struct func_args2 *p  = (struct func_args2 *) args;
 schnaps_real rhoL = p->argD1;
 schnaps_real rhoR = p->argD2;
 schnaps_real R1L = p->argD3;
 schnaps_real R1R = p->argD4;
 return (rhoR/rhoL) - m1_riemann_li(1, R1L, x) * m1_riemann_li(2, x, R1R);
}



schnaps_real m1_riemann_li(
   const int i12,
   const schnaps_real R1m,
   const schnaps_real R1p)
{
 if (R1p < R1m) return m1_riemann_shock(i12, R1m, R1p);
 else return m1_riemann_rarefaction(i12, R1m, R1p) ;
}
 

 
schnaps_real m1_riemann_shock(
   const int i12, 
   const schnaps_real R1m,
   const schnaps_real R1p)
{
 schnaps_real R2p;
 schnaps_real R2m;
 schnaps_real a;
 schnaps_real b;
 schnaps_real c;
 schnaps_real d;
 R2p = m1_get_r2(R1p);
 R2m = m1_get_r2(R1m);
 a = R2p - R1p * R1p;
 b = - R2p - R2m + 2 * R1p * R1m;
 c = R2m - R1m * R1m;
 d = b * b - 4 * a * c;
 if (i12 == 1) return ( - b + sqrt(d)) / (2 * a);
 if (i12 == 2) return ( - b - sqrt(d)) / (2 * a);
}
 
 
 
schnaps_real m1_riemann_rarefaction(
   const int i12,
   const schnaps_real R1m,
   const schnaps_real R1p)
{
 return exp(m1_riemann_get_invariant_gsl(i12, R1m, R1p));
}



schnaps_real m1_riemann_get_r1_gsl(
   const int i12, 
   const schnaps_real xi,
   const schnaps_real tol) 
{
 int status = GSL_CONTINUE;
 int iter = 0;
 int max_iter = 100;
 schnaps_real x_lo;
 schnaps_real x_hi;
 schnaps_real r = 0;
 schnaps_real dx = 1e-4;
 schnaps_real lower;
 schnaps_real upper;
 schnaps_real rootVP = 0.66;

 /* Setup intervals for root searching */
 if (i12 == 1) {
  if (xi > 0) {
   x_lo = 0.65;
   x_hi = 1 - dx;
  } else {
   x_lo = -1 + dx;
   x_hi = 0.69;
  }
 }
 if (i12 == 2) {
  if (xi > 0) {
   x_lo = - 0.68;
   x_hi = 1 - dx;
  } else {
   x_lo = -1 + dx;
   x_hi = -0.65;
  }
 }
  
 struct func_args args = {i12, xi};
 const gsl_root_fsolver_type *T;
 gsl_root_fsolver *s;
 gsl_function F;
 F.function = &m1_riemann_to_min_gsl;
 F.params = &args;
 T = gsl_root_fsolver_brent;
 s = gsl_root_fsolver_alloc(T);
 gsl_root_fsolver_set(s, &F, x_lo, x_hi);
 while (status == GSL_CONTINUE && iter < max_iter) {
  iter++;
  status = gsl_root_fsolver_iterate(s);
  r = gsl_root_fsolver_root(s);
  x_lo = gsl_root_fsolver_x_lower(s);
  x_hi = gsl_root_fsolver_x_upper(s);
  status = gsl_root_test_interval(x_lo, x_hi, 0, tol);
 }
 gsl_root_fsolver_free(s); 
 return r;	
};
 
 
 
schnaps_real m1_riemann_invariant_to_int_gsl(
   const schnaps_real R1,
   void *args) 
{
 int i12 = *(int *)args;
 return 1.0 / (m1_get_eigen_i(i12, R1) - R1);
}



schnaps_real m1_riemann_get_invariant_gsl(
   int i12, 
   const schnaps_real R1m,
   const schnaps_real R1p) 
{
 gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
 schnaps_real result, error;

 gsl_function F;
 F.function = &m1_riemann_invariant_to_int_gsl;
 F.params = &i12;
 size_t nevals;
 gsl_integration_qags (&F, R1m, R1p, 1e-5, 1e-5, 1000,w, &result, &error); 
 // gsl_integration_qng(&F, R1m, R1p, 1e-5, 1e-5, &result, &error, &nevals); 
 // printf ("result     = % .18f\n", result);
 // printf ("estimated error = % .18f\n", error);
 // printf ("intervals    = %u\n", w->size);
 gsl_integration_workspace_free(w);
 return result;
};



schnaps_real m1_riemann_to_min_gsl(const schnaps_real x, void *args) {
 struct func_args *p = (struct func_args *) args;
 int i12 = p->argI1;
 schnaps_real xi = p->argD1;
 return xi - m1_get_eigen_i(i12, x);
}
#endif
 


#pragma start_foo
void m1_imposed_gaussian_density(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 schnaps_real w_src_watt = 0.01;
 schnaps_real dt_mc = 0.1;
 /*Case of isotrop intensity*/

 w[0] = gaussian_density(x,t) + 10;
 w[1] = 0;
 w[2] = 0;
 w[3] = 0;
}



void m1_imposed_test_average(
    const schnaps_real *x,
    const schnaps_real t,
    schnaps_real *w)
{ 
 w[0] = 0;
 schnaps_real a = 0;
 schnaps_real b = 0.125;
 int j = 0;
 
 if (x[0] >= 0.5 )
   w[0] = 1;
 
 // if (x[0] >= 3*0.125 && x[0] <= 4*0.125 )
   // w[0] = 2;
  
 // if (x[0] >= 6*0.125 && x[0] <= 7*0.125 )
   // w[0] = 3;
  
 // for (int i = 0; i < 5; i++) {
  // j=2*i+2;

  // if (x[0] >= j*0.125 && x[0] <= (j+1)*0.125 ) {
   
   // printf("x =(%f,%f,%f)\n",x[0],x[1],x[2]);
   // printf("j= %d,  a=%f, b=%f\n",j, j*b, (j+1)*b);
  // }
 // }
}



void m1_init_test_average(
    schnaps_real *x,
    schnaps_real *w)
{
 m1_imposed_test_average(x, 0, w);
}


void m1_init_gaussian_density(
    schnaps_real *x, 
    schnaps_real *w)
{
 m1_imposed_gaussian_density(x, 0, w);
}



#ifdef _WITH_GSL
void m1_imposed_riemann_density(
     const schnaps_real *x,
     const schnaps_real t,
     schnaps_real *w)
{
 
 schnaps_real uL[2] = {0.7,0.4};
 schnaps_real uR[2] = {0.8,0.3};

 
 if (t ==0 ) {
  if (x[0] < 0) {
   w[0] = uL[0];
   w[1] = uL[1];
   w[2] = 0;
   w[3] = 0;
   } else { 
    w[0] = uR[0];
    w[1] = uR[1];
    w[2] = 0;
    w[3] = 0;
   }
   return;
 } else {
  schnaps_real u[2] = {0, 0};
  m1_exact_riemann_solver_gsl(x[0] / t, uL, uR, u);
  w[0] = u[0];
  w[1] = u[1];
  w[2] = 0;
  w[3] = 0;
 }
 return;
}
#endif


void m1_discrete_compute_fk_jac(
    schnaps_real xk[3], 
    schnaps_real w[4],
    schnaps_real fk[3],
    schnaps_real jac[9])
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 
 /* gk entries */
 schnaps_real Ix = 0.0;
 schnaps_real Iy = 0.0;
 schnaps_real Iz = 0.0;
 schnaps_real wd = 0.0;

 /* Jacobian entries */
 schnaps_real DIxDbx = 0.0;
 schnaps_real DIxDby = 0.0;
 schnaps_real DIxDbz = 0.0;
 schnaps_real DIyDbx = 0.0;
 schnaps_real DIyDby = 0.0;
 schnaps_real DIyDbz = 0.0;
 schnaps_real DIzDbx = 0.0;
 schnaps_real DIzDby = 0.0;
 schnaps_real DIzDbz = 0.0;
 
 /* Common elements */
 schnaps_real fiwi;
 schnaps_real fiwivx;
 schnaps_real fiwivy;
 schnaps_real fiwivz;
 for (int i = 0; i < ad->nb_v; i++) {
            
  fiwi = ad->w_tab[i] * exp(+ xk[0] * ad->v_tab[3 * i + 0] 
                + xk[1] * ad->v_tab[3 * i + 1]
                + xk[2] * ad->v_tab[3 * i + 2]);

  fiwivx = fiwi * ad->v_tab[3 * i + 0];
  fiwivy = fiwi * ad->v_tab[3 * i + 1];
  fiwivz = fiwi * ad->v_tab[3 * i + 2];

  wd += fiwi;
  Ix += fiwivx;
  Iy += fiwivy;
  Iz += fiwivz;

  DIxDbx += fiwivx * ad->v_tab[3 * i + 0];
  DIxDby += fiwivx * ad->v_tab[3 * i + 1];
  DIxDbz += fiwivx * ad->v_tab[3 * i + 2];

  DIyDbx += fiwivy * ad->v_tab[3 * i + 0];
  DIyDby += fiwivy * ad->v_tab[3 * i + 1];
  DIyDbz += fiwivy * ad->v_tab[3 * i + 2];

  DIzDbx += fiwivz * ad->v_tab[3 * i + 0];
  DIzDby += fiwivz * ad->v_tab[3 * i + 1];
  DIzDbz += fiwivz * ad->v_tab[3 * i + 2];
 }
  
 jac[0] = (DIxDbx * wd - Ix * Ix) / (wd * wd);
 jac[1] = (DIxDby * wd - Ix * Iy) / (wd * wd);
 jac[2] = (DIxDbz * wd - Ix * Iz) / (wd * wd);

 jac[3] = (DIyDbx * wd - Iy * Ix) / (wd * wd);
 jac[4] = (DIyDby * wd - Iy * Iy) / (wd * wd);
 jac[5] = (DIyDbz * wd - Iy * Iz) / (wd * wd);

 jac[6] = (DIzDbx * wd - Iz * Ix) / (wd * wd);
 jac[7] = (DIzDby * wd - Iz * Iy) / (wd * wd);
 jac[8] = (DIzDbz * wd - Iz * Iz) / (wd * wd);
 
 
 
 
 // printf("Jacobian\n");
 // printf("[ %f,  %f  %f ]\n",jac[0],jac[1],jac[2]);
 // printf("[ %f,  %f  %f ]\n",jac[3],jac[4],jac[5]);
 // printf("[ %f,  %f  %f ]\n",jac[6],jac[7],jac[8]);
  
 fk[0] = - (Ix/wd - w[1]/w[0]);
 fk[1] = - (Iy/wd - w[2]/w[0]);
 fk[2] = - (Iz/wd - w[3]/w[0]);
 
 // printf("Fk vector\n");
 // printf("[ %f,  %f  %f ]\n",fk[0],fk[1],fk[2]);
}


void m1_discrete_newton(
    schnaps_real w[4],
    schnaps_real xk[3])
{
 schnaps_real norm_fk;
 schnaps_real jac[9];
 schnaps_real fk[3];
 schnaps_real delta_x[3];
 int i = 0;
 int sigma[3] = {0, 1, 2};
 
 while(i < _NewtonIter) {
  m1_discrete_compute_fk_jac(xk, w, fk, jac);
  PLU_Square(3, jac, sigma);
  PLU_Solve(3, jac, sigma, fk, delta_x);
  norm_fk = sqrt(+ fk[0] * fk[0]
          + fk[1] * fk[1]
          + fk[2] * fk[2]);
  // printf("Fk norm\n");
  // printf("||fk|| = %.15f \n",norm_fk);
  xk[0] += delta_x[0]; 
  xk[1] += delta_x[1]; 
  xk[2] += delta_x[2];
  if (norm_fk <= _NewtonTol) {
   // printf("Newton iter : %d\n",i);
   return;
  }
  i=i+1;
  // printf("b vector\n");
  // printf("[ %f,  %f  %f ]\n",xk[0],xk[1],xk[2]);
 }
}



void m1_num_flux_kinetic(
    schnaps_real * wL,
    schnaps_real * wR,
    schnaps_real * vnorm,
    schnaps_real * flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 /* 4*Pi Simplification done */
 schnaps_real vin;
 schnaps_real IL_norm;
 schnaps_real IR_norm;
 schnaps_real bL_norm;
 schnaps_real bR_norm;
 schnaps_real bLx; 
 schnaps_real bLy; 
 schnaps_real bLz; 
 schnaps_real bRx; 
 schnaps_real bRy; 
 schnaps_real bRz; 
 schnaps_real aL;
 schnaps_real aR;
 schnaps_real fiR;
 schnaps_real fiL;
 schnaps_real vix;
 schnaps_real viy;
 schnaps_real viz;
 schnaps_real wi;
 schnaps_real wifivin;
 schnaps_real r1_tol = 1e-8; /*only works at 1e-5*/
 int iter = 0;
 
 flux[0] = 0;
 flux[1] = 0;
 flux[2] = 0;
 flux[3] = 0;

 /* Left cell */
 IL_norm = sqrt(wL[1] * wL[1] + wL[2] * wL[2] + wL[3] * wL[3]);
 
 // bL_norm = m1_get_b_newton(IL_norm / wL[0], &iter);
 // bL_norm = m1_get_b_brendt_gsl(IL_norm / wL[0]);
 
#ifdef _WITH_GSL
 bL_norm = m1_get_b_mixed_gsl(IL_norm / wL[0]);
#else
 bL_norm = m1_get_b_newton(IL_norm / wL[0], &iter);
#endif
 
 
 if (bL_norm <= 10e-12) {
  bLx = bL_norm;
  bLy = bL_norm;
  bLz = bL_norm;
  schnaps_real t1 = bL_norm * bL_norm;
  /* Use 5-th Order taylor for a*/
  aL = wL[0] * (0.1e1 - (0.1e1 / 0.6e1) * bL_norm * bL_norm
    + (0.7e1 / 360e0) * t1 * t1);
 } else {
  bLx = bL_norm * wL[1] / IL_norm;
  bLy = bL_norm * wL[2] / IL_norm;
  bLz = bL_norm * wL[3] / IL_norm;
  aL = wL[0] * bL_norm / sinh(bL_norm);
 }
 
 /* Right cell */
 IR_norm = sqrt(wR[1] * wR[1] + wR[2] * wR[2] + wR[3] * wR[3]);
 // bR_norm = m1_get_b_newton(IR_norm / wR[0], &iter);
 // bR_norm = m1_get_b_brendt_gsl(IR_norm / wR[0]);
 
#ifdef _WITH_GSL
 bR_norm = m1_get_b_mixed_gsl(IR_norm / wR[0]);
#else
 bR_norm = m1_get_b_newton(IR_norm / wR[0], &iter);
#endif

 
 if (bR_norm <= 10e-12) {
  bRx = bR_norm;
  bRy = bR_norm;
  bRz = bR_norm;
  /* Use 5-th Order taylor for a*/
  schnaps_real t1 = bR_norm * bR_norm;
  aR = wR[0] * (0.1e1 - (0.1e1 / 0.6e1) * bR_norm * bR_norm
    + (0.7e1 / 360e0) * t1 * t1) ;
 } else {
  bRx = bR_norm * wR[1] / IR_norm;
  bRy = bR_norm * wR[2] / IR_norm;
  bRz = bR_norm * wR[3] / IR_norm;
  aR = wR[0] * bR_norm / sinh(bR_norm);
 }

 /* Compute numerical integrals for kinetic flux */ 
 for (int i = 0; i < ad->nb_v; i++) {
  
   vix = ad->v_tab[3 * i];
   viy = ad->v_tab[3 * i + 1];
   viz = ad->v_tab[3 * i + 2];
   wi = ad->w_tab[i];
  
   vin = vix * vnorm[0] + viy * vnorm[1] + viz * vnorm[2];
  /* L->R use left reconstructed pdf */
  if (vin >= 0) {
   fiL = aL * exp(bLx * vix + bLy * viy + bLz * viz);
   wifivin = wi * fiL * vin; 
     
   flux[0] += wifivin; 
   flux[1] += wifivin * vix; 
   flux[2] += wifivin * viy; 
   flux[3] += wifivin * viz; 
   
  /* R->L use right reconstructed pdf */
  } else {
   fiR = aR * exp(bRx * vix + bRy * viy + bRz * viz);
   wifivin = wi * fiR * vin; 
   
   flux[0] += wifivin; 
   flux[1] += wifivin * vix; 
   flux[2] += wifivin * viy; 
   flux[3] += wifivin * viz; 
  }
 }
}


void m1_num_flux_rusanov(
    schnaps_real * wL,
    schnaps_real * wR,
    schnaps_real * vnorm,
    schnaps_real * flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 /* 4*Pi Simplification done */
 schnaps_real IL_norm;
 schnaps_real IR_norm;
 schnaps_real bL_norm;
 schnaps_real bR_norm;
 schnaps_real RL;
 schnaps_real RR;
 schnaps_real r1_tol = 1e-8; /*only works at 1e-5*/
 schnaps_real FL[3];
 schnaps_real FR[3]; 
 int iter = 0;
 
 flux[0] = 0;
 flux[1] = 0;
 flux[2] = 0;
 flux[3] = 0;

 /* Left cell */
 IL_norm = sqrt(wL[1] * wL[1] + wL[2] * wL[2] + wL[3] * wL[3]);
 
 // bL_norm = m1_get_b_newton(IL_norm / wL[0], &iter);
 // bL_norm = m1_get_b_brendt_gsl(IL_norm / wL[0]);
 //printf("wL = %f,IL= %f\n",wL[0], IL_norm);
#ifdef _WITH_GSL
 bL_norm = m1_get_b_mixed_gsl(IL_norm / wL[0]) + 1.e-10;
#else
 bL_norm = m1_get_b_newton(IL_norm / wL[0], &iter);
#endif
 
 RL = IL_norm / (wL[0] + 1.e-10);
 IL_norm = IL_norm + 1.e-10;
 
 //printf("wR = %f,IR= %f\n",wR[0], IR_norm);
 /* Right cell */
 IR_norm = sqrt(wR[1] * wR[1] + wR[2] * wR[2] + wR[3] * wR[3]);
 // bR_norm = m1_get_b_newton(IR_norm / wR[0], &iter);
 // bR_norm = m1_get_b_brendt_gsl(IR_norm / wR[0]);
 
#ifdef _WITH_GSL
 bR_norm = m1_get_b_mixed_gsl(IR_norm / wR[0]) + 1.e-10;
#else
 bR_norm = m1_get_b_newton(IR_norm / wR[0], &iter);
#endif

 RR = IR_norm / (wR[0] + 1.e-10);
 IR_norm = IR_norm + 1.e-10;

 flux[0] = 0.5 * (wL[1] + wR[1] - (wR[0] - wL[0])) * vnorm[0] 
		  + 0.5 * (wL[2] + wR[2] - (wR[0] - wL[0])) * vnorm[1] 
		  + 0.5 * (wL[3] + wR[3] - (wR[0] - wL[0])) * vnorm[2]; 
 
 FL[0] = ((1-3*RL/bL_norm) * (wL[1]*wL[1])/(IL_norm*IL_norm) + RL/bL_norm) * wL[0];
 FL[1] =  (1-3*RL/bL_norm) * (wL[2]*wL[1])/(IL_norm*IL_norm) * wL[0];
 FL[2] =  (1-3*RL/bL_norm) * (wL[3]*wL[1])/(IL_norm*IL_norm) * wL[0];
 
 FR[0] = ((1-3*RR/bR_norm) * (wR[1]*wR[1])/(IR_norm*IR_norm) + RR/bR_norm) * wR[0];
 FR[1] =  (1-3*RR/bR_norm) * (wR[2]*wR[1])/(IR_norm*IR_norm) * wR[0];
 FR[2] =  (1-3*RR/bR_norm) * (wR[3]*wR[1])/(IR_norm*IR_norm) * wR[0];
 
 flux[1] = 0.5 * (FL[0] + FR[0] - (wR[1] - wL[1])) * vnorm[0] 
		  + 0.5 * (FL[1] + FR[1] - (wR[1] - wL[1])) * vnorm[1] 
		  + 0.5 * (FL[2] + FR[2] - (wR[1] - wL[1])) * vnorm[2]; 
		  
 FL[0] =  (1-3*RL/bL_norm) * (wL[1]*wL[2])/(IL_norm*IL_norm) * wL[0];
 FL[1] = ((1-3*RL/bL_norm) * (wL[2]*wL[2])/(IL_norm*IL_norm) + RL/bL_norm ) * wL[0];
 FL[2] =  (1-3*RL/bL_norm) * (wL[3]*wL[2])/(IL_norm*IL_norm) * wL[0];
 
 FR[0] =  (1-3*RR/bR_norm) * (wR[1]*wR[2])/(IR_norm*IR_norm) * wR[0];
 FR[1] = ((1-3*RR/bR_norm) * (wR[2]*wR[2])/(IR_norm*IR_norm) + RR/bR_norm) * wR[0];
 FR[2] =  (1-3*RR/bR_norm) * (wR[3]*wR[2])/(IR_norm*IR_norm) * wR[0];
 
 flux[2] = 0.5 * (FL[0] + FR[0] - (wR[2] - wL[2])) * vnorm[0] 
		  + 0.5 * (FL[1] + FR[1] - (wR[2] - wL[2])) * vnorm[1] 
		  + 0.5 * (FL[2] + FR[2] - (wR[2] - wL[2])) * vnorm[2]; 
		  
 FL[0] =  (1-3*RL/bL_norm) * (wL[1]*wL[3])/(IL_norm*IL_norm) * wL[0];
 FL[1] =  (1-3*RL/bL_norm) * (wL[2]*wL[3])/(IL_norm*IL_norm) * wL[0];
 FL[2] =  ((1-3*RL/bL_norm) * (wL[3]*wL[3])/(IL_norm*IL_norm) + RL/bL_norm ) * wL[0];
 
 FR[0] =  (1-3*RR/bR_norm) * (wR[1]*wR[3])/(IR_norm*IR_norm) * wR[0];
 FR[1] =  (1-3*RR/bR_norm) * (wR[2]*wR[3])/(IR_norm*IR_norm) * wR[0];
 FR[2] =  ((1-3*RR/bR_norm) * (wR[3]*wR[3])/(IR_norm*IR_norm) + RR/bR_norm ) * wR[0];
 
 flux[3] = 0.5 * (FL[0] + FR[0] - (wR[3] - wL[3])) * vnorm[0] 
		  + 0.5 * (FL[1] + FR[1] - (wR[3] - wL[3])) * vnorm[1] 
		  + 0.5 * (FL[2] + FR[2] - (wR[3] - wL[3])) * vnorm[2]; 
 //printf("vnorm %f, %f, %f \n",vnorm[1],vnorm[2],vnorm[3]);
 //printf("flux %f, %f, %f \n", flux[1],flux[2],flux[3]);
}


void randfrom(schnaps_real min, schnaps_real max,schnaps_real *result) 
{
  schnaps_real range = (max - min); 
  schnaps_real div = RAND_MAX / range;
  *result = min + (rand() / div);
}

void randn(schnaps_real *result) 
{
	schnaps_real u1;
	schnaps_real u2;
	
  randfrom(0., 1., &u1);
  randfrom(0., 1., &u2);
  *result = sqrt(-2. * log(u1)) * cos(2. * M_PI * u2);
}


void m1_boundary_flux(
    schnaps_real *x,
    schnaps_real t,
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 schnaps_real vin;
 schnaps_real IL_norm;
 schnaps_real bL_norm;
 schnaps_real bLx; 
 schnaps_real bLy; 
 schnaps_real bLz; 
 schnaps_real bL_dot_vn;
 schnaps_real b_spec_x; 
 schnaps_real b_spec_y; 
 schnaps_real b_spec_z; 
 schnaps_real aL;
 schnaps_real fiL;
 schnaps_real fiR;
 schnaps_real vix;
 schnaps_real viy;
 schnaps_real viz;
 schnaps_real wi;
 schnaps_real rdx;
 schnaps_real rdy;
 schnaps_real rdz;
 schnaps_real rd_norm;
 schnaps_real b_rand_x; 
 schnaps_real b_rand_y; 
 schnaps_real b_rand_z; 
 schnaps_real wifivin = 0;
 schnaps_real int_vn = 0;
 schnaps_real r1_tol = 1e-5; /*only works at 1e-5*/
 schnaps_real flux_out[4] = {0, 0, 0, 0};
 schnaps_real flux_spec[4] = {0, 0, 0, 0};
 schnaps_real flux_rand[4] = {0, 0, 0, 0};
 schnaps_real flux_diff[4] = {0, 0, 0, 0};
 int iter = 0;
 
 
 // flux[0] = 0;
 // flux[1] = 0;
 // flux[2] = 0;
 // flux[3] = 0;

 
 // schnaps_real wR[4] = {0, 0, 0, 0};
 // m1_imposed_gaussian_density(x, t, wR);
 // m1_num_flux_kinetic(wL, wR, vnorm, flux);
 // return;
 
 /* Left cell */
 IL_norm = sqrt(wL[1] * wL[1] + wL[2] * wL[2] + wL[3] * wL[3]);
 
 // bL_norm = m1_get_b_newton(IL_norm / wL[0], &iter);
 // bL_norm = m1_get_b_brendt_gsl(IL_norm / wL[0]);
 
#ifdef _WITH_GSL
 bL_norm = m1_get_b_mixed_gsl(IL_norm / wL[0]);
#else
 bL_norm = m1_get_b_newton(IL_norm / wL[0], &iter);
#endif
 

 if (bL_norm <= 10e-12) {
  bLx = bL_norm;
  bLy = bL_norm;
  bLz = bL_norm;
  schnaps_real t1 = bL_norm * bL_norm;
  /* Use 5-th Order taylor for a*/
  aL = wL[0] * (0.1e1 - (0.1e1 / 0.6e1) * t1
    + (0.7e1 / 360e0) * t1 * t1);
 } else {
  bLx = bL_norm * wL[1] / IL_norm;
  bLy = bL_norm * wL[2] / IL_norm;
  bLz = bL_norm * wL[3] / IL_norm;
  aL = wL[0] * bL_norm / sinh(bL_norm);
 }

 /* Specular direction : Method 1 */
 bL_dot_vn = bLx * vnorm[0] + bLy * vnorm[1] + bLz * vnorm[2]; 
 b_spec_x = bLx - 2.0 * vnorm[0] * bL_dot_vn;
 b_spec_y = bLy - 2.0 * vnorm[1] * bL_dot_vn;
 b_spec_z = bLz - 2.0 * vnorm[2] * bL_dot_vn;
 
 
 randn( &rdx);
 randn( &rdy); 
 randn( &rdz);
 rd_norm = sqrt(rdx * rdx + rdy * rdy + rdz * rdz); 
 b_rand_x = bL_norm * rdx / rd_norm; 
 b_rand_y = bL_norm * rdy / rd_norm; 
 b_rand_z = bL_norm * rdz / rd_norm; 
 
 
 if ((b_rand_x * vnorm[0] + b_rand_y * vnorm[1] + b_rand_z * vnorm[2]) > 0){
  b_rand_x = - b_rand_x; 
  b_rand_y = - b_rand_y; 
  b_rand_z = - b_rand_z; 
 } 
 
 /* Compute Numerical Integrals */
 for (int i = 0; i < ad->nb_v; i++) {
  
   vix = ad->v_tab[3 * i];
   viy = ad->v_tab[3 * i + 1];
   viz = ad->v_tab[3 * i + 2];
   wi = ad->w_tab[i];
  
   vin = vix * vnorm[0] + viy * vnorm[1] + viz * vnorm[2];
 
  /* L->R use left reconstructed pdf */
  if (vin >= 0) {
   fiL = aL * exp(bLx * vix + bLy * viy + bLz * viz);
   wifivin = wi * fiL * vin; 

   flux_out[0] += wifivin; 
   flux_out[1] += wifivin * vix; 
   flux_out[2] += wifivin * viy; 
   flux_out[3] += wifivin * viz; 
   /* Used for diffuse condition */
   int_vn += wi * vin;
   
  /* R->L use right reconstructed pdf */
  } else {
   /* Specular Part */
   fiR = aL * exp(b_spec_x * vix + b_spec_y * viy + b_spec_z * viz);
   wifivin = wi * fiR * vin; 
         
   flux_spec[0] += wifivin; 
   flux_spec[1] += wifivin * vix; 
   flux_spec[2] += wifivin * viy; 
   flux_spec[3] += wifivin * viz; 
   
   fiR = aL * exp(b_rand_x * vix + b_rand_y * viy + b_rand_z * viz);
   wifivin = wi * fiR * vin; 
         
   flux_rand[0] += wifivin; 
   flux_rand[1] += wifivin * vix; 
   flux_rand[2] += wifivin * viy; 
   flux_rand[3] += wifivin * viz;     
  }
 }
 
 /*Diffusive Part*/
 for (int i = 0; i < ad->nb_v; i++) {
      
  vix = ad->v_tab[3 * i];
  viy = ad->v_tab[3 * i + 1];
  viz = ad->v_tab[3 * i + 2];
  wi = ad->w_tab[i];
  
  vin = vix * vnorm[0] + viy * vnorm[1] + viz * vnorm[2];
  
  fiR = (flux_out[0] / int_vn);
  wifivin = wi * fiR * vin;
  
  if (vin < 0) {
   flux_diff[0] += wifivin; 
   flux_diff[1] += wifivin * vix; 
   flux_diff[2] += wifivin * viy; 
   flux_diff[3] += wifivin * viz; 
  }
 }
 
 schnaps_real beta = 0;
  
 // printf("Flux[0] out/in = (%f, %f) \n", flux_out[0], beta*flux_spec[0] + (1-beta) * flux_diff[0]);
 // printf("Flux[1] out/in = (%f, %f) \n", flux_out[1], beta*flux_spec[1] + (1-beta) * flux_diff[1]);
 // printf("Flux[2] out/in = (%f, %f) \n", flux_out[2], beta*flux_spec[2] + (1-beta) * flux_diff[2]);
 // printf("Flux[3] out/in = (%f, %f) \n", flux_out[3], beta*flux_spec[3] + (1-beta) * flux_diff[3]);
  
 flux[0] = flux_out[0] + beta * flux_spec[0] + (1-beta) * flux_diff[0];
 flux[1] = flux_out[1] + beta * flux_spec[1] + (1-beta) * flux_diff[1];
 flux[2] = flux_out[2] + beta * flux_spec[2] + (1-beta) * flux_diff[2];
 flux[3] = flux_out[3] + beta * flux_spec[3] + (1-beta) * flux_diff[3];
 
 //flux[0] = flux_out[0] + beta * flux_spec[0] + (1-beta) * flux_rand[0];
 //flux[1] = flux_out[1] + beta * flux_spec[1] + (1-beta) * flux_rand[1];
 //flux[2] = flux_out[2] + beta * flux_spec[2] + (1-beta) * flux_rand[2];
 //flux[3] = flux_out[3] + beta * flux_spec[3] + (1-beta) * flux_rand[3];
 
 return;
}




void m1_num_flux_kinetic_discrete(
    schnaps_real * wL,
    schnaps_real * wR,
    schnaps_real * vnorm,
    schnaps_real * flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 /* 4*Pi Simplification done */
 schnaps_real vn;
 schnaps_real norm_IL;
 schnaps_real norm_IR;
 schnaps_real b_norm_L;
 schnaps_real b_norm_R;
 schnaps_real bL[3] = {0, 0, 0};
 schnaps_real bR[3] = {0, 0, 0};
 schnaps_real aL;
 schnaps_real aR;
 schnaps_real fR;
 schnaps_real fL;

 schnaps_real r1_tol = 1e-5; /*only works at 1e-5*/
 
 flux[0] = 0;
 flux[1] = 0;
 flux[2] = 0;
 flux[3] = 0;

 /* Left cell */
 m1_discrete_newton(wL, bL);
 b_norm_L = sqrt(bL[0] * bL[0] + bL[1] * bL[1] + bL[2] * bL[2]);
 // printf("b_norm_L = (%f, %f, %f)\n",bL[0], bL[1], bL[2]);
 if (b_norm_L <= r1_tol) {
  aL = wL[0]; 
  bL[0] = 0;
  bL[1] = 0;
  bL[2] = 0;
 } else {
  aL = wL[0] * b_norm_L / sinh(b_norm_L);
 }
 
 /* Right cell */
 m1_discrete_newton(wR, bR);
 b_norm_R = sqrt(bR[0] * bR[0] + bR[1] * bR[1] + bR[2] * bR[2]);

 if (b_norm_R <= r1_tol) {
  aR = wR[0]; 
  bR[0] = 0;
  bR[1] = 0;
  bR[2] = 0;
 } else {
  aR = wR[0] * b_norm_R / sinh(b_norm_R);
 }

 /* Compute numerical integrals for kinetic flux */ 
 for (int i = 0; i < ad->nb_v; i++) {
   vn = + ad->v_tab[3 * i + 0] * vnorm[0]
      + ad->v_tab[3 * i + 1] * vnorm[1]
      + ad->v_tab[3 * i + 2] * vnorm[2];
 
  /* Outflow on the L->R interface. Use left reconstructed pdf */
  if (vn >= 0) {
   fL = aL * exp(+ bL[0] * ad->v_tab[3 * i + 0]
          + bL[1] * ad->v_tab[3 * i + 1]
          + bL[2] * ad->v_tab[3 * i + 2]);
           
   flux[0] += ad->w_tab[i] * fL * vn; 
   flux[1] += ad->w_tab[i] * fL * vn * ad->v_tab[3 * i + 0]; 
   flux[2] += ad->w_tab[i] * fL * vn * ad->v_tab[3 * i + 1]; 
   flux[3] += ad->w_tab[i] * fL * vn * ad->v_tab[3 * i + 2]; 
   
  /* Inflow on the L->R interface. Use right reconstruced pdf */
  } else {
   fR = aR * exp(+ bR[0] * ad->v_tab[3 * i + 0]
          + bR[1] * ad->v_tab[3 * i + 1]
          + bR[2] * ad->v_tab[3 * i + 2]);
          
   flux[0] += ad->w_tab[i] * fR * vn; 
   flux[1] += ad->w_tab[i] * fR * vn * ad->v_tab[3 * i + 0]; 
   flux[2] += ad->w_tab[i] * fR * vn * ad->v_tab[3 * i + 1]; 
   flux[3] += ad->w_tab[i] * fR * vn * ad->v_tab[3 * i + 2]; 
  }
 }
 return;
}
 
 
 
void m1_boundary_flux_discrete(
    schnaps_real *x,
    schnaps_real t,
    schnaps_real *wL,
    schnaps_real *vnorm,
    schnaps_real *flux)
{
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 schnaps_real vn;
 schnaps_real norm_IL;
 schnaps_real b_norm;
 schnaps_real bL[3] = {0, 0, 0}; 
 schnaps_real bvn;
 schnaps_real b_spec_x; 
 schnaps_real b_spec_y; 
 schnaps_real b_spec_z; 
 schnaps_real a;
 schnaps_real f;
 schnaps_real r1_tol = 1e-5; /*only works at 1e-5*/
 
 int iter = 0;
 
 schnaps_real flux_out[4] = {0, 0, 0, 0};
 schnaps_real flux_spec[4] = {0, 0, 0, 0};
 schnaps_real flux_diff[4] = {0, 0, 0, 0};
 
 flux[0] = 0;
 flux[1] = 0;
 flux[2] = 0;
 flux[3] = 0;
 
 // schnaps_real flx, fly, flz;
 // schnaps_real wR[4] = {0, 0, 0, 0};
 // m1_imposed_gaussian_density(x, t, wR);
 // m1_num_flux_kinetic(wL, wR, vnorm, flux);
 // return;
 int diff = 1;
 /* Left cell */
 m1_discrete_newton(wL, bL);
 b_norm = sqrt(bL[0] * bL[0] + bL[1] * bL[1] + bL[2] * bL[2]);
 // printf("bL_norm = (%f, %f, %f)\n",bL[0], bL[1], bL[2]);
 if (b_norm <= r1_tol) {
  a = wL[0]; 
  bL[0] = 0;
  bL[1] = 0;
  bL[2] = 0;
 } else {
  // printf("NON NULL \n");
  // diff = 0;
  a = wL[0] * b_norm / sinh(b_norm);
 }
 
 /* Specular direction */
 bvn = bL[0] * vnorm[0] + bL[1] * vnorm[1] + bL[2] * vnorm[2]; 
 b_spec_x = bL[0] - 2.0 * vnorm[0] * bvn;
 b_spec_y = bL[1] - 2.0 * vnorm[1] * bvn;
 b_spec_z = bL[2] - 2.0 * vnorm[2] * bvn;
  
 /* Compute Numerical Integrals */
  for (int i = 0; i < ad->nb_v; i++) {
   
   vn = + ad->v_tab[3 * i + 0] * vnorm[0]
      + ad->v_tab[3 * i + 1] * vnorm[1]
      + ad->v_tab[3 * i + 2] * vnorm[2];
   
   if (vn >= 0) {
    
    f = a * exp(+ bL[0] * ad->v_tab[3 * i + 0]
          + bL[1] * ad->v_tab[3 * i + 1]
          + bL[2] * ad->v_tab[3 * i + 2]);
          
    flux_out[0] += ad->w_tab[i] * f * vn; 
    flux_out[1] += ad->w_tab[i] * f * vn * ad->v_tab[3 * i + 0]; 
    flux_out[2] += ad->w_tab[i] * f * vn * ad->v_tab[3 * i + 1]; 
    flux_out[3] += ad->w_tab[i] * f * vn * ad->v_tab[3 * i + 2]; 
    
   } else {

    f = a * exp(+ b_spec_x * ad->v_tab[3 * i + 0]
          + b_spec_y * ad->v_tab[3 * i + 1]
          + b_spec_z * ad->v_tab[3 * i + 2]);

    flux_spec[0] += ad->w_tab[i] * f * vn; 
    flux_spec[1] += ad->w_tab[i] * f * vn * ad->v_tab[3 * i + 0]; 
    flux_spec[2] += ad->w_tab[i] * f * vn * ad->v_tab[3 * i + 1]; 
    flux_spec[3] += ad->w_tab[i] * f * vn * ad->v_tab[3 * i + 2]; 
    
    /*Diffusive Part*/

    flux_diff[0] += ad->w_tab[i] * wL[0] * vn; 
    flux_diff[1] += ad->w_tab[i] * wL[0] * vn * ad->v_tab[3 * i + 0]; 
    flux_diff[2] += ad->w_tab[i] * wL[0] * vn * ad->v_tab[3 * i + 1]; 
    flux_diff[3] += ad->w_tab[i] * wL[0] * vn * ad->v_tab[3 * i + 2]; 
   }
  }
  
  schnaps_real beta = 1;
  
  printf("Flux[0] out/in = (%f, %f) \n", flux_out[0], beta*flux_spec[0] + (1-beta) * flux_diff[0]);
  printf("Flux[1] out/in = (%f, %f) \n", flux_out[1], beta*flux_spec[1] + (1-beta) * flux_diff[1]);
  printf("Flux[2] out/in = (%f, %f) \n", flux_out[2], beta*flux_spec[2] + (1-beta) * flux_diff[2]);
  printf("Flux[3] out/in = (%f, %f) \n", flux_out[3], beta*flux_spec[3] + (1-beta) * flux_diff[3]);
  
  flux[0] = flux_out[0] + beta * flux_spec[0] + (1-beta) * flux_diff[0];
  flux[1] = flux_out[1] + beta * flux_spec[1] + (1-beta) * flux_diff[1];
  flux[2] = flux_out[2] + beta * flux_spec[2] + (1-beta) * flux_diff[2];
  flux[3] = flux_out[3] + beta * flux_spec[3] + (1-beta) * flux_diff[3];

 return;
} 


/* M1 diagnosis routines */


void m1_compute_rcp(
    Simulation *simu,
    schnaps_real *rcp_density) 
{
 /*
	 * Compute the spatial integral within the detector
  * 
	 */
 __constant acoustic_data *ad = &schnaps_acoustic_data;
 schnaps_real density; 
 schnaps_real intensity_x; 
 schnaps_real intensity_y;
 schnaps_real intensity_z;
 schnaps_real v_g = 0;
 schnaps_real v_r = 0;
 schnaps_real rcp_value = 0;
 schnaps_real wpg;
 schnaps_real det;
 schnaps_real xphy[3];
 schnaps_real xpgref[3];
 schnaps_real dtau[3][3];
 schnaps_real codtau[3][3];
 
 int npg = 0;
 int imem_density = 0;
 int imem_intensity_x = 0;
 int imem_intensity_y = 0;
 int imem_intensity_z = 0;
 int imem_intensity_norm = 0;
 int imem = 0;
 
 /* macromesh */
 for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
  field *f = simu->fd + ie;
  
  /* gauss point */
  npg = NPG(f->deg, f->raf);
  for (int ipg = 0; ipg < npg; ipg++) {
   
   /* get memory index */
   imem_density = f->varindex(f->deg, 
                 f->raf,
                 f->model.m,
                 ipg,
                 0);
   
   imem_intensity_x = f->varindex(f->deg, 
                   f->raf,
                   f->model.m,
                   ipg,
                   1);
                   
   imem_intensity_y = f->varindex(f->deg, 
                   f->raf,
                   f->model.m,
                   ipg,
                   2);
                   
   imem_intensity_z = f->varindex(f->deg, 
                   f->raf,
                   f->model.m,
                   ipg,
                   3);
                   
   imem_intensity_norm = f->varindex(f->deg, 
                    f->raf,
                    f->model.m,
                    ipg,
                    4);
   /* pick value in wn array */  
   density = f->wn[imem_density]; 
   intensity_x = f->wn[imem_intensity_x];
   intensity_y = f->wn[imem_intensity_y];
   intensity_z = f->wn[imem_intensity_z];
   
   /* get GP coordinates in ref FE */
   ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
  
   /* get GE physical coordinates */
   schnaps_ref2phy(f->physnode,
           xpgref,
           0, 
           -1,
           xphy,
           dtau,
           codtau,
           NULL,
           NULL); 
          
   det = dot_product(dtau[0], codtau[0]);
      
   /* compute receptors values */
   v_r += wpg * det * sn_rcp_sphere_mask(xphy, ad->rcp_data);
   // v_g += wpg * det * gaussian_density(xphy,0);
   rcp_value += (density * wpg * det 
          * sn_rcp_sphere_mask(xphy, ad->rcp_data));
   
   /* compute and store norm of intensity vector */
   f->wn[imem_intensity_norm] = sqrt(+ intensity_x*intensity_x 
                     + intensity_y*intensity_y
                     + intensity_z*intensity_z);             
  }
 }
 *rcp_density = rcp_value/v_r;
 return;
}

#pragma end_foo
#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include "maxwell.h"

#include <dlfcn.h>

// void (*NumFlux)(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux) = NULL; 
// void (*BoundaryFlux)(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vn, schnaps_real *flux) = NULL;
// void (*Source)(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source) = NULL;
// void (*InitData)(schnaps_real *x, schnaps_real *w) = NULL;
// void (*ImposedData)(const schnaps_real *x, const schnaps_real t, schnaps_real *w) = NULL;


fluxptr numflux(const char *name)
{
  if(strcmp(name, "VecTransNumFlux2d") == 0)
    return &VecTransNumFlux2d;

  if(strcmp(name, "TransNumFlux2d") == 0)
    return &TransNumFlux2d;

  if(strcmp(name, "Maxwell3DNumFluxClean_upwind") == 0)
    return &Maxwell3DNumFluxClean_upwind;

  printf("Numerical flux %s not found!\n", name);
  assert(false);
  return 0;
}

bfluxptr bflux(const char *name)
{
  if(strcmp(name, "TransBoundaryFlux2d") == 0)
    return &TransBoundaryFlux2d;

  if(strcmp(name, "Maxwell3DBoundaryFlux_upwind") == 0)
    return &Maxwell3DBoundaryFlux_upwind;

  if(strcmp(name, "Maxwell3DBoundaryFluxClean_upwind") == 0)
    return &Maxwell3DBoundaryFluxClean_upwind;

  printf("Boundary flux %s not found!\n", name);
  assert(false);
  return 0;
}

initdataptr initdata(const char *name)
{
  if(strcmp(name, "TransInitData2d") == 0)
    return &TransInitData2d;

  if(strcmp(name, "Maxwell3DInitData") == 0)
    return &Maxwell3DInitData;

  printf("Init data %s not found!\n", name);
  assert(false);
  return 0;
}

imposeddataptr imposeddata(const char *name)
{
  if(strcmp(name, "TransImposedData2d") == 0)
    return &TransImposedData2d;

  if(strcmp(name, "TestTransImposedData2d") == 0)
    return &TestTransImposedData2d;

  if(strcmp(name, "Maxwell3DImposedData") == 0)
    return &Maxwell3DImposedData;

  printf("Imposed data %s not found!\n", name);
  assert(false);
  return 0;
}

#pragma start_opencl
void TransNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real* vnorm, schnaps_real* flux)
{
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
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
void TransNumFlux2d(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux)
{
  const schnaps_real transport_v2d[] = {sqrt(0.5), sqrt(0.5), 0};
  //const real transport_v2d[] = {1,0, 0};
  schnaps_real vn
    = transport_v2d[0] * vnorm[0]
    + transport_v2d[1] * vnorm[1]
    + transport_v2d[2] * vnorm[2];
  schnaps_real vnp = vn > 0 ? vn : 0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  /* if (fabs(vnorm[2])>1e-6) { */
  /*   printf("vnds %lf %lf %lf \n", vnorm[0], vnorm[1], vnorm[2]); */
  /* } */
  // verify that 2d computations are actually
  // activated
  //assert(fabs(vnorm[2]) < 1e-8);
}
#pragma end_opencl

#pragma start_opencl
void VecTransNumFlux2d(__private schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux)
{
  schnaps_real vn = sqrt(0.5) * vnorm[0] + sqrt(0.5) * vnorm[1];
  schnaps_real vnp = vn > 0.0 ? vn : 0.0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  vn = -vn;
  vnp = vn > 0.0 ? vn : 0.0;
  vnm = vn - vnp;
  flux[1] = vnp * wL[1] + vnm * wR[1];
}
#pragma end_opencl

#pragma start_opencl
void TransBoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm, schnaps_real *flux)
{
  schnaps_real wR[1];
  TransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
}
#pragma end_opencl

#pragma start_opencl
void TransBoundaryFlux2d(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm, schnaps_real *flux)
{
  schnaps_real wR[1];
  TransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
  //printf("imposed trans wR=%f x=%f %f %f\n",wR[0],x[0],x[1],x[2]);
}
#pragma end_opencl

// m = 2 test-case
#pragma start_opencl
void VecTransImposedData2d(const schnaps_real *x, const schnaps_real t, schnaps_real* w)
{
  schnaps_real vx = sqrt(0.5) * x[0] + sqrt(0.5) * x[1];
  schnaps_real xx = vx - t;
  w[0] = xx * xx;
  w[0] = sin(xx);
  xx = vx + t;
  w[1] = xx * xx;
  w[1] = sin(xx);
  /* w[0] = 1000; */
  /* w[1] = 2000; */
}
#pragma end_opencl

#pragma start_opencl
void VecTransBoundaryFlux2d(schnaps_real *x, schnaps_real t,
			    schnaps_real *wL, schnaps_real *vnorm,
			    schnaps_real *flux)
{
  schnaps_real wR[2];
  VecTransImposedData2d(x, t, wR);
  VecTransNumFlux2d(wL, wR, vnorm, flux);
}
#pragma end_opencl

void TransInitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TransImposedData(x, t, w);
}

void TransInitData2d(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TransImposedData2d(x, t, w);
}

void VecTransInitData2d(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  VecTransImposedData2d(x, t, w);
}

#pragma start_opencl
void TransImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w)
{
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
  schnaps_real vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];
  schnaps_real xx = vx - t;
  w[0] = cos(xx);
}
#pragma end_opencl

#pragma start_opencl
void TransImposedData2d(const schnaps_real *x, const schnaps_real t, schnaps_real *w)
{
  schnaps_real vx  = sqrt(0.5) * x[0] + sqrt(0.5) * x[1];
  schnaps_real xx = vx - t;
  w[0] = sin(3*xx);
}
#pragma end_opencl

#pragma start_opencl
void TestTransBoundaryFlux(schnaps_real *x, schnaps_real t,
			       schnaps_real *wL, schnaps_real* vnorm,
			       schnaps_real* flux) {
  schnaps_real wR[1];
  TestTransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
}
#pragma end_opencl

void TestTransBoundaryFlux2d(schnaps_real *x, schnaps_real t, schnaps_real *wL,
				 schnaps_real* vnorm, schnaps_real* flux) {
  schnaps_real wR[1];
  TestTransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
}

void TestTransInitData(schnaps_real *x, schnaps_real *w)
{
  schnaps_real t = 0;
  TestTransImposedData(x, t, w);
}

void TestTransInitData2d(schnaps_real *x, schnaps_real *w)
{
  schnaps_real t = 0;
  TestTransImposedData2d(x, t, w);
}

#pragma start_opencl
void TestTransImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  const schnaps_real sqrt_third =  sqrt(1.0/3.0);
  const schnaps_real transport_v[] = {sqrt_third, sqrt_third, sqrt_third};
  //const real transport_v[] = {1,0,0};
  schnaps_real vx
    = transport_v[0] * x[0]
    + transport_v[1] * x[1]
    + transport_v[2] * x[2];
  schnaps_real xx = vx - t;
  //w[0] = 1;
  w[0] = sin(3*xx);
  w[0] = xx * xx;
  //w[0]=0;
}
#pragma end_opencl

void TestTransImposedData2d(const schnaps_real *x, const schnaps_real t, schnaps_real *w) {
  const schnaps_real transport_v2d[] = {sqrt(0.5), sqrt(0.5), 0};
  schnaps_real vx
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  schnaps_real xx = vx - t;
  w[0] = xx * xx;
  w[0] = sin(3*xx);
}

void set_global_m(int m0)
{
  m = m0;
}

// Set the parameters for the Vlasov equation that are stored in the
// global space of model.h
void set_vlasov_params(Model *mod)
{
  m = mod->m;
  assert(m > 0);
  vlasov_mx = mod->vlasov_mx;
  vlasov_my = mod->vlasov_my;
  vlasov_mz = mod->vlasov_mz;
  assert(m == vlasov_mx * vlasov_my * vlasov_mz);
  vlasov_vmax = mod->vlasov_vmax;
}

void vlaTransInitData2d(schnaps_real *x, schnaps_real *w)
{
  schnaps_real t = 0;
  vlaTransImposedData2d(x, t, w);
}

// Return the component of the vlasov velocity with index id.
schnaps_real vlasov_vel(const int id, const int md, schnaps_real vlasov_vmax)
{
  int mid = md / 2;
  schnaps_real dv = vlasov_vmax / mid;
  return (id - mid) * dv;
}

void vlaTransNumFlux2d(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux)
{
  // 3m multiplies
  for(int ix = 0; ix < vlasov_mx; ++ix) {
    schnaps_real vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);

    for(int iy = 0; iy < vlasov_my; ++iy) {
      schnaps_real vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);

      schnaps_real vn = vx * vnorm[0]	+ vy * vnorm[1];
      schnaps_real vnp = vn > 0 ? vn : 0;
      schnaps_real vnm = vn - vnp;

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      flux[im] = vnp * wL[im] + vnm * wR[im];
    }
  }
  // Verify that 2d computations are actually activated
  assert(fabs(vnorm[2]) < _SMALL);
}

void vlaTransBoundaryFlux2d(schnaps_real *x, schnaps_real t,
			    schnaps_real *wL, schnaps_real *vnorm,
			    schnaps_real* flux)
{
  schnaps_real wR[m];
  vlaTransImposedData2d(x, t, wR);
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

// compact support bump (C-infinity, but not analytic):
schnaps_real compact_bump(schnaps_real r)
{
  if(fabs(r) >= 0.5)
    return 0;
  return exp(-1.0 / (1.0 - 4.0 * r * r));
}

schnaps_real icgaussian(schnaps_real r, schnaps_real sigma)
{
  schnaps_real PI = 4.0 * atan(1.0);
  return exp(-r * r / sigma) / (sigma * sqrt(2.0 * PI));
}

// 6th-degree polynomial with compact support
schnaps_real compact_poly6(schnaps_real r)
{
  if (fabs(2 * r) > 1)
    return 0;
  schnaps_real rrm1 = 2 * r - 1;
  schnaps_real rrp1 = 2 * r + 1;
  return -35.0 / 16.0 * rrm1 * rrm1 * rrm1 * rrp1 * rrp1 * rrp1;
}

void vlaTransImposedData2d(const schnaps_real *x, const schnaps_real t, schnaps_real *w)
{
  //real PI = 4.0 * atan(1.0);
  //real s2pi = sqrt(2.0 * PI);
  //real xval = 1.0;
  schnaps_real sigma = 0.1;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    schnaps_real vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    schnaps_real px = x[0] - vx * t;

    for(int iy = 0; iy < vlasov_my; ++iy) {
      schnaps_real vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      schnaps_real py = x[1] - vy * t;

      schnaps_real r = sqrt(px * px + py * py);
      //real pi = 4.0 * atan(1.0);
      schnaps_real pr = icgaussian(r, sigma);

      schnaps_real vr = sqrt(vx * vx + vy * vy);
      schnaps_real pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemcracs2014_imposed_data(const schnaps_real *x, const schnaps_real t, schnaps_real *w)
{
  schnaps_real sigma = 1.0;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    schnaps_real vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    schnaps_real px = x[0] - vx * t;

    for(int iy = 0; iy < vlasov_my; ++iy) {
      schnaps_real vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      schnaps_real py = x[1] - vy * t;

      schnaps_real r = sqrt(px * px + py * py);
      schnaps_real pr = compact_bump(r);
      //real pr = compact_poly6(r);

      schnaps_real vr = sqrt(vx * vx + vy * vy);
      schnaps_real pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemcracs2014a_imposed_data(const schnaps_real *x, const schnaps_real t, schnaps_real *w)
{
  schnaps_real sigma = 1.0;

  for(int ix = 0; ix < vlasov_mx; ++ix) {
    schnaps_real vx = vlasov_vel(ix, vlasov_mx, vlasov_vmax);
    schnaps_real px = x[0] - vx * t;

    for(int iy = 0; iy < vlasov_my; ++iy) {
      schnaps_real vy = vlasov_vel(iy, vlasov_my, vlasov_vmax);
      schnaps_real py = x[1] - vy * t;

      schnaps_real r = sqrt(px * px + py * py);
      //real pr = compact_bump(r);
      schnaps_real pr = compact_poly6(r);

      schnaps_real vr = sqrt(vx * vx + vy * vy);
      schnaps_real pvr = icgaussian(vr, sigma);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      w[im] = pr * pvr;
    }
  }
}

void cemracs2014_TransInitData(schnaps_real *x, schnaps_real *w)
{
  schnaps_real t = 0;
  cemcracs2014_imposed_data(x, t, w);
}

void cemracs2014a_TransInitData(schnaps_real *x, schnaps_real *w)
{
  schnaps_real t = 0;
  cemcracs2014a_imposed_data(x, t, w);
}

void cemracs2014_TransBoundaryFlux(schnaps_real *x, schnaps_real t,
				   schnaps_real *wL, schnaps_real *vnorm,
				   schnaps_real *flux)
{
  schnaps_real wR[m];
  for(unsigned int i = 0; i < m; ++i)
    wR[i] = 0;
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

void OneSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source)
{
  for(int i = 0; i < m; ++i) {
    source[i] = 1.0;
  }
}

void schnaps_model_load(char* model_source, Model *m){
  
  void (*foo)(schnaps_real *x, schnaps_real *w);
  
  
  int stat =  system("cp ../test/model_trans.cl model_trans.c");
  stat =  system("gcc -fpic -shared model_trans.c -o ./libmodel_trans.so");
  //stat = system("rm ../test/model_trans.c");

  void *flib = dlopen("./libmodel_trans.so", RTLD_GLOBAL | RTLD_LAZY);
  
  m->NumFlux = dlsym(flib,"NumFlux"); 
  m->BoundaryFlux = dlsym(flib,"BoundaryFlux");
  m->ImposedData = dlsym(flib,"ImposedData");
  m->InitData = dlsym(flib,"InitData");
  m->Source = NULL;

  
  // foo = dlsym(flib,"InitData");
  // printf("InitData = %p \n", m->InitData);
  // printf("ImposedData = %p \n", m->ImposedData);

  // schnaps_real x[3] = {1, 1, 1};
  // schnaps_real w[3] = {0, 0, 0};

  // m->InitData(x, w);
  // printf("x=%f %f %f w=%f %f, %f \n", x[0], x[1], x[2], w[0], w[1], w[2]);
 
}

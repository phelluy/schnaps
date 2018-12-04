#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "maxwell.h"

void Coil2DImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[])
{
  schnaps_real r = x[0] * x[0] + x[1] * x[1];
  w[0] = 0;
  w[1] = 0;
  w[2] = r > 1 ? 0 : 1;
  w[3] = 0;
  w[4] = 0;
  w[5] = 0;
  w[6] = 0;
}

void coil_pre_dtfields(void *simu);

void coil_pre_dtfields(void *simu){
  Simulation * simu2=simu;
  AccumulateParticles(simu,simu2->w);
}


void Coil2DSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source)
{
  // w: (Ex, Ey, Hz, Hz, \lambda, rho, Jx, Jy)
  
  // FIXME add documentation

  static int icall = 0;
  
  const schnaps_real khi = 1.0;
  source[0] = -w[4];
  source[1] = -w[5];
  source[2] = 0;
  source[3] = 0;//instead of khi * w[6]: we want div E =0
  source[4] = 0;
  source[5] = 0;
  source[6] = 0;

  icall++;
  //printf("source call %d w=%f\n",icall,w[4]);


  
}




void Coil2DBoundaryFlux(schnaps_real x[3], schnaps_real t, schnaps_real wL[], schnaps_real *vnorm,
			schnaps_real *flux)
{
  schnaps_real wR[7];
  Coil2DImposedData(x, t, wR);
  Maxwell2DNumFlux_upwind(wL, wR, vnorm, flux);
}

void Coil2DInitData(schnaps_real x[3], schnaps_real w[])
{
  schnaps_real t = 0;
  Coil2DImposedData(x, t, w);
}

int TestCoil2D(void)
{
  bool test = true;

  char *mshname =  "../test/testmacromesh.msh";
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,mshname);
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);


  // test gmsh file reading
  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");
  Detect2DMacroMesh(&mesh);
  assert(mesh.is2d);
  BuildConnectivity(&mesh);
  //PrintMacroMesh(&m);

  Model model;

  model.m = 7; // num of conservative variables

  model.NumFlux = Maxwell2DNumFlux_upwind;
  model.BoundaryFlux = Coil2DBoundaryFlux;
  model.InitData = Coil2DInitData;
  model.Source = Coil2DSource;
  //model.Source = NULL;
  model.ImposedData = Coil2DImposedData;
    
  int deg[]={2, 2, 0};
  int raf[]={4, 4, 1};

  CheckMacroMesh(&mesh, deg, raf);

  
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.pre_dtfields = coil_pre_dtfields; // must be DEFINED after init

  printf("cfl param=%f \n", simu.hmin);
   
  
  PIC pic;
  simu.pic = &pic;

  InitPIC(&pic,100); 
  CreateCoil2DParticles(&pic, &mesh);
  PlotParticles(&pic, &mesh);

  // time evolution
  schnaps_real tmax = 0.5;
  simu.cfl=0.2;
  simu.vmax = 1;
  RK2(&simu, tmax);

  // for a good comparison with the exact solution
  // we have to cancel the data used for the source term

  for(int ie = 0; ie < simu.macromesh.nbelems; ie++){
    field *f = simu.fd + ie;
    int offset = ie * f->wsize;
    //real *wf = f->wn;
    schnaps_real *wf = simu.w + offset;
    //assert(1==3); xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++){
      int iv = 4;
      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      wf[imem] = 0;
      iv = 5;
      imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      wf[imem]=0;
      iv = 6;
      imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      wf[imem]=0;
    } 
  }
 

  
  // Save the results and the error
  PlotFields(2, false, &simu, NULL, "dgvisu.msh");
  PlotFields(2, true, &simu, "error", "dgerror.msh");

  DisplaySimulation(&simu);

  schnaps_real dd = L2error(&simu);
  schnaps_real tolerance = 0.3;
  test = test && (dd < tolerance);
  printf("L2 error: %f\n", dd);

  FreeMacroMesh(&mesh);

  return test;
}

int main(void) {
  int resu = TestCoil2D();
  if (resu) 
    printf("Coil2D test OK!\n");
  else 
    printf("Coil2D failed !\n");
  return !resu;
}

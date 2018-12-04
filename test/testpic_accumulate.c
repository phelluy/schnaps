#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int TestPICAccumulate(void);

int main(void) {
  // Unit tests
  int resu=TestPICAccumulate();
  if (resu) printf("PIC accumulate test OK !\n");
  else printf("PIC accumulate  test failed !\n");
  return !resu;
} 

void Maxwell2DConstInitData(schnaps_real x[3], schnaps_real w[]) {
  w[0]=0;
  w[1]=0;
  w[2]=1;
  w[3]=0;
  w[4]=0;
  w[5]=0;
  w[6]=0;
}


// some unit tests of the macromesh code
int TestPICAccumulate(void)
{
  bool test = true;

  char *mshname =  "../test/testmacromesh.msh";
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,mshname);
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  int deg4[]={4, 4, 4};
  int raf4[]={1, 1, 1};

  CheckMacroMesh(&mesh, deg4, raf4);

  // test gmsh file reading
  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");
  BuildConnectivity(&mesh);
  CheckMacroMesh(&mesh, deg4, raf4);
  //PrintMacroMesh(&m);

  PIC pic;

  InitPIC(&pic,1); 
  CreateParticles(&pic,&mesh);
  PlotParticles(&pic,&mesh);

  Model model;

  model.m = 7; // num of conservative variables

  /* f.model.NumFlux = Maxwell2DNumFlux; */
  /* f.model.BoundaryFlux = Maxwell2DBoundaryFlux; */
  model.InitData = Maxwell2DConstInitData;
  /* f.model.ImposedData = Maxwell2DImposedData; */
    
  int deg[]={1, 1, 0};
  int raf[]={1, 1, 1};

  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  simu.pic = &pic;

  // place the particle 
  pic.xv[0]=0;
  pic.xv[1]=0;
  pic.xv[2]=0.5;
  schnaps_real xref[3];
  pic.cell_id[0]=NumElemFromPoint(&mesh,pic.xv,xref);
  pic.xv[0]=xref[0];  
  pic.xv[1]=xref[1];  
  pic.xv[2]=xref[2];  
  pic.xv[3]=1;
  pic.xv[4]=0;
  pic.xv[5]=0;

  pic.weight = 1;

  PlotParticles(&pic,&mesh);


  AccumulateParticles(&simu, simu.w);

  int ie=2;
  int ipg=2;
  int iv=4;

  field *f = simu.fd + ie;

  int imem=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);

  printf("w=%f wex=%f\n",f->wn[imem],1/1.96);
  test = test && (fabs(f->wn[imem]-1/1.96) < 1e-8);


  Displayfield(f);

  FreeMacroMesh(&mesh);




  return test;
}

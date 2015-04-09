#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  // Unit tests
  int resu=TestPIC();
  if (resu) printf("PIC test OK !\n");
  else printf("PIC test failed !\n");
  return !resu;
} 

void Maxwell2DConstInitData(double x[3], double w[]) {
  w[0]=1;
  w[1]=0;
  w[2]=0;
  w[3]=0;
}


// some unit tests of the macromesh code
int TestPIC(void)
{
  MacroMesh m;

  bool test=true;

  int param[]={4, 4, 4, 1, 1, 1, 0};
  
  field f;
  // test gmsh file reading
  ReadMacroMesh(&(f.macromesh), "test/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));
  CheckMacroMesh(&(f.macromesh), param);
  //PrintMacroMesh(&m);

  PIC pic;

  InitPIC(&pic,10); 
  CreateParticles(&pic,&(f.macromesh));
  PlotParticles(&pic,&(f.macromesh));

  f.model.m = 4; // num of conservative variables

  /* f.model.NumFlux = Maxwell2DNumFlux; */
  /* f.model.BoundaryFlux = Maxwell2DBoundaryFlux; */
  f.model.InitData = Maxwell2DConstInitData;
  /* f.model.ImposedData = Maxwell2DImposedData; */
  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  Initfield(&f);

  pic.dt=0.01;
  for(int iter=0;iter<100;iter++){
    PushParticles(&f,&pic);
  }
   PlotParticles(&pic,&(f.macromesh));
 

  return test;
}

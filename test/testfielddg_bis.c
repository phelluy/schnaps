#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>


int TestfieldDG_bis(void){

  int test = true;

  Model model;
  
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;
  model.Source = NULL;

  int deg[]={4, 4, 4};
  int raf[]={1, 1, 1};
  /* int deg[]={1, 1, 1}; */
  /* int raf[]={1, 1, 1}; */
 
  MacroMesh mesh;
  //ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  ReadMacroMesh(&mesh,"../test/testcube2.msh");
  BuildConnectivity(&mesh);

  /* real A[3][3] = {{10,2 , 0}, {0, 1, -0.1}, {0, 0.1,1}}; */
  /* real x0[3] = {1, 2, 3}; */
  /* AffineMapMacroMesh(&mesh,A,x0); */

  CheckMacroMesh(&mesh, deg, raf);

  //PrintMacroMesh(&mesh);

  Simulation simu;
  EmptySimulation(&simu);
  
  InitSimulation(&simu, &mesh, deg, raf, &model);

  simu.tnow = 0;

  DtFields_bis(&simu, simu.w, simu.dtw);
  
  for(int i = 0; i < simu.wsize; i++){
    printf("i=%d dtw=%.8e\n",i,simu.dtw[i]);
  }
  //DisplaySimulation(&simu);

  //PlotFields(0, false, &simu, NULL, "visu_spu.msh");
  //PlotFields(0, true, &simu, "error", "error.msh");

  //Test the time derivative with the exact solution
  schnaps_real test2 = 0;
  for(int i = 0;
      i < model.m * mesh.nbelems * NPG(deg,raf);
      i++){
    schnaps_real errloc = fabs(4 * simu.w[i] - pow(simu.dtw[i], 2));
    //real errloc = pow(simu.dtw[i], 2.);
   //errloc = fabs(pow(simu.dtw[i], 2));
    test2 += errloc * errloc;
    test = test && errloc < 1e-2;
    //printf("i=%d err=%f \n",i,4 * w[i] - pow(dtw[i], 2));
    //assert(test);
  }

  printf("error=%.10e\n",sqrt(test2/ (mesh.nbelems * NPG(deg,raf)) ));

  FreeMacroMesh(&mesh);
  
  return test;
};

int main(void) {
  // Unit tests
  int resu = TestfieldDG_bis();
  if (resu) printf("field DG bis test OK !\n");
  else printf("field DG bis test failed !\n");
  return !resu;
} 

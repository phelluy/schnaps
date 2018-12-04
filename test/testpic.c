#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int TestPIC(void);


int main(void) {
  // Unit tests
  int resu=TestPIC();
  if (resu) printf("PIC test OK !\n");
  else printf("PIC test failed !\n");
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
int TestPIC(void)
{
  bool test = true;

  char *mshname =  "../test/testmacromesh.msh";
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,mshname);
  Detect2DMacroMesh(&mesh);
  BuildConnectivity(&mesh);

  int deg[]={4, 4, 4};
  int raf[]={1, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);

  // test gmsh file reading
  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");
  BuildConnectivity(&mesh);
  CheckMacroMesh(&mesh, deg, raf);
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
    
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  // place the particle at (0,1,0) and v=(1,0,0)
  pic.xv[0]=0;
  pic.xv[1]=1;
  pic.xv[2]=0.5;
  schnaps_real xref[3];
  pic.cell_id[0]=NumElemFromPoint(&mesh,pic.xv,xref);
  pic.xv[0]=xref[0];  
  pic.xv[1]=xref[1];  
  pic.xv[2]=xref[2];  
  pic.xv[3]=1;
  pic.xv[4]=0;
  pic.xv[5]=0;

  simu.pic = &pic;

  schnaps_real final_pos_phy[3]={0,-1,0.5};
  schnaps_real final_pos[3];
  int final_cell=NumElemFromPoint(&mesh,
				  final_pos_phy,final_pos);

  pic.dt=0.001;
  for(int iter=0;iter<3141;iter++){
    PushParticles(&simu,&pic);
  }
  PlotParticles(&pic,&mesh);
 

  printf("Dist=%f\n",Dist(pic.xv,final_pos));

  test = test && (Dist(pic.xv,final_pos) < 1e-3);

  test = test && (final_cell == pic.cell_id[0]);

  AccumulateParticles(&simu, simu.w);

  FreeMacroMesh(&mesh);




  return test;
}

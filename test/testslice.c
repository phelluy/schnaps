#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test/test.h"
#include "quantities_vp.h"
#include "gyro.h"

int TestSlice(void);

int main(void) {
  
  // unit tests
    
  int resu=TestSlice();
	 
  if (resu) printf("slice test OK !\n");
  else printf("slice test failed !\n");

  return !resu;
} 


int TestSlice(void) { 

  bool test=true;

  // read the gmsh file
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../geo/cylindre.msh");

  mesh.period[2]=2;
  BuildConnectivity(&mesh);

  int vec=1;
  
    
  int deg[]={1, 1, 1};
  int raf[]={1, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);

  Simulation simu;
  EmptySimulation(&simu);

  
  Model model;

  //extern KineticData  schnaps_kinetic_data;
  
  KineticData *kd = &schnaps_kinetic_data;
  // nbelemv = 10
  // deg_v = 4
  InitKineticData(&schnaps_kinetic_data,1,2);
  kd->solve_quasineutrality = true;
  
  printf("_MV=%d\n",kd->mv);
  printf("_INDEX_MAX=%d\n",kd->index_max);
  printf("_INDEX_MAX_KIN=%d\n",kd->index_max_kin);

  model.m= kd->index_max; // num of conservative variables
  model.NumFlux=GyroUpwindNumFlux;
  //model.NumFlux=GyroZeroNumFlux;
  //model.NumFlux=NULL;
  model.BoundaryFlux=GyroBoundaryFlux;
  model.InitData=GyroInitData;
  model.ImposedData=GyroImposedData;
  model.Source = NULL;
  //model.Source = GyroSource;


  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=kd->index_phi;

  InitContinuousSolver(&ps,&simu,1,nb_var,listvar);

  test= test && (ps.slice_size == 8);


  return test; 

};

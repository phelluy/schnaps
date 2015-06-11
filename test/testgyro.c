#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "gyro.h"



int main(void) {
  
  // unit tests
    
  int resu=TestGyro();
	 
  if (resu) printf("gyro test OK !\n");
  else printf("gyro test failed !\n");

  return !resu;
} 


int TestGyro(void) {

  bool test=true;

  field f;
  init_empty_field(&f);

  int vec=1;
  
  f.model.m=_INDEX_MAX; // num of conservative variables
  f.model.NumFlux=Gyro_Lagrangian_NumFlux;
  //f.model.NumFlux=NULL;
  f.model.BoundaryFlux=Gyro_Lagrangian_BoundaryFlux;
  f.model.InitData=GyroInitData;
  f.model.ImposedData=GyroImposedData;
  //f.model.Source = NULL;
  f.varindex=GenericVarindex;
    
  f.interp.interp_param[0]= f.model.m;//_MV;  // _M
  f.interp.interp_param[1]=1;  // x direction degree
  f.interp.interp_param[2]=1;  // y direction degree
  f.interp.interp_param[3]=1;  // z direction degree
  f.interp.interp_param[4]=1;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement
  // read the gmsh file
  //ReadMacroMesh(&(f.macromesh),"geo/cube.msh");
  ReadMacroMesh(&(f.macromesh),"geo/cylindre.msh");
  // try to detect a 2d mesh
  //bool is1d=Detect1DMacroMesh(&(f.macromesh));
  //assert(is1d);

  // mesh preparation
  f.macromesh.period[2]=2;
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  Initfield(&f);
  f.model.Source = GyroSource;
  //f.update_before_rk=NULL;
  f.update_after_rk=NULL; 
  f.vmax = _VMAX; // maximal wave speed 
  //f.macromesh.is1d=true;
  //f.is1d=true;

  // prudence...
  CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

  printf("cfl param =%f\n",f.hmin);

  // time derivative
  //dtfield(&f);
  //Displayfield(&f);
  //assert(1==2);
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  f.model.cfl=0.2;
  real dt = 0;
  RK2(&f,0.1, dt);
 
  // save the results and the error
  Plotfield(0,(1==0),&f,"sol","dgvisu.msh");
  //Plotfield(0,(1==1),&f,"error","dgerror.msh");

  double dd=L2error(&f);
  //double dd_l2_vel =GyroL2VelError(&f)
  //double dd_Kinetic=L2_Kinetic_error(&f);
  
  //printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  printf("tnow is  %lf\n",f.tnow);
  Velocity_distribution_plot(f.wn);
  //test= test && (dd<3e-4);


  return test;

};

#include "solverpoisson.h"
#include "geometry.h"
#include "quantities_vp.h"


void SolvePoisson(field *f,double * w,int type_bc, double bc_l, double bc_r){

  // for the moment, works only for the 1d case
  assert(f->is1d);

  // assembly of the rigidity matrix

  Skyline sky;

  // number of equation of the Poisson solver
  // = number of nodes in the mesh
  int degx=f->interp.interp_param[1];
  int nelx=f->interp.interp_param[4];
  double xmin=0;
  double xmax=1;  // TO DO: compute the maximal x coordinate
  int neq=degx*nelx+1;
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m=f->model.m;

  InitSkyline(&sky,neq);

  // compute the profile of the matrix
  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	int ino=iloc + ie * degx;
	int jno=jloc + ie * degx;
	SwitchOn(&sky,ino,jno);
      }
    }
  }

  AllocateSkyline(&sky);

  // local matrix (assuming identical finite elements)
  double aloc[degx+1][degx+1];
  for(int iloc=0;iloc<degx+1;iloc++){
    for(int jloc=0;jloc<degx+1;jloc++){
      aloc[iloc][jloc]=0;
    }
  }
  for(int ipg=0;ipg<degx+1;ipg++){
    double omega=wglop(degx,ipg);
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	double dxi=dlag(degx,iloc,ipg);
	double dxj=dlag(degx,jloc,ipg);
	aloc[iloc][jloc]+=dxi*dxj*omega*nelx;
      }
    }
  }

  // assembly of the matrix
  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	int ino=iloc + ie * degx;
	int jno=jloc + ie * degx;
	double val = aloc[iloc][jloc];
	SetSkyline(&sky,ino,jno,val);
      }
    }
  }


  // dirichlet boundary condition at the first and last location
  if(type_bc == _Dirichlet_Poisson_BC){
    SetSkyline(&sky,0,0,1e20);
    SetSkyline(&sky,neq-1,neq-1,1e20);
  }

  //DisplaySkyline(&sky);

  FactoLU(&sky);

    // source assembly 
  double source[neq];
  for(int i=0;i<neq;i++){
    source[i]=0;
  }

  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      double omega=wglop(degx,iloc);
      int ino=iloc + ie * degx;  
      int imem=f->varindex(f->interp_param,0,iloc+ie*(degx+1),_INDEX_RHO);
      double charge=w[imem];     
      source[ino]+= charge*omega/nelx;
    }
  }

  // Apply dirichlet Boundary condition
  if(type_bc == _Dirichlet_Poisson_BC){
    source[0]=1e20*bc_l;
    source[neq-1]=1e20*bc_r;
  }
  
  double sol[neq];
  SolveSkyline(&sky,source,sol);


  // now put the solution at the right place
  for(int ie=0;ie<nelx;ie++){
     for(int ipg=0;ipg<degx+1;ipg++){
       // position in the continuous vector
       int ino=ipg + ie * degx;
       // position in the DG vector
       int imem=f->varindex(f->interp_param,0,ipg+ie*(degx+1),_INDEX_PHI);
       w[imem]=sol[ino];
     }
  }
	
  FreeSkyline(&sky);


  Compute_electric_field(f,w);

}

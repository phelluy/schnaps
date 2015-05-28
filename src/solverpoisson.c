#include "solverpoisson.h"
#include "geometry.h"
#include "quantities_vp.h"



int CompareFatNode(const void* a,const void* b){

  FatNode* fna = (FatNode*) a;
  FatNode* fnb = (FatNode*) b;

  int r = fna->x_int[0]-fnb->x_int[0];
  if (r==0)
    r = fna->x_int[1]-fnb->x_int[1];
  if (r==0)
    r = fna->x_int[2]-fnb->x_int[2];
  return r;

}

int BuildFatNodeList(field* f,FatNode* fn_list){

  int nb_dg_nodes =  NPG(f->interp_param+1) * f->macromesh.nbelems;
  
  fn_list = malloc(nb_dg_nodes * sizeof(FatNode));
  assert(fn_list);

  int big_int = 1 << 28; // 2**28 = 268 435 456


  int ino=0;
  real* xmin=f->macromesh.xmin;
  real* xmax=f->macromesh.xmax;
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    
    for(int ipg = 0; ipg < NPG(f->interp_param + 1); ipg++) {
      real xpg[3];
      real xref[3];
      ref_pg_vol(f->interp_param + 1, ipg, xref, NULL, NULL);
      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, NULL,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      fn_list[ino].x[0]=xpg[0];
      fn_list[ino].x[1]=xpg[1];
      fn_list[ino].x[2]=xpg[2];

      fn_list[ino].x_int[0]=(int) (xpg[0]-xmin[0])/(xmax[0]-xmin[0]) * big_int;
      fn_list[ino].x_int[1]=(int) (xpg[0]-xmin[1])/(xmax[0]-xmin[1]) * big_int;
      fn_list[ino].x_int[2]=(int) (xpg[0]-xmin[2])/(xmax[0]-xmin[2]) * big_int;
      
      fn_list[ino].dg_index = ino;

      ino++;
    }
  }

  assert(ino == nb_dg_nodes);
  return nb_dg_nodes;

}

void InitPoissonSolver(PoissonSolver* ps, field* fd,int charge_index){

  ps->fd = fd;
  ps->charge_index = charge_index;

  ps->nb_dg_nodes =  NPG(fd->interp_param+1) * fd->macromesh.nbelems;

  // first step: paste the nodes of the DG mesh
  BuildFatNodeList(fd,ps->fn_list);

  


}

void SolvePoisson1D(field *f,real * w,int type_bc, real bc_l, real bc_r){


  real charge_average;
  charge_average=0;

  if(type_bc == _Periodic_Poisson_BC){
    charge_average=Computation_charge_average(f,w);
    bc_l=0;
    bc_r=0;
  }
  else {
    charge_average=0;
  }
  
  // for the moment, works only for the 1d case
  assert(f->macromesh.is1d);

  // assembly of the rigidity matrix

  Skyline sky;

  // number of equation of the Poisson solver
  // = number of nodes in the mesh
  int degx=f->interp.interp_param[1];
  int nelx=f->interp.interp_param[4];
  real xmin=0;
  real xmax=1;  // TO DO: compute the maximal x coordinate
  int neq=degx*nelx+1;
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m = f->model.m;

  InitSkyline(&sky, neq);

  // compute the profile of the matrix
  for(int ie = 0; ie < nelx; ie++){
    for(int iloc = 0; iloc < degx + 1; iloc++){
      for(int jloc = 0; jloc < degx + 1; jloc++){
	int ino = iloc + ie * degx;
	int jno = jloc + ie * degx;
	SwitchOn(&sky, ino, jno);
      }
    }
  }

  AllocateSkyline(&sky);

  // local matrix (assuming identical finite elements)
  real aloc[degx+1][degx+1];
  for(int iloc=0;iloc<degx+1;iloc++){
    for(int jloc=0;jloc<degx+1;jloc++){
      aloc[iloc][jloc]=0;
    }
  }
  for(int ipg=0;ipg<degx+1;ipg++){
    real omega=wglop(degx,ipg);
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	real dxi=dlag(degx,iloc,ipg);
	real dxj=dlag(degx,jloc,ipg);
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
	real val = aloc[iloc][jloc];
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
  real source[neq];
  for(int i=0;i<neq;i++){
    source[i]=0;
  }

  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      real omega=wglop(degx,iloc);
      int ino=iloc + ie * degx;  
      int imem=f->varindex(f->interp_param,0,iloc+ie*(degx+1),_INDEX_RHO);
      real charge=w[imem];          
      source[ino]+= (charge-charge_average)*omega/nelx;
    }
  }

  // Apply dirichlet Boundary condition
  
  source[0]=1e20*bc_l;
  source[neq-1]=1e20*bc_r;
  
  real sol[neq];
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

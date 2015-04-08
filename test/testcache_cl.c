#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int Testfield(void){
  int test = true;

  double cfl = 0.5;
  int degx = 2;
  int nraf = 4;
  double vmax = 1.0;
  int mx = 2;
  int my = 4;

  field f;
  f.varindex = GenericVarindex;
  f.model.vlasov_mz = 1;
  f.model.cfl = cfl;
  f.model.NumFlux = vlaTransNumFlux2d;
  f.model.vlasov_mx = mx;
  f.model.vlasov_my = my;
  f.model.vlasov_vmax = vmax;
  f.model.m = f.model.vlasov_mx * f.model.vlasov_my * f.model.vlasov_mz;
  f.model.BoundaryFlux = vlaTransBoundaryFlux2d;
  f.model.InitData = vlaTransInitData2d;
  f.model.ImposedData = vlaTransImposedData2d;

  f.interp.interp_param[0] = f.model.m; // _M
  f.interp.interp_param[1] = degx; // x direction degree
  f.interp.interp_param[2] = degx; // y direction degree
  f.interp.interp_param[3] = 1; // z direction degree
  f.interp.interp_param[4] = nraf; // x direction refinement
  f.interp.interp_param[5] = nraf; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  set_vlasov_params(&(f.model));

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));

  Initfield(&f);
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  int raf[3]={nraf,nraf,1};
  int deg[3]={degx,degx,1};
  int ic[3]={1,0,3};
  int ix[3]={1,2,0};

  int ipg,ic0[3],ix0[3];

  for(int ii=0;ii<3;ii++){
    ic0[ii]=ic[ii];
    ix0[ii]=ix[ii];
  }

  xyz_to_ipg(raf,deg,ic,ix,&ipg);
  ipg_to_xyz(raf,deg,ic,ix,&ipg);

  for(int ii=0;ii<3;ii++){
    test = test && (ic0[ii]==ic[ii]);
    test = test && (ix0[ii]==ix[ii]);
  }

  for(int ii=0;ii<3;ii++){
    ic[ii]=raf[ii]-1;
    ix[ii]=deg[ii];
  }
  

  int worksize;  // number of gauss points in the macrocell

  xyz_to_ipg(raf,deg,ic,ix,&worksize);
  worksize++;

  printf("worksize=%d\n",worksize);

  int groupsize; // number of gauss points in a subcell

  for(int ii=0;ii<3;ii++){
    ic[ii]=0;
    ix[ii]=deg[ii];
  }

  xyz_to_ipg(raf,deg,ic,ix,&groupsize);
  groupsize++;
  printf("groupsize=%d\n",groupsize);

  int cache_size_in=(deg[0]+3)*(deg[1]+3)*(deg[2]+3)*f.model.m;
  printf("cache_size_in=%d\n",cache_size_in);
  int cache_size_out=(deg[0]+1)*(deg[1]+1)*(deg[2]+1)*f.model.m;
  printf("cache_size_out=%d\n",cache_size_out);
    
  int prec=4;
  printf("memory usage=%d\n",(cache_size_out+cache_size_in)*prec);
 
  int group_id=4;
  assert(group_id < worksize/groupsize);
  // kernel simulation
  for(int local_id=0;  local_id < groupsize;local_id++){
    int global_id=local_id+group_id*groupsize;

    // prefetching
    int nfetch=(f.model.m * cache_size_in) / cache_size_out;
    
    int ifetch;
    for(ifetch=0;ifetch<nfetch;ifetch++){
      // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
      // (local_id,ifetch)-> imem2 (varindex_local)
      
      // wnloc[imem2]=wn[imem];
    }
    // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
    // (local_id,ifetch)-> imem2 (varindex_local)
    
    // wnloc[imem2]=wn[imem]; if imem2 < cache_size_in;     
    

    // barrier

    // work on wnloc and dtwnloc ..........


    // barrier

    // back to global mem
    for(ifetch=0;ifetch<f.model.m;ifetch++){
      // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
      // (local_id,ifetch)-> imem2 (varindex_local2)
      // dtwn[imem]=dtwnloc[imem2]
    }
    
    
  }

  
  return test;
}

int main(void) {
  // Unit tests
  int resu = Testfield();
  if (resu)
    printf("field test OK !\n");
  else 
    printf("field test failed !\n");
  return !resu;
} 

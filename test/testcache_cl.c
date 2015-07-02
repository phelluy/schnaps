#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestCache(void){
  int test = true;

  real cfl = 0.5;
  int degx = 1;
  int nraf = 1;
  int nrafx = 2;
  real vmax = 1.0;
  int mx = 2;
  int my = 1;

  field f;
  init_empty_field(&f);
  
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
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = nrafx; // x direction refinement
  f.interp.interp_param[5] = nraf; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  set_vlasov_params(&(f.model));

  ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  // Try to detect a 2d mesh
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);  
  BuildConnectivity(&(f.macromesh));

  Initfield(&f);
  //CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  int raf[3]={nrafx,nraf,1};
  int deg[3]={degx,degx,0};
  int ic[3]={0,0,0};
  int ix[3]={0,1,0};

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
 
  int group_id=0;
  assert(group_id < worksize/groupsize);
  // kernel simulation
  for(int local_id=0;  local_id < groupsize;local_id++){
    int global_id=local_id+group_id*groupsize;

    // get the central cell indices
    int ic0[3],ix0[3];
    ipg_to_xyz(raf,deg,ic0,ix0,&global_id);

    // prefetching
    int nfetch=(f.model.m * cache_size_in) / cache_size_out;
    printf("nfetch=%d\n",nfetch);

    
    int ifetch;
    for(ifetch=0;ifetch<nfetch;ifetch++){

      // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
      int cache_in=local_id+ifetch*groupsize;
      int iw = local_id % f.model.m;
      int ipg_in= cache_in / nfetch ;

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
  int resu = TestCache();
  if (resu)
    printf("field test OK !\n");
  else 
    printf("field test failed !\n");
  return !resu;
} 

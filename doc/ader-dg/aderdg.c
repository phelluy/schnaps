#include "aderdg.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

void InitADERDG(ADERDG* adg,int nbelems,double xmin,double xmax){

  adg->nbelems=nbelems;
  adg->nbfaces=nbelems+1;

  adg->xmin=xmin;
  adg->xmax=xmax;

  adg->face = malloc(sizeof(double) * adg->nbfaces);
  assert(adg->face);
  
  adg->dx=(xmax-xmin)/nbelems;

  for(int i=0;i<adg->nbfaces;i++){
    adg->face[i]=xmin+i*adg->dx;
  }

  adg->dt = adg->dx;  // to do put the velocity
  
  adg->ncfl=1;

  adg->cell_level = malloc(sizeof(int) * adg->nbelems);
  assert(adg->cell_level);
  for(int i=0;i<adg->nbelems;i++){
    adg->cell_level[i]=0;
  }

  adg->face_level = malloc(sizeof(int) * adg->nbfaces);
  assert(adg->face_level);
  for(int i=0;i<adg->nbfaces;i++){
    adg->face_level[i]=0;
  }

  adg->wnow = malloc(sizeof(double) * adg->nbelems * _M * (_D+1));
  assert(adg->wnow);

  for(int ie=0;ie<adg->nbelems;ie++){
    for(int j=0;j<_D+1;j++){
      double h=adg->face[ie+1]-adg->face[ie];
      double x=adg->face[ie]+ h * gauss_lob_point[gauss_lob_offset[_D]+j];
      ExactSol(x,0,adg->wnow + ( ie * (_D+1) + j) * _M );
    }
  }

  adg->wnext = malloc(sizeof(double) * adg->nbelems * _M * (_D+1));
  assert(adg->wnext);
  
  memcpy(adg->wnext,adg->wnow,sizeof(double) * adg->nbelems * _M * (_D+1));
    



  adg->face_level = malloc(sizeof(int) * adg->nbfaces);
  assert(adg->face_level);
  
  

}

void ExactSol(double x,double t,double* w){

  for(int iv=0;iv<_M;iv++){
    w[iv]=sin(2*M_PI*(x-t));
  }

}



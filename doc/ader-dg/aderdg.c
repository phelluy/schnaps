#include "aderdg.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int main(void){

  ADERDG adg;

  InitADERDG(&adg,10,-1,1);

  Plot(&adg);

  return 0;

}

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

void Plot(ADERDG* adg)
{

  FILE * gnufile;
  gnufile = fopen("adgplot.dat", "w" );

  for(int ie=0;ie<adg->nbelems;ie++){
    for(int j=0;j<_D+1;j++){
      double wex[_M];
      double h=adg->face[ie+1]-adg->face[ie];
      double x=adg->face[ie]+ h * gauss_lob_point[gauss_lob_offset[_D]+j];
      ExactSol(x,adg->tnow,wex);
      fprintf(gnufile,"%f ",x);
      for (int iv=0;iv<_M;++iv){
	fprintf(gnufile,"%f ",adg->wnow[(ie * (_D+1) + j) * _M + iv]);	
      }
      for (int iv=0;iv<_M;++iv){
	fprintf(gnufile,"%f ",wex[iv]);	
      }
      fprintf(gnufile,"\n");
    }
  }

  fclose(gnufile);
 

}

void ExactSol(double x,double t,double* w){

  w[0]=sin(2*M_PI*(x-t));
  w[1]=sin(2*M_PI*(x+t));

}


void Predictor(ADERDG* adg,double h,double* wcell_init,
	       double s,double* wcell)
{

  double w0[_D+1],w[_D+1];

  

  for(int iv=0; iv < _M; ++iv){
    for(int ipg=0;ipg <_D+1;++ipg){
      w0[ipg]=wcell_init[iv+_M*ipg];
    }

    double t=s/h*velocity[iv];

    switch(_D) {

    case 1 :
      w[0] = (1 + t) * w0[0] - t * w0[1];
      w[1] = t * w0[0] + (1 - t) * w0[1];
      break;

    case 2 :
      w[0] = (1 + 3 * t + 2 * t * t) * w0[0] 
	+ (-4 * t * t - 4 * t) * w0[1] + (2 * t * t + t) * w0[2];
      w[1] = (2 * t * t + t) * w0[0] 
	+ (1 - 4 * t * t) * w0[1] + (2 * t * t - t) * w0[2];
      w[2] = (2 * t * t - t) * w0[0] 
	+ (-4 * t * t + 4 * t) * w0[1] + (1 + 2 * t * t - 3 * t) * w0[2];
      break;

    case 3 :
      w[0] = 0.1000000000e1 * (0.1e1 + 0.6e1 * t + 0.10e2 * t * t + 0.5e1 * pow(t, 0.3e1)) * w0[0] - 0.1809016994e1 * (0.1065247584e2 * t + 0.4472135954e1 + 0.618033988e1 * t * t) * t * w0[1] + 0.6909830056e0 * (0.2065247584e2 * t + 0.4472135954e1 + 0.1618033988e2 * t * t) * t * w0[2] - 0.1000000000e1 * t * (0.1e1 + 0.5e1 * t + 0.5e1 * t * t) * w0[3];
      w[1] = 0.5000000000e0 * t * (0.1170820393e2 * t + 0.3236067977e1 + 0.10e2 * t * t) * w0[0] - 0.3618033987e0 * (-0.2763932023e1 + 0.2763932023e2 * t * t + 0.3090169942e2 * pow(t, 0.3e1)) * w0[1] + 0.1000000000e1 * t * (0.1118033988e2 * t * t + 0.5e1 * t - 0.2236067977e1) * w0[2] - 0.5000000000e0 * t * (0.1708203931e1 * t - 0.1236067977e1 + 0.10e2 * t * t) * w0[3];
      w[2] = -0.5000000000e0 * t * (0.1708203931e1 * t + 0.1236067977e1 - 0.10e2 * t * t) * w0[0] - 0.1000000000e1 * t * (0.1118033988e2 * t * t - 0.5e1 * t - 0.2236067977e1) * w0[1] + 0.1381966011e0 * (0.7236067977e1 - 0.7236067977e2 * t * t + 0.8090169942e2 * pow(t, 0.3e1)) * w0[2] + 0.5000000000e0 * t * (0.1170820393e2 * t - 0.3236067977e1 - 0.10e2 * t * t) * w0[3];
      w[3] = 0.1000000000e1 * t * (0.1e1 - 0.5e1 * t + 0.5e1 * t * t) * w0[0] - 0.1809016993e1 * t * (0.1708203931e1 + 0.618033988e1 * t * t - 0.788854382e1 * t) * w0[1] + 0.6909830058e0 * t * (0.1170820393e2 + 0.1618033988e2 * t * t - 0.2788854382e2 * t) * w0[2] - 0.1000000000e1 * (-0.1e1 + 0.6e1 * t - 0.10e2 * t * t + 0.5e1 * pow(t, 0.3e1)) * w0[3];
      break;

    case 4 :
      w[0] = 0.9999999983e0 * (0.1e1 + 0.10e2 * t + 0.14e2 * pow(t, 0.4e1) + 0.35e2 * pow(t, 0.3e1) + 0.30e2 * t * t) * w0[0] + 0.6756502479e1 * (-0.2e1 - 0.8417424305e1 * t - 0.1125227292e2 * t * t - 0.4834848610e1 * pow(t, 0.3e1)) * t * w0[1] + 0.5333333325e1 * t * (0.1e1 + 0.8e1 * t + 0.7e1 * pow(t, 0.3e1) + 0.14e2 * t * t) * w0[2] - 0.1410164176e1 * (0.2e1 + 0.1758257570e2 * t + 0.3874772708e2 * t * t + 0.2316515139e2 * pow(t, 0.3e1)) * t * w0[3] + 0.9999999983e0 * t * (0.1e1 + 0.9e1 * t + 0.14e2 * pow(t, 0.3e1) + 0.21e2 * t * t) * w0[4];
      w[1] = 0.7142857131e-1 * t * (0.3474772708e2 + 0.2012340896e3 * t + 0.3546242389e3 * t * t + 0.196e3 * pow(t, 0.3e1)) * w0[0] + 0.1378878057e0 * (0.725227292e1 - 0.1692197014e3 * t * t - 0.387731045e3 * pow(t, 0.3e1) - 0.2369075819e3 * pow(t, 0.4e1)) * w0[1] + 0.7619047607e0 * t * (0.14e2 * t + 0.49e2 * pow(t, 0.3e1) - 0.4582575695e1 + 0.6415605973e2 * t * t) * w0[2] - 0.2014520251e0 * (-0.7582575695e1 + 0.1158257570e2 * t + 0.1592340896e3 * t * t + 0.1621560597e3 * pow(t, 0.3e1)) * t * w0[3] + 0.7142857131e-1 * t * (-0.725227292e1 + 0.876591040e1 * t + 0.1586242389e3 * t * t + 0.196e3 * pow(t, 0.3e1)) * w0[4];
      w[2] = 0.2499999996e0 * t * (-0.3e1 - 0.6e1 * t + 0.56e2 * pow(t, 0.3e1) + 0.28e2 * t * t) * w0[0] + 0.1689125620e1 * (0.1582575695e1 + 0.4834848610e1 * t - 0.633030278e1 * t * t - 0.1933939444e2 * pow(t, 0.3e1)) * t * w0[1] + 0.3333333329e0 * (0.3e1 + 0.112e3 * pow(t, 0.4e1) - 0.40e2 * t * t) * w0[2] - 0.3525410440e0 * (0.7582575695e1 - 0.2316515139e2 * t - 0.3033030278e2 * t * t + 0.9266060556e2 * pow(t, 0.3e1)) * t * w0[3] + 0.2499999996e0 * t * (0.3e1 - 0.6e1 * t + 0.56e2 * pow(t, 0.3e1) - 0.28e2 * t * t) * w0[4];
      w[3] = -0.7142857131e-1 * t * (-0.725227292e1 - 0.876591040e1 * t + 0.1586242389e3 * t * t - 0.196e3 * pow(t, 0.3e1)) * w0[0] + 0.9652146398e0 * (-0.1582575695e1 - 0.2417424305e1 * t + 0.3323408960e2 * t * t - 0.3384394027e2 * pow(t, 0.3e1)) * t * w0[1] - 0.7619047607e0 * t * (-0.14e2 * t - 0.49e2 * pow(t, 0.3e1) - 0.4582575695e1 + 0.6415605973e2 * t * t) * w0[2] - 0.2877886073e-1 * (-0.3474772708e2 + 0.8107802986e3 * t * t - 0.1857731045e4 * pow(t, 0.3e1) + 0.1135092418e4 * pow(t, 0.4e1)) * w0[3] - 0.7142857131e-1 * t * (0.3474772708e2 - 0.2012340896e3 * t + 0.3546242389e3 * t * t - 0.196e3 * pow(t, 0.3e1)) * w0[4];
      w[4] = 0.9999999983e0 * t * (-0.1e1 - 0.21e2 * t * t + 0.9e1 * t + 0.14e2 * pow(t, 0.3e1)) * w0[0] - 0.1166666665e1 * t * (-0.2417424305e1 + 0.2125227292e2 * t - 0.4683484861e2 * t * t + 0.28e2 * pow(t, 0.3e1)) * w0[1] + 0.5333333325e1 * t * (-0.1e1 + 0.8e1 * t + 0.7e1 * pow(t, 0.3e1) - 0.14e2 * t * t) * w0[2] + 0.1166666665e1 * t * (0.1158257570e2 - 0.4874772708e2 * t + 0.6516515139e2 * t * t - 0.28e2 * pow(t, 0.3e1)) * w0[3] + 0.9999999983e0 * (0.1e1 - 0.10e2 * t + 0.30e2 * t * t - 0.35e2 * pow(t, 0.3e1) + 0.14e2 * pow(t, 0.4e1)) * w0[4];
      break;

    default :
      assert(1==2);

    }

    for(int ipg=0;ipg <_D+1;++ipg){
      wcell[iv+_M*ipg]=w[ipg];
    }
   


    
 
  }
}


//return glop weight i
double wglop(int deg,int i) {
  return gauss_lob_weight[gauss_lob_offset[deg] + i];
}

// returns glop i
double glop(int deg,int i){
  return gauss_lob_point[i+gauss_lob_offset[deg]];
}

// return the 1d derivative of lagrange polynomial ib at glop ipg
double dlag(int deg,int ib,int ipg) 
{
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg]+ib*(deg+1)+ipg];
}




#include "macromesh.h"
#include "geometry.h"
#include "global.h"
#include "interpolation.h"
#include "test.h"
//#include "model.h"
//#include "field.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>


int main(void) {
  
  // unit tests
    
  int resu=TestInterpolation();
	 

  if (resu) printf("Interpolation test OK !\n");
  else printf("Interpolation test failed !\n");

  return !resu;
} 








int TestInterpolation(void){
  
  int test= (1==1);

  assert(_DEGX < 5);
  assert(_DEGY < 5);
  assert(_DEGZ < 5);

  // reference element
  double physnode[20][3];
  physnode[0][0] = 0.2;
    physnode[0][1] = 0;
    physnode[0][2] = 0;
    physnode[1][0] = 1;
    physnode[1][1] = 0;
    physnode[1][2] = 0;
    physnode[2][0] = 1;
    physnode[2][1] = 1;
    physnode[2][2] = 0;
    physnode[3][0] = 0;
    physnode[3][1] = 1;
    physnode[3][2] = 0;
    physnode[4][0] = 0;
    physnode[4][1] = 0;
    physnode[4][2] = 1;
    physnode[5][0] = 1;
    physnode[5][1] = 0;
    physnode[5][2] = 1;
    physnode[6][0] = 1;
    physnode[6][1] = 1;
    physnode[6][2] = 1;
    physnode[7][0] = 0;
    physnode[7][1] = 1;
    physnode[7][2] = 1;
    physnode[8][0] = 0.1e1 / 0.2e1;
    physnode[8][1] = 0;
    physnode[8][2] = 0;
    physnode[9][0] = 0;
    physnode[9][1] = 0.1e1 / 0.2e1;
    physnode[9][2] = 0;
    physnode[10][0] = 0;
    physnode[10][1] = 0;
    physnode[10][2] = 0.1e1 / 0.2e1;
    physnode[11][0] = 1;
    physnode[11][1] = 0.1e1 / 0.2e1;
    physnode[11][2] = 0;
    physnode[12][0] = 1;
    physnode[12][1] = 0;
    physnode[12][2] = 0.1e1 / 0.2e1;
    physnode[13][0] = 0.1e1 / 0.2e1;
    physnode[13][1] = 1;
    physnode[13][2] = 0;
    physnode[14][0] = 1;
    physnode[14][1] = 1;
    physnode[14][2] = 0.1e1 / 0.2e1;
    physnode[15][0] = 0;
    physnode[15][1] = 1;
    physnode[15][2] = 0.1e1 / 0.2e1;
    physnode[16][0] = 0.1e1 / 0.2e1;
    physnode[16][1] = 0;
    physnode[16][2] = 1;
    physnode[17][0] = 0;
    physnode[17][1] = 0.1e1 / 0.2e1;
    physnode[17][2] = 1;
    physnode[18][0] = 1;
    physnode[18][1] = 0.1e1 / 0.2e1;
    physnode[18][2] = 1;
    physnode[19][0] = 0.1e1 / 0.2e1;
    physnode[19][1] = 1;
    physnode[19][2] = 1;

    // transform the ref element with an affine map
    // for more generality

    for(int i=0;i<20;i++){
      AffineMap(&(physnode[i][0]));
    }


  // store the interpolation properties
  // 3 orders
  // 3 refinement parameters
  // 1 integer for converting face pg id into a volume pg id
  int deg[7];
  int nraf[3];


  // test that the ref_ipg function
  // is compatible with ref_pg_vol
  for(int d=1;d<5;d++){
    printf("Degree=%d\n",d);
    
    deg[0]=d;
    deg[1]=d;
    deg[2]=d;
    deg[3]=1;
    deg[4]=1;
    deg[5]=1;
    
    nraf[0]=deg[3];
    nraf[1]=deg[4];
    nraf[2]=deg[5];

    double xref1[3],xref2[3],xphy[3];

    for(int ipg=0;ipg<NPG(deg);ipg++){
      double wpg;
      ref_pg_vol(deg,ipg,xref1,&wpg);
      Ref2Phy(physnode,
             xref1,
             0,
             -1,
             xphy,
             0,
             0,
             0,
             0);
      Phy2Ref(physnode,xphy,xref2);
      test=test &&(ipg==ref_ipg(deg,xref2));
      printf("ipg=%d ipg2=%d\n",ipg,ref_ipg(deg,xref2));
      assert(test);
    }
  }

#define _MAXDEGTEST 5

  // test that the ref_pg_face function
  // is compatible with ref_pg_vol
  for(int d=1;d<_MAXDEGTEST;d++){
    printf("Degree=%d\n",d);
    
    deg[0]=d;
    deg[1]=d;
    deg[2]=d;
    deg[3]=1;
    deg[4]=1;
    deg[5]=1;
    
    nraf[0]=deg[3];
    nraf[1]=deg[4];
    nraf[2]=deg[5];

    double xref1[3],xref2[3],xphy[3];
    for(int ifa=0;ifa<6;ifa++){
      for(int ipgf=0;ipgf<NPGF(deg,ifa);ipgf++){
	double wpg;
	ref_pg_face(deg,ifa,ipgf,xref1,&wpg);
	int ipg=deg[6];
	Ref2Phy(physnode,
		xref1,
		0,
		-1,
		xphy,
		0,
		0,
		0,
		0);
	Phy2Ref(physnode,xphy,xref2);
	test=test &&(ipg==ref_ipg(deg,xref2));
	assert(test);
      }
    }
  }


  // test green formula for Gauss-Lobatto points


  for(int d0=1;d0<_MAXDEGTEST;d0++){
  for(int d1=1;d1<_MAXDEGTEST;d1++){
  for(int d2=1;d2<_MAXDEGTEST;d2++){
    printf("Degree=%d %d %d\n",d0,d1,d2);

    deg[0]=d0;
    deg[1]=d1;
    deg[2]=d2;
    deg[3]=3;
    deg[4]=1;
    deg[5]=2;

    nraf[0]=deg[3];
    nraf[1]=deg[4];
    nraf[2]=deg[5];

    int npg;      

    // check that integration by parts works with
    // two arbitrary polynomials

    npg = NPG(deg);
    double f[npg];
    double g[npg];
    double xref[3],omega;


    // Define test functions f and g
    for(int ipg=0;ipg<npg;ipg++){
      ref_pg_vol(deg, ipg, xref, &omega);
      double xphy[3];
      double dtau[3][3],codtau[3][3];
      Ref2Phy(physnode,xref,NULL,-1,
              xphy,dtau,codtau,NULL,NULL);
      f[ipg] = pow(xphy[1],(double) d1);
      //printf("%d %f\n",ipg,f[ipg]);
      g[ipg] =   pow(xphy[0],(double) d0);//xphy[0] ;
    }

    // Computation of the two integrants of green's formula
    // with null condition on the boundary
    // k corresponds to x, y or z
    int k = 1;
    double int_dfg = 0;
    double int_fdg = 0;
    // Loop on  Gauss-Lobatto points for the computation
    // of derivative of f against g
    for(int ipg=0;ipg<npg;ipg++){
      double dkf = 0;
      double dphiref[3];
      double dphi[3];
      //printf(" ipg= %d \n", ipg);
      for(int l=0;l<npg;l++){
	grad_psi_pg(deg,l,ipg,dphiref);
	ref_pg_vol(deg, ipg, xref, &omega);
        double xphy[3];
        double dtau[3][3],codtau[3][3];
        Ref2Phy(physnode,xref,dphiref,-1,
              xphy,dtau,codtau,dphi,NULL);
	//printf("ipg= %d xphy=%f\n", l,xphy[0]);
        double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	  dtau[0][2]*codtau[0][2];
	dkf += f[l] * dphi[k]/det;
		
      }
      ref_pg_vol(deg, ipg, xref, &omega);
      double xphy[3];
      double dtau[3][3],codtau[3][3];
      Ref2Phy(physnode,xref,NULL,-1,
              xphy,dtau,codtau,NULL,NULL);
      //printf(" x=%f f=%f dkf= %f\n",xphy[0],f[ipg],dkf);
      double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	dtau[0][2]*codtau[0][2];
      //printf("det=%f\n",det);
      int_dfg += dkf * g[ipg] * omega * det;
    }
    // Loop on  Gauss-Lobatto points for the computation
    // of f against derivative of g
    double sum_wpg = 0;
    for(int ipg=0;ipg<npg;ipg++){
      double dkg = 0;
      double dphi[3],dphiref[3];
      for(int l=0;l<npg;l++){
	grad_psi_pg(deg,l,ipg,dphiref);
	ref_pg_vol(deg, ipg, xref, &omega);
        double xphy[3];
        double dtau[3][3],codtau[3][3];
        Ref2Phy(physnode,xref,dphiref,-1,
              xphy,dtau,codtau,dphi,NULL);
	//printf("ipg= %d xphy=%f\n", l,xphy[0]);
        double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	  dtau[0][2]*codtau[0][2];
	dkg += dphi[k] * g[l]/det;
      }
      ref_pg_vol(deg, ipg, xref, &omega);
      double xphy[3];
      double dtau[3][3],codtau[3][3];
      Ref2Phy(physnode,xref,NULL,-1,
              xphy,dtau,codtau,NULL,NULL);
      //printf("  dkf = %.2e dkf_ex = %.2e\n", dkf, pow(xref[0],d-1)*d);
        double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	  dtau[0][2]*codtau[0][2];
      int_fdg += f[ipg] * dkg * omega *det;
      sum_wpg+=omega;
    }

    printf("  derivative of f against g = %f\n", int_dfg);
    printf("  f against derivative of g = %f\n", int_fdg);
 

    //assert(1==2);
    double int_fgn = 0;
    // add the boundary term: integral on the faces
    // of f * g * normal  ds
    for(int ifa=0;ifa<6;ifa++){
      double xphy[3],vnds[3];
      double dtau[3][3],codtau[3][3];
      int npg_f = NPGF(deg,ifa);
      for(int ipgf=0;ipgf<npg_f;ipgf++){
	ref_pg_face(deg,ifa,ipgf,xref,&omega);
 	
	Ref2Phy(physnode,xref,NULL,ifa,
                xphy,dtau,codtau,NULL,vnds);
	int_fgn += f[deg[6]] * g[deg[6]] * vnds[k] * omega;
      }
    }

    test = (test && fabs(int_dfg+int_fdg-int_fgn)<1e-8);

    printf("  int_fgn = %.2e sum of ints = %.2e sum_wpg=%f \n", int_fgn,int_dfg+int_fdg-int_fgn,sum_wpg);
    
  }}}


  return test;

}

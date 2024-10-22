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

int TestInterpolation(void);

int main(void) {
  // unit tests
  int resu=TestInterpolation();
  if (resu) printf("Interpolation test OK !\n");
  else printf("Interpolation test failed !\n");
  return !resu;
} 

int TestInterpolation(void){
  int test = true;
  // reference element
  schnaps_real physnode[20][3];
  physnode[0][0] = 0.;
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
  //real A[3][3]={{1,2,1},{0,-1,4},{7,8,-5}};
  //real x0[3]={0,0,0};
  schnaps_real A[3][3]={{1,0,0},{0,1,0},{0,0,1}};
  schnaps_real x0[3]={0,0,0};
  for(int i = 0; i < 20; i++){
    AffineMap(&(physnode[i][0]),A,x0);
  }

  // store the interpolation properties
  // 3 orders
  // 3 refinement parameters
  // 1 integer for converting face pg id into a volume pg id
  int deg[7];
  int nraf[3];

  // test that the ref_ipg function
  // is compatible with ref_pg_vol

  deg[0]=1;
  deg[1]=1;
  deg[2]=1;
  
  nraf[0]=1;
  nraf[1]=2;
  nraf[2]=1;

  schnaps_real xref[3], xin[3], wpg;

  for(int ipg=0; ipg < NPG_CG(deg, nraf); ipg++){
    ref_pg_vol_CG(deg, nraf, ipg, xref, &wpg, xin);

    printf("ipg=%d xref=%f %f %f xin=%f %f %f wpg=%f \n",ipg,xref[0],xref[1],xref[2],
	   xin[0],xin[1],xin[2],
	   wpg);
  }

  assert(NPG_CG(deg, nraf) == 12);

  assert(fabs(xref[1]-1) < _SMALL);
  
  for(int d = 1; d < 5; d++){
    printf("Degree=%d\n", d);
    
    deg[0] = d;
    deg[1] = d;
    deg[2] = d;
    deg[3] = 1;
    deg[4] = 2;
    deg[5] = 1;
    
    nraf[0] = deg[3];
    nraf[1] = deg[4];
    nraf[2] = deg[5];

    schnaps_real xref1[3], xref2[3], xphy[3], xref_in[3];

    for(int ipg = 0; ipg < NPG(deg,nraf); ipg++){
      schnaps_real wpg;
      ref_pg_vol(deg, nraf, ipg, xref1, &wpg, xref_in);
      //printf("xref_in %f %f %f \n",xref_in[0],xref_in[1],xref_in[2]);
      schnaps_ref2phy(physnode,
	      xref_in,
	      0,
	      -1,
	      xphy,
	      0,
	      0,
	      0,
	      0);
   
      schnaps_phy2ref(physnode, xphy, xref2);
      test = test && (ipg == ref_ipg(deg, nraf, xref2));
      //printf("ipg=%d ipg2=%d\n",ipg,ref_ipg(deg,xref2));
      assert(test);
    }

    for(int ipg = 0; ipg < NPG_CG(deg,nraf); ipg++){
      schnaps_real wpg;
      ref_pg_vol_CG(deg, nraf, ipg, xref1, &wpg, xref_in);
      //printf("xref_in %f %f %f \n",xref_in[0],xref_in[1],xref_in[2]);
      schnaps_ref2phy(physnode,
	      xref_in,
	      0,
	      -1,
	      xphy,
	      0,
	      0,
	      0,
	      0);
   
      schnaps_phy2ref(physnode, xphy, xref2);
      test = test && (ipg == ref_ipg_CG(deg, nraf, xref2));
      //printf("ipg=%d ipg2=%d\n",ipg,ref_ipg(deg,xref2));
      assert(test);
    }
  }

  //}

  


#define _MAXDEGTEST 4

  // test that ipg_to_xyz and xyz_to_ipg are inverses of each other

#define _DMAX 5

  for(int d0=0;d0<_DMAX;d0++){
    for(int d1=0;d1<_DMAX;d1++){
      for(int d2=0;d2<_DMAX;d2++){
	for(int r0=0;r0<_DMAX;r0++){
	  for(int r1=0;r1<_DMAX;r1++){
	    for(int r2=0;r2<_DMAX;r2++){
	      int deg[6]={d0,d1,d2,r0,r1,r2};
	      int* raf=deg+3;
	      for (int ipg=0;ipg<NPG(deg,raf);ipg++){
		int ix[3],ic[3];
		int ipgcopy=ipg;
		//printf("avant %d %d %d npg=%d \n",d0,d1,d2,NPG(deg));
		ipg_to_xyz(raf,deg,ic,ix,&ipgcopy);
		//printf("apres %d %d %d \n",d0,d1,d2);
		xyz_to_ipg(raf,deg,ic,ix,&ipgcopy);
		assert(ipg==ipgcopy);
	      }
	    }
	  }
	}
      }
    }
  }

    
    

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

    schnaps_real xref1[3],xref2[3],xphy[3],xref_in[3];
    for(int ifa=0;ifa<6;ifa++){
      for(int ipgf=0;ipgf<NPGF(deg, nraf, ifa);ipgf++){
	schnaps_real wpg;
	int ipg = ref_pg_face(deg, nraf,ifa,ipgf,xref1,&wpg,NULL);
	ref_pg_vol(deg, nraf,ipg,xref1,&wpg,xref_in);
	schnaps_ref2phy(physnode,
		xref_in,
		0,
		-1,
		xphy,
		0,
		0,
		0,
		0);
	schnaps_phy2ref(physnode,xphy,xref2);
	test=test &&(ipg==ref_ipg(deg,nraf,xref2));
	assert(test);
      }
    }
  } // end degree loop


  // test that the ref_pg_face_CG function
  // is compatible with ref_pg_vol_CG
  for(int d=1;d<_MAXDEGTEST;d++){
    printf("CG test ref_pg_face Degree=%d\n",d);
    
    deg[0]=d;
    deg[1]=d;
    deg[2]=d;
    
    nraf[0]=1;
    nraf[1]=2;
    nraf[2]=1;

    schnaps_real xref1[3],xref2[3],xphy[3],xref_in[3];
    for(int ifa=0;ifa<6;ifa++){
      for(int ipgf=0;ipgf<NPGF_CG(deg, nraf, ifa);ipgf++){
	schnaps_real wpg;
	int ipg = ref_pg_face_CG(deg, nraf,ifa,ipgf,xref1,&wpg,NULL);
	ref_pg_vol_CG(deg, nraf,ipg,xref1,&wpg,xref_in);
	schnaps_ref2phy(physnode,
		xref_in,
		0,
		-1,
		xphy,
		0,
		0,
		0,
		0);
	schnaps_phy2ref(physnode,xphy,xref2);
	test=test && (ipg == ref_ipg_CG(deg, nraf, xref2));
	printf("ipgf=%d ipg=%d ref_ipg=%d\n\n",ipgf,ipg,ref_ipg_CG(deg, nraf, xref2));
	assert(test);
      }
    }
  } // end degree loop


  //assert(1==2);

  // test green formula for Gauss-Lobatto points


  for(int d0=1;d0<_MAXDEGTEST;d0++){
    for(int d1=1;d1<_MAXDEGTEST;d1++){
      for(int d2=1;d2<_MAXDEGTEST;d2++){
	printf("Degree=%d %d %d\n",d0,d1,d2);

	deg[0]=d0;
	deg[1]=d1;
	deg[2]=d2;
	deg[3]=3;
	deg[4]=2;
	deg[5]=1;

	nraf[0]=deg[3];
	nraf[1]=deg[4];
	nraf[2]=deg[5];

	int npg;      

	// check that integration by parts works with
	// two arbitrary polynomials

	npg = NPG(deg,nraf);
	schnaps_real f[npg];
	schnaps_real g[npg];
	schnaps_real xref[3],omega;


	// Define test functions f and g
	for(int ipg = 0; ipg < npg; ipg++){
	  ref_pg_vol(deg, nraf, ipg, xref, &omega,NULL);
	  schnaps_real xphy[3];
	  schnaps_real dtau[3][3], codtau[3][3];
	  schnaps_ref2phy(physnode, xref, NULL, -1,
		  xphy, dtau, codtau, NULL, NULL);
	  //f[ipg] = pow(xphy[1],(real) d1);
	  f[ipg] = pow(xref[1], (schnaps_real) d1);
	  //printf("%d %f\n",ipg,f[ipg]);
	  //g[ipg] = pow(xphy[0],(real) d0)*pow(xphy[1],(real) d1);//xphy[0] ;
	  g[ipg] = pow(xref[0], (schnaps_real) d0) * pow(xref[1], 2.0 * d1);//xphy[0] ;
	}

	// Computation of the two integrants of green's formula
	// k corresponds to x, y or z
	int k = 2;
	schnaps_real int_dfg = 0;
	schnaps_real int_fdg = 0;
	// Loop on  Gauss-Lobatto points for the computation
	// of derivative of f against g
	for(int ipg = 0; ipg < npg; ipg++){
	  schnaps_real dkf = 0;
	  schnaps_real dphiref[3];
	  schnaps_real dphi[3];
	  //printf(" ipg= %d \n", ipg);
	  for(int l = 0; l < npg; l++){
	    grad_psi_pg(deg, nraf, l, ipg, dphiref);
	    ref_pg_vol(deg, nraf, ipg, xref, &omega, NULL);
	    schnaps_real xphy[3];
	    schnaps_real dtau[3][3], codtau[3][3];
	    schnaps_ref2phy(physnode, xref, dphiref, -1,
		    xphy, dtau, codtau, dphi, NULL);
	    //printf("ipg= %d xphy=%f\n", l,xphy[0]);
	    schnaps_real det 
	      = dtau[0][0] * codtau[0][0]
	      + dtau[0][1] * codtau[0][1]
	      + dtau[0][2] * codtau[0][2];
	    dkf += f[l] * dphi[k] / det;
		
	  }
	  ref_pg_vol(deg, nraf, ipg, xref, &omega,NULL);
	  schnaps_real xphy[3];
	  schnaps_real dtau[3][3], codtau[3][3];
	  schnaps_ref2phy(physnode, xref, NULL, -1,
		  xphy, dtau, codtau, NULL, NULL);
	  //printf(" x=%f f=%f dkf= %f\n",xphy[0],f[ipg],dkf);
	  schnaps_real det
	    = dtau[0][0] * codtau[0][0]
	    + dtau[0][1] * codtau[0][1]
	    + dtau[0][2] * codtau[0][2];
	  //printf("det=%f\n",det);
	  int_dfg += dkf * g[ipg] * omega * det;
	}
	// Loop on  Gauss-Lobatto points for the computation
	// of f against derivative of g
	schnaps_real sum_wpg = 0;
	for(int ipg = 0;ipg < npg; ipg++){
	  schnaps_real dkg = 0;
	  schnaps_real dphi[3], dphiref[3];
	  for(int l = 0; l < npg; l++){
	    grad_psi_pg(deg, nraf, l, ipg, dphiref);
	    ref_pg_vol(deg, nraf, ipg, xref, &omega,NULL);
	    schnaps_real xphy[3];
	    schnaps_real dtau[3][3], codtau[3][3];
	    schnaps_ref2phy(physnode, xref, dphiref, -1,
		    xphy, dtau, codtau, dphi, NULL);
	    //printf("ipg= %d xphy=%f\n", l,xphy[0]);
	    schnaps_real det
	      = dtau[0][0] * codtau[0][0]
	      + dtau[0][1] * codtau[0][1]
	      + dtau[0][2] * codtau[0][2];
	    dkg += dphi[k] * g[l] / det;
	  }
	  ref_pg_vol(deg, nraf, ipg, xref, &omega,NULL);
	  schnaps_real xphy[3];
	  schnaps_real dtau[3][3], codtau[3][3];
	  schnaps_ref2phy(physnode, xref, NULL, -1,
		  xphy, dtau, codtau, NULL, NULL);
	  //printf("  dkf = %.2e dkf_ex = %.2e\n", dkf, pow(xref[0],d-1)*d);
	  schnaps_real det
	    = dtau[0][0] * codtau[0][0]
	    + dtau[0][1] * codtau[0][1]
	    + dtau[0][2] * codtau[0][2];
	  int_fdg += f[ipg] * dkg * omega *det;
	  sum_wpg += omega;
	}

	printf("  derivative of f against g = %f\n", int_dfg);
	printf("  f against derivative of g = %f\n", int_fdg);
 
	schnaps_real int_fgn = 0;
	// add the boundary term: integral on the faces
	// of f * g * normal  ds
	for(int ifa=0;ifa<6;ifa++){
	  schnaps_real xphy[3], vnds[3];
	  schnaps_real dtau[3][3], codtau[3][3];
	  int npg_f = NPGF(deg, nraf, ifa);
	  for(int ipgf = 0; ipgf < npg_f; ipgf++){
	    int ninvol = ref_pg_face(deg, nraf, ifa, ipgf, xref, &omega, NULL);
 	
	    schnaps_ref2phy(physnode, xref, NULL, ifa,
		    xphy, dtau, codtau, NULL, vnds);
	    int_fgn += f[ninvol] * g[ninvol] * vnds[k] * omega;
	  }
	}

	test = (test && fabs(int_dfg+int_fdg-int_fgn)<_SMALL);
	printf("  int_fgn = %.2e sum of ints = %.2e sum_wpg=%f \n", 
	       int_fgn, int_dfg + int_fdg - int_fgn, sum_wpg);
      }
    }
  }
  return test;
}

#include "macromesh.h"
#include "geometry.h"
#include "interpolation.h"
#include "test.h"
#include "model.h"
#include "field.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

// some unit tests of the macromesh code
int TestMacroMesh(void){

  int test=1;
  MacroMesh m;
  
  ReadMacroMesh(&m,"../geo/testmacromesh.msh");
  BuildConnectivity(&m);
  PrintMacroMesh(&m);

  test = (m.nbelems == 5);
  test =  (test && m.nbnodes == 50);

  return test;

}



int TestGeometry(void){

  int test=1;
  MacroMesh mc;

  ReadMacroMesh(&mc,"../geo/testgeometry.msh");
  BuildConnectivity(&mc);
  PrintMacroMesh(&mc);

  // test the geometric transformation

  double xref[3]={0.1,0.3,0.7};
  double v[3]={0.1,0.3,0.7};
  double physnode[20*3];

  double xphy[3],dtau[9];

  for(int inoloc=0;inoloc<20;inoloc++){
    int ino=mc.elem2node[0*20+inoloc];
    physnode[3*inoloc+0]=mc.node[3*ino+0]; //x
    physnode[3*inoloc+1]=mc.node[3*ino+1]; //y
    physnode[3*inoloc+2]=mc.node[3*ino+2]; //z
  }

  Ref2Phy(physnode,xref,0,-1,xphy,dtau,0,0,0);

  printf("xphy= %f %f %f \n",xphy[0],xphy[1],xphy[2]);

  Phy2Ref(physnode,xphy,xref);

  printf("xref= %f %f %f \n",xref[0],xref[1],xref[2]);

  v[0]-=xref[0];
  v[1]-=xref[1];
  v[2]-=xref[2];

  double d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  test = (d < 1e-12);

  return test;

}

int TestInterpolation(void){
  
  int test= (1==1);

  // test green formula for Gauss-Lobatto points

  // store the interpolation properties
  // 3 orders
  // 3 refinement parameters
  // 1 integer for converting face pg id into a volume pg id
  int deg[7];
  int nraf[3];

  for(int d=1;d<5;d++){
    printf("Degree=%d\n",d);

    deg[0]=d;
    deg[1]=d;
    deg[2]=d;
    deg[3]=3;
    deg[4]=1;
    deg[5]=2;

    nraf[0]=deg[3];
    nraf[1]=deg[4];
    nraf[2]=deg[5];

    int npg;

    // Reference element
    double physnode[20*3];
    physnode[0*3+0] = 0;
    physnode[0*3+1] = 0;
    physnode[0*3+2] = 0;
    physnode[1*3+0] = 1;
    physnode[1*3+1] = 0;
    physnode[1*3+2] = 0;
    physnode[2*3+0] = 1;
    physnode[2*3+1] = 1;
    physnode[2*3+2] = 0;
    physnode[3*3+0] = 0;
    physnode[3*3+1] = 1;
    physnode[3*3+2] = 0;
    physnode[4*3+0] = 0;
    physnode[4*3+1] = 0;
    physnode[4*3+2] = 1;
    physnode[5*3+0] = 1;
    physnode[5*3+1] = 0;
    physnode[5*3+2] = 1;
    physnode[6*3+0] = 1;
    physnode[6*3+1] = 1;
    physnode[6*3+2] = 1;
    physnode[7*3+0] = 0;
    physnode[7*3+1] = 1;
    physnode[7*3+2] = 1;
    physnode[8*3+0] = 0.1e1 / 0.2e1;
    physnode[8*3+1] = 0;
    physnode[8*3+2] = 0;
    physnode[9*3+0] = 0;
    physnode[9*3+1] = 0.1e1 / 0.2e1;
    physnode[9*3+2] = 0;
    physnode[10*3+0] = 0;
    physnode[10*3+1] = 0;
    physnode[10*3+2] = 0.1e1 / 0.2e1;
    physnode[11*3+0] = 1;
    physnode[11*3+1] = 0.1e1 / 0.2e1;
    physnode[11*3+2] = 0;
    physnode[12*3+0] = 1;
    physnode[12*3+1] = 0;
    physnode[12*3+2] = 0.1e1 / 0.2e1;
    physnode[13*3+0] = 0.1e1 / 0.2e1;
    physnode[13*3+1] = 1;
    physnode[13*3+2] = 0;
    physnode[14*3+0] = 1;
    physnode[14*3+1] = 1;
    physnode[14*3+2] = 0.1e1 / 0.2e1;
    physnode[15*3+0] = 0;
    physnode[15*3+1] = 1;
    physnode[15*3+2] = 0.1e1 / 0.2e1;
    physnode[16*3+0] = 0.1e1 / 0.2e1;
    physnode[16*3+1] = 0;
    physnode[16*3+2] = 1;
    physnode[17*3+0] = 0;
    physnode[17*3+1] = 0.1e1 / 0.2e1;
    physnode[17*3+2] = 1;
    physnode[18*3+0] = 1;
    physnode[18*3+1] = 0.1e1 / 0.2e1;
    physnode[18*3+2] = 1;
    physnode[19*3+0] = 0.1e1 / 0.2e1;
    physnode[19*3+1] = 1;
    physnode[19*3+2] = 1;


    // transform the ref element with an affine map
    // for more generality

    for(int i=0;i<20;i++){
      AffineMap(&(physnode[i*3]));
    }
      

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
      double dtau[9],codtau[9];
      Ref2Phy(physnode,xref,NULL,0,
              xphy,dtau,codtau,NULL,NULL);
      f[ipg] = pow(xphy[1],(double) d);
      //printf("%d %f\n",ipg,f[ipg]);
      g[ipg] =   pow(xphy[0],(double) d-1);//xphy[0] ;
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
        double dtau[9],codtau[9];
        Ref2Phy(physnode,xref,dphiref,0,
              xphy,dtau,codtau,dphi,NULL);
	//printf("ipg= %d xphy=%f\n", l,xphy[0]);
        double det=dtau[0]*codtau[0]+dtau[1]*codtau[1]+dtau[2]*codtau[2];
	dkf += f[l] * dphi[k]/det;
		
      }
      ref_pg_vol(deg, ipg, xref, &omega);
      double xphy[3];
      double dtau[9],codtau[9];
      Ref2Phy(physnode,xref,NULL,0,
              xphy,dtau,codtau,NULL,NULL);
      printf(" x=%f f=%f dkf= %f\n",xphy[0],f[ipg],dkf);
      double det=dtau[0]*codtau[0]+dtau[1]*codtau[1]+dtau[2]*codtau[2];
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
        double dtau[9],codtau[9];
        Ref2Phy(physnode,xref,dphiref,0,
              xphy,dtau,codtau,dphi,NULL);
	//printf("ipg= %d xphy=%f\n", l,xphy[0]);
        double det=dtau[0]*codtau[0]+dtau[1]*codtau[1]+dtau[2]*codtau[2];
	dkg += dphi[k] * g[l]/det;
      }
      ref_pg_vol(deg, ipg, xref, &omega);
      double xphy[3];
      double dtau[9],codtau[9];
      Ref2Phy(physnode,xref,NULL,0,
              xphy,dtau,codtau,NULL,NULL);
      //printf("  dkf = %.2e dkf_ex = %.2e\n", dkf, pow(xref[0],d-1)*d);
      double det=dtau[0]*codtau[0]+dtau[1]*codtau[1]+dtau[2]*codtau[2];
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
      double dtau[9],codtau[9];
      int npg_f = NPGF(deg,ifa);
      for(int ipgf=0;ipgf<npg_f;ipgf++){
	ref_pg_face(deg,ifa,ipgf,xref,&omega);
        //printf("ifa=%d ipgf=%d xref=%f %f %f omega=%f \n",ifa,ipgf,
	//xref[0],xref[1],xref[2],omega);

	// if (ifa == 2 && d == 1 && ipgf == 4) {
	//   printf("ipgf=%d , ipgv=%d\n",ipgf,deg[6]);
	//   assert(1==2);
	// }
	
	Ref2Phy(physnode,xref,NULL,ifa,
                xphy,dtau,codtau,NULL,vnds);
	// printf("dtau = \n");
	// for(int ii=0;ii<3;ii++){
	//   for(int jj=0;jj<3;jj++){
        //     printf("%f ", dtau[3*ii+jj]);
	//   }
        //   printf("\n ");
	// }
	// double fpg = 0;
	// double gpg = 0;
	// double psi;
	// for(int l=0;l<npg;l++){
	//   hexa_subcell_elem_interp->psi_ref(deg,l,xref,&psi,NULL);
	//   fpg += f[l] * psi;
	//   gpg += g[l] * psi;
	// }
	//	fpg = pow(xref[0],d);
	//      f[ipg] = 1;
	//      g[ipg] = xref[1];
	//gpg = 2 * fpg;
	//int_fgn += fpg * gpg * vnds[k] * omega;
	// the call to has computed the vol pg index from ifa and ipgf
	// the result is stored in deg[6]
	//printf("ifa=%d ipgf=%d ipgv=%d omega=%f\n",ifa,ipgf,deg[6],omega);
	//printf("vnds=%f %f %f\n",vnds[0],vnds[1],vnds[2]);
	int_fgn += f[deg[6]] * g[deg[6]] * vnds[k] * omega;
      }
    }

    test = (test && fabs(int_dfg+int_fdg-int_fgn)<1e-8);

    printf("  int_fgn = %.2e sum of ints = %.2e sum_wpg=%f \n", int_fgn,int_dfg+int_fdg-int_fgn,sum_wpg);
    
  }


  return test;

}

void AffineMap(double* x){

  double A[3][3]={1,2,1,0,-1,4,7,8,-5};
  //double A[3][3]={1,0,0,0,2,0,0,0,1};
  double x0[3]={10,0,-4};
  //double x0[3]={0,0,0};

  double newx[3];

  for(int i=0;i<3;i++){
    newx[i]=x0[i];
    for(int j=0;j<3;j++){
      newx[i]+=A[i][j]*x[j];
    }
  }
  x[0]=newx[0];
  x[1]=newx[1];
  x[2]=newx[2];

}

int TestModel(void){

  int test= (1==1);
  // creation of a simple transport model
  Model tr;
  tr.m=1; // only one conservative variable
  tr.NumFlux=TransportNumFlux;
  tr.BoundaryFlux=TransportBoundaryFlux;
  tr.InitData=TransportInitData;
  tr.ImposedData=TransportImposedData;

  double wL[tr.m];
  double wR[tr.m];
  double flux1[tr.m],flux2[tr.m];

  double x[3]={1,1,2};
  double t=0;
  double vn[3]={1/sqrt(3),1/sqrt(3),-1/sqrt(3)};

  tr.InitData(x,wR);
  tr.NumFlux(wL,wR,vn,flux1);
  printf("f %f \n",flux1[0]);
  tr.BoundaryFlux(x,t,wL,vn,flux2);
  printf("f %f \n",flux2[0]);

  double err=fabs(flux2[0]-flux1[0]);

  test=(err < 1e-8);

  return test;

};

int TestField(void){

  int test = (1==1);

  Field f;
  f.model.m=1; // only one conservative variable
  f.model.NumFlux=TransportNumFlux;
  f.model.BoundaryFlux=TransportBoundaryFlux;
  f.model.InitData=TransportInitData;
  f.model.ImposedData=TransportImposedData;
  f.varindex=GenericVarindex;

  ReadMacroMesh(&(f.macromesh),"../geo/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));

  InitField(&f);

  DisplayField(&f,"testvisufield.msh");
  
  return test;


};

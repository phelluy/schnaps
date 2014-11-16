#pragma OPENCL EXTENSION cl_khr_fp64: enable

//#define _FPREC
// for single precision computations
//#define _FPREC f

//! Gauss LObatto Points (GLOP) up to order 4
__constant double gauss_lob_point[] = {
  0.5,
  0,
  1,
  0,
  0.5,
  1,
  0,
  0.276393202250021030359082633127,
  0.723606797749978969640917366873,
  1,
  0,
  0.172673164646011428100853771877, 
  0.5, 
  0.827326835353988571899146228123,
  1
};

//! GLOP weights up to order 4
__constant double gauss_lob_weight[] = {
  1,
  0.5,
  0.5,
  0.166666666666666666666666666667,
  0.666666666666666666666666666668,
  0.166666666666666666666666666667,
  0.0833333333333333333333333333333,
  0.416666666666666666666666666666,
  0.416666666666666666666666666666,
  0.0833333333333333333333333333333,
  0.05,
  0.272222222222222222222222222223,
  0.355555555555555555555555555556,
  0.272222222222222222222222222219,
  0.05
};

//! indirection for finding the GLOP
//! data for a given degree in the previous arrays
__constant int gauss_lob_offset[] = {0, 1, 3, 6, 10};


__constant double gauss_lob_dpsi[] = {
  0.00000000000000000000000000000000000000000000000000000000000,
  -1.,
  -1.,
  1.,
  1.,
  -3.,
  -1.,
  1.,
  4.,
  0.,
  -4.,
  -1.,
  1.,
  3.,
  -6,
  -1.61803398874989484820458683436,
  .618033988749894848204586834362,
  -1,
  8.09016994374947424102293417177,
  0,
  -2.23606797749978969640917366872,
  3.09016994374947424102293417184,
  -3.09016994374947424102293417182,
  2.23606797749978969640917366872,
  0,
  -8.09016994374947424102293417177,
  1,
  -.618033988749894848204586834362,
  1.61803398874989484820458683436,
  6,
  -10,
  -2.48198050606196571569743868436,
  .75,
  -.518019493938034284302561315632,
  1,
  13.5130049774484800076860550594,
  0,
  -2.67316915539090667050969419631,
  1.52752523165194666886268239794,
  -2.82032835588485332564727827404,
  -5.33333333333333333333333333336,
  3.49148624377587810025755976667,
  0,
  -3.49148624377587810025755976662,
  5.33333333333333333333333333336,
  2.82032835588485332564727827399,
  -1.52752523165194666886268239791,
  2.67316915539090667050969419635,
  0,
  -13.5130049774484800076860550594,
  -1,
  .518019493938034284302561315631,
  -.75,
  2.48198050606196571569743868437,
  10
};

//! indirection for finding the GLOP
//! data for a given degree in the previous arrays
__constant int gauss_lob_dpsi_offset[] = {0, 1, 5, 14, 30};


double dlag(int deg,int ib,int ipg);
// return the 1d derivative of lagrange polynomial ib at glop ipg
double dlag(int deg,int ib,int ipg){

  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg]+ib*(deg+1)+ipg];

}

int varindex(__constant int* param, int elem, int ipg, int iv);

// memory location of w : component iv, macrocell elem and
// gauss point id in the macrocell ipg 
int varindex(__constant int* param, int elem, int ipg, int iv){

  int npg= (param[1]+1)*(param[2]+1)*(param[3]+1)
    *param[4]*param[5]*param[6];

  return iv + param[0] * ( ipg + npg * elem);

}

void NumFlux(double wL[],double wR[],double* vnorm,double* flux);

void NumFlux(double wL[],double wR[],double* vnorm,double* flux){
  

  double s2=0.707106781186547524400844362105;
  double vn =
    s2 * vnorm[0] +
    s2 * vnorm[1];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];


};




//!  \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
double wglop(int deg,int i);
double wglop(int deg,int i){
  return gauss_lob_weight[gauss_lob_offset[deg]+i];
}

void get_dtau(double x,double y,double z,
	      __constant double* physnode,double dtau[][3]);



// compute the volume terms on one macrocell
__kernel
void DGVolume(
	      __constant int* param,        // interp param
	      __constant int* ie,            // macrocel index
	      __constant double* physnode,  // macrocell nodes
              __global double* wn,       // field values
	      __global double* dtwn){       // time derivative
  
  const int m = param[0];
  const int deg[3]={param[1],
			 param[2],
			 param[3]};
  const int npg[3] = {deg[0]+1,
			   deg[1]+1,
			   deg[2]+1};
  const int nraf[3]={param[4],
			  param[5],
			  param[6]};  
  
  //const int nnpg=(param[1]+1)*(param[2]+1)*(param[3]+1) *
  // (param[4])*(param[5])*(param[6]);
  
  // subcell id
  int icell=get_group_id(0);
  int ic[3];
  ic[0]=icell%nraf[0];
  ic[1]=(icell/nraf[0])%nraf[1];
  ic[2]=icell/nraf[0]/nraf[1];
  
  
  // gauss point id where
  // we compute the jacobian
  int ipg=get_local_id(0);
  int p[3];
  p[0]=ipg%npg[0];
  p[1]=(ipg/npg[0])%npg[1];
  p[2]=ipg/npg[0]/npg[1];
  
  // ref coordinates
  double hx=1/(double) nraf[0];
  double hy=1/(double) nraf[1];
  double hz=1/(double) nraf[2];

  int offset[3];

  offset[0]=gauss_lob_offset[deg[0]]+p[0];
  offset[1]=gauss_lob_offset[deg[1]]+p[1];
  offset[2]=gauss_lob_offset[deg[2]]+p[2];

  double x=hx*(ic[0]+gauss_lob_point[offset[0]]);
  double y=hy*(ic[1]+gauss_lob_point[offset[1]]);
  double z=hz*(ic[2]+gauss_lob_point[offset[2]]);

  double wpg=hx*hy*hz*gauss_lob_weight[offset[0]]*
    gauss_lob_weight[offset[1]]*
    gauss_lob_weight[offset[2]];

  double dtau[3][3];
  get_dtau(x,y,z,physnode,dtau);
  
  double codtau[3][3];

  codtau[0][0] = dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
  codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
  codtau[0][2] = dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
  codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
  codtau[1][1] = dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
  codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
  codtau[2][0] = dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
  codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
  codtau[2][2] = dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];

  /* double det=dtau[0][0]*dtau[1][1]*dtau[2][2]-dtau[0][0]*dtau[1][2]*dtau[2][1]-dtau[1][0]*dtau[0][1]*dtau[2][2]+ */
  /*   dtau[1][0]*dtau[0][2]*dtau[2][1]+dtau[2][0]*dtau[0][1]*dtau[1][2]-dtau[2][0]*dtau[0][2]*dtau[1][1]; */

  /* printf("det=%f\n",det); */

#define _M 1  /// to do let schnaps specify m !!!!!!!!

  double wL[_M];
  for(int iv=0; iv < m; iv++){
    // gauss point id in the macrocell
    int ipgL=npg[0]*npg[1]*npg[2]*icell+p[0]+npg[0]*(p[1]+npg[1]*p[2]);
    int imemL=varindex(param,*ie,ipgL,iv);
    //int imemL= iv + m * ( get_global_id(0) + nnpg * *ie);
    wL[iv] = wn[imemL]; 
  }


  //for(int dim0 = 0; dim0 < 3; dim0++){
  for(int dim0 = 0; dim0 < 2; dim0++){
  int q[3]={p[0],p[1],p[2]};
    for(int iq = 0; iq < npg[dim0]; iq++){
      q[dim0]=(p[dim0]+iq)%npg[dim0];
      double dphiref[3]={0,0,0};
      dphiref[dim0]=dlag(deg[dim0],q[dim0],p[dim0])*nraf[dim0];
      double dphi[3];
      for(int ii=0;ii<3;ii++){
	dphi[ii]=0;
	for(int jj=0;jj<3;jj++){
	  dphi[ii]+=codtau[ii][jj]*dphiref[jj];
	}
      }
      double flux[_M];
      NumFlux(wL,wL,dphi,flux); // to do: let schnaps gives fluxnum

      for(int iv=0; iv < m; iv++){
        int ipgR=npg[0]*npg[1]*npg[2]*icell+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
	int imemR=varindex(param,*ie,ipgR,iv); // to do !
	int ipgL=npg[0]*npg[1]*npg[2]*icell+p[0]+npg[0]*(p[1]+npg[1]*p[2]);
	int imemL=varindex(param,*ie,ipgL,iv);
	if (imemL==7 && iq==iq){
	  printf("gpu ipgL=%d ipgR=%d iq=%d flux=%f\n",ipgL,ipgR,iq,flux[0]);
	  printf("gpu dim0=%d q=%d %d %d\n",dim0,q[0],q[1],q[2]);
	}

    	dtwn[imemR]+=flux[iv]*wpg;
      }


    }
  }

}

// apply division by the mass matrix on one macrocell
__kernel
void DGMass(
	    __constant int* param,        // interp param
            __constant int* ie,            // macrocel index
            __constant double* physnode,  // macrocell nodes
            __global double* dtwn){       // time derivative
  
  int ipg=get_global_id(0);
  int npg=(param[1]+1)*(param[2]+1)*(param[3]+1) *
         (param[4])*(param[5])*(param[6]);

  double dtau[3][3],x,y,z,wpg;
  //ref_pg_vol(param+1,ipg,xpgref,&wpg,NULL);
  int ix = ipg % (param[1] + 1);
  ipg/=(param[1] + 1);

  int iy = ipg % (param[2] + 1);
  ipg/=(param[2] + 1);

  int iz = ipg % (param[3] + 1);
  ipg/=(param[3] + 1);

  int ncx= ipg % param[4];
  double hx=1/(double) param[4];
  ipg/=param[4];

  int ncy= ipg % param[5];
  double hy=1/(double) param[5];
  ipg/=param[5];

  int ncz= ipg;
  double hz=1/(double) param[6];

  int offset[3];

  offset[0]=gauss_lob_offset[param[1]]+ix;
  offset[1]=gauss_lob_offset[param[2]]+iy;
  offset[2]=gauss_lob_offset[param[3]]+iz;

  x=hx*(ncx+gauss_lob_point[offset[0]]);
  y=hy*(ncy+gauss_lob_point[offset[1]]);
  z=hz*(ncz+gauss_lob_point[offset[2]]);


  wpg=hx*hy*hz*gauss_lob_weight[offset[0]]*
    gauss_lob_weight[offset[1]]*
    gauss_lob_weight[offset[2]];


  // end of ref_pg_vol
  //////////////////////////////////////////////

  //Ref2Phy(physnode, // phys. nodes
  //        xpgref,  // xref
  //        NULL,-1, // dpsiref,ifa
  //        NULL,dtau,  // xphy,dtau
  //        codtau,NULL,NULL); // codtau,dpsi,vnds

   get_dtau(x,y,z,physnode,dtau);

   //codtau[0][0] = dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
   //codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
   //codtau[0][2] = dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
   //codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
   //codtau[1][1] = dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
   //codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
   //codtau[2][0] = dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
   //codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
   //codtau[2][2] = dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
   //double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
   // dtau[0][2]*codtau[0][2];
   
  // end of Ref2Phy
  //////////////////////////////////////////////////////////

  double det=dtau[0][0]*dtau[1][1]*dtau[2][2]-dtau[0][0]*dtau[1][2]*dtau[2][1]-dtau[1][0]*dtau[0][1]*dtau[2][2]+
    dtau[1][0]*dtau[0][2]*dtau[2][1]+dtau[2][0]*dtau[0][1]*dtau[1][2]-dtau[2][0]*dtau[0][2]*dtau[1][1];

  for(int iv=0;iv<param[0];iv++){
    // varindex
    int imem=iv + param[0] * ( get_global_id(0) + npg * *ie);
    // end of varindex
    /////////////////////////////////////
    //printf("imem=%d dtw=%f\n",imem,dtwn[imem]);
    //printf("det=%f wpg=%f imem=%d h=%f %f %f\n",det,wpg,imem,hx,hy,hz);

    dtwn[imem]/=(wpg*det);
  }

}


void get_dtau(double x,double y,double z,__constant double* p,double dtau[][3]){

 // gradient of the shape functions and value (4th component)
  // of the shape functions
  /* double gradphi[20][3]; */
  /* //double x,y,z; */
  /* // this fills the values of gradphi */
  /* gradphi[0][0] = (-1 + z) * (-1 + y) * (2 * y + 2 * z + 4 * x - 3); */
  /* gradphi[0][1] = (-1 + z) * (-1 + x) * (2 * x + 2 * z - 3 + 4 * y); */
  /* gradphi[0][2] = (-1 + y) * (-1 + x) * (2 * x + 2 * y - 3 + 4 * z); */
  /* gradphi[1][0] = (-1 + z) * (-1 + y) * (-2 * y - 2 * z - 1 + 4 * x); */
  /* gradphi[1][1] = x * (-1 + z) * (-2 * z + 1 - 4 * y + 2 * x); */
  /* gradphi[1][2] = x * (-1 + y) * (-2 * y + 1 - 4 * z + 2 * x); */
  /* gradphi[2][0] = -y * (-1 + z) * (2 * y - 2 * z - 3 + 4 * x); */
  /* gradphi[2][1] = -x * (-1 + z) * (-2 * z - 3 + 4 * y + 2 * x); */
  /* gradphi[2][2] = -x * y * (2 * y - 4 * z - 1 + 2 * x); */
  /* gradphi[3][0] = -y * (-1 + z) * (-2 * y + 2 * z - 1 + 4 * x); */
  /* gradphi[3][1] = -(-1 + z) * (-1 + x) * (2 * x + 2 * z + 1 - 4 * y); */
  /* gradphi[3][2] = -y * (-1 + x) * (2 * x - 1 + 4 * z - 2 * y); */
  /* gradphi[4][0] = -z * (-1 + y) * (2 * y - 2 * z - 1 + 4 * x); */
  /* gradphi[4][1] = -z * (-1 + x) * (2 * x - 2 * z - 1 + 4 * y); */
  /* gradphi[4][2] = -(-1 + y) * (-1 + x) * (2 * x - 4 * z + 2 * y + 1); */
  /* gradphi[5][0] = -z * (-1 + y) * (-2 * y + 2 * z + 4 * x - 3); */
  /* gradphi[5][1] = -x * z * (2 * z - 4 * y - 1 + 2 * x); */
  /* gradphi[5][2] = -x * (-1 + y) * (-2 * y - 3 + 4 * z + 2 * x); */
  /* gradphi[6][0] = y * z * (2 * y + 2 * z + 4 * x - 5); */
  /* gradphi[6][1] = x * z * (2 * z + 4 * y - 5 + 2 * x); */
  /* gradphi[6][2] = x * y * (2 * y + 4 * z - 5 + 2 * x); */
  /* gradphi[7][0] = y * z * (-2 * y - 2 * z + 1 + 4 * x); */
  /* gradphi[7][1] = z * (-1 + x) * (2 * x - 2 * z + 3 - 4 * y); */
  /* gradphi[7][2] = y * (-1 + x) * (2 * x - 2 * y + 3 - 4 * z); */
  /* gradphi[8][0] = -4 * (-1 + z) * (-1 + y) * (2 * x - 1); */
  /* gradphi[8][1] = -4 * x * (-1 + z) * (-1 + x); */
  /* gradphi[8][2] = -4 * x * (-1 + y) * (-1 + x); */
  /* gradphi[9][0] = -4 * y * (-1 + z) * (-1 + y); */
  /* gradphi[9][1] = -4 * (-1 + z) * (2 * y - 1) * (-1 + x); */
  /* gradphi[9][2] = -4 * y * (-1 + y) * (-1 + x); */
  /* gradphi[10][0] = -4 * z * (-1 + z) * (-1 + y); */
  /* gradphi[10][1] = -4 * z * (-1 + z) * (-1 + x); */
  /* gradphi[10][2] = -4 * (2 * z - 1) * (-1 + y) * (-1 + x); */
  /* gradphi[11][0] = 4 * y * (-1 + z) * (-1 + y); */
  /* gradphi[11][1] = 4 * x * (-1 + z) * (2 * y - 1); */
  /* gradphi[11][2] = 4 * x * y * (-1 + y); */
  /* gradphi[12][0] = 4 * z * (-1 + z) * (-1 + y); */
  /* gradphi[12][1] = 4 * x * z * (-1 + z); */
  /* gradphi[12][2] = 4 * x * (2 * z - 1) * (-1 + y); */
  /* gradphi[13][0] = 4 * y * (-1 + z) * (2 * x - 1); */
  /* gradphi[13][1] = 4 * x * (-1 + z) * (-1 + x); */
  /* gradphi[13][2] = 4 * x * y * (-1 + x); */
  /* gradphi[14][0] = -4 * y * z * (-1 + z); */
  /* gradphi[14][1] = -4 * x * z * (-1 + z); */
  /* gradphi[14][2] = -4 * x * y * (2 * z - 1); */
  /* gradphi[15][0] = 4 * y * z * (-1 + z); */
  /* gradphi[15][1] = 4 * z * (-1 + z) * (-1 + x); */
  /* gradphi[15][2] = 4 * y * (2 * z - 1) * (-1 + x); */
  /* gradphi[16][0] = 4 * z * (-1 + y) * (2 * x - 1); */
  /* gradphi[16][1] = 4 * x * z * (-1 + x); */
  /* gradphi[16][2] = 4 * x * (-1 + y) * (-1 + x); */
  /* gradphi[17][0] = 4 * y * z * (-1 + y); */
  /* gradphi[17][1] = 4 * z * (2 * y - 1) * (-1 + x); */
  /* gradphi[17][2] = 4 * y * (-1 + y) * (-1 + x); */
  /* gradphi[18][0] = -4 * y * z * (-1 + y); */
  /* gradphi[18][1] = -4 * x * z * (2 * y - 1); */
  /* gradphi[18][2] = -4 * x * y * (-1 + y); */
  /* gradphi[19][0] = -4 * y * z * (2 * x - 1); */
  /* gradphi[19][1] = -4 * x * z * (-1 + x); */
  /* gradphi[19][2] = -4 * x * y * (-1 + x); */
  /* for(int ii=0;ii<3;ii++){ */
  /*   for(int jj=0;jj<3;jj++){ */
  /*     dtau[ii][jj]=0; */
  /*   } */
  /*   for(int i=0;i<20;i++){ */
  /*     //printf("xyzphy=%f %f %f \n",physnode[3*i+0],physnode[3*i+1],physnode[3*i+2]); */
  /*     for(int jj=0;jj<3;jj++){ */
  /*       dtau[ii][jj]+=physnode[3*i+ii]*gradphi[i][jj];; */
  /*     } */
  /*   } */
  /* } */

dtau[0][0]=2*(-1+z)*(-1+y)*(-1+x)*p[0]+2*x*(-1+z)*(-1+y)*p[3]-2*x*y*(-1+z)*p[6]-2*y*(-1+z)*(-1+x)*p[9]-2*z*(-1+y)*(-1+x)*p[12]-2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]+2*y*z*(-1+x)*p[21]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[0]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[3]-y*(-1+z)*(2*y-2*z-3+2*x)*p[6]-y*(-1+z)*(2*x+2*z+1-2*y)*p[9]-z*(-1+y)*(2*x+2*y-2*z+1)*p[12]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[15]+y*z*(2*y+2*z-5+2*x)*p[18]+y*z*(2*x-2*z+3-2*y)*p[21]-4*(-1+z)*(-1+y)*(-1+x)*p[24]-4*x*(-1+z)*(-1+y)*p[24]-4*y*(-1+z)*(-1+y)*p[27]-4*z*(-1+z)*(-1+y)*p[30]+4*y*(-1+z)*(-1+y)*p[33]+4*z*(-1+z)*(-1+y)*p[36]+4*y*(-1+z)*(-1+x)*p[39]+4*x*y*(-1+z)*p[39]-4*y*z*(-1+z)*p[42]+4*y*z*(-1+z)*p[45]+4*z*(-1+y)*(-1+x)*p[48]+4*x*z*(-1+y)*p[48]+4*y*z*(-1+y)*p[51]-4*y*z*(-1+y)*p[54]-4*y*z*(-1+x)*p[57]-4*x*y*z*p[57];

dtau[0][1]=2*(-1+z)*(-1+y)*(-1+x)*p[0]-2*x*(-1+z)*(-1+y)*p[3]-2*x*y*(-1+z)*p[6]+2*y*(-1+z)*(-1+x)*p[9]-2*z*(-1+y)*(-1+x)*p[12]+2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]-2*y*z*(-1+x)*p[21]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[0]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[3]-x*(-1+z)*(2*y-2*z-3+2*x)*p[6]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[9]-z*(-1+x)*(2*x+2*y-2*z+1)*p[12]-x*z*(-2*y+2*z-3+2*x)*p[15]+x*z*(2*y+2*z-5+2*x)*p[18]+z*(-1+x)*(2*x-2*z+3-2*y)*p[21]-4*x*(-1+z)*(-1+x)*p[24]-4*(-1+z)*(-1+y)*(-1+x)*p[27]-4*y*(-1+z)*(-1+x)*p[27]-4*z*(-1+z)*(-1+x)*p[30]+4*x*(-1+z)*(-1+y)*p[33]+4*x*y*(-1+z)*p[33]+4*x*z*(-1+z)*p[36]+4*x*(-1+z)*(-1+x)*p[39]-4*x*z*(-1+z)*p[42]+4*z*(-1+z)*(-1+x)*p[45]+4*x*z*(-1+x)*p[48]+4*z*(-1+y)*(-1+x)*p[51]+4*y*z*(-1+x)*p[51]-4*x*z*(-1+y)*p[54]-4*x*y*z*p[54]-4*x*z*(-1+x)*p[57];

dtau[0][2]=2*(-1+z)*(-1+y)*(-1+x)*p[0]-2*x*(-1+z)*(-1+y)*p[3]+2*x*y*(-1+z)*p[6]-2*y*(-1+z)*(-1+x)*p[9]+2*z*(-1+y)*(-1+x)*p[12]-2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]-2*y*z*(-1+x)*p[21]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[0]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[3]-x*y*(2*y-2*z-3+2*x)*p[6]-y*(-1+x)*(2*x+2*z+1-2*y)*p[9]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[12]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[15]+x*y*(2*y+2*z-5+2*x)*p[18]+y*(-1+x)*(2*x-2*z+3-2*y)*p[21]-4*x*(-1+y)*(-1+x)*p[24]-4*y*(-1+y)*(-1+x)*p[27]-4*(-1+z)*(-1+y)*(-1+x)*p[30]-4*z*(-1+y)*(-1+x)*p[30]+4*x*y*(-1+y)*p[33]+4*x*(-1+z)*(-1+y)*p[36]+4*x*z*(-1+y)*p[36]+4*x*y*(-1+x)*p[39]-4*x*y*(-1+z)*p[42]-4*x*y*z*p[42]+4*y*(-1+z)*(-1+x)*p[45]+4*y*z*(-1+x)*p[45]+4*x*(-1+y)*(-1+x)*p[48]+4*y*(-1+y)*(-1+x)*p[51]-4*x*y*(-1+y)*p[54]-4*x*y*(-1+x)*p[57];

dtau[1][0]=2*(-1+z)*(-1+y)*(-1+x)*p[1]+2*x*(-1+z)*(-1+y)*p[4]-2*x*y*(-1+z)*p[7]-2*y*(-1+z)*(-1+x)*p[10]-2*z*(-1+y)*(-1+x)*p[13]-2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]+2*y*z*(-1+x)*p[22]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[1]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[4]-y*(-1+z)*(2*y-2*z-3+2*x)*p[7]-y*(-1+z)*(2*x+2*z+1-2*y)*p[10]-z*(-1+y)*(2*x+2*y-2*z+1)*p[13]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[16]+y*z*(2*y+2*z-5+2*x)*p[19]+y*z*(2*x-2*z+3-2*y)*p[22]-4*(-1+z)*(-1+y)*(-1+x)*p[25]-4*x*(-1+z)*(-1+y)*p[25]-4*y*(-1+z)*(-1+y)*p[28]-4*z*(-1+z)*(-1+y)*p[31]+4*y*(-1+z)*(-1+y)*p[34]+4*z*(-1+z)*(-1+y)*p[37]+4*y*(-1+z)*(-1+x)*p[40]+4*x*y*(-1+z)*p[40]-4*y*z*(-1+z)*p[43]+4*y*z*(-1+z)*p[46]+4*z*(-1+y)*(-1+x)*p[49]+4*x*z*(-1+y)*p[49]+4*y*z*(-1+y)*p[52]-4*y*z*(-1+y)*p[55]-4*y*z*(-1+x)*p[58]-4*x*y*z*p[58];

dtau[1][1]=2*(-1+z)*(-1+y)*(-1+x)*p[1]-2*x*(-1+z)*(-1+y)*p[4]-2*x*y*(-1+z)*p[7]+2*y*(-1+z)*(-1+x)*p[10]-2*z*(-1+y)*(-1+x)*p[13]+2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]-2*y*z*(-1+x)*p[22]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[1]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[4]-x*(-1+z)*(2*y-2*z-3+2*x)*p[7]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[10]-z*(-1+x)*(2*x+2*y-2*z+1)*p[13]-x*z*(-2*y+2*z-3+2*x)*p[16]+x*z*(2*y+2*z-5+2*x)*p[19]+z*(-1+x)*(2*x-2*z+3-2*y)*p[22]-4*x*(-1+z)*(-1+x)*p[25]-4*(-1+z)*(-1+y)*(-1+x)*p[28]-4*y*(-1+z)*(-1+x)*p[28]-4*z*(-1+z)*(-1+x)*p[31]+4*x*(-1+z)*(-1+y)*p[34]+4*x*y*(-1+z)*p[34]+4*x*z*(-1+z)*p[37]+4*x*(-1+z)*(-1+x)*p[40]-4*x*z*(-1+z)*p[43]+4*z*(-1+z)*(-1+x)*p[46]+4*x*z*(-1+x)*p[49]+4*z*(-1+y)*(-1+x)*p[52]+4*y*z*(-1+x)*p[52]-4*x*z*(-1+y)*p[55]-4*x*y*z*p[55]-4*x*z*(-1+x)*p[58];

dtau[1][2]=2*(-1+z)*(-1+y)*(-1+x)*p[1]-2*x*(-1+z)*(-1+y)*p[4]+2*x*y*(-1+z)*p[7]-2*y*(-1+z)*(-1+x)*p[10]+2*z*(-1+y)*(-1+x)*p[13]-2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]-2*y*z*(-1+x)*p[22]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[1]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[4]-x*y*(2*y-2*z-3+2*x)*p[7]-y*(-1+x)*(2*x+2*z+1-2*y)*p[10]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[13]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[16]+x*y*(2*y+2*z-5+2*x)*p[19]+y*(-1+x)*(2*x-2*z+3-2*y)*p[22]-4*x*(-1+y)*(-1+x)*p[25]-4*y*(-1+y)*(-1+x)*p[28]-4*(-1+z)*(-1+y)*(-1+x)*p[31]-4*z*(-1+y)*(-1+x)*p[31]+4*x*y*(-1+y)*p[34]+4*x*(-1+z)*(-1+y)*p[37]+4*x*z*(-1+y)*p[37]+4*x*y*(-1+x)*p[40]-4*x*y*(-1+z)*p[43]-4*x*y*z*p[43]+4*y*(-1+z)*(-1+x)*p[46]+4*y*z*(-1+x)*p[46]+4*x*(-1+y)*(-1+x)*p[49]+4*y*(-1+y)*(-1+x)*p[52]-4*x*y*(-1+y)*p[55]-4*x*y*(-1+x)*p[58];

dtau[2][0]=2*(-1+z)*(-1+y)*(-1+x)*p[2]+2*x*(-1+z)*(-1+y)*p[5]-2*x*y*(-1+z)*p[8]-2*y*(-1+z)*(-1+x)*p[11]-4*y*(-1+z)*(-1+y)*p[29]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[2]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[5]-y*(-1+z)*(2*y-2*z-3+2*x)*p[8]-y*(-1+z)*(2*x+2*z+1-2*y)*p[11]-z*(-1+y)*(2*x+2*y-2*z+1)*p[14]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[17]+y*z*(2*y+2*z-5+2*x)*p[20]+y*z*(2*x-2*z+3-2*y)*p[23]-4*(-1+z)*(-1+y)*(-1+x)*p[26]-4*x*(-1+z)*(-1+y)*p[26]-2*z*(-1+y)*(-1+x)*p[14]-2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]+2*y*z*(-1+x)*p[23]-4*z*(-1+z)*(-1+y)*p[32]+4*y*(-1+z)*(-1+y)*p[35]+4*z*(-1+z)*(-1+y)*p[38]+4*y*(-1+z)*(-1+x)*p[41]+4*x*y*(-1+z)*p[41]-4*y*z*(-1+z)*p[44]+4*y*z*(-1+z)*p[47]+4*z*(-1+y)*(-1+x)*p[50]+4*x*z*(-1+y)*p[50]+4*y*z*(-1+y)*p[53]-4*y*z*(-1+y)*p[56]-4*y*z*(-1+x)*p[59]-4*x*y*z*p[59];

dtau[2][1]=2*(-1+z)*(-1+y)*(-1+x)*p[2]-2*x*(-1+z)*(-1+y)*p[5]-2*x*y*(-1+z)*p[8]+2*y*(-1+z)*(-1+x)*p[11]-2*z*(-1+y)*(-1+x)*p[14]+2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]-2*y*z*(-1+x)*p[23]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[2]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[5]-x*(-1+z)*(2*y-2*z-3+2*x)*p[8]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[11]-z*(-1+x)*(2*x+2*y-2*z+1)*p[14]-x*z*(-2*y+2*z-3+2*x)*p[17]+x*z*(2*y+2*z-5+2*x)*p[20]+z*(-1+x)*(2*x-2*z+3-2*y)*p[23]-4*x*(-1+z)*(-1+x)*p[26]-4*(-1+z)*(-1+y)*(-1+x)*p[29]-4*y*(-1+z)*(-1+x)*p[29]-4*z*(-1+z)*(-1+x)*p[32]+4*x*(-1+z)*(-1+y)*p[35]+4*x*y*(-1+z)*p[35]+4*x*z*(-1+z)*p[38]+4*x*(-1+z)*(-1+x)*p[41]-4*x*z*(-1+z)*p[44]+4*z*(-1+z)*(-1+x)*p[47]+4*x*z*(-1+x)*p[50]+4*z*(-1+y)*(-1+x)*p[53]+4*y*z*(-1+x)*p[53]-4*x*z*(-1+y)*p[56]-4*x*y*z*p[56]-4*x*z*(-1+x)*p[59];

dtau[2][2]=2*(-1+z)*(-1+y)*(-1+x)*p[2]-2*x*(-1+z)*(-1+y)*p[5]+2*x*y*(-1+z)*p[8]-2*y*(-1+z)*(-1+x)*p[11]+2*z*(-1+y)*(-1+x)*p[14]-2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]-2*y*z*(-1+x)*p[23]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[2]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[5]-x*y*(2*y-2*z-3+2*x)*p[8]-y*(-1+x)*(2*x+2*z+1-2*y)*p[11]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[14]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[17]+x*y*(2*y+2*z-5+2*x)*p[20]+y*(-1+x)*(2*x-2*z+3-2*y)*p[23]-4*x*(-1+y)*(-1+x)*p[26]-4*y*(-1+y)*(-1+x)*p[29]-4*(-1+z)*(-1+y)*(-1+x)*p[32]-4*z*(-1+y)*(-1+x)*p[32]+4*x*y*(-1+y)*p[35]+4*x*(-1+z)*(-1+y)*p[38]+4*x*z*(-1+y)*p[38]+4*x*y*(-1+x)*p[41]-4*x*y*(-1+z)*p[44]-4*x*y*z*p[44]+4*y*(-1+z)*(-1+x)*p[47]+4*y*z*(-1+x)*p[47]+4*x*(-1+y)*(-1+x)*p[50]+4*y*(-1+y)*(-1+x)*p[53]-4*x*y*(-1+y)*p[56]-4*x*y*(-1+x)*p[59];

}

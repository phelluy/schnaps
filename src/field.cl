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


//!  \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
double wglop(int deg,int i){
  return gauss_lob_weight[gauss_lob_offset[deg]+i];
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

   printf("param=%d %d %d %d %d %d %d ie=%d dtw=%f\n",param[0],
	  param[1],param[2],param[3],param[4],param[5],param[6],
	  *ie,dtwn[0]);

  for(int i=0;i<20;i++){
    //printf("ie=%d physnode[%d]=%f %f %f\n",*ie,i,physnode[3*i+0],physnode[3*i+1],physnode[3*i+2]);
  }

  double dtau[3][3],codtau[3][3],x,y,z,wpg;
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

  offset[0]=gauss_lob_offset[param[0]]+ix;
  offset[1]=gauss_lob_offset[param[1]]+iy;
  offset[2]=gauss_lob_offset[param[2]]+iz;

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

  // gradient of the shape functions and value (4th component)
  // of the shape functions
  double gradphi[20][3];
  // this fills the values of gradphi
  gradphi[0][0] = (-1 + z) * (-1 + y) * (2 * y + 2 * z + 4 * x - 3);
  gradphi[0][1] = (-1 + z) * (-1 + x) * (2 * x + 2 * z - 3 + 4 * y);
  gradphi[0][2] = (-1 + y) * (-1 + x) * (2 * x + 2 * y - 3 + 4 * z);
  gradphi[1][0] = (-1 + z) * (-1 + y) * (-2 * y - 2 * z - 1 + 4 * x);
  gradphi[1][1] = x * (-1 + z) * (-2 * z + 1 - 4 * y + 2 * x);
  gradphi[1][2] = x * (-1 + y) * (-2 * y + 1 - 4 * z + 2 * x);
  gradphi[2][0] = -y * (-1 + z) * (2 * y - 2 * z - 3 + 4 * x);
  gradphi[2][1] = -x * (-1 + z) * (-2 * z - 3 + 4 * y + 2 * x);
  gradphi[2][2] = -x * y * (2 * y - 4 * z - 1 + 2 * x);
  gradphi[3][0] = -y * (-1 + z) * (-2 * y + 2 * z - 1 + 4 * x);
  gradphi[3][1] = -(-1 + z) * (-1 + x) * (2 * x + 2 * z + 1 - 4 * y);
  gradphi[3][2] = -y * (-1 + x) * (2 * x - 1 + 4 * z - 2 * y);
  gradphi[4][0] = -z * (-1 + y) * (2 * y - 2 * z - 1 + 4 * x);
  gradphi[4][1] = -z * (-1 + x) * (2 * x - 2 * z - 1 + 4 * y);
  gradphi[4][2] = -(-1 + y) * (-1 + x) * (2 * x - 4 * z + 2 * y + 1);
  gradphi[5][0] = -z * (-1 + y) * (-2 * y + 2 * z + 4 * x - 3);
  gradphi[5][1] = -x * z * (2 * z - 4 * y - 1 + 2 * x);
  gradphi[5][2] = -x * (-1 + y) * (-2 * y - 3 + 4 * z + 2 * x);
  gradphi[6][0] = y * z * (2 * y + 2 * z + 4 * x - 5);
  gradphi[6][1] = x * z * (2 * z + 4 * y - 5 + 2 * x);
  gradphi[6][2] = x * y * (2 * y + 4 * z - 5 + 2 * x);
  gradphi[7][0] = y * z * (-2 * y - 2 * z + 1 + 4 * x);
  gradphi[7][1] = z * (-1 + x) * (2 * x - 2 * z + 3 - 4 * y);
  gradphi[7][2] = y * (-1 + x) * (2 * x - 2 * y + 3 - 4 * z);
  gradphi[8][0] = -4 * (-1 + z) * (-1 + y) * (2 * x - 1);
  gradphi[8][1] = -4 * x * (-1 + z) * (-1 + x);
  gradphi[8][2] = -4 * x * (-1 + y) * (-1 + x);
  gradphi[9][0] = -4 * y * (-1 + z) * (-1 + y);
  gradphi[9][1] = -4 * (-1 + z) * (2 * y - 1) * (-1 + x);
  gradphi[9][2] = -4 * y * (-1 + y) * (-1 + x);
  gradphi[10][0] = -4 * z * (-1 + z) * (-1 + y);
  gradphi[10][1] = -4 * z * (-1 + z) * (-1 + x);
  gradphi[10][2] = -4 * (2 * z - 1) * (-1 + y) * (-1 + x);
  gradphi[11][0] = 4 * y * (-1 + z) * (-1 + y);
  gradphi[11][1] = 4 * x * (-1 + z) * (2 * y - 1);
  gradphi[11][2] = 4 * x * y * (-1 + y);
  gradphi[12][0] = 4 * z * (-1 + z) * (-1 + y);
  gradphi[12][1] = 4 * x * z * (-1 + z);
  gradphi[12][2] = 4 * x * (2 * z - 1) * (-1 + y);
  gradphi[13][0] = 4 * y * (-1 + z) * (2 * x - 1);
  gradphi[13][1] = 4 * x * (-1 + z) * (-1 + x);
  gradphi[13][2] = 4 * x * y * (-1 + x);
  gradphi[14][0] = -4 * y * z * (-1 + z);
  gradphi[14][1] = -4 * x * z * (-1 + z);
  gradphi[14][2] = -4 * x * y * (2 * z - 1);
  gradphi[15][0] = 4 * y * z * (-1 + z);
  gradphi[15][1] = 4 * z * (-1 + z) * (-1 + x);
  gradphi[15][2] = 4 * y * (2 * z - 1) * (-1 + x);
  gradphi[16][0] = 4 * z * (-1 + y) * (2 * x - 1);
  gradphi[16][1] = 4 * x * z * (-1 + x);
  gradphi[16][2] = 4 * x * (-1 + y) * (-1 + x);
  gradphi[17][0] = 4 * y * z * (-1 + y);
  gradphi[17][1] = 4 * z * (2 * y - 1) * (-1 + x);
  gradphi[17][2] = 4 * y * (-1 + y) * (-1 + x);
  gradphi[18][0] = -4 * y * z * (-1 + y);
  gradphi[18][1] = -4 * x * z * (2 * y - 1);
  gradphi[18][2] = -4 * x * y * (-1 + y);
  gradphi[19][0] = -4 * y * z * (2 * x - 1);
  gradphi[19][1] = -4 * x * z * (-1 + x);
  gradphi[19][2] = -4 * x * y * (-1 + x);
  for(int ii=0;ii<3;ii++){
    for(int jj=0;jj<3;jj++){
      dtau[ii][jj]=0;
    }
    for(int i=0;i<20;i++){
      //printf("xyzphy=%f %f %f \n",physnode[3*i+0],physnode[3*i+1],physnode[3*i+2]);
      for(int jj=0;jj<3;jj++){
        dtau[ii][jj]+=physnode[3*i+ii]*gradphi[i][jj];;
      }
    }
  }
  codtau[0][0] = dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
  codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
  codtau[0][2] = dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
  codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
  codtau[1][1] = dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
  codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
  codtau[2][0] = dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
  codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
  codtau[2][2] = dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];

  // end of Ref2Phy
  //////////////////////////////////////////////////////////

  double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
    dtau[0][2]*codtau[0][2];

  //printf("det=%f\n",det);
  for(int iv=0;iv<param[0];iv++){
    // varindex
    int imem=iv + param[0] * ( get_global_id(0) + npg * *ie);
    // end of varindex
    /////////////////////////////////////
    printf("imem=%d dtw=%f\n",imem,dtwn[imem]);

    dtwn[imem]/=(wpg*det);
  }

}

#include "macromesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "geometry.h"
#include "interpolation.h"
#include <math.h>


void ReadMacroMesh(MacroMesh* m,char* filename){

  FILE* f=NULL;

  char* line=NULL;
  size_t linesize=0;
  //size_t ret;

  printf("Read mesh file %s\n",filename);

  f=fopen(filename,"r");
  assert(f != NULL);

  do {
    getline(&line,&linesize,f);
  }
  while(strcmp(line,"$Nodes\n") != 0);   

  // read the nodes data
  getline(&line,&linesize,f);
  m->nbnodes=atoi(line);
  printf("nbnodes=%d\n",m->nbnodes);

  m->node=malloc(3 * m->nbnodes * sizeof(double));
  assert(m->node);

  for(int i=0;i<m->nbnodes;i++){    
    getdelim(&line,&linesize,(int) ' ',f); // node number
    getdelim(&line,&linesize,(int) ' ',f); // x
    m->node[3*i+0]=atof(line);
    getdelim(&line,&linesize,(int) ' ',f); // y
    m->node[3*i+1]=atof(line);
    getline(&line,&linesize,f); // z (end of the line)
    m->node[3*i+2]=atof(line);
    /* printf("Node %d x=%f y=%f z=%f\n",i, */
    /* 	   m->node[3*i+0], */
    /* 	   m->node[3*i+1], */
    /* 	   m->node[3*i+2]); */
  }    
  getline(&line,&linesize,f);
  //printf("%s",line);
  // check that we have reached the end of nodes
  assert(strcmp(line,"$EndNodes\n") == 0);
  
  // Now read all the elements of the mesh
  do {
    getline(&line,&linesize,f);
  }
  while(strcmp(line,"$Elements\n") != 0);   

  // size of the gmsh elems list 
  getline(&line,&linesize,f);
  int nball=atoi(line);
  // allocate to a too big size
  m->elem2node=malloc(20 * sizeof(int) * nball);
  assert(m->elem2node);

  // now count only the H20 elems (code=17)
  m->nbelems=0;
  int countnode=0;
  for(int i=0;i<nball;i++){
    getdelim(&line,&linesize,(int) ' ',f); // elem number
    getdelim(&line,&linesize,(int) ' ',f); // elem type
    int elemtype=atoi(line);
    if (elemtype != 17) {
      getline(&line,&linesize,f);
    }
    else {
      m->nbelems++;
      getdelim(&line,&linesize,(int) ' ',f); //useless code
      getdelim(&line,&linesize,(int) ' ',f); //useless code
      getdelim(&line,&linesize,(int) ' ',f); //useless code
      for(int j=0;j<19;j++){
	getdelim(&line,&linesize,(int) ' ',f);
	//printf("%d ",atoi(line));
	m->elem2node[countnode]=atoi(line)-1;
	countnode++;
      }
      getline(&line,&linesize,f);
      //printf("%d\n",atoi(line));
      m->elem2node[countnode]=atoi(line)-1;
      countnode++;
    }
  }	
  getline(&line,&linesize,f);
  //printf("%s",line);
  // check that we have reached the end of nodes
  assert(strcmp(line,"$EndElements\n") == 0);
  printf("nbelems=%d\n",m->nbelems);
  //m->elem2node=realloc(m->elem2node,20 * sizeof(int) * m->nbelems);
  //assert(m->elem2node);

  m->elem2elem=NULL;

}


void AffineMap(double* x){

  //double A[3][3]={{1,2,1},{0,-1,4},{7,8,-5}};
  double A[3][3]={{0,-1,0},{-1,0,0},{0,0,-1}};
  //double A[3][3]={1,0,0,0,2,0,0,0,1};
  double x0[3]={0,0,0};
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


void AffineMapMacroMesh(MacroMesh* m){

    for(int ino=0;ino<m->nbnodes;ino++){
      AffineMap(&(m->node[ino*3]));
    }
}



// display macromesh data on standard output
void PrintMacroMesh(MacroMesh* m){
  printf("Print macromesh...\n");
  int start=1;
  printf("nbnodes=%d\n",m->nbnodes);
  for(int i=0;i<m->nbnodes;i++){
    printf("node %d x=%f y=%f z=%f\n",i+start,
	   m->node[3*i+0],
	   m->node[3*i+1],
	   m->node[3*i+2]);
  }
  printf("nbelems=%d\n",m->nbelems);
  for(int i=0;i<m->nbelems;i++){
    printf("elem %d -> ",i+start);
    for(int j=0;j<20;j++){
      printf("%d ",m->elem2node[20*i+j]+start);
    }
    printf("\n");
  }

  for(int i=0;i<m->nbelems;i++){
    printf("elem %d voisins: ",i+start);
    for(int j=0;j<6;j++){
      printf("%d ",m->elem2elem[6*i+j]+start);
    }
    printf("\n");
  }
}

// build others connectivity arrays
void BuildConnectivity(MacroMesh* m){

  printf("Build connectivity...\n");

  assert(m->elem2elem==NULL);
  
  // build a list of faces
  // each face is made of four corners of 
  // the hexaedron mesh
  Face4Sort* face;
  Face4Sort* f;

  face=malloc(6*sizeof(Face4Sort)*m->nbelems);

  assert(face);

  int face2locnode[6][4]={
    {0,1,5,4},
    {1,2,6,5},
    {2,3,7,6},
    {0,4,7,3},
    {5,6,7,4},
    {0,3,2,1},
};


  for(int ie=0;ie<m->nbelems;ie++){
    for(int ifa=0;ifa<6;ifa++){
      f=face+ifa+6*ie;
      for(int ino=0;ino<4;ino++){
	f->node[ino]=m->elem2node[face2locnode[ifa][ino]+20*ie];
      }
      f->left=ie;
      f->locfaceleft=ifa;
      f->right=-1;
      f->locfaceright=-1;
      OrderFace4Sort(f);
      /* printf("elem=%d ifa=%d left=%d nodes %d %d %d %d\n",ie,ifa, */
      /* 	     f->left,f->node[0], */
      /* 	     f->node[1],f->node[2],f->node[3]); */
    }
  }

  // now sort the list of faces
  qsort(face,6*m->nbelems,sizeof(Face4Sort),CompareFace4Sort);
  // check
  /* for(int ie=0;ie<m->nbelems;ie++){ */
  /*   for(int ifa=0;ifa<6;ifa++){ */
  /*     f=face+ifa+6*ie; */
  /*     printf("left=%d right=%d, nodes %d %d %d %d\n", */
  /* 	     f->left,f->right,f->node[0], */
  /* 	     f->node[1],f->node[2],f->node[3]); */
  /*   } */
  /* } */

  //assert(1==2);

  // allocate element connectivity array
  m->elem2elem=malloc(6 * m->nbelems * sizeof(int));
  for(int i=0;i<6 * m->nbelems;i++){
    m->elem2elem[i]=-1;
  }

  // now, two successive equal faces
  // correspond to two neighbours in the 
  // element list
  for(int ifa=0;ifa<6*m->nbelems-1;ifa++){
    Face4Sort* f1=face+ifa;
    Face4Sort* f2=face+ifa+1;
    if (CompareFace4Sort(f1,f2)==0){
      int ie1=f1->left;
      int if1=f1->locfaceleft;
      int ie2=f2->left;
      int if2=f2->locfaceleft;
	m->elem2elem[if1+6*ie1]=ie2;
	m->elem2elem[if2+6*ie2]=ie1;
    }
  }

  free(face);

  // check
  /* for(int ie=0;ie<m->nbelems;ie++){ */
  /*   for(int ifa=0;ifa<6;ifa++){ */
  /*     printf("elem=%d face=%d, voisin=%d\n", */
  /* 	     ie,ifa,m->elem2elem[ifa+6*ie]); */
  /*   } */
  /* } */

  
}

// compare two integers
int CompareInt(const void* a,const void* b){
  return(*(int*)a - *(int*)b);
}

// sort the nodes list of the face
void OrderFace4Sort(Face4Sort* f){
  qsort ( f->node, 4, sizeof(int), CompareInt);
}

// compare two ordered four-corner faces
// lexicographical order
int CompareFace4Sort(const void* a,const void* b){
  Face4Sort* f1= (Face4Sort*)a; 
  Face4Sort* f2= (Face4Sort*)b;

  int r= f1->node[0]-f2->node[0];
  if (r==0) r=f1->node[1]-f2->node[1];
  if (r==0) r=f1->node[2]-f2->node[2];
  if (r==0) r=f1->node[3]-f2->node[3];
  return r;
  
};

void CheckMacroMesh(MacroMesh* m){

  double dtau[9],codtau[9];
  double physnode[20*3];

  double face_centers[6][3]={
    {0.5,0.0,0.5},
    {1.0,0.5,0.5},
    {0.5,1.0,0.5},
    {0.0,0.5,0.5},
    {0.5,0.5,1.0},
    {0.5,0.5,0.0},
  };

  /* double refnormal[6][3]={{0,-1,0},{1,0,0}, */
  /* 			  {0,1,0},{-1,0,0}, */
  /* 			  {0,0,1},{0,0,-1}}; */


  for(int ie=0;ie<m->nbelems;ie++){
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=m->elem2node[20*ie+inoloc];
      physnode[inoloc*3+0]=m->node[3*ino+0];
      physnode[inoloc*3+1]=m->node[3*ino+1];
      physnode[inoloc*3+2]=m->node[3*ino+2];
    }
    
    // test that the ref_ipg function
    // is compatible with ref_pg_vol
    int param[7]={_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
    for(int ipg=0;ipg<NPG(param);ipg++){
      double xref1[3],xref2[3],xphy[3];
      double wpg;
      ref_pg_vol(param,ipg,xref1,&wpg);
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
      assert(ipg==ref_ipg(param,xref2));
    }

    // middle of the element
    double xrefm[3]={0.5,0.5,0.5},xphym[3];
    Ref2Phy(physnode,
    	    xrefm,
    	    0,-1, // dphiref,ifa
    	    xphym,NULL, // xphy dtau
    	    NULL,NULL,NULL); // codtau,dphi,vnds

    for(int ifa=0;ifa<6;ifa++){
      // middle of the face
      double xphyfa[3],vnds[3];
      Ref2Phy(physnode,
	      face_centers[ifa],
	      0,ifa, // dphiref,ifa
	      xphyfa,dtau,
	      codtau,NULL,vnds); // codtau,dphi,vnds
      double det=codtau[0]*dtau[0]+codtau[1]*dtau[1]+codtau[2]*dtau[2];

      // check volume  orientation
      assert(det >0);
      
      double vec[3]={xphyfa[0]-xphym[0],
		     xphyfa[1]-xphym[1],xphyfa[2]-xphym[2]};
      
      // check face orientation
      assert(vnds[0]*vec[0]+vnds[1]*vec[1]+vnds[2]*vec[2] > 0);

      // check compatibility between face and volume numbering
      //if (m->elem2elem[6*ie+ifa]>=0){
        for(int ipgf=0;ipgf<NPGF(param,ifa);ipgf++){
          double xpgref[3],wpg;
          // get the coordinates of the Gauss point
          ref_pg_face(param,ifa,ipgf,xpgref,&wpg);
          // recover the volume gauss point from
          // the face index
          int ipgv=param[6];
          double xpgref2[3],wpg2;
          ref_pg_vol(param,ipgv,xpgref2,&wpg2);
          assert(Dist(xpgref,xpgref2)<1e-8);
        }
        //}


    }


  }
  



};



  


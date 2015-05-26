#include "pic.h"
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

void InitPIC(PIC* pic,int n){

  pic->nbparts=n;
  pic->xv=malloc(sizeof(real) * 6 * n);
  assert(pic->xv);
  pic->cell_id=malloc(sizeof(int)  * n);
  pic->old_cell_id=malloc(sizeof(int)  * n);
  assert(pic->cell_id);
}

void FreePIC(PIC* pic){

  pic->nbparts=0;
  
  if (pic->xv != NULL) free(pic->xv );
  if (pic->cell_id != NULL) free(pic->cell_id);
  if (pic->old_cell_id != NULL) free(pic->old_cell_id);
}

real corput(int n,int k1,int k2){
  real corput=0;
  real s=1;
  while(n > 0){
    s/=k1;
    corput+=(k2*n%k1)%k1*s;
    n/=k1;
  }
  return corput;
}

void BoxMuller3d(real *xx,int* k1, int* k2)
{
  real x1, x2, rsq;
  
  static int iset=0;
  static int n[4]={0,0,0,0};
  static int n2=0;
  static real gset;
  
  if (iset == 0)
    {
      do {
	x1=2.0*corput(n[0]++,k1[0],k2[0]) -1.0;
	x2=2.0*corput(n[1]++,k1[1],k2[1]) -1.0;
	/* x1=2.0*rand()/RAND_MAX -1.0; */
	/* x2=2.0*rand()/RAND_MAX -1.0; */
	rsq=x1*x1+x2*x2;
      } while (rsq >= 1.0 || rsq == 0);
      rsq = sqrt(-2.0*log(rsq)/rsq);
      xx[2]=x1*rsq;
      iset=1;
      gset = x2*rsq;
    } else
    {
      iset=0;
      xx[2] = gset;
    } 
  do {
    x1 = 2.0 * corput(n[2]++,k1[2],k2[2]) - 1.0;
    x2 = 2.0 * corput(n[3]++,k1[3],k2[3]) - 1.0;
    /* x1 = 2.0 * rand()/RAND_MAX - 1.0; */
    /* x2 = 2.0 * rand()/RAND_MAX - 1.0; */
    rsq = x1*x1 + x2*x2;
  } while ( rsq >= 1.0 || rsq == 0);
  
  rsq = sqrt( -2.0 * log(rsq)/rsq );
  xx[0] = x1 * rsq;
  xx[1] = x2 * rsq;

}


void CreateParticles(PIC* pic,MacroMesh *m){

  int k1[7]={3,5,7,11,13,17,19};
  int k2[7]={2,3,5,7,11,13,17};

  int n=0;
  int np=0;
  
  real vt=1;

  real xp[3];
  real vp[3];

  //srand(time(NULL));
  //int r = rand();


  while(np < pic->nbparts){

    for(int idim=0;idim<3;idim++){
      real r=corput(n,k1[idim],k2[idim]);
      //real r=rand();
      //r/=RAND_MAX;
      xp[idim]=m->xmin[idim]+r*
	(m->xmax[idim]-m->xmin[idim]);
    }

    real xref[3];

    int num_elem=NumElemFromPoint(m,xp,xref);
    //printf("numelem=%d %f %f %f\n",num_elem,xp[0],xp[1],xp[2]);
    
    if (num_elem != -1) {

      /* pic->xv[np*6+0]=xp[0]; */
      /* pic->xv[np*6+1]=xp[1]; */
      /* pic->xv[np*6+2]=xp[2]; */
      pic->xv[np*6+0]=xref[0];
      pic->xv[np*6+1]=xref[1];
      pic->xv[np*6+2]=xref[2];

      real vp[3]={1,0,0};
      BoxMuller3d(vp,k1+3,k2+3);

	/*      for(int idim=0;idim<3;idim++){
	vp[idim]=vt*(corput(n,k1[idim+3],k2[idim+3])-0.5);
	vp[idim]=0;  //!!!!!!!!!!!!!!!!!!
	}*/

      pic->xv[np*6+3]=vp[0];
      pic->xv[np*6+4]=vp[1];
      pic->xv[np*6+5]=vp[2];

      pic->cell_id[np]=num_elem;
      pic->old_cell_id[np]=num_elem;

      np++;
    }

    n++;

    
  }
  

}

void AccumulateParticles(PIC* pic,field *f){
  
  for(int i=0;i<pic->nbparts;i++) {
    
    int ie=pic->old_cell_id[i];
    
    if (ie < 0) goto nexti;
 
    int npg=NPG(f->interp_param + 1);
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    real dtau[3][3], codtau[3][3];
    Ref2Phy(physnode, // phys. nodes
	    pic->xv + 6 * i, // xref
	    NULL, -1, // dpsiref, ifa
	    NULL, dtau, // xphy, dtau
	    codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);

    for(int ib=0;ib < npg;ib++){
      real wpg;
      ref_pg_vol(f->interp_param + 1, ib, NULL, &wpg, NULL);
      printf("det=%f wpg=%f \n", det, wpg);
      wpg *= det;
      real psi;
      psi_ref(f->interp_param+1,ib,pic->xv + 6*i,&psi,NULL);

      int iv=6;  // rho index
      int imem = f->varindex(f->interp_param, ie, ib, iv);
      f->wn[imem]+= psi/wpg;
 
      iv=4;  // j1 index
      imem = f->varindex(f->interp_param, ie, ib, iv);
      f->wn[imem]+= pic->xv[6*i+3]*psi/wpg;

      iv=5;  // j2 index
      imem = f->varindex(f->interp_param, ie, ib, iv);
      f->wn[imem]+= pic->xv[6*i+4]*psi/wpg;
    }
      
    
    
    
  nexti:  
    assert(1==1);
  }
  
}



void PlotParticles(PIC* pic,MacroMesh *m){
  
  FILE * gmshfile;
  FILE * gnufile;
  gmshfile = fopen("partplot.msh", "w" );
  gnufile = fopen("partplot.dat", "w" );

  float x,y,z,vx,vy,vz;
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(real));
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", pic->nbparts);
  /* fic << "$MeshFormat"<<endl; */
  /* fic << "2 0 8" << endl; */
  /* fic << "$EndMeshFormat"<<endl; */
  /* fic << "$Nodes" << endl; */
  
  /* fic << NbPart <<endl;       */
  /* cout << "NbPartFinal " << NbPart << endl; */

  for(int i=0;i<pic->nbparts;i++) {

    int ie=pic->old_cell_id[i];
    
    
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = m->elem2node[20*ie+inoloc];
      physnode[inoloc][0] = m->node[3 * ino + 0];
      physnode[inoloc][1] = m->node[3 * ino + 1];
      physnode[inoloc][2] = m->node[3 * ino + 2];
    }
   
    real xphy[3];
    Ref2Phy(physnode, // phys. nodes
	    &(pic->xv[6*i]), // xref
	    NULL, -1, // dpsiref, ifa
	    xphy, NULL, // xphy, dtau
	    NULL, NULL, NULL); // codtau, dpsi, vnds

    
    /* x=pic->xv[6*i+0]; */
    /* y=pic->xv[6*i+1]; */
    /* z=pic->xv[6*i+2]; */

    x=xphy[0];
    y=xphy[1];
    z=xphy[2];

    vx=pic->xv[6*i+3];
    vy=pic->xv[6*i+4];
    vz=pic->xv[6*i+5];
    fprintf(gmshfile,"%d %f %f %f \n",i+1,x,y,z);
    /* fic << i+1 << " "<<x<<" "<<y<<" "<<0<<endl;	 */
    fprintf(gnufile,"%f %f %f %f %f %f  \n",x,y,z,vx,vy,vz);
  }
  
  fprintf(gmshfile, "$EndNodes\n");
  //fic << "$EndNodes"<<endl;

  fclose(gmshfile);
  fclose(gnufile);

}

void PushParticles(field *f,PIC* pic){

  for(int i=0;i<pic->nbparts;i++) {
    

    // jacobian of tau at the particle
    real physnode[20][3];
    int ie=pic->cell_id[i];
    if (ie >=0) {
      real w[f->model.m];
      InterpField(f,pic->cell_id[i],&(pic->xv[6*i]),w);

      for(int inoloc = 0; inoloc < 20; inoloc++) {
	int ino = f->macromesh.elem2node[20*ie+inoloc];
	physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
	physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
	physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
      }

      real dtau[3][3],codtau[3][3];
      
      real xphy[3];
      Ref2Phy(physnode, // phys. nodes
	      &(pic->xv[6*i]), // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real det = dot_product(dtau[0], codtau[0]);
      
      
      //printf("w=%f %f %f %f\n",w[0],w[1],w[2],w[3]);
      
      // 2D motion
      
      real vref[3];
      pic->xv[6*i+3+0] +=pic->dt * (w[0]+w[2]*pic->xv[6*i+4]);
      pic->xv[6*i+3+1] +=pic->dt * (w[1]-w[2]*pic->xv[6*i+3]);
      pic->xv[6*i+3+2] +=0;


      for(int ii=0;ii<3;ii++){
	vref[ii]=0;
	for(int jj=0;jj<3;jj++){
	  vref[ii] += codtau[jj][ii] * pic->xv[6*i+3+jj] / det ;
	}
      }

      // TO DO: check if the particle is in a new element
      pic->xv[6*i+0]+=pic->dt * vref[0];
      pic->xv[6*i+1]+=pic->dt * vref[1];
      pic->xv[6*i+2]+=pic->dt * vref[2];

      bool is_out = (pic->xv[6*i+0] < 0 || pic->xv[6*i+0] > 1) ||
	(pic->xv[6*i+1] < 0 || pic->xv[6*i+1] > 1)  ||
	(pic->xv[6*i+2] < 0 || pic->xv[6*i+2] > 1);
      
      //is_out = false;

      if (is_out) {
	Ref2Phy(physnode, // phys. nodes
		&(pic->xv[6*i]), // xref
		NULL, -1, // dpsiref, ifa
		xphy, NULL, // xphy, dtau
		NULL, NULL, NULL); // codtau, dpsi, vnds
	int old=pic->cell_id[i];
	printf("oldref=%f %f %f \n",pic->xv[6*i+0],
	       pic->xv[6*i+1],pic->xv[6*i+2]);
	pic->cell_id[i]= NumElemFromPoint(&(f->macromesh),
					  xphy,
					  &(pic->xv[6*i]));
	if (pic->cell_id[i] != -1) pic->old_cell_id[i]=pic->cell_id[i]; 
	printf("newref=%f %f %f \n",pic->xv[6*i+0],
	       pic->xv[6*i+1],pic->xv[6*i+2]);
	printf("change elem: %d -> %d \n",old,pic->cell_id[i]);
      }

    

    } // if ie >= 0
  }
  

}

  


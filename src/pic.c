#include "pic.h"
#include <assert.h>
#include <time.h>
#include <stdlib.h>

void InitPIC(PIC* pic,int n){

  pic->nbparts=n;
  pic->xv=malloc(sizeof(double) * 6 * n);
  assert(pic->xv);
  pic->cell_id=malloc(sizeof(int)  * n);
  assert(pic->cell_id);
}

void FreePIC(PIC* pic){

  pic->nbparts=0;
  
  if (pic->xv != NULL) free(pic->xv );
  if (pic->cell_id != NULL) free(pic->cell_id);
}

double corput(int n,int k1,int k2){
  double corput=0;
  double s=1;
  while(n > 0){
    s/=k1;
    corput+=(k2*n%k1)%k1*s;
    n/=k1;
  }
  return corput;
}


void CreateParticles(PIC* pic,MacroMesh *m){

  const int k1[6]={3,5,7,11,13,17};
  const int k2[6]={2,3,5,7,11,13};

  int n=0;
  int np=0;
  
  double vt=1;

  double xp[3];
  double vp[3];

  //srand(time(NULL));
int r = rand();


  while(np < pic->nbparts){

    for(int idim=0;idim<3;idim++){
      double r=corput(n,k1[idim],k2[idim]);
      //double r=rand();
      //r/=RAND_MAX;
      xp[idim]=m->xmin[idim]+r*
	(m->xmax[idim]-m->xmin[idim]);
    }

    int num_elem=NumElemFromPoint(m,xp,NULL);
    //printf("numelem=%d %f %f %f\n",num_elem,xp[0],xp[1],xp[2]);
    
    if (num_elem != -1) {

      pic->xv[np*6+0]=xp[0];
      pic->xv[np*6+1]=xp[1];
      pic->xv[np*6+2]=xp[2];

      for(int idim=0;idim<3;idim++){
	vp[idim]=vt*(corput(n,k1[idim+3],k2[idim+3])-0.5);
      }

      pic->xv[np*6+3]=vp[0];
      pic->xv[np*6+4]=vp[1];
      pic->xv[np*6+5]=vp[2];
      
      pic->cell_id[np]=num_elem;

      np++;
    }

    n++;

    
  }


}


void PlotParticles(PIC* pic){
  
  FILE * gmshfile;
  FILE * gnufile;
  gmshfile = fopen("partplot.msh", "w" );
  gnufile = fopen("partplot.dat", "w" );

  float x,y,z,vx,vy,vz;
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(double));
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", pic->nbparts);
  /* fic << "$MeshFormat"<<endl; */
  /* fic << "2 0 8" << endl; */
  /* fic << "$EndMeshFormat"<<endl; */
  /* fic << "$Nodes" << endl; */
  
  /* fic << NbPart <<endl;       */
  /* cout << "NbPartFinal " << NbPart << endl; */

  for(int i=0;i<pic->nbparts;i++) {

    x=pic->xv[6*i+0];
    y=pic->xv[6*i+1];
    z=pic->xv[6*i+2];
    vx=pic->xv[6*i+3];
    vy=pic->xv[6*i+4];
    vz=pic->xv[6*i+5];
    fprintf(gmshfile,"%d %f %f %f \n",i+1,x,y,z);
    /* fic << i+1 << " "<<x<<" "<<y<<" "<<0<<endl;	 */
    fprintf(gnufile,"%f %f %f \n",x,y,z);
  }
  
  fprintf(gmshfile, "$EndNodes\n");
  //fic << "$EndNodes"<<endl;

  fclose(gmshfile);
  fclose(gnufile);

}
  


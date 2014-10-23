#include "field.h"
#include "geometry.h"
#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "global.h"
#include <math.h>

// param[0] = M
// param[1] = deg x
// param[2] = deg y
// param[3] = deg z
// param[4] = raf x
// param[5] = raf y
// param[6] = raf z
int GenericVarindex(int* param, int elem, int ipg, int iv){

  int npg= (param[1]+1)*(param[2]+1)*(param[3]+1)
    *param[4]*param[5]*param[6];

  return iv + param[0] * ( ipg + npg * elem);

}

void InitField(Field* f){

  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  double w[f->model.m];
  double xpg[3];
  double xref[3],omega;
  double dtau[3][3];
  double physnode[20][3];

  f->is2d = false;

  // a copy for avoiding too much "->"
  for(int ip=0;ip<8;ip++){
    f->interp_param[ip]=f->interp.interp_param[ip];
  }

  int nmem=f->model.m * f->macromesh.nbelems * 
    NPG(f->interp_param+1);
  printf("allocate %d doubles\n",nmem);
  f->wn=malloc(nmem * sizeof(double));
  assert(f->wn);	       
  f->wnp1=malloc(nmem * sizeof(double));
  assert(f->wnp1);	       
  f->dtwn=malloc(nmem * sizeof(double));
  assert(f->dtwn);	       

  f->tnow=0;

  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      ref_pg_vol(f->interp_param+1, ipg, xref, &omega,NULL);
      Ref2Phy(physnode,
	      xref,
	      0,-1, // dphiref,ifa
              xpg,dtau,  
	      NULL,NULL,NULL); // codtau,dphi,vnds
      // check the reverse transform at all the GLOPS
      double xref2[3];
      Phy2Ref(physnode,xpg,xref2);
      assert(Dist(xref,xref2) < 1e-8);
      
      f->model.InitData(xpg,w);
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	f->wn[imem]=w[iv];
      }
    }
  }

  // compute cfl parameter min_i vol_i/surf_i
  f->hmin=1e10;

  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    double codtau[3][3],dtau[3][3];
    double vol=0,surf=0;
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }
    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],wpg;
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      NULL,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	dtau[0][2]*codtau[0][2];
      vol+=wpg*det;  
    }
    for(int ifa=0;ifa<6;ifa++){
      // loop on the faces
      for(int ipgf=0;ipgf<NPGF(f->interp_param+1,ifa);ipgf++){
	double xpgref[3],wpg;
	//double xpgref2[3],wpg2;
	// get the coordinates of the Gauss point
	ref_pg_face(f->interp_param+1,ifa,ipgf,xpgref,&wpg,NULL);
	double vnds[3];
	Ref2Phy(physnode,
		xpgref,
		NULL,ifa, // dpsiref,ifa
		NULL,dtau,
		codtau,NULL,vnds); // codtau,dpsi,vnds
	surf+=sqrt(vnds[0]*vnds[0]+vnds[1]*vnds[1]+vnds[2]*vnds[2])*wpg;
      }
    }    
    f->hmin = f->hmin < vol/surf ? f->hmin : vol/surf;

  }

  // now take into account the polynomial degree
  int maxd=f->interp_param[1];
  maxd = maxd > f->interp_param[2] ? maxd : f->interp_param[2];
  maxd = maxd > f->interp_param[3] ? maxd : f->interp_param[3];
  
  f->hmin/=(maxd+1);

  printf("hmin=%f\n",f->hmin);

};

// display the field on screen
void DisplayField(Field* f){
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
 
  printf("Display field...\n");
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    printf("elem %d\n",ie);
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xref[3],xphy[3],wpg;
      ref_pg_vol(f->interp_param+1,ipg,xref,&wpg,NULL);

      printf("Gauss point %d %f %f %f \n",ipg,xref[0],xref[1],xref[2]);
    printf("dtw= ");
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	printf("%f ",f->dtwn[imem]);
      }
      printf("\n");
    printf("w= ");
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	printf("%f ",f->wn[imem]);
      }
      printf("\n");
    }
  }



};



// save the results in the gmsh format
// typplot: index of the plotted variable
// int compare == true -> compare with the exact value
void PlotField(int typplot,int compare,Field* f,char* filename){

  const int hexa64ref[3*64]={
    0,0,3,
    3,0,3,
    3,3,3,
    0,3,3,
    0,0,0,3,0,0,3,3,0,0,3,0,
    1,0,3,2,0,3,0,1,3,0,2,3,0,0,2,0,0,1,3,1,3,3,2,3,
    3,0,2,3,0,1,2,3,3,1,3,3,3,3,2,3,3,1,0,3,2,0,3,1,
    1,0,0,2,0,0,0,1,0,0,2,0,3,1,0,3,2,0,2,3,0,1,3,0,
    1,1,3,1,2,3,2,2,3,2,1,3,1,0,2,2,0,2,2,0,1,1,0,1,
    0,1,2,0,1,1,0,2,1,0,2,2,3,1,2,3,2,2,3,2,1,3,1,1,
    2,3,2,1,3,2,1,3,1,2,3,1,1,1,0,2,1,0,2,2,0,1,2,0,
    1,1,2,2,1,2,2,2,2,1,2,2,1,1,1,2,1,1,2,2,1,1,2,1};

  int* elem2nodes = f->macromesh.elem2node;
  double* node = f->macromesh.node;


  FILE * gmshfile;
  gmshfile = fopen( filename, "w" );

  // data plots
  int mw = f->model.m;
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  int nraf[3]={f->interp_param[4],f->interp_param[5],f->interp_param[6]};
  int deg[3]={f->interp_param[1],f->interp_param[2],f->interp_param[3]};
  // refinement size in each direction
  double hh[3]={1./nraf[0],1./nraf[1],1./nraf[2]};

  int npgv = NPG(f->interp_param+1);
  int nnodes = 20;

  double Xn[nnodes][3];
  double Xr[3];
  double Xphy[3];

  // header
  fprintf(gmshfile,"$MeshFormat\n2.2 0 %d\n",(int) sizeof(double));
  //int one=1;
  //fwrite((char*) &one,sizeof(int),1,gmshfile);
  fprintf(gmshfile,"$EndMeshFormat\n$Nodes\n%d\n",
	  f->macromesh.nbelems*nraf[0]*nraf[1]*nraf[2]*64);

  int nb_plotnodes=f->macromesh.nbelems*nraf[0]*nraf[1]*nraf[2]*64;
  double* value=malloc(nb_plotnodes*sizeof(double));
  assert(value);
  int nodecount=0;
  // nodes
  for(int i=0;i<f->macromesh.nbelems;i++){
    // get the nodes of element L
    for(int ino=0;ino<nnodes;ino++){
      int numnoe=elem2nodes[nnodes*i+ino];
      for(int ii=0;ii<3;ii++){
        Xn[ino][ii]=node[3*numnoe+ii];
      }
    }
    // loop on the macro elem subcells
    int icL[3];
    // loop on the subcells
    for(icL[0]=0;icL[0]<nraf[0];icL[0]++){
      for(icL[1]=0;icL[1]<nraf[1];icL[1]++){
	for(icL[2]=0;icL[2]<nraf[2];icL[2]++){
	  // get the left subcell id
	  int ncL=icL[0]+nraf[0]*(icL[1]+nraf[1]*icL[2]);
	  // first glop index in the subcell
	  //int offsetL=(deg[0]+1)*(deg[1]+1)*(deg[2]+1)*ncL;
	  
	  for(int ino=0;ino<64;ino++){
	    int nnoe=64*i+ino+1;
	    Xr[0]=(double) (hexa64ref[3*ino+0]) / 3;
	    Xr[1]=(double) (hexa64ref[3*ino+1]) / 3;
	    Xr[2]=(double) (hexa64ref[3*ino+2]) / 3;
	    
	    Xr[0] = icL[0]*hh[0]+ Xr[0] * hh[0];
	    Xr[1] = icL[1]*hh[1]+ Xr[1] * hh[1];
	    Xr[2] = icL[2]*hh[2]+ Xr[2] * hh[2];
	    
	    Ref2Phy(Xn,
		    Xr,
		    NULL,
		    -1,
		    Xphy,
		    NULL,
		    NULL,
		    NULL,
		    NULL);

	    double Xplot[3];
	    Xplot[0]=Xphy[0];
	    Xplot[1]=Xphy[1];
	    Xplot[2]=Xphy[2];

	    value[nodecount]=0;
	    for(int ib=0;ib<npgv;ib++){
	      double psi;
	      psi_ref_subcell(f->interp_param+1,icL, ib, Xr, &psi, NULL);
	      
	      int vi = f->varindex(f->interp_param, i, ib, typplot);
	      value[nodecount] += psi * f->wn[vi];
	    }

	    // compare with an
	    // exact solution
	    if (compare){
	      double wex[f->model.m];
	      f->model.ImposedData(Xphy,f->tnow,wex);
	      value[nodecount] -= wex[typplot];
	    }
	    nodecount++;


	    // fwrite((char*) &nnoe,sizeof(int),1,gmshfile);
	    // fwrite((char*) &(Xplot[0]),sizeof(double),1,gmshfile);
	    // fwrite((char*) &(Xplot[1]),sizeof(double),1,gmshfile);
	    // fwrite((char*) &(Xplot[2]),sizeof(double),1,gmshfile);
	    fprintf(gmshfile,"%d %f %f %f\n",nnoe,Xplot[0],Xplot[1],Xplot[2]);
	    
	  }
	}
      }
    }
  }

  fprintf(gmshfile,"$EndNodes\n");



  // elements
  fprintf(gmshfile,"$Elements\n");
  fprintf(gmshfile,"%d\n",f->macromesh.nbelems*nraf[0]*nraf[1]*nraf[2]);



  int elm_type=92;
  //int num_elm_follow=f->macromesh.nbelems;
  int num_tags=0;

  // fwrite((char*) &elm_type,sizeof(int),1,gmshfile);
  // fwrite((char*) &num_elm_follow,sizeof(int),1,gmshfile);
  // fwrite((char*) &num_tags,sizeof(int),1,gmshfile);

  for(int i=0;i<f->macromesh.nbelems;i++){

    // loop on the macro elem subcells
    int icL[3];
    // loop on the subcells
    for(icL[0]=0;icL[0]<nraf[0];icL[0]++){
      for(icL[1]=0;icL[1]<nraf[1];icL[1]++){
	for(icL[2]=0;icL[2]<nraf[2];icL[2]++){
	  // get the subcell id
	  int ncL=icL[0]+nraf[0]*(icL[1]+nraf[1]*icL[2]);
	  // first glop index in the subcell
	  //int offsetL=(deg[0]+1)*(deg[1]+1)*(deg[2]+1)*ncL;

	  // global subcell id
	  int numelem=ncL+i*nraf[0]*nraf[1]*nraf[2]+1;
	  //fwrite((char*) &numelem,sizeof(int),1,gmshfile);
	  fprintf(gmshfile,"%d ",numelem);
	  fprintf(gmshfile,"%d ",elm_type);
	  fprintf(gmshfile,"%d ",num_tags);
    
    
	  for(int ii=0;ii<64;ii++){
	    int numnoe=64*(i*nraf[0]*nraf[1]*nraf[2]+ncL) + ii +1;
	    //fwrite((char*) &numnoe,sizeof(int),1,gmshfile);
	    fprintf(gmshfile,"%d ",numnoe);
	  }
	  fprintf(gmshfile,"\n");
	}
      }
    }
  }
  
  fprintf(gmshfile,"$EndElements\n");

  fclose(gmshfile);
  assert(1==2);
  // now display data
    
    fprintf(gmshfile,"$NodeData\n");
    fprintf(gmshfile,"1\n");
    fprintf(gmshfile,"\"Field %d\"\n",typplot);

    double t = 0;

    fprintf(gmshfile,"1\n%f\n3\n0\n1\n",t);

    fprintf(gmshfile,"%d\n",64*f->macromesh.nbelems);

    
    for(int i=0;i<f->macromesh.nbelems;i++){
      for(int ino=0;ino<20;ino++){
	int numnoe=elem2nodes[nnodes*i+ino];
	for(int ii=0;ii<3;ii++){
	  Xn[ino][ii]=node[3*numnoe+ii];
	}
      }
      
      // data at the eight nodes
      for(int ii=0;ii<64;ii++){
	int nodenumber=64*i + ii +1;
	
	Xr[0]=(double) (hexa64ref[3*ii+0]) / 3;
	Xr[1]=(double) (hexa64ref[3*ii+1]) / 3;
	Xr[2]=(double) (hexa64ref[3*ii+2]) / 3;
	
	Ref2Phy(Xn,
		Xr,
		NULL,
		-1,
		Xphy,
		NULL,
		NULL,
		NULL,
		NULL);
	

	double value=0;
	for(int ib=0;ib<npgv;ib++){
	  double psi;
	  psi_ref(f->interp_param+1, ib, Xr, &psi, NULL);
	  
	  int vi = f->varindex(f->interp_param, i, ib, typplot);
	  value += psi * f->wn[vi];
	}

	// compare with an
	// exact solution
        if (compare){
          double wex[f->model.m];
          f->model.ImposedData(Xphy,f->tnow,wex);
          value -= wex[typplot];
        }


	//fwrite(const void *ptr, size_t size_of_elements,
	// size_t number_of_elements, FILE *a_file);
	//fwrite((char*) &nodenumber, sizeof(int),1,gmshfile);
	//fwrite((char*) &value, sizeof(double),1,gmshfile);
	//fprintf(gmshfile,"%d %f\n",nodenumber,value);
	fprintf(gmshfile,"%d %f\n",nodenumber,value);
      }

    }

    fprintf(gmshfile,"\n$EndNodeData\n");

    
    fclose(gmshfile);
  

}

// inter-subcell fluxes
void DGSubCellInterface(Field* f){

  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }


    int nraf[3]={f->interp_param[4],f->interp_param[5],f->interp_param[6]};
    int deg[3]={f->interp_param[1],f->interp_param[2],f->interp_param[3]};

    int icL[3];
    // loop on the subcells
    for(icL[0]=0;icL[0]<nraf[0];icL[0]++){
      for(icL[1]=0;icL[1]<nraf[1];icL[1]++){
	for(icL[2]=0;icL[2]<nraf[2];icL[2]++){
	  // get the left subcell id
	  int ncL=icL[0]+nraf[0]*(icL[1]+nraf[1]*icL[2]);
	  // first glop index in the subcell
	  int offsetL=(deg[0]+1)*(deg[1]+1)*(deg[2]+1)*ncL;
	  // sweeping subcell faces in the three directions
	  for(int dim0=0;dim0<3;dim0++){
	    // compute the subface flux only
	    // if we do not touch the subcell boundary
	    // along the current direction dim0
	    if (icL[dim0] != nraf[dim0]-1) {
	      int icR[3]={icL[0],icL[1],icL[2]};
	      // the right cell index corresponds to 
	      // an increment in the dim0 direction
	      icR[dim0]++;
	      int ncR=icR[0]+nraf[0]*(icR[1]+nraf[1]*icR[2]);
	      int offsetR=(deg[0]+1)*(deg[1]+1)*(deg[2]+1)*ncR;
	      // now loop on the left glops of the subface
	      int iL[3],dim1=(dim0+1)%3,dim2=(dim0+2)%3;
	      iL[dim0]=deg[dim0];
	      for(iL[dim1]=0;iL[dim1]<deg[dim1]+1;iL[dim1]++){
		for(iL[dim2]=0;iL[dim2]<deg[dim2]+1;iL[dim2]++){
		  // find the right and left glops volume indices
		  int iR[3]={iL[0],iL[1],iL[2]};
		  iR[dim0]=0;
		  int ipgL=offsetL+iL[0]+(deg[0]+1)*(iL[1]+(deg[1]+1)*iL[2]);
		  int ipgR=offsetR+iR[0]+(deg[0]+1)*(iR[1]+(deg[1]+1)*iR[2]);
		  //printf("ipgL=%d ipgR=%d\n",ipgL,ipgR);
		  double xref[3],wpg;		  
		  ref_pg_vol(f->interp_param+1,ipgL,xref,&wpg,NULL);
		  // mapping from the ref glop to the physical glop
		  double dtau[3][3],codtau[3][3];
		  Ref2Phy(physnode,
			  xref,
			  NULL,  // dphiref
			  -1,    // ifa                                 
			  NULL,  // xphy  
			  dtau,
			  codtau,
			  NULL,  // dphi
			  NULL);  // vnds       
		  // we compute ourself the normal vector because
		  // we have to take into account the subcell surface
		  double vnds[3];
		  double h1h2=1./nraf[dim1]/nraf[dim2];
		  vnds[0]=codtau[0][dim0]*h1h2;
		  vnds[1]=codtau[1][dim0]*h1h2;
		  vnds[2]=codtau[2][dim0]*h1h2;
		  // numerical flux from the left and right state and
		  // normal vector
		  double wL[f->model.m],wR[f->model.m],flux[f->model.m];
		  for(int iv=0;iv<f->model.m;iv++){
		    int imem=f->varindex(f->interp_param,ie,ipgL,iv);
		    wL[iv]=f->wn[imem];
		    imem=f->varindex(f->interp_param,ie,ipgR,iv);
		    wR[iv]=f->wn[imem];
		  }
		  f->model.NumFlux(wL,wR,vnds,flux);
		  /* printf("vnds %f %f %f flux %f wpg %f\n", */
		  /* 	 vnds[0],vnds[1],vnds[2], */
		  /* 	 flux[0],wpg); */
		  // subcell ref surface glop weight
		  wpg=wglop(deg[dim1],iL[dim1])*wglop(deg[dim2],iL[dim2]);
		  // finally distribute the flux on the two sides
		  for(int iv=0;iv<f->model.m;iv++){
		    int imem=f->varindex(f->interp_param,ie,ipgL,iv);
		    f->dtwn[imem]-=flux[iv]*wpg;
		  }
		  for(int iv=0;iv<f->model.m;iv++){
		    int imem=f->varindex(f->interp_param,ie,ipgR,iv);
		    f->dtwn[imem]+=flux[iv]*wpg;
		  }

		  
		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop
	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

  } // macro elem loop

}




// compute the Discontinuous Galerkin inter-macrocells boundary terms
void DGMacroCellInterface(Field* f){

 // init to zero the time derivative
  int sizew=0;
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	f->dtwn[imem]=0;
        sizew++;
      }
    }
  }
  assert(sizew==f->macromesh.nbelems * f->model.m * NPG(f->interp_param+1));

  // assembly of the surface terms
  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the 6 faces
    // or four faces for 2d computations
    int nbfa=6;
    if (f->is2d) nbfa=4;
    for(int ifa=0;ifa<nbfa;ifa++){
      // get the right elem or the boundary id
      int ieR=f->macromesh.elem2elem[6*ie+ifa];
      double physnodeR[20][3];
      if (ieR >= 0) {
      	for(int inoloc=0;inoloc<20;inoloc++){
      	  int ino=f->macromesh.elem2node[20*ieR+inoloc];
      	  physnodeR[inoloc][0]=f->macromesh.node[3*ino+0];
      	  physnodeR[inoloc][1]=f->macromesh.node[3*ino+1];
      	  physnodeR[inoloc][2]=f->macromesh.node[3*ino+2];
      	}
      }
      
      // loop on the glops (numerical integration)
      // of the face ifa
      for(int ipgf=0;ipgf<NPGF(f->interp_param+1,ifa);ipgf++){
  	double xpgref[3],xpgref_in[3],wpg;
  	//double xpgref2[3],wpg2;
  	// get the coordinates of the Gauss point
	// and coordinates of a point slightly inside the
	// opposite element in xref_in
  	ref_pg_face(f->interp_param+1,ifa,ipgf,xpgref,&wpg,xpgref_in);

  	// recover the volume gauss point from
  	// the face index
  	int ipg=f->interp_param[7];
  	// get the left value of w at the gauss point
  	double wL[f->model.m],wR[f->model.m];
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(f->interp_param,ie,ipg,iv);
  	  wL[iv]=f->wn[imem];
  	}
  	// the basis functions is also the gauss point index
  	int ib=ipg;
  	// normal vector at gauss point ipg
  	double dtau[3][3],codtau[3][3],xpg[3];
  	double vnds[3];
  	Ref2Phy(physnode,
  		xpgref,
  		NULL,ifa, // dpsiref,ifa
  		xpg,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds
  	double flux[f->model.m];
  	if (ieR >=0) {  // the right element exists
  	  // find the corresponding point in the right elem
	  double xpg_in[3];
  	Ref2Phy(physnode,
  		xpgref_in,
  		NULL,ifa, // dpsiref,ifa
  		xpg_in,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds
  	  double xref[3];
	  Phy2Ref(physnodeR,xpg_in,xref);
  	  int ipgR=ref_ipg(f->interp_param+1,xref);
	  double xpgR[3],xrefR[3],wpgR;
	  ref_pg_vol(f->interp_param+1, ipgR, xrefR, &wpgR,NULL);
	  Ref2Phy(physnodeR,
		  xrefR,
		  NULL,-1, // dphiref,ifa
		  xpgR,NULL,  
		  NULL,NULL,NULL); // codtau,dphi,vnds
	  assert(Dist(xpgR,xpg)<1e-10);
  	  for(int iv=0;iv<f->model.m;iv++){
  	    int imem=f->varindex(f->interp_param,ieR,ipgR,iv);
  	    wR[iv]=f->wn[imem];
  	  }
  	  // int_dL F(wL,wR,grad phi_ib )
  	  f->model.NumFlux(wL,wR,vnds,flux);
  	}
  	else { //the right element does not exist
  	  f->model.BoundaryFlux(xpg,f->tnow,wL,vnds,flux);
  	}
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(f->interp_param,ie,ib,iv);
  	  f->dtwn[imem]-=flux[iv]*wpg;
  	}
	
      }

    }
  }



}


// compute the Discontinuous Galerkin volume terms
void DGVolume(Field* f){

  // assembly of the volume terms
  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // mass matrix
    double masspg[NPG(f->interp_param+1)];
    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],wpg;
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);

      // get the value of w at the gauss point
      double w[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // loop on the basis functions
      for(int ib=0;ib<NPG(f->interp_param+1);ib++){
	// gradient of psi_ib at gauss point ipg
	double dpsiref[3],dpsi[3];
	double dtau[3][3],codtau[3][3];//,xpg[3];
	grad_psi_pg(f->interp_param+1,ib,ipg,dpsiref);
	Ref2Phy(physnode, // phys. nodes
		xpgref,  // xref
		dpsiref,-1, // dpsiref,ifa
		NULL,dtau,  // xphy,dtau
		codtau,dpsi,NULL); // codtau,dpsi,vnds
	// remember the diagonal mass term
	if (ib == ipg){
	  double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	    dtau[0][2]*codtau[0][2];
	  masspg[ipg]=wpg*det;
	}
	// int_L F(w,w,grad phi_ib )
	double flux[f->model.m];
	f->model.NumFlux(w,w,dpsi,flux);
	for(int iv=0;iv<f->model.m;iv++){
	  int imem=f->varindex(f->interp_param,ie,ib,iv);
	  f->dtwn[imem]+=flux[iv]*wpg;
	}
      }
    }

    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      // apply the inverse of the diagonal mass matrix
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	(f->dtwn[imem])/=masspg[ipg];
      }
    }

  }

  

}


// apply the Discontinuous Galerkin approximation for computing
// the time derivative of the field
void dtField(Field* f){

  // interpolation params
  // warning: this is ugly, but the last
  // parameter is used for computing the volume
  // GLOP index from the face GLOP index...
  // ugly too: the first parameter is not used by all
  // utilities. we have sometimes to jump over : pass param+1
  // instead of param...
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};


  // init to zero the time derivative
  int sizew=0;
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	f->dtwn[imem]=0;
        sizew++;
      }
    }
  }
  assert(sizew==f->macromesh.nbelems * f->model.m * NPG(f->interp_param+1));

  // assembly of the surface terms
  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the 6 faces
    // or four faces for 2d computations
    int nbfa=6;
    if (f->is2d) nbfa=4;
    for(int ifa=0;ifa<nbfa;ifa++){
      // get the right elem or the boundary id
      int ieR=f->macromesh.elem2elem[6*ie+ifa];
      double physnodeR[20][3];
      if (ieR >= 0) {
      	for(int inoloc=0;inoloc<20;inoloc++){
      	  int ino=f->macromesh.elem2node[20*ieR+inoloc];
      	  physnodeR[inoloc][0]=f->macromesh.node[3*ino+0];
      	  physnodeR[inoloc][1]=f->macromesh.node[3*ino+1];
      	  physnodeR[inoloc][2]=f->macromesh.node[3*ino+2];
      	}
      }
      
      // loop on the glops (numerical integration)
      // of the face ifa
      for(int ipgf=0;ipgf<NPGF(f->interp_param+1,ifa);ipgf++){
  	double xpgref[3],xpgref_in[3],wpg;
  	//double xpgref2[3],wpg2;
  	// get the coordinates of the Gauss point
	// and coordinates of a point slightly inside the
	// opposite element in xref_in
  	ref_pg_face(f->interp_param+1,ifa,ipgf,xpgref,&wpg,xpgref_in);

  	// recover the volume gauss point from
  	// the face index
  	int ipg=f->interp_param[7];
  	// get the left value of w at the gauss point
  	double wL[f->model.m],wR[f->model.m];
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(f->interp_param,ie,ipg,iv);
  	  wL[iv]=f->wn[imem];
  	}
  	// the basis functions is also the gauss point index
  	int ib=ipg;
  	// normal vector at gauss point ipg
  	double dtau[3][3],codtau[3][3],xpg[3];
  	double vnds[3];
  	Ref2Phy(physnode,
  		xpgref,
  		NULL,ifa, // dpsiref,ifa
  		xpg,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds
  	double flux[f->model.m];
  	if (ieR >=0) {  // the right element exists
  	  // find the corresponding point in the right elem
	  double xpg_in[3];
  	Ref2Phy(physnode,
  		xpgref_in,
  		NULL,ifa, // dpsiref,ifa
  		xpg_in,dtau,
  		codtau,NULL,vnds); // codtau,dpsi,vnds
  	  double xref[3];
	  Phy2Ref(physnodeR,xpg_in,xref);
  	  int ipgR=ref_ipg(f->interp_param+1,xref);
	  double xpgR[3],xrefR[3],wpgR;
	  ref_pg_vol(f->interp_param+1, ipgR, xrefR, &wpgR,NULL);
	  Ref2Phy(physnodeR,
		  xrefR,
		  NULL,-1, // dphiref,ifa
		  xpgR,NULL,  
		  NULL,NULL,NULL); // codtau,dphi,vnds
	  assert(Dist(xpgR,xpg)<1e-10);
  	  for(int iv=0;iv<f->model.m;iv++){
  	    int imem=f->varindex(f->interp_param,ieR,ipgR,iv);
  	    wR[iv]=f->wn[imem];
  	  }
  	  // int_dL F(wL,wR,grad phi_ib )
  	  f->model.NumFlux(wL,wR,vnds,flux);
  	}
  	else { //the right element does not exist
  	  f->model.BoundaryFlux(xpg,f->tnow,wL,vnds,flux);
  	}
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(f->interp_param,ie,ib,iv);
  	  f->dtwn[imem]-=flux[iv]*wpg;
  	}
	
      }

    }
  }

  // assembly of the volume terms
  // loop on the elements
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // mass matrix
    double masspg[NPG(f->interp_param+1)];
    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],wpg;
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);

      // get the value of w at the gauss point
      double w[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // loop on the basis functions
      for(int ib=0;ib<NPG(f->interp_param+1);ib++){
	// gradient of psi_ib at gauss point ipg
	double dpsiref[3],dpsi[3];
	double dtau[3][3],codtau[3][3];//,xpg[3];
	grad_psi_pg(f->interp_param+1,ib,ipg,dpsiref);
	Ref2Phy(physnode, // phys. nodes
		xpgref,  // xref
		dpsiref,-1, // dpsiref,ifa
		NULL,dtau,  // xphy,dtau
		codtau,dpsi,NULL); // codtau,dpsi,vnds
	// remember the diagonal mass term
	if (ib == ipg){
	  double det=dtau[0][0]*codtau[0][0]+dtau[0][1]*codtau[0][1]+
	    dtau[0][2]*codtau[0][2];
	  masspg[ipg]=wpg*det;
	}
	// int_L F(w,w,grad phi_ib )
	double flux[f->model.m];
	f->model.NumFlux(w,w,dpsi,flux);
	for(int iv=0;iv<f->model.m;iv++){
	  int imem=f->varindex(f->interp_param,ie,ib,iv);
	  f->dtwn[imem]+=flux[iv]*wpg;
	}
      }
    }

    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      // apply the inverse of the diagonal mass matrix
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	(f->dtwn[imem])/=masspg[ipg];
      }
    }

  }
    
};

// time integration by a second order Runge-Kutta algorithm 
void RK2(Field* f,double tmax){

  double vmax=1; // to be changed for another model !!!!!!!!!
  double cfl=0.05;

  double dt = cfl * f->hmin / vmax;
  int itermax=tmax/dt+1;
  int freq=itermax/10;

  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  int sizew=f->macromesh.nbelems * f->model.m * NPG(f->interp_param+1);
 
  int iter=0;

  while(f->tnow<tmax){
    if (iter%freq==0)
      printf("t=%f iter=%d/%d dt=%f\n",f->tnow,iter,itermax,dt);
    // predictor
    dtField(f);
    for(int iw=0;iw<sizew;iw++){
      f->wnp1[iw]=f->wn[iw]+ dt/2 * f->dtwn[iw]; 
    }
    //exchange the field pointers 
    double *temp;
    temp=f->wnp1;
    f->wnp1=f->wn;
    f->wn=temp;

    // corrector
    f->tnow+=dt/2;
    dtField(f);
    for(int iw=0;iw<sizew;iw++){
      f->wnp1[iw]+=dt*f->dtwn[iw];
    }
    f->tnow+=dt/2;
    iter++;
    //exchange the field pointers 
    temp=f->wnp1;
    f->wnp1=f->wn;
    f->wn=temp;

  }
  printf("t=%f iter=%d/%d dt=%f\n",f->tnow,iter,itermax,dt);

}
// time integration by a second order Runge-Kutta algorithm
//  with memory copy instead of pointers exchange 
void RK2Copy(Field* f,double tmax){

  double vmax=1; // to be changed for another model !!!!!!!!!
  double cfl=0.05;

  double dt = cfl * f->hmin / vmax;

  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  int sizew=f->macromesh.nbelems * f->model.m * NPG(f->interp_param+1);
 
  int iter=0;

  while(f->tnow<tmax){
    printf("t=%f iter=%d dt=%f\n",f->tnow,iter,dt);
    // predictor
    dtField(f);
    for(int iw=0;iw<sizew;iw++){
      f->wnp1[iw]=f->wn[iw]+ dt/2 * f->dtwn[iw]; 
    }
    //exchange the field pointers 
    for(int iw=0;iw<sizew;iw++){
      double temp=f->wn[iw];
      f->wn[iw]=f->wnp1[iw];
      f->wnp1[iw]=temp;
    }
    // corrector
    f->tnow+=dt/2;
    dtField(f);
    for(int iw=0;iw<sizew;iw++){
      f->wnp1[iw]+=dt*f->dtwn[iw];
    }
    f->tnow+=dt/2;
    iter++;
    //exchange the field pointers 
    for(int iw=0;iw<sizew;iw++){
      f->wn[iw]=f->wnp1[iw];
    }
 
  }
}

// compute the normalized L2 distance with the imposed data
double L2error(Field* f){

  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  double error=0;
  double moy=0; // mean value
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],xphy[3],wpg;
      double dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
		xpgref,  // xref
		NULL,-1, // dpsiref,ifa
		xphy,dtau,  // xphy,dtau
		codtau,NULL,NULL); // codtau,dpsi,vnds
      double det=dtau[0][0]*codtau[0][0]+
        dtau[0][1]*codtau[0][1]+dtau[0][2]*codtau[0][2]; 
      double w[f->model.m],wex[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      f->model.ImposedData(xphy,f->tnow,wex);
      for(int iv=0;iv<f->model.m;iv++){
        error+=pow(w[iv]-wex[iv],2)*wpg*det;
        moy+=pow(w[iv],2)*wpg*det;
      }
    }
  }
  return sqrt(error)/sqrt(moy);
}

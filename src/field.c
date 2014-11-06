#include "field.h"
#include "geometry.h"
#include "interpolation.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "global.h"
#include <math.h>
#include <pthread.h>

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
      double dtau[3][3];
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
      double codtau[3][3],dtau[3][3];
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      NULL,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2];
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
	double codtau[3][3],dtau[3][3];
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

  // now take into account the polynomial degree and the refinement
  int maxd=f->interp_param[1];
  maxd = maxd > f->interp_param[2] ? maxd : f->interp_param[2];
  maxd = maxd > f->interp_param[3] ? maxd : f->interp_param[3];
  
  f->hmin/=((maxd+1)*f->interp_param[4]);

  printf("hmin=%f\n",f->hmin);

};

// display the field on screen
void DisplayField(Field* f){
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
 
  printf("Display field...\n");
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    printf("elem %d\n",ie);
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xref[3],wpg;
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
  //int param[8]={f->model.m,_DEGX,_DEGY,_DEGZ,_RAFX,_RAFY,_RAFZ,0};
  int nraf[3]={f->interp_param[4],f->interp_param[5],f->interp_param[6]};
  // refinement size in each direction
  double hh[3]={1./nraf[0],1./nraf[1],1./nraf[2]};

  int npgv = NPG(f->interp_param+1);
  int nnodes = 20;

  double physnode[nnodes][3];
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
        physnode[ino][ii]=node[3*numnoe+ii];
      }
    }
    // loop on the macro elem subcells
    int icL[3];
    // loop on the subcells
    for(icL[0]=0;icL[0]<nraf[0];icL[0]++){
      for(icL[1]=0;icL[1]<nraf[1];icL[1]++){
	for(icL[2]=0;icL[2]<nraf[2];icL[2]++){
	  // get the left subcell id
	  // first glop index in the subcell
	  //int offsetL=(deg[0]+1)*(deg[1]+1)*(deg[2]+1)*ncL;
	  
	  for(int ino=0;ino<64;ino++){
	    Xr[0]=(double) (hexa64ref[3*ino+0]) / 3;
	    Xr[1]=(double) (hexa64ref[3*ino+1]) / 3;
	    Xr[2]=(double) (hexa64ref[3*ino+2]) / 3;
	    
	    Xr[0] = icL[0]*hh[0]+ Xr[0] * hh[0];
	    Xr[1] = icL[1]*hh[1]+ Xr[1] * hh[1];
	    Xr[2] = icL[2]*hh[2]+ Xr[2] * hh[2];

	    for(int ii=0;ii<3;ii++){
	      assert(Xr[ii]<1+1e-10 && Xr[ii]>-1e-10);
	    }
	    
	    Ref2Phy(physnode,
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
	    double testpsi=0;
	    for(int ib=0;ib<npgv;ib++){
	      double psi;
	      psi_ref_subcell(f->interp_param+1,icL, ib, Xr, &psi, NULL);
	      testpsi+=psi;
	      int vi = f->varindex(f->interp_param, i, ib, typplot);
	      value[nodecount] += psi * f->wn[vi];
	    }
	    assert(fabs(testpsi-1)<1e-10);

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
	    fprintf(gmshfile,"%d %f %f %f\n",nodecount,Xplot[0],Xplot[1],Xplot[2]);
	    
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

  // now display data
    
  fprintf(gmshfile,"$NodeData\n");
  fprintf(gmshfile,"1\n");
  fprintf(gmshfile,"\"Field %d\"\n",typplot);
  
  double t = 0;

  fprintf(gmshfile,"1\n%f\n3\n0\n1\n",t);

  fprintf(gmshfile,"%d\n",nb_plotnodes);

    
  for(int ino=0;ino<nb_plotnodes;ino++){

    //fwrite(const void *ptr, size_t size_of_elements,
    // size_t number_of_elements, FILE *a_file);
    //fwrite((char*) &nodenumber, sizeof(int),1,gmshfile);
    //fwrite((char*) &value, sizeof(double),1,gmshfile);
    //fprintf(gmshfile,"%d %f\n",nodenumber,value);
    fprintf(gmshfile,"%d %f\n",ino+1,value[ino]);      

  }
    
  /* for(int i=0;i<f->macromesh.nbelems;i++){ */
  /*   for(int ino=0;ino<20;ino++){ */
  /* 	int numnoe=elem2nodes[nnodes*i+ino]; */
  /* 	for(int ii=0;ii<3;ii++){ */
  /* 	  physnode[ino][ii]=node[3*numnoe+ii]; */
  /* 	} */
  /*   } */
      
  /*   // data at the eight nodes */
  /*   for(int ii=0;ii<64;ii++){ */
  /* 	int nodenumber=64*i + ii +1; */
	
  /* 	Xr[0]=(double) (hexa64ref[3*ii+0]) / 3; */
  /* 	Xr[1]=(double) (hexa64ref[3*ii+1]) / 3; */
  /* 	Xr[2]=(double) (hexa64ref[3*ii+2]) / 3; */
	
  /* 	Ref2Phy(physnode, */
  /* 		Xr, */
  /* 		NULL, */
  /* 		-1, */
  /* 		Xphy, */
  /* 		NULL, */
  /* 		NULL, */
  /* 		NULL, */
  /* 		NULL); */
	

  /* 	double value=0; */
  /* 	for(int ib=0;ib<npgv;ib++){ */
  /* 	  double psi; */
  /* 	  psi_ref(f->interp_param+1, ib, Xr, &psi, NULL); */
	  
  /* 	  int vi = f->varindex(f->interp_param, i, ib, typplot); */
  /* 	  value += psi * f->wn[vi]; */
  /* 	} */

  /* 	// compare with an */
  /* 	// exact solution */
  /*     if (compare){ */
  /*       double wex[f->model.m]; */
  /*       f->model.ImposedData(Xphy,f->tnow,wex); */
  /*       value -= wex[typplot]; */
  /*     } */


  /* 	//fwrite(const void *ptr, size_t size_of_elements, */
  /* 	// size_t number_of_elements, FILE *a_file); */
  /* 	//fwrite((char*) &nodenumber, sizeof(int),1,gmshfile); */
  /* 	//fwrite((char*) &value, sizeof(double),1,gmshfile); */
  /* 	//fprintf(gmshfile,"%d %f\n",nodenumber,value); */
  /* 	fprintf(gmshfile,"%d %f\n",nodenumber,value); */
  /*   } */

  /* } */

  fprintf(gmshfile,"\n$EndNodeData\n");

    
  fclose(gmshfile);
 
  free(value);

}

// inter-subcell fluxes
void* DGSubCellInterface(void* mc){

  MacroCell* mcell = (MacroCell*) mc;

  Field* f= mcell->field;

  // loop on the elements
  for (int ie=mcell->first_cell;ie<mcell->last_cell_p1;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }


    const int nraf[3]={f->interp_param[4],
		       f->interp_param[5],
		       f->interp_param[6]};
    const int deg[3]={f->interp_param[1],
		      f->interp_param[2],
		      f->interp_param[3]};
    const int npg[3] = {deg[0]+1,
			deg[1]+1,
			deg[2]+1};
    const int m = f->model.m;


    // loop on the subcells
    //#pragma omp parallel for collapse(3)
    for(int icL0 = 0; icL0 < nraf[0]; icL0++){
      for(int icL1 = 0; icL1 < nraf[1]; icL1++){
	for(int icL2 = 0; icL2 < nraf[2]; icL2++){

	  int icL[3]={icL0,icL1,icL2};

	  // get the left subcell id
	  int ncL=icL[0]+nraf[0]*(icL[1]+nraf[1]*icL[2]);
	  // first glop index in the subcell
	  int offsetL=npg[0]*npg[1]*npg[2]*ncL;

	  // sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++){
	    // compute the subface flux only
	    // if we do not touch the subcell boundary
	    // along the current direction dim0
	    if (icL[dim0] != nraf[dim0]-1) {
	      int icR[3]={icL[0],icL[1],icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR=icR[0]+nraf[0]*(icR[1]+nraf[1]*icR[2]);
	      int offsetR=npg[0]*npg[1]*npg[2]*ncR;

	      // FIXME: write only write to L-values (and do both
	      // faces) to parallelise better.

	      const int altdim1[3]={1,0,0};
	      const int altdim2[3]={2,2,1};

	      // now loop on the left glops of the subface
	      //int dim1=(dim0+1)%3, dim2=(dim0+2)%3;
	      int dim1=altdim1[dim0], dim2=altdim2[dim0];
	      int iL[3];
	      iL[dim0] = deg[dim0];
	      for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++){
		for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++){
		  // find the right and left glops volume indices
		  
		  int iR[3] = {iL[0],iL[1],iL[2]};
		  iR[dim0] = 0;

		  int ipgL=offsetL+iL[0]+(deg[0]+1)*(iL[1]+(deg[1]+1)*iL[2]);
		  int ipgR=offsetR+iR[0]+(deg[0]+1)*(iR[1]+(deg[1]+1)*iR[2]);
		  //printf("ipgL=%d ipgR=%d\n",ipgL,ipgR);

		  // Compute the normal vector for integrating on the
		  // face
		  double vnds[3];
		  {
		    double xref[3], wpg3;		  
		    ref_pg_vol(f->interp_param+1,ipgL,xref,&wpg3,NULL);
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
		    // we compute ourself the normal vector because we
		    // have to take into account the subcell surface
		    
		    double h1h2=1./nraf[dim1]/nraf[dim2];
		    vnds[0] = codtau[0][dim0]*h1h2;
		    vnds[1] = codtau[1][dim0]*h1h2;
		    vnds[2] = codtau[2][dim0]*h1h2;
		  }

		  // numerical flux from the left and right state and
		  // normal vector
		  double wL[m],wR[m],flux[m];
		  for(int iv=0; iv < m; iv++){
		    int imemL=f->varindex(f->interp_param,ie,ipgL,iv);
		    int imemR=f->varindex(f->interp_param,ie,ipgR,iv);
		    wL[iv] = f->wn[imemL];
		    wR[iv] = f->wn[imemR];
		  }
		  f->model.NumFlux(wL,wR,vnds,flux);

		  // subcell ref surface glop weight
		  double wpg
		    = wglop(deg[dim1],iL[dim1]) 
		    * wglop(deg[dim2],iL[dim2]);

		  /* printf("vnds %f %f %f flux %f wpg %f\n", */
		  /* 	 vnds[0],vnds[1],vnds[2], */
		  /* 	 flux[0],wpg); */
		  
		  // finally distribute the flux on the two sides
		  for(int iv=0; iv < m; iv++){
		    int imemL = f->varindex(f->interp_param, ie, ipgL, iv);
		    int imemR = f->varindex(f->interp_param, ie, ipgR, iv);
		    f->dtwn[imemL] -= flux[iv] * wpg;
		    f->dtwn[imemR] += flux[iv] * wpg;
		  }
		  
		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop
	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

  } // macro elem loop

  return NULL;

}




// compute the Discontinuous Galerkin inter-macrocells boundary terms
void* DGMacroCellInterface(void* mc){

  MacroCell* mcell = (MacroCell*) mc;

  Field* f= mcell->field;

  int iparam[8];
  for(int ip=0;ip<8;ip++) iparam[ip]=f->interp_param[ip];
    
  // init to zero the time derivative
  for (int ie=mcell->first_cell;ie<mcell->last_cell_p1;ie++){
    for(int ipg=0;ipg<NPG(iparam+1);ipg++){
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(iparam,ie,ipg,iv);
	f->dtwn[imem]=0;
      }
    }
  }
  //assert(sizew==f->macromesh.nbelems * f->model.m * NPG(iparam+1));

  // assembly of the surface terms
  // loop on the elements
  for (int ie=mcell->first_cell;ie<mcell->last_cell_p1;ie++){
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
	//      for(int ipgf=0;ipgf<NPGF(iparam+1,ifa);ipgf++){ // FIXME?

  	double xpgref[3],xpgref_in[3],wpg;
  	//double xpgref2[3],wpg2;
  	// get the coordinates of the Gauss point
	// and coordinates of a point slightly inside the
	// opposite element in xref_in
  	ref_pg_face(iparam+1,ifa,ipgf,xpgref,&wpg,xpgref_in);

  	// recover the volume gauss point from
  	// the face index
  	int ipg=iparam[7];
  	// get the left value of w at the gauss point
  	double wL[f->model.m],wR[f->model.m];
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(iparam,ie,ipg,iv);
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
  	  int ipgR=ref_ipg(iparam+1,xref);
	  double xpgR[3],xrefR[3],wpgR;
	  ref_pg_vol(iparam+1, ipgR, xrefR, &wpgR,NULL);
	  Ref2Phy(physnodeR,
		  xrefR,
		  NULL,-1, // dphiref,ifa
		  xpgR,NULL,  
		  NULL,NULL,NULL); // codtau,dphi,vnds

	  assert(Dist(xpgR,xpg)<1e-10);
  	  for(int iv=0;iv<f->model.m;iv++){
  	    int imem=f->varindex(iparam,ieR,ipgR,iv);
  	    wR[iv]=f->wn[imem];
  	  }
  	  // int_dL F(wL,wR,grad phi_ib )
  	  f->model.NumFlux(wL,wR,vnds,flux);

   
  	}
  	else { //the right element does not exist
  	  f->model.BoundaryFlux(xpg,f->tnow,wL,vnds,flux);
  	}
  	for(int iv=0;iv<f->model.m;iv++){
  	  int imem=f->varindex(iparam,ie,ib,iv);
  	  f->dtwn[imem]-=flux[iv]*wpg;
  	}
	
      }

    }
  }

  return NULL;
}

// apply division by the mass matrix
void* DGMass(void* mc){

  MacroCell* mcell = (MacroCell*) mc;

  Field* f= mcell->field;

  // loop on the elements
  for (int ie=mcell->first_cell;ie<mcell->last_cell_p1;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){

      double dtau[3][3],codtau[3][3],xpgref[3],wpg;
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      NULL,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0]*codtau[0][0]
	+ dtau[0][1]*codtau[0][1]
	+ dtau[0][2]*codtau[0][2];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	f->dtwn[imem]/=(wpg*det);
      }
    }
  }

  return NULL;


}

// compute the Discontinuous Galerkin volume terms
// fast version
void* DGVolume(void* mc){

  MacroCell* mcell = (MacroCell*) mc;

  Field* f= mcell->field;

  // loop on the elements
  for (int ie=mcell->first_cell;ie<mcell->last_cell_p1;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    const int m = f->model.m;
    const int deg[3]={f->interp_param[1],
		      f->interp_param[2],
		      f->interp_param[3]};
    const int npg[3] = {deg[0]+1,
			deg[1]+1,
			deg[2]+1};
    const int nraf[3]={f->interp_param[4],
		       f->interp_param[5],
		       f->interp_param[6]};

    const unsigned int sc_npg=npg[0]*npg[1]*npg[2];

    int f_interp_param[8]= {f->interp_param[0],
			    f->interp_param[1],
			    f->interp_param[2],
			    f->interp_param[3],
			    f->interp_param[4],
			    f->interp_param[5],
			    f->interp_param[6],
			    f->interp_param[7]};

    // loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++){
      for(int icL1 = 0; icL1 < nraf[1]; icL1++){
	for(int icL2 = 0; icL2 < nraf[2]; icL2++){

	  int icL[3] = {icL0,icL1,icL2};
	  // get the L subcell id
	  int ncL=icL[0]+nraf[0]*(icL[1]+nraf[1]*icL[2]);
	  // first glop index in the subcell
	  int offsetL=npg[0]*npg[1]*npg[2]*ncL;

	  // compute all of the xref for the subcell
	  double *xref0 = malloc(sc_npg * sizeof(double));
	  double *xref1 = malloc(sc_npg * sizeof(double));
	  double *xref2 = malloc(sc_npg * sizeof(double));
	  double *omega = malloc(sc_npg * sizeof(double));
	  int *imems = malloc(m * sc_npg * sizeof(double*));
	  int pos=0;
	  for(unsigned int p=0; p < sc_npg; ++p) {
	    double xref[3];
	    double tomega;

	    ref_pg_vol(f->interp_param+1,offsetL+p,xref,&tomega,NULL);
	    xref0[p] = xref[0];
	    xref1[p] = xref[1];
	    xref2[p] = xref[2];
	    omega[p] = tomega;
	    
	    for(int im=0; im < m; ++im) {
	      imems[pos++] = f->varindex(f_interp_param,ie,offsetL+p,im);
	    }
	  }

	  // loop in the "cross" in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++){
	    // point p at which we compute the flux
    
	    for(int p0 = 0; p0 < npg[0]; p0++){
	      for(int p1 = 0; p1 < npg[1]; p1++){
		for(int p2 = 0; p2 < npg[2]; p2++){
		  double wL[m],flux[m];
		  int p[3]={p0,p1,p2};
		  int ipgL=offsetL+p[0]+npg[0]*(p[1]+npg[1]*p[2]);
		  for(int iv=0; iv < m; iv++){
		    ///int imemL=f->varindex(f_interp_param,ie,ipgL,iv);
		    wL[iv] = f->wn[imems[ipgL-offsetL+iv]];

		    //wL[iv] = f->wn[imemL];
		  }
		  int q[3]={p[0],p[1],p[2]};
		  // loop on the direction dim0 on the "cross"
		  for(int iq = 0; iq < npg[dim0]; iq++){
		    q[dim0]=(p[dim0]+iq)%npg[dim0];
		    double dphiref[3]={0,0,0};
		    // compute grad phi_q at glop p
		    dphiref[dim0]=dlag(deg[dim0],q[dim0],p[dim0])*nraf[dim0];

		    double xrefL[3]={xref0[ipgL-offsetL],
		    		     xref1[ipgL-offsetL],
		    		     xref2[ipgL-offsetL]};
		    double wpgL=omega[ipgL-offsetL];
		    /* double xrefL[3], wpgL; */
		    /* ref_pg_vol(f->interp_param+1,ipgL,xrefL,&wpgL,NULL); */

		    // mapping from the ref glop to the physical glop
		    double dtau[3][3],codtau[3][3],dphiL[3];
		    Ref2Phy(physnode,
			    xrefL,
			    dphiref,  // dphiref
			    -1,    // ifa                                 
			    NULL,  // xphy  
			    dtau,
			    codtau,
			    dphiL,  // dphi
			    NULL);  // vnds   
    
		    f->model.NumFlux(wL,wL,dphiL,flux);

		    int ipgR=offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		    for(int iv=0; iv < m; iv++){
		      //int imemR=f->varindex(f_interp_param,ie,ipgR,iv);
		      f->dtwn[imems[ipgR-offsetL+iv]]+=flux[iv]*wpgL;
		    }
		  } // iq
		} // p2
	      } // p1
	    } // p0
	    
	  } // dim loop
	  
	  free(omega);
	  free(xref0);
	  free(xref1);
	  free(xref2);
	  free(imems);

	} // icl2
      } //icl1
    } // icl0
  }

  return NULL;


}

// compute the Discontinuous Galerkin volume terms
// slow version
void DGVolumeSlow(Field* f){

  // assembly of the volume terms
  // loop on the elements
  //#pragma omp parallel for
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
	  double det
	    = dtau[0][0] * codtau[0][0]
	    + dtau[0][1] * codtau[0][1]
	    + dtau[0][2] * codtau[0][2];
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


  MacroCell mcell[f->macromesh.nbelems];

  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    mcell[ie].field=f;
    mcell[ie].first_cell=ie;
    mcell[ie].last_cell_p1=ie+1;
  }

#ifdef _WITH_PTHREAD

  // we will have only one flying thread per 
  // macrocell
  pthread_t tmcell[f->macromesh.nbelems];
  int status;

  // launch a thread for each macro cell
  // computation of the inter subcell fluxes
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    status=pthread_create (&(tmcell[ie]),   // thread
			   NULL,                   // default attributes
			   DGMacroCellInterface,    // called function
			   (void*) (mcell+ie));  // function params
    assert(status==0);
    //DGMacroCellInterface((void*) (mcell+ie));
  }
  // wait the end of the threads before next step
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    pthread_join(tmcell[ie], NULL);
  }

  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    status=pthread_create (&(tmcell[ie]),   // thread
			   NULL,                   // default attributes
			   DGSubCellInterface,    // called function
			   (void*) (mcell+ie));  // function params
    assert(status==0);
    //DGSubCellInterface((void*) (mcell+ie));
  }
  // wait the end of the threads before next step
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    pthread_join(tmcell[ie], NULL);
  }

  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    status=pthread_create (&(tmcell[ie]),   // thread
			   NULL,                   // default attributes
			   DGVolume,    // called function
			   (void*) (mcell+ie));  // function params
    assert(status==0);
    //DGVolume((void*) (mcell+ie));
  }
  // wait the end of the threads before next step
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    pthread_join(tmcell[ie], NULL);
  }

  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    status=pthread_create (&(tmcell[ie]),   // thread
			   NULL,                   // default attributes
			   DGMass,    // called function
			   (void*) (mcell+ie));  // function params
    assert(status==0);
    //DGMass((void*) (mcell+ie));
  }
  // wait the end of the threads before next step
  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    pthread_join(tmcell[ie], NULL);
  }

#else


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
  for(int ie=0; ie < f->macromesh.nbelems; ++ie) {
    DGMacroCellInterface((void*) (mcell+ie));
    DGSubCellInterface((void*) (mcell+ie));
    DGVolume((void*) (mcell+ie));
    DGMass((void*) (mcell+ie));
  }
  
#endif


}

// apply the Discontinuous Galerkin approximation for computing
// the time derivative of the field
void dtFieldSlow(Field* f){

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
  //#pragma omp parallel for
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
  //#pragma omp parallel for
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
  //#pragma omp parallel for
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
	  double det
	    = dtau[0][0] * codtau[0][0]
	    + dtau[0][1] * codtau[0][1]
	    + dtau[0][2] * codtau[0][2];
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
  //#pragma omp parallel for
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
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
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

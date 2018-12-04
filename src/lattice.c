#include "lattice.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"

/**************************************************************************************/
void Lattice_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
 LatticeData * ld=&schnaps_lattice_data;
  for(int i = 0;i < ld->index_max_q+1;i++){
    schnaps_real vn=0;
    for(int dim = 0;dim < ld->d; dim++){
      vn += ld->q_tab[i][dim]*vnorm[dim];
    }
    schnaps_real vnp = vn>0 ? vn : 0;
    schnaps_real vnm = vn-vnp;
    flux[i] = vnp * wL[i] + vnm * wR[i];
    //printf("i :%i \t Flux :%e \n",i, flux[i]); 
  }
  flux[ld->index_rho]=0; 
  flux[ld->index_ux]=0; 
  flux[ld->index_uy]=0;
  flux[ld->index_uz]=0;
  flux[ld->index_temp]=0; 
  flux[ld->index_p]=0; 
}
/**************************************************************************************/
void Lattice_OneNodeNumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  LatticeData * ld=&schnaps_lattice_data;
  int inode=ld->current_node_index;
  schnaps_real vn=0;
  for(int dim = 0;dim < ld->d; dim++){
    vn += ld->q_tab[inode][dim]*vnorm[dim];
  }
  schnaps_real vnp = vn>0 ? vn : 0;
  schnaps_real vnm = vn-vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
}
/**************************************************************************************/
void Compute_distribution_eq(Simulation *simu, schnaps_real * w_eq) {
  LatticeData * ld=&schnaps_lattice_data;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      	int irho=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_rho);
        int iux=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_ux);
        int iuy=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uy);
        int itemp=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_temp);
        schnaps_real rho = f->wn[irho];
        schnaps_real temp = f->wn[itemp];
        schnaps_real  ux = f->wn[iux];
        schnaps_real uy = f->wn[iuy];
        schnaps_real uz= 0.0;
        schnaps_real p = 0.0;
        //schnaps_real uz = f->wn[iuz];
        //schnaps_real u2=  ux*ux + uy*uy;
        for(int iv=0;iv<ld->index_max_q+1;iv++){
          int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
          //schnaps_real uv = ux * ld->q_tab[iv][0] + uy * ld->q_tab[iv][1];
          //schnaps_real temp_val=ld->w_tab[iv]*rho*(1.0+uv+0.5*uv*uv-0.5*t*u2);	  
          schnaps_real temp_val=ld->feq(iv,ld, rho,ux,uy,uz,temp,p);
          w_eq[ikin]= temp_val;
        }
    } 
  }
}
/**************************************************************************************/
void Lattice_Dummy_InitData(schnaps_real x[3],schnaps_real w[]){
  w[0]=0.0;
}
/**************************************************************************************/
schnaps_real feq_isothermal_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p){
    LatticeData *ld= lattice;
    schnaps_real u2= (ux * ux + uy * uy)/temp;
    schnaps_real uv= (ux *ld->q_tab[i_node][0] + uy * ld->q_tab[i_node][1])/temp;
    schnaps_real feq= ld->w_tab[i_node] * rho * (1.0+uv+0.5 * (uv * uv- u2));
    return feq;
}
/**************************************************************************************/
schnaps_real feq_isothermal_linearwave_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p){
    LatticeData *ld= lattice;
    schnaps_real uv= (ux *ld->q_tab[i_node][0] + ld->q_tab[i_node][1]* uy)/temp;
    schnaps_real feq= ld->w_tab[i_node]* rho * (1.0+uv);
    return feq;
}
/**************************************************************************************/
void Compute_moments(Simulation *simu) {
  LatticeData * ld=&schnaps_lattice_data;
  //
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    //
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      	int irho=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_rho);
        int iux=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_ux);
        int iuy=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uy);
        int iuz=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uz);
        int it=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_temp);
      f->wn[irho]=0.0;
      f->wn[iux]=0.0;
      f->wn[iuy]=0.0;
      f->wn[iuz]=0.0;
      for(int iv=0;iv<ld->index_max_q+1;iv++){
        int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
        //
        f->wn[irho]+=f->wn[ikin];
        f->wn[iux]+=f->wn[ikin]*ld->q_tab[iv][0];
        f->wn[iuy]+=f->wn[ikin]*ld->q_tab[iv][1];
      };
	    f->wn[iux]=f->wn[iux]/f->wn[irho];
	    f->wn[iuy]=f->wn[iuy]/f->wn[irho];
    }
  }
}
/**************************************************************************************/

void Compute_relaxation(Simulation *simu, schnaps_real * w_eq) {
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real nu=0;
  nu=simu->dt/(ld->tau+0.5*simu->dt);
  //
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
  	for(int iv=0;iv<ld->index_max_q+1;iv++){
      int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      f->wn[ikin]=f->wn[ikin]-nu*(f->wn[ikin]-w_eq[ikin]);
    };
    };
  };
}
/*************************************************************************************/
void LatticeThetaTimeScheme(Simulation *simu,Model *model_advec, schnaps_real tmax, schnaps_real dt){
  MatrixStorage ms=SKYLINE;
  LatticeData *ld=&schnaps_lattice_data;
  int nb_nodes=ld->q; // velocity nodes on the lattice
  int ig_glob=0,ig_node=0;
  field *f_glob,*f_node;
  // 
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};  
  int nb_ipg = NPG(deg, nraf);
  //
  Simulation *simu_advec;
  simu_advec=malloc(sizeof(Simulation));
  EmptySimulation(simu_advec);
  InitSimulation(simu_advec, &(simu->macromesh), deg, nraf, model_advec);
  simu_advec->vmax=simu->vmax;
  simu_advec->cfl=0.0;
  simu_advec->nb_diags=0;
  simu_advec->pre_dtfields=NULL;
  simu_advec->post_dtfields=NULL;
  simu_advec->update_after_rk=NULL;
  //
  LinearSolver solver_imp[nb_nodes]; 
  LinearSolver solver_exp[nb_nodes];
  schnaps_real *res_advec=  calloc(simu_advec->wsize,sizeof(schnaps_real*)); 
  //
  schnaps_real theta=0.5; // implicit explicit splitting parameter 
  simu->dt=dt;
  int itermax=(int) (tmax/dt)+1;
  printf("Called with tmax=%f dt=%f Nb iterations:%i \n",tmax,dt,itermax);
  simu->itermax_rk=itermax;
  int freq = (1 >= simu->itermax_rk / 10)? 1 : simu->itermax_rk / 10;
  simu->tnow=0.0;
  simu_advec->dt=dt;
  simu_advec->itermax_rk=itermax;
  simu_advec->tnow=0.0;
  for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
    simu->fd[ie].tnow=simu->tnow;
    simu->fd[ie].tnow=simu_advec->tnow;
  }
  // Diagnostics (this should be elsewhere, some timetraces  module ?
  int size_diags = simu->nb_diags * itermax;
  if(simu->nb_diags != 0) {
     simu->Diagnostics = malloc(size_diags * sizeof(schnaps_real));
  };
  //
  time_t t_start,t_end; // time measurements for op factorization
  t_start=time(NULL);
  //  Solvers Init/Assembly
  printf("Sparse Linear Solvers init");
  for (int isim=0;isim < nb_nodes;isim++){
    ld->current_node_index=isim;
    //
    InitImplicitLinearSolver(simu_advec, &solver_imp[isim],ms);
    InitImplicitLinearSolver(simu_advec, &solver_exp[isim],ms);
    //
  }
  // End Operators init
  //
  // Time loop start
  for (int iter=0;iter < itermax;iter++){
    //
    if (iter%freq == 0){
      printf(" iter %i/%i t=%f\n",iter,itermax,simu->tnow);
    }
    //
    simu->iter_time_rk=iter;
    //
    if (iter==0){
      t_start=time(NULL);
    }
    //
    if (simu->pre_dtfields !=NULL){
      simu->pre_dtfields(simu);
    };
    // now loop on velocity nodes
    for(int isim=0;isim < nb_nodes;isim++){
      ld->current_node_index=isim;
      // dispatch main w to per node w's
      for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
        f_glob=simu->fd+ie;
        for(int ipg=0; ipg < nb_ipg;ipg++){
          f_node= simu_advec->fd+ie;
          ig_glob=f_glob->varindex(f_glob->deg,f_glob->raf,f_glob->model.m,ipg,isim);
          ig_node=f_node->varindex(f_node->deg,f_node->raf,1,ipg,0);
          f_node->wn[ig_node]= f_glob->wn[ig_glob];
        //
        }; // ipg end loop glops
      }; // ie end loop macroelements
      // end of data dispatch
      // Solvers assembly and factorization if necessary
      solver_imp[isim].rhs_is_assembly=false;
      solver_exp[isim].rhs_is_assembly=false;
      if (iter==0){
          solver_imp[isim].mat_is_assembly=false;
          solver_exp[isim].mat_is_assembly=false;
      }
      else
      {
          solver_imp[isim].mat_is_assembly=true;
          solver_exp[isim].mat_is_assembly=true;
      };
      //
      simu_advec->tnow=simu->tnow;
      for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
        simu_advec->fd[ie].tnow=simu_advec->tnow;
      }
      AssemblyImplicitLinearSolver(simu_advec, &solver_exp[isim],-(1.0-theta),simu_advec->dt);
      simu_advec->tnow=simu->tnow+simu->dt;
      for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
        simu_advec->fd[ie].tnow=simu_advec->tnow;
      }
      AssemblyImplicitLinearSolver(simu_advec, &solver_imp[isim],theta,simu_advec->dt);
      // compute residual
      MatVect(&solver_exp[isim], simu_advec->w, res_advec);
      //
      for(int i=0;i<solver_imp[isim].neq;i++){
        solver_imp[isim].rhs[i]=-solver_exp[isim].rhs[i]+solver_imp[isim].rhs[i]+res_advec[i];
      }
      //
      solver_imp[isim].solver_type=LU;
      Advanced_SolveLinearSolver(&solver_imp[isim],simu_advec);
      //
      for(int i=0;i<solver_imp[isim].neq;i++){
        simu_advec->w[i]=solver_imp[isim].sol[i];
      }
      // collect nodes ws to main w
      for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
        f_glob=simu->fd+ie;
        for(int ipg=0; ipg < nb_ipg;ipg++){
          f_node= simu_advec->fd+ie;
          ig_glob=f_glob->varindex(f_glob->deg,f_glob->raf,f_glob->model.m,ipg,isim);
          ig_node=f_node->varindex(f_node->deg,f_node->raf,1,ipg,0);
          f_glob->wn[ig_glob]= f_node->wn[ig_node]; 
        }; // ipg end loop glops
      }; // ie end loop macroelements
    }; // isim end loop on velocity node 
    // post advec ops
    if (simu->post_dtfields !=NULL){
      simu->post_dtfields(simu);
    }
    if (simu->update_after_rk !=NULL){
      simu->update_after_rk(simu,simu->w);
    }
    if (iter==0){
      t_end= time(NULL);
      printf("First step duration %ld\n", t_end-t_start);
    }
    simu->tnow += simu->dt;
    simu_advec->tnow=simu->tnow;
    for(int ie=0; ie < simu->macromesh.nbelems; ++ie){
      simu->fd[ie].tnow=simu->tnow;
      simu_advec->fd[ie].tnow=simu->tnow;
    }
    //
  }; // iter end of time loop 
  //
  if (res_advec !=NULL){
    free(res_advec);
  }
  if (simu_advec!= NULL){
    free(simu_advec);
  }
  //
}
/***************************************************************************************************************/
void PlotExtFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename) {
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
    simu->fd[0].raf[1],
    simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
    simu->fd[0].deg[1],
    simu->fd[0].deg[2]};
  const int npg[3] = {deg[0] + 1,
    deg[1] + 1,
    deg[2] + 1};
  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];
  //
  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(schnaps_real));
  //
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  //int npgv = NPG(deg, nraf);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        }
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

        schnaps_real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};
        schnaps_real testpsi = 0;
        ////////////////////////////////////////
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projection to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
            // loop on basis function indexed by glop
              int ix[3];
              int ib;
              for (ix[0]=0;ix[0]< npg[0]; ix[0]++){
                for (ix[1]=0;ix[1]< npg[1]; ix[1]++){
                  for (ix[2]=0;ix[2]< npg[2]; ix[2]++){
                    xyz_to_ipg(nraf,deg,jcL,ix,&ib);
                    schnaps_real psi;
                    psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
                    testpsi += psi;
                    value[nodecount] += psi * w_in[ib];
                  };  // end ix[2]
                };  // end ix[1]
              };  // end ix[0]
        }; //end loop neighbour subcell 2
        };//end loop neighbour subcell 1
        }; //end loop neighbour subcell 0
        assert(fabs(testpsi-1) < _SMALL);
      ///////////////////////////////////////
      // Compare with an exact solution
      nodecount++;
      fprintf(gmshfile, "%d %f %f %f\n", nodecount,
        Xplot[0], Xplot[1], Xplot[2]);
      } // end loop Hex64 fe nodes
    } // end loop subcell 2
    } // end loop subcell 1
    } // end loop subcell 0
  } // end loop macrocells

  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  fprintf(gmshfile, "%d\n",
  simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);
  int elm_type = 92;
  int num_tags = 0;
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
  	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    // Get the subcell id
    int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
    // Global subcell id
    int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
    fprintf(gmshfile, "%d ", numelem);
    fprintf(gmshfile, "%d ", elm_type);
    fprintf(gmshfile, "%d ", num_tags);
    for(int ii = 0; ii < 64; ii++) {
      int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
      fprintf(gmshfile, "%d ", numnoe);
    }
    fprintf(gmshfile, "\n");
    }
    }
    }
  } 
  fprintf(gmshfile, "$EndElements\n");
  // Now display data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field : anonymous\n");
  else
    fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);
  //
  for(int ino = 0; ino < nb_plotnodes; ino++) {
    fprintf(gmshfile, "%d %f\n", ino + 1, value[ino]);
  }
  fprintf(gmshfile, "\n$EndNodeData\n");
  //
  fclose(gmshfile);
  free(value);
}
void PlotExtFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename){
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  int one = 1;
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // loop on basis function indexed by glop
            int ix[3];
            int ib;
            for (ix[0]=0;ix[0]< npg[0]; ix[0]++){
              for (ix[1]=0;ix[1]< npg[1]; ix[1]++){
                for (ix[2]=0;ix[2]< npg[2]; ix[2]++){
                  xyz_to_ipg(nraf,deg,jcL,ix,&ib);
                  schnaps_real psi;
                  psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
                  testpsi += psi;
                  value[nodecount] += psi * w_in[ib];
                };  // end ix[2]
              };  // end ix[1]
            };  // end ix[0]
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        nodecount++;
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "field : anonymous");
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);

  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
  free(value);
} 
/////////////////////////////////////
void PlotExtScalarFieldBinMultitime(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename,int create_file, schnaps_real t,int istep){
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  if (create_file==1){
  gmshfile = fopen(filename, "w" );
  }
  else
  {
  gmshfile = fopen(filename, "a" );
  }
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  if (create_file==1){
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  }
  int one = 1;
  if (create_file==1){
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  }
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  if (create_file==1){
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  }
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // loop on basis function indexed by glop
            int ix[3];
            int ib;
            for (ix[0]=0;ix[0]< npg[0]; ix[0]++){
              for (ix[1]=0;ix[1]< npg[1]; ix[1]++){
                for (ix[2]=0;ix[2]< npg[2]; ix[2]++){
                  xyz_to_ipg(nraf,deg,jcL,ix,&ib);
                  schnaps_real psi;
                  psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
                  testpsi += psi;
                  value[nodecount] += psi * w_in[ib];
                };  // end ix[2]
              };  // end ix[1]
            };  // end ix[0]
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        nodecount++;
        if (create_file==1){
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
        }
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  if (create_file==1){
  fprintf(gmshfile, "$EndNodes\n");
  }
  // Elements
  if (create_file==1){
    fprintf(gmshfile, "$Elements\n");
    int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
    fprintf(gmshfile, "%d\n",nb_elements);
    int elm_type = 92; //
    int num_tags = 0;
    int num_elm_follow= nb_elements;
    int elem_header[3] = {elm_type, num_elm_follow, num_tags};
    fwrite(elem_header, sizeof(int), 3, gmshfile);
    for(int i = 0; i < simu->macromesh.nbelems; i++) {
        // Loop on the subcells
        int icL[3];
        for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
        for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
        for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
        // Get the subcell id
        int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
        // Global subcell id
        int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
        int elm_data_size=1+num_tags+64; 
        int elem_data[elm_data_size];
        elem_data[0]= numelem;
          for(int ii = 0; ii < 64; ii++) {
            int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
            elem_data[ii+1]= numnoe;
          };
        fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
        };
        };
        };
    };
    fprintf(gmshfile, "$\nEndElements\n");
  }
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  int number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field:anonymous \n");
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  int number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  int number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  int number_of_components=1;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  fclose(gmshfile);
  free(value);
} 
void PlotVecFieldsBinSparseMultitime(int typplot[3], int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep){

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  if (create_file==1){
  gmshfile = fopen(filename, "w" );
  //printf("Creating file\n");
  }
  else
  {
  //printf("Appending data\n");
  gmshfile = fopen(filename, "a" );
  }
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];
  bool is_valid_field[3];
  int maxfield_index=simu->fd->model.m;
  for (int idim=0;idim < 3;idim++){
    is_valid_field[idim]=((typplot[idim] < maxfield_index) && (typplot[idim] > -1));
    //printf("idim %i typplot %i isvalid %i\n",idim,typplot[idim],is_valid_field[idim]);
  }
  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  if (create_file==1){
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  }
  int one = 1;
  if (create_file==1){
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  }
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  if (create_file==1){
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  }
  //
  schnaps_real *value = calloc(3*nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // loop on basis function indexed by glop
            int ix[3];
            int ib;
            for (ix[0]=0;ix[0]< npg[0]; ix[0]++){
              for (ix[1]=0;ix[1]< npg[1]; ix[1]++){
                for (ix[2]=0;ix[2]< npg[2]; ix[2]++){
                  xyz_to_ipg(nraf,deg,jcL,ix,&ib);
                  schnaps_real psi;
                  psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
                  testpsi += psi;
                  for (int idim=0;idim < 3; idim++){
                    if (is_valid_field[idim]){
                    int vi = f->varindex(f->deg, f->raf, f->model.m, ib, typplot[idim]);
                    value[3*nodecount+idim] += psi * f->wn[vi];
                    }
                    else{
                    value[3*nodecount+idim] += 0.0;
                    }
                  } // end idim
                  }; // end ix[2]
                }; // end ix[1]
              }; // end ix[0]
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        // Compare with an exact solution
        if (compare) {
          schnaps_real wex[f->model.m];
          f->model.ImposedData(Xphy, f->tnow, wex);
          for (int idim=0;idim<3;idim++){
          value[3*nodecount+idim] -= wex[typplot[idim]];
          }
        }
        nodecount++;
        if (create_file==1){
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
        }
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  if (create_file==1){
  fprintf(gmshfile, "$EndNodes\n");
  }
  // Elements
  if (create_file==1){
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  }
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  int number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL){
    int idim=0;
    fprintf(gmshfile, "\"field %d\"\n", typplot[idim]);
  }
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  int number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  int number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  int number_of_components=3;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    for (int idim=0; idim < 3;idim++){
    fwrite(&value[3*(ino-1)+idim],sizeof(double),1,gmshfile);
    }
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  fclose(gmshfile);
  free(value);
} 
//
void PlotFieldsBinSparseMultitime(int typplot, int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep){

  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  if (create_file==1){
  gmshfile = fopen(filename, "w" );
  //printf("Creating file\n");
  }
  else
  {
  //printf("Appending data\n");
  gmshfile = fopen(filename, "a" );
  }
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  if (create_file==1){
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  }
  int one = 1;
  if (create_file==1){
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  }
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  if (create_file==1){
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  }
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // loop on basis function indexed by glop
            int ix[3];
            int ib;
            for (ix[0]=0;ix[0]< npg[0]; ix[0]++){
              for (ix[1]=0;ix[1]< npg[1]; ix[1]++){
                for (ix[2]=0;ix[2]< npg[2]; ix[2]++){
                  xyz_to_ipg(nraf,deg,jcL,ix,&ib);
                  schnaps_real psi;
                  psi_ref_subcell(f->deg, f->raf, icL, ib, Xr, &psi, NULL);
                  testpsi += psi;
                  int vi = f->varindex(f->deg, f->raf, f->model.m, ib, typplot);
                  value[nodecount] += psi * f->wn[vi];
                };  // end ix[2]
              };  // end ix[1]
            };  // end ix[0]
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        // Compare with an exact solution
        if (compare) {
          schnaps_real wex[f->model.m];
          f->model.ImposedData(Xphy, f->tnow, wex);
          value[nodecount] -= wex[typplot];
        }
        nodecount++;
        if (create_file==1){
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
        }
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  if (create_file==1){
  fprintf(gmshfile, "$EndNodes\n");
  }
  // Elements
  if (create_file==1){
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  }
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  int number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", typplot);
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  int number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  int number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  int number_of_components=1;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  fclose(gmshfile);
  free(value);
} 
//**************************************************************************************************************************//
// lattice diagnostics (1D time traces) 
void Store_Lattice_diags(Simulation *simu){
  LatticeData *ld=&schnaps_lattice_data;
  schnaps_real tnow=simu->tnow;
  int iter= simu->iter_time_rk;
  int maxiter=simu->itermax_rk;
  int nb_diags=simu->nb_diags;
  //
  schnaps_real diag_vals[nb_diags];
  for (int irank=0;irank < nb_diags;irank++){
    diag_vals[irank]=0.0;
  };
  // actual collection
  ld->collect_diags(simu,diag_vals);
  //
  for (int irank=0;irank < nb_diags;irank++){
    simu->Diagnostics[iter* nb_diags+irank]= diag_vals[irank];
  };
}
void Dump_Lattice_Diagnostics(Simulation *simu,char simtag[3]){
  FILE *diagfile;
  assert((simu->Diagnostics != NULL));
  int nb_diags=simu->nb_diags;
  schnaps_real dt=simu->dt;
  char filename[sizeof("lbm_diag_TAG_raf000_cfl0.000.dat")];
  int raf=simu->fd[0].raf[0];
  schnaps_real cfl=simu->cfl;
  sprintf(filename,"lbm_diag_%s_raf%03d_cfl%1.3f.dat",simtag,raf,cfl);
  diagfile = fopen(filename,"w");
  for (int it=0; it < simu->itermax_rk;it++){
    schnaps_real tnow = dt * (schnaps_real) it;
    fprintf(diagfile,"%f\t",tnow);
    for (int irank=0;irank< simu->nb_diags;irank++){
      fprintf(diagfile,"%f\t",simu->Diagnostics[it*nb_diags+irank]);
    }
    fprintf(diagfile,"\n");
  };
  fclose(diagfile);
}
//
/**************************************************************************************/
void Compute_and_dump_Vorticity_2d(Simulation *simu,char *filename, int create_file, schnaps_real t, int istep){
    LatticeData *ld=&schnaps_lattice_data;
    field *f = simu->fd;
    int ifield_ux=ld->index_ux;
    int ifield_uy=ld->index_uy;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *temp_gradient;
    schnaps_real *vort;
    vort= calloc(nb_dof,sizeof(schnaps_real)); 
    temp_gradient= calloc(wdsize,sizeof(schnaps_real));
    //
    Compute_derivative(simu,temp_gradient,ifield_uy);
    //
    for (int i=0; i< nb_dof; i++){
      vort[i]= temp_gradient[i];
    };
    Compute_derivative(simu,temp_gradient,ifield_ux);
    schnaps_real max_vort=-10000.0;
    for (int i=0; i< nb_dof; i++){
      vort[i]= vort[i]-temp_gradient[i+nb_dof];
      if (vort[i] > max_vort){
        max_vort=vort[i];
      }
    };
    printf("Max vort %f\n",max_vort); 
    //
    if (temp_gradient != NULL){
      free(temp_gradient);
    }
    //
    //PlotGenericFieldAsciiSparse(simu,vort,"vorticity",filename);
    //PlotGenericFieldBinSparse(simu,vort,"vorticity",filename);
    PlotExtScalarFieldBinMultitime(simu, vort,"vorticity",filename,create_file,t,istep);
    int typplot[3]= {ld->index_ux,ld->index_uy,ld->index_uz};
    PlotVecFieldsBinSparseMultitime(typplot,false,simu,"u",filename,0,t,istep);
    //
    if (vort != NULL){
      free(vort);
    }
}
//

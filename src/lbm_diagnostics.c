#include "lbm_diagnostics.h"
/***************************************************************************************************************/
void LBM_PlotExtFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename) {
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
void LBM_PlotExtFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename){
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
void LBM_PlotExtScalarFieldBinMultitime(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename,int create_file, schnaps_real t,int istep){
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
        }; // end jcL[2]
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
void LBM_PlotVecFieldsBinSparseMultitime(int typplot[3], int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep){

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
  // check existence of fields, out of bounds field values lead to zero padding (useful for 1D/2D)
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
            }; // end jcL[2]
          }; // end jcL[1]
        };// end jcL[0]
        assert(fabs(testpsi-1) < _SMALL);
        // Compare with an exact solution
        if (compare) {
          schnaps_real wex[f->model.m];
          f->model.ImposedData(Xphy, f->tnow, wex);
          for (int idim=0;idim<3;idim++){
          if (is_valid_field[idim]){
          value[3*nodecount+idim] -= wex[typplot[idim]];
          }
          else{
          value[3*nodecount+idim] -= 0.0;
          }
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
  // add the norm as scalar field
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  number_of_string_tags=1;
  fprintf(gmshfile, "%i\n",number_of_string_tags);
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field %d\"\n", *typplot);
  else
  fprintf(gmshfile, "\"field: %s\"\n", "norm");
  //
  number_of_real_tags=1;
  fprintf(gmshfile,"%i\n",number_of_real_tags);
  fprintf(gmshfile,"%f\n",t);
  number_of_int_tags=3;
  fprintf(gmshfile,"%i\n",number_of_int_tags);
  fprintf(gmshfile,"%i\n",istep);
  number_of_components=1;
  fprintf(gmshfile,"%i\n",number_of_components);
  fprintf(gmshfile, "%i\n", nb_plotnodes);
  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    //
    schnaps_real tmpval = 0.0;
    for (int idim=0;idim<3;idim++){
      tmpval += value[3*(ino-1)+idim] * value[3*(ino-1)+idim];
    }
    tmpval = sqrt(tmpval);
    //
    fwrite(&tmpval,sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");
  //
  fclose(gmshfile);
  free(value);
} 
//
void LBM_PlotFieldsBinSparseMultitime(int typplot, int compare, Simulation* simu, char *fieldname,char *filename, int create_file, schnaps_real t, int istep){

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
        int jcL[3]={0,0,0};
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projection to neighbouring subcells (!! CAUTION )
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
            }; //end jcL[2]
          }; //end jcL[1]
        }; //end jcL[0]
        //
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
      } // end for ino
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
void LBM_Store_Lattice_diags(LBMSimulation *lbsimu){
  schnaps_real tnow=lbsimu->tnow;
  int iter= lbsimu->iter_time;
  int maxiter=lbsimu->itermax;
  int nb_diags_macro=lbsimu->macro_simu.nb_diags;
  int nb_diags_micro=lbsimu->micro_simu.nb_diags;
  //
  schnaps_real macro_diag_vals[nb_diags_macro];
  schnaps_real micro_diag_vals[nb_diags_micro];
  for (int irank=0;irank < nb_diags_macro;irank++){
    macro_diag_vals[irank]=0.0;
  };
  for (int irank=0;irank < nb_diags_micro;irank++){
    micro_diag_vals[irank]=0.0;
  };
  // actual collection
  if (lbsimu->collect_diags){
  lbsimu->collect_diags(lbsimu,macro_diag_vals,micro_diag_vals);
  }
  //
  for (int irank=0;irank < nb_diags_macro;irank++){
    lbsimu->macro_simu.Diagnostics[iter* nb_diags_macro+irank]= macro_diag_vals[irank];
  };
  for (int irank=0;irank < nb_diags_micro;irank++){
    lbsimu->micro_simu.Diagnostics[iter* nb_diags_micro+irank]= micro_diag_vals[irank];
  };
}
void LBM_Dump_Lattice_Diagnostics(LBMSimulation *lbsimu,char simtag[3]){
  FILE *diagfile_macro;
  FILE *diagfile_micro;
  int nb_diags_macro=lbsimu->macro_simu.nb_diags;
  int nb_diags_micro=lbsimu->micro_simu.nb_diags;
  schnaps_real dt=lbsimu->micro_simu.dt;
  int itermax=lbsimu->micro_simu.itermax_rk;
  printf("micro nb iter %i\n",itermax);
  char filename_macro[sizeof("lbm_diag_macro_TAG_raf000_cfl0.000.dat")];
  char filename_micro[sizeof("lbm_diag_micro_TAG_raf000_cfl0.000.dat")];
  int raf=lbsimu->macro_simu.fd[0].raf[0];
  schnaps_real cfl=lbsimu->micro_simu.cfl;
  sprintf(filename_macro,"lbm_diag_macro_%s_raf%03d_cfl%1.3f.dat",simtag,raf,cfl);
  sprintf(filename_micro,"lbm_diag_micro_%s_raf%03d_cfl%1.3f.dat",simtag,raf,cfl);
  if (nb_diags_macro > 0){
  diagfile_macro = fopen(filename_macro,"w");
  
  for (int it=0; it < lbsimu->micro_simu.itermax_rk;it++){
    schnaps_real tnow = dt * (schnaps_real) it;
    fprintf(diagfile_macro,"%f\t",tnow);
    for (int irank=0;irank< nb_diags_macro;irank++){
      fprintf(diagfile_macro,"%f\t",lbsimu->macro_simu.Diagnostics[it*nb_diags_macro+irank]);
    }
    fprintf(diagfile_macro,"\n");
  };
  fclose(diagfile_macro);
  }
  //
  if (nb_diags_micro >0){
  diagfile_micro = fopen(filename_micro,"w");
  for (int it=0; it < lbsimu->micro_simu.itermax_rk;it++){
    schnaps_real tnow = dt * (schnaps_real) it;
    fprintf(diagfile_micro,"%f\t",tnow);
    for (int irank=0;irank< nb_diags_micro;irank++){
      fprintf(diagfile_micro,"%f\t",lbsimu->micro_simu.Diagnostics[it*nb_diags_micro+irank]);
    }
    fprintf(diagfile_micro,"\n");
  };
  fclose(diagfile_micro);
  }
}
//
/**************************************************************************************/
/*void LBM_Compute_and_dump_Vorticity_2d(LBMSimulation *lbsimu,char *filename, int create_file, schnaps_real t, int istep){*/
/*    LatticeData *ld=&schnaps_lattice_data;*/
/*    field *f = simu->fd;*/
/*    int ifield_ux=ld->index_ux;*/
/*    int ifield_uy=ld->index_uy;*/
/*    int nb_dof=simu->wsize/f->model.m;*/
/*    int wdsize= 3 * nb_dof;*/
/*    schnaps_real *temp_gradient;*/
/*    schnaps_real *vort;*/
/*    vort= calloc(nb_dof,sizeof(schnaps_real)); */
/*    temp_gradient= calloc(wdsize,sizeof(schnaps_real));*/
/*    //*/
/*    Compute_derivative(simu,temp_gradient,ifield_uy);*/
/*    //*/
/*    for (int i=0; i< nb_dof; i++){*/
/*      vort[i]= temp_gradient[i];*/
/*    };*/
/*    Compute_derivative(simu,temp_gradient,ifield_ux);*/
/*    schnaps_real max_vort=-10000.0;*/
/*    for (int i=0; i< nb_dof; i++){*/
/*      vort[i]= vort[i]-temp_gradient[i+nb_dof];*/
/*      if (vort[i] > max_vort){*/
/*        max_vort=vort[i];*/
/*      }*/
/*    };*/
/*    printf("Max vort %f\n",max_vort); */
/*    //*/
/*    if (temp_gradient != NULL){*/
/*      free(temp_gradient);*/
/*    }*/
/*    //*/
/*    //PlotGenericFieldAsciiSparse(simu,vort,"vorticity",filename);*/
/*    //PlotGenericFieldBinSparse(simu,vort,"vorticity",filename);*/
/*    LBM_PlotExtScalarFieldBinMultitime(simu, vort,"vorticity",filename,create_file,t,istep);*/
/*    int typplot[3]= {ld->index_ux,ld->index_uy,ld->index_uz};*/
/*    LBM_PlotVecFieldsBinSparseMultitime(typplot,false,simu,"u",filename,0,t,istep);*/
/*    //*/
/*    if (vort != NULL){*/
/*      free(vort);*/
/*    }*/
/*}*/
//

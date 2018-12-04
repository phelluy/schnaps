#include "io.h"
#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>

#define _DEBUG_XDMF_STRUCT_NEW

void schnaps_simulation_xdmf_plot_fields_xml_structured_glops(Simulation *simu, char *filename){

  //
  bool plot_2D = simu->macromesh.is2d;
	bool is_file_writer=true;

    FILE * xdmffile;
    xdmffile = fopen(filename, "w" );
    fprintf(xdmffile,"<?xml version=\"1.0\" ?>\n");
    fprintf(xdmffile,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xdmffile,"<Xdmf Version=\"2.0\" xmlns:xi=\"[http://www.w3.org/2001/XInclude]\">\n");
    fprintf(xdmffile,"<Domain>\n");
    fprintf(xdmffile,"<Grid Name=\"TheMesh\" GridType=\"Collection\">\n");
   fprintf(xdmffile,"<Time Value=\"%f\"/>\n",simu->tnow);

  for(int i = 0; i < simu->macromesh.nbelems; i++){
  	field *f = NULL;
 		f = simu->fd + i;
   
    plot_2D = plot_2D && (0 == f->deg[2]);
    plot_2D = plot_2D && (1 == f->raf[2]);

		int nraf[3];
		int deg[3];
		for (int k=0; k < 3 ; k++){
			nraf[k] = f->raf[k];
			deg[k] = f->deg[k];
		}
	  int nb_points = (nraf[0] * deg[0] +1) * (nraf[1] * deg[1] +1) * (nraf[2] *deg[2] +1) ;


      fprintf(xdmffile,"<Grid Name=\"TestMesh%d\">\n",i);

      if (plot_2D){
      fprintf(xdmffile,"<Topology TopologyType = \"2DSMesh\"  Dimensions =\"%d %d\">\n",
        (nraf[0] * deg[0]+1),
        (nraf[1] * deg[1]+1));

      }else {
      fprintf(xdmffile,"<Topology TopologyType = \"3DSMesh\"  Dimensions =\"%d %d %d\">\n",
        (nraf[0] * deg[0]+1),
        (nraf[1] * deg[1]+1),
        (nraf[2] * deg[2]+1));
      }
      fprintf(xdmffile,"</Topology>\n");
      if (plot_2D){
        fprintf(xdmffile,"<Geometry Type=\"XY\" >\n");
        fprintf(xdmffile,"<DataItem Name =\"xyvalues\" Format=\"xml\" NumberType=\"Float\" Dimensions=\"%d 2\">\n",nb_points);
      } else{
        fprintf(xdmffile,"<Geometry Type=\"XYZ\" >\n");
        fprintf(xdmffile,"<DataItem Name =\"xyzvalues\" Format=\"xml\" NumberType=\"Float\" Dimensions=\"%d 3\">\n",nb_points);
      }

    int ic[3];
    int ix[3];
    int ixmax[3];
    int icount=0;

    int value_size = f->model.m * nb_points;
    schnaps_real *value = malloc(f->model.m * nb_points * sizeof(schnaps_real));

    for (ic[0] = 0; ic[0] < nraf[0]; ic[0]++){
        if (ic[0] < nraf[0]-1){ ixmax[0] = deg[0];} 
          else{ixmax[0] = deg[0]+1;}
    for (ix[0] = 0;ix[0] < ixmax[0]; ix[0]++){
    for (ic[1] = 0; ic[1] < nraf[1]; ic[1]++){
        if (ic[1] < nraf[1]-1){ ixmax[1] = deg[1];} 
          else{ixmax[1] = deg[1]+1;}
    for (ix[1] = 0;ix[1] < ixmax[1]; ix[1]++){
    for (ic[2] = 0; ic[2] < nraf[2]; ic[2]++){
        if (ic[2] < nraf[2]-1){ ixmax[2] = deg[2];} 
          else{ixmax[2] = deg[2]+1;}
    for (ix[2] = 0;ix[2] < ixmax[2]; ix[2]++){
        
        schnaps_real xpgref[3],xpg[3];
        int ipg=0;
        xyz_to_ipg(f->raf,f->deg,ic,ix,&ipg);
        ref_pg_vol(f->deg,f->raf,ipg,xpgref,NULL,NULL);
        schnaps_ref2phy(simu->fd[i].physnode, xpgref, NULL, -1, xpg, NULL,  NULL, NULL, NULL);
        
        double Xplot[3] ={ (double)xpg[0],  (double)xpg[1],  (double)xpg[2]};
          if (plot_2D){
          fprintf(xdmffile,"%f %f\n",Xplot[0],Xplot[1]);
          } else{
          fprintf(xdmffile,"%f %f %f\n",Xplot[0],Xplot[1],Xplot[2]);
          }
        
          for (int iv = 0; iv < f->model.m; iv++){
            int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
            value[iv * nb_points + icount] = f->wn[imem];
          }
        
        icount++;
      } // end ix[2]
      } // end ic[2]
      } // end ix[1]
    } // end ic[1]
    } // end ix[0]
    } // end ic[0]

  
	  fprintf(xdmffile,"</DataItem>\n"); 
    fprintf(xdmffile,"</Geometry>\n");

 
		  for (int iv = 0; iv < f->model.m; iv++){

		  fprintf(xdmffile, "<Attribute Name=\"m%d\" center=\"node\" >\n",iv);

      if (plot_2D){
        fprintf(xdmffile,"<DataItem Format=\"xml\" NumberType=\"Float\" Dimensions=\"%d %d\">\n",
        (nraf[0] * deg[0]+1),
        (nraf[1] * deg[1]+1));
      } else{
        fprintf(xdmffile,"<DataItem Format=\"xml\" NumberType=\"Float\" Dimensions=\"%d %d %d\">\n",
        (nraf[0] * deg[0]+1),
        (nraf[1] * deg[1]+1),
        (nraf[2] * deg[2]+1));
      }
		  int icount= 0;
      for (ic[0] = 0; ic[0] < nraf[0]; ic[0]++){
          if (ic[0] < nraf[0]-1){ ixmax[0] = deg[0];} 
            else{ixmax[0] = deg[0]+1;}
      for (ix[0] = 0;ix[0] < ixmax[0]; ix[0]++){
      for (ic[1] = 0; ic[1] < nraf[1]; ic[1]++){
          if (ic[1] < nraf[1]-1){ ixmax[1] = deg[1];} 
            else{ixmax[1] = deg[1]+1;}
      for (ix[1] = 0;ix[1] < ixmax[1]; ix[1]++){
      for (ic[2] = 0; ic[2] < nraf[2]; ic[2]++){
          if (ic[2] < nraf[2]-1){ ixmax[2] = deg[2];} 
            else{ixmax[2] = deg[2]+1;}
      for (ix[2] = 0;ix[2] < ixmax[2]; ix[2]++){
        int ipg=0;
        xyz_to_ipg(f->raf,f->deg,ic,ix,&ipg);
			  int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);

        fprintf(xdmffile,"%f\n",value[iv * nb_points + icount]);

			  icount++;

        } // end ix[2]
        } // end ic[2]
        } // end ix[1]
      } // end ic[1]
      } // end ix[0]
      } // end ic[0]
	
      fprintf(xdmffile,"</DataItem>\n");
		  fprintf(xdmffile, "</Attribute>\n"); 
		  //
		  } // end iv

      //
	    fprintf(xdmffile,"</Grid>\n");

      free(value);
      
  } // end i (macroelement)
  
  fprintf(xdmffile,"</Grid>\n");
  fprintf(xdmffile,"</Domain>\n");
  
  
  fprintf(xdmffile,"</Xdmf>\n");
  fclose(xdmffile);	
  
  
}



#ifdef _WITH_HDF5
void schnaps_simulation_dump_field_state(Simulation *simu,char* filename_pref){
	// hdf5 specific data
  hid_t       file_id =0;
  hid_t       dset_id =0;
  hid_t       dspace_id =0;
  hid_t       adspace_id = 0;
  hid_t       attr_id  = 0;
  hsize_t     dims_data[1]  = {0};
  herr_t      status  = 0;
  //

  char filename[strlen(filename_pref) + sizeof("_000_000000.hdf5")];
  sprintf(filename,"%s_%06d.hdf5",filename_pref,simu->iter_time_rk);
  file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  //
  hsize_t dims1[1] = {1};
  schnaps_real rdata[1] = {0};
  int idata[1] = {0};
  dspace_id = H5Screate_simple(1, dims1, NULL);
  // tnow
  rdata[0] = simu->tnow;
#ifdef _DOUBLE_PRECISION
  dset_id = H5Dcreate2(file_id,"tnow", H5T_NATIVE_DOUBLE, dspace_id, 
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,rdata);
#else
  dset_id = H5Dcreate2(file_id,"tnow", H5T_NATIVE_FLOAT, dspace_id, 
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,rdata);
#endif
  status = H5Dclose(dset_id);
  //
  // iter
  idata[0]= simu->iter_time_rk;
  dset_id = H5Dcreate2(file_id,"iter", H5T_NATIVE_INT, dspace_id, 
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,idata);
  status = H5Dclose(dset_id);
  //
  status = H5Sclose(dspace_id);
    //

  for (int ie = 0; ie < simu->macromesh.nbelems; ie++){

    field *f = simu->fd + ie;
    dims_data[0] = f->wsize;

    schnaps_real* buff = NULL;
    buff = f->wn;
    assert(buff);

    // add some attributes (light metadata) to dataset
    hsize_t adims[1];
    int attr_data1[1] = {0};

    //
    // Micro fields
    //
    char dset_name[sizeof("data_000000")];
    sprintf(dset_name,"data_%06d",ie);
    //
    dspace_id = H5Screate_simple(1, dims_data, NULL);
    //
#ifdef _DOUBLE_PRECISION
    dset_id = H5Dcreate2(file_id,dset_name, H5T_NATIVE_LDOUBLE, dspace_id, 
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,buff);
#else
    dset_id = H5Dcreate2(file_id,dset_name, H5T_NATIVE_FLOAT, dspace_id, 
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buff);
#endif

    // wsize
    adims[0] = 1;
    adspace_id = H5Screate_simple(1, adims, NULL);
    attr_id = H5Acreate2 (dset_id, "wsize", H5T_NATIVE_INT, adspace_id, 
                              H5P_DEFAULT, H5P_DEFAULT);
    attr_data1[0] = f->wsize;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, attr_data1);
    status = H5Aclose(attr_id);
    // m
    attr_id = H5Acreate2 (dset_id, "m", H5T_NATIVE_INT, adspace_id, 
                              H5P_DEFAULT, H5P_DEFAULT);
    attr_data1[0] = f->model.m;
    status = H5Awrite(attr_id, H5T_NATIVE_INT, attr_data1);
    status = H5Aclose(attr_id);
    // 
    status = H5Sclose(adspace_id);
    // new adspace_id for deg and raf
    adims[0] = 3;
    adspace_id = H5Screate_simple(1, adims, NULL);
    // deg
    attr_id = H5Acreate2 (dset_id, "deg", H5T_NATIVE_INT, adspace_id, 
                              H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_INT, f->deg);
    status = H5Aclose(attr_id);
    // raf
    attr_id = H5Acreate2 (dset_id, "raf", H5T_NATIVE_INT, adspace_id, 
                              H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attr_id, H5T_NATIVE_INT, f->raf);
    status = H5Aclose(attr_id);
    //
    status = H5Sclose(adspace_id);
    // end of macro dataset attributes writing
    status = H5Dclose(dset_id);
    status = H5Sclose(dspace_id);

  } // end ie
  status = H5Fclose(file_id);
}


void schnaps_simulation_load_field_state(Simulation *simu,char* filename_pref){

  printf("Loading simulation field data at t=%f iter=%d\n",simu->tnow,simu->iter_time_rk);
  //
	// hdf5 specific data
  hid_t       file_id;
  hid_t       dset_id;
  hsize_t     dims_mac_data[1];
  hsize_t     dims_mic_data[1];
  herr_t      status;
  //

  char filename[strlen(filename_pref) + sizeof("_000_000000.hdf5")];
  sprintf(filename,"%s_%06d.hdf5",filename_pref,simu->iter_time_rk);
  printf("simu : \n",simu->iter_time_rk);
  file_id = H5Fopen(filename,H5F_ACC_RDONLY,H5P_DEFAULT);
  assert(file_id > -1); // abort if file does not ex
  // read global attributes
  schnaps_real tnow;
  dset_id = H5Dopen(file_id,"tnow", H5P_DEFAULT);
#ifdef _DOUBLE_PRECISION
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&tnow);
#else
  status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&tnow);
#endif
  printf("Simu tnow set by code to %f - loaded from file %f\n",simu->tnow,tnow);
  status = H5Dclose(dset_id);
  int niter;
  dset_id = H5Dopen(file_id,"iter", H5P_DEFAULT);
  status = H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&niter);
  printf("Simu iter set by code to %d - loaded from file %d\n",simu->iter_time_rk,niter);
  status = H5Dclose(dset_id);
    //

  for (int ie = 0; ie < simu->macromesh.nbelems; ie++){

    field *f = simu->fd + ie;
    schnaps_real* buff = NULL;
    //
    buff = f->wn;

    //
    char dset_name[sizeof("data_000000")];
    sprintf(dset_name,"data_%06d",ie);
    //
    dset_id = H5Dopen(file_id,dset_name, H5P_DEFAULT);
#ifdef _DOUBLE_PRECISION
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,buff);
#else
    status = H5Dread(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buff);
#endif
    status = H5Dclose(dset_id);
    } // end ie

    // close file
  status = H5Fclose(file_id);
};
#endif


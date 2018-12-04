#include "dg_implicit.h"
#include "lobatto.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// gcc fe_implicit.c skyline.c -O -lm



int main(void) {

  //TestPLU();

  double xmin = -2;
  double xmax = 2;

  // cfl computed from the element size dx
  double cfl_dx = 1;
  // the actual cfl is computed from the smallest
  // distance bewteen two Gauss-Lobatto points
  double cfl = cfl_dx * (glop(DEG, 1) - glop(DEG, 0));
  //double cfl = cfl_dx;

  printf("cfl=%f\n",cfl);

  //dcmplx c = (1. + I) / 2;
  dcmplx c =1./2;

  double tmax = 0.4;

  static galerkin gal;

  gal.t = 0;
  gal_construct(&gal, xmin, xmax, cfl, tmax, c);
  int iter=0;

  int methode = 0;
  
  while(gal.t < tmax){

    dcmplx tnow = gal.t;

    gal.t += 2* creal(gal.dt);
    // correction last time step
    if (gal.t > tmax) {
      double dt = 2* creal(gal.dt) - gal.t + tmax;
      gal.dt = dt * c;
      printf("last step\n");
      gal.t = tmax;
    }
   
    methode = 2;

    // méthode 1 S(1/2) C S(1/2)
   if (methode == 1) {
    
      iter++;
      printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,
	     creal(gal.dt),cimag(gal.dt));

      gal.dt = conj(gal.dt);
      tnow += gal.dt;
      gal_step_nonlin(&gal, tnow);


      //projection
      for(int ino = 0; ino < NB_NODES; ino++){
	dcmplx wloc[M];
	for(int iv = 0; iv < M; iv++) wloc[iv] =gal.wnm1[vindex(ino,iv)];
	bgk_relax(wloc, gal.dt/2);
	for(int iv = 0; iv < M; iv++) gal.wnm1[vindex(ino,iv)]=wloc[iv];
      }

      iter++;
      printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,
	     creal(gal.dt),cimag(gal.dt));
      gal.dt = conj(gal.dt);
      tnow += gal.dt;
      gal_step_nonlin(&gal, tnow);
    }
   // méthode 2: S(1/2)S(1/2) C(1)
    else if (methode == 2){
      iter++;
      printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,
	     creal(gal.dt),cimag(gal.dt));
      
      gal.dt = conj(gal.dt);
      tnow += gal.dt;
      gal_step_nonlin(&gal, tnow);
      iter++;
      printf("iter=%d t=%f dt=(%f,%f)\n",iter,gal.t,
	     creal(gal.dt),cimag(gal.dt));
      
      gal.dt = conj(gal.dt);
      tnow += gal.dt;
      gal_step_nonlin(&gal, tnow);
      for(int ino = 0; ino < NB_NODES; ino++){
	dcmplx wloc[M];
	for(int iv = 0; iv < M; iv++) wloc[iv] =gal.wnm1[vindex(ino,iv)];
	bgk_relax(wloc, gal.dt*2);
	//bgk_project(wloc);
	for(int iv = 0; iv < M; iv++) gal.wnm1[vindex(ino,iv)]=wloc[iv];
      }

    } else {
      assert(1==4);
    }

  }

  //last projection
  if (methode == 2) {
    for(int ino = 0; ino < NB_NODES; ino++){
      dcmplx wloc[M];
      for(int iv = 0; iv < M; iv++) wloc[iv] =gal.wnm1[vindex(ino,iv)];
      bgk_project(wloc);
      for(int iv = 0; iv < M; iv++) gal.wnm1[vindex(ino,iv)]=wloc[iv];
      for(int iv = 0; iv < M; iv++) gal.wn[vindex(ino,iv)]=wloc[iv];
    }
  }

  
  gal_plot(&gal);

  for(int iv = 0; iv < M; iv++)
    printf("erreur L2 var %d=%e\n", iv,creal(gal_L2_error(&gal, iv)));

}


void source(double x, double t, dcmplx* w, dcmplx* source){

  //source[0] = 2*VMAX*x;
  source[0] = 0;
  if (M ==2) {
    //source[1] = -2*VMAX*x;
    source[1] = 0;
  }
  if (M == 3) {
    //source[1] = -2*VMAX*x;
    source[1] = 0;
    source[2] = 0;
  }

}


void bgk_project(dcmplx *w){

  dcmplx r = w[0] + w[1] + w[2];
  dcmplx q = w[2] - w[0];
  dcmplx u = q / r;

  dcmplx weq[3] = {q * (u - 1) / 2 + CSON * CSON * r / 2,
		   r * (1 - u * u - CSON * CSON),
		   q * (u + 1) / 2 + CSON * CSON * r / 2};

  w[0] = weq[0];
  w[1] = weq[1];
  w[2] = weq[2];

}


int vindex(int ipg, int iv)
{
  //return iv + ipg * M;
  return iv * NB_NODES + ipg;
}

int connec(int ie, int iloc)
{
  return (DEG + 1) * ie + iloc;
}


void gal_construct(galerkin *gal,
		   double xmin,
		   double xmax,
		   double cfl,
		   double tmax,
		   dcmplx smul)
{

  gal->xmin = xmin;
  gal->xmax = xmax;


  gal->cfl = cfl;
  gal->tmax = tmax;

  gal->smul = smul;

  gal->dx = (xmax - xmin) / NB_ELEMS;

  if (VMAX != 0) {
    gal->dt  = cfl * gal->dx / VMAX * smul;
  } else {
    gal->dt  = cfl * gal->dx * smul;
  }

  printf("dx=%f dt=(%f,%f)\n",gal->dx,creal(gal->dt),cimag(gal->dt));

  // nodes
  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int iloc = 0; iloc <= DEG; iloc++){
      double xi = glop(DEG, iloc);
      double x = gal->xmin + (ie + xi) * gal->dx;
      int ino = connec(ie, iloc);
      gal->xnode[ino] = x;
    }
  }

  // initial condition

  for(int ino = 0; ino < NB_NODES; ino++){
    double x = gal->xnode[ino];
    double t=0;
    dcmplx w[M];
    solexacte(x, t, w);
    for(int iv = 0; iv < M; iv++){
      int imem = vindex(ino,iv);
      gal->wn[imem] = w[iv];
      gal->wnm1[imem] = w[iv];
    }
  }

}

double gal_L2_error(galerkin *gal, int numvar)
{

  double gal_err = 0;
  double sav_err = 0;
  double vol = 0;
  
  FILE * ru_file = NULL;
  ru_file = fopen("fe_sav.bin", "rb" );
  dcmplx *rho_sav = NULL;
  dcmplx *u_sav = NULL;
  dcmplx *w_sav = NULL;
  double *x_sav = NULL;
  int deg_sav = 0;
  int nbelems_sav = 0;
  int nbnodes_sav = 0;

  if (ru_file){
    size_t ret;
    ret = fread(&nbelems_sav, sizeof(int), 1, ru_file);
    ret = fread(&deg_sav, sizeof(int), 1, ru_file);
    nbnodes_sav = (deg_sav + 1) * nbelems_sav;
    printf("deg=%d nbel=%d nbno=%d\n",deg_sav,nbelems_sav,nbnodes_sav),
      rho_sav = malloc(nbnodes_sav * sizeof(dcmplx));
    w_sav = malloc(nbnodes_sav * M * sizeof(dcmplx));
    u_sav = malloc(nbnodes_sav * sizeof(dcmplx));
    x_sav = malloc(nbnodes_sav * sizeof(double));
    for(int ino = 0; ino < nbnodes_sav; ino++){
      ret = fread(x_sav + ino, sizeof(double), 1, ru_file);
      for(int iv=0; iv < M; iv++) {
	ret = fread(w_sav + iv * nbnodes_sav + ino, sizeof(dcmplx), 1, ru_file);
      }
      ret = fread(rho_sav + ino, sizeof(dcmplx), 1, ru_file);
      ret = fread(u_sav + ino, sizeof(dcmplx), 1, ru_file);
      /* printf("ino=%d x=%f w=%f %f %f r=%f u=%f\n",ino, */
      /* 	     x_sav[ino], */
      /* 	     creal(w_sav[0 * nbnodes_sav + ino]), */
      /* 	     creal(w_sav[1 * nbnodes_sav + ino]), */
      /* 	     creal(w_sav[2 * nbnodes_sav + ino]), */
      /* 	     creal(rho_sav[ino]), */
      /* 	     creal(u_sav[ino])); */
    }
  } else {
    printf("no saved data for comparison\n");
  }
  
  for(int ie = 0; ie < NB_ELEMS; ie++){
    for(int kloc = 0; kloc <= DEG; kloc++){   // k: quadrature point
      int kno = connec(ie, kloc);
      double x = gal->xnode[kno];
      dcmplx w[M];
      dcmplx wex[M];
      solexacte(x, gal->t, wex);
      //printf("erreur temps=%f\n",gal->t);
      for(int iv = 0; iv < M; iv++){
	int kw = vindex(kno,iv);
	w[iv] =  gal->wn[kw];
      }

      dcmplx ws[M];
      if (ru_file) {
	interpolate(x_sav, w_sav,
		    deg_sav, nbnodes_sav,
		    x, ws);
      }
      
      vol += wglop(DEG, kloc) * gal->dx;
      
      double cerr = cabs(w[numvar]-wex[numvar]);
      double cerr_sav = 1e10;
      if (ru_file) cerr_sav = cabs(w[numvar]-ws[numvar]);
      //printf("%f %f\n",cabs(w[numvar]-wex[numvar]),vol);
      
      gal_err += cerr * cerr * wglop(DEG, kloc) * gal->dx;
      sav_err += cerr_sav * cerr_sav * wglop(DEG, kloc) * gal->dx;
    }
  }
  gal_err = sqrt(gal_err / vol);
  if (ru_file) {
    sav_err = sqrt(sav_err / vol);
    printf("err_sav=%e\n",sav_err);
  }

  if (ru_file){
    free(rho_sav);
    free(w_sav);
    free(u_sav);
    free(x_sav);
    fclose(ru_file);
  }


  
  return gal_err;
}

void lagrange(double *p, double *subdiv,
	      int deg, int ii, double x) {
  *p = 1;
  const int npg = deg + 1;
  for(int j = 0; j < npg; j++) {
    if (j != ii) {
      *p *= (x - subdiv[j]) / (subdiv[ii] - subdiv[j]);
    }
  }
}

// interpolate another function from a different mesh
// xi: interp points
// w_sav: stored data
// deg_sav, nbnodes_sav: quantities of saved data
// x: point where is computed the interp values
// ws: computed interpolation
void interpolate(double *xi, dcmplx *wi,
		 int deg, int nbnodes,
		 double x, dcmplx *wloc){

  double xmin = xi[0];
  double xmax = xi[nbnodes - 1];
  int nbelems = nbnodes / (deg + 1);
  assert(nbnodes % (deg + 1) == 0);

  assert(x >= xmin && x <= xmax);
  if (x == xmax) x = x - 1e-10 * (xmax-xmin);

  
  int numelem = (x-xmin)/(xmax-xmin) * nbelems;

  for(int iv = 0; iv < M; iv++) wloc[iv] = 0;
  
  for(int i = 0; i < deg + 1; i++){
    double val;
    double *sub = xi + numelem * (deg + 1);
    lagrange(&val, sub,
	     deg, i, x);
    //printf("x=%f ie=%d i=%d xi=%f val=%f\n",x,numelem,i,sub[i],val);
    int ino = numelem * (deg + 1) + i;
    for(int iv = 0; iv < M; iv++)
      wloc[iv] += wi[iv * nbnodes + ino] * val;
  }

  //printf("w=%f %f %f\n\n",creal(wloc[0]),creal(wloc[1]),creal(wloc[2]));


}


// compute and invert the local implicit matrix
void loc_assembly(galerkin *gal, dcmplx *mat, double vit){
  
  int n = DEG + 1;

  dcmplx a[n * n];
  for(int i = 0; i < n * n; i++) a[i] = 0;

  for(int i = 0; i <= DEG; i++) a[n * i + i] = 1; 

  dcmplx z = vit * gal->dt / gal->dx * THETA;
  double wg = wglop(DEG, 0);
  
  if (vit > 0) {
    a[0] += z / wg;
  } else if (vit < 0) {
    a[n * n -1] -= z / wg;
  } else{
    //assert(vit != 0);
  }
  
  for(int i = 0; i <= DEG; i++){
    for(int j = 0; j <= DEG; j++){
      a[n * i + j] += z * dlag(DEG, j, i);
    }
  }

  InvertSquare(n, a, mat); 

  // test
  /* for(int i = 0; i <= DEG; i++){ */
  /*   for(int j = 0; j <= DEG; j++){ */
  /*     dcmplx v = 0; */
  /*     for(int k =  0; k <= DEG; k++){ */
  /* 	v += a[n * i + k] * mat[n * k + j]; */
  /*     } */
  /*     v=a[n * i + j]; */
  /*     printf("i=%d j=%d aij=(%f,%f)\n",i,j, */
  /* 	     creal(v),cimag(v)); */
  /*   } */
  /* } */
  //assert(1==3);


}

// compute and invert the local implicit matrix
void loc_assembly_nonlin(galerkin *gal, dcmplx *mat, double vit){
  
  int n = M * (DEG + 1);

  dcmplx a[n * n];
  for(int i = 0; i < n * n; i++) a[i] = 0;

  for(int i = 0; i < n; i++) a[n * i + i] = 1; 

  dcmplx z = vit * gal->dt / gal->dx * THETA;
  double wg = wglop(DEG, 0);

  int iv;
  
  if (vit > 0) {
    int iloc = 0;
    iv = 2;
    int i = iloc * M + iv;  
    a[n * i + i] += z / wg;
  } else if (vit < 0) {
    int iloc = DEG;
    iv = 0;
    int i = iloc * M + iv;  
    a[n * i + i] -= z / wg;
  } else{
    iv = 1;
    z = 0;
    //assert(vit != 0);
  }
  
  for(int iloc = 0; iloc <= DEG; iloc++){
    int i = iloc * M + iv;  
    for(int jloc = 0; jloc <= DEG; jloc++){
      int j = jloc * M + iv;  
      a[n * i + j] += z * dlag(DEG, jloc, iloc);
    }
  }

  InvertSquare(n, a, mat); 

  // test
  /* for(int i = 0; i <= DEG; i++){ */
  /*   for(int j = 0; j <= DEG; j++){ */
  /*     dcmplx v = 0; */
  /*     for(int k =  0; k <= DEG; k++){ */
  /* 	v += a[n * i + k] * mat[n * k + j]; */
  /*     } */
  /*     v=a[n * i + j]; */
  /*     printf("i=%d j=%d aij=(%f,%f)\n",i,j, */
  /* 	     creal(v),cimag(v)); */
  /*   } */
  /* } */
  //assert(1==3);


}


void gal_step(galerkin *gal, dcmplx tnow)
{


  int n = DEG + 1;

  dcmplx mat[n * n];
  // negative velocity

  for(int ivel = 1; ivel >= -1; ivel--){
    int iv = ivel + 1;
    loc_assembly(gal, mat, ivel * VMAX);

    for(int ieref = 0; ieref < NB_ELEMS ; ieref++) {
      if (ivel == 0) break;
      int ie = ieref;
      if (ivel < 0) ie = NB_ELEMS -1 - ieref;

      int pm = (ivel + 1) / 2; // 0 if negative 1 if positive 
    
      // last glop in cell
      int ino = connec(ie, DEG * (1 - pm));
      int iw = vindex(ino, iv);
      dcmplx wvnm1[M];
      dcmplx wvn[M];
      dcmplx res[DEG + 1];
      // boundary condition or upwind data
      if (ieref > 0) {
	wvnm1[iv] = gal->wnm1[vindex(connec(ie - ivel, pm * DEG),iv)];
	wvn[iv] = gal->wn[vindex(connec(ie - ivel, pm * DEG),iv)];
      } else {
	solexacte(gal->xnode[ino],tnow,wvnm1);
	solexacte(gal->xnode[ino],tnow + gal->dt,wvn);
	//printf("wR=(%f,%f)\n",creal(wv[iv]),cimag(wv[iv]));
      }

      // explicit part
      for(int iloc = 0; iloc <= DEG; iloc++)
	res[iloc] = 0; 

      res[DEG * (1- pm)] += (gal->wnm1[iw] - wvnm1[iv]) 
	* VMAX * gal->dt * (THETA - 1) / gal->dx / wglop(DEG, 0);

      for(int iloc = 0; iloc <= DEG; iloc++) {
	for(int jloc = 0; jloc <= DEG; jloc++) {
	  res[iloc] += ivel * VMAX * gal->dt * (THETA - 1) / gal->dx *
	    dlag(DEG, jloc, iloc) * gal->wnm1[iw - DEG * (1 - pm) + jloc];
	}
      }

      gal->wnm1[iw] -= wvn[iv] *  (-VMAX) *
	gal->dt * THETA / gal->dx / wglop(DEG, DEG);
      // local implicit scheme
      iw -= DEG * (1 - pm); // return to first glop for matrix product
      for(int iloc = 0; iloc <= DEG; iloc++){
	gal->wn[iw + iloc] = 0;
	for(int jloc = 0; jloc <= DEG; jloc++)
	  gal->wn[iw + iloc] += mat[n * iloc + jloc] *
	    (gal->wnm1[iw + jloc] + res[jloc]) ;
      }

    }

  }

  // update
  for(int iw = 0; iw < WSIZE; iw++){
    gal->wnm1[iw] = gal->wn[iw];
  }

  // for(int ino = 0; ino < NB_NODES; ino++){ */
  //  dcmplx wloc[M]; 
  //  for(int iv = 0; iv < M; iv++) wloc[iv] =gal->wn[vindex(ino,iv)];  
  //  bgk_project(wloc); 
  //   for(int iv = 0; iv < M; iv++) gal->wn[vindex(ino,iv)]=wloc[iv];  
  //}


  // update
  for(int iw = 0; iw < WSIZE; iw++) {
    gal->wnm1[iw] = gal->wn[iw];
  }
}


void gal_step_nonlin(galerkin *gal, dcmplx tnow)
{


  int n = DEG + 1;

  dcmplx mat[n * n];

  for(int ivel = -1; ivel <= 1; ivel++){
    int iv = ivel + 1;
    loc_assembly(gal, mat, ivel * VMAX);

    for(int ieref = 0; ieref < NB_ELEMS ; ieref++) {
      if (ivel == 0) break;
      int ie = ieref;
      if (ivel < 0) ie = NB_ELEMS -1 - ieref;

      int pm = (ivel + 1) / 2; // 0 if negative 1 if positive 
    
      // last glop in cell
      int ino = connec(ie, DEG * (1 - pm));
      //int iw = vindex(ino, iv);
      dcmplx wvnm1[M];
      dcmplx wvn[M];
      //dcmplx res[DEG + 1];
      // boundary condition or upwind data
      if (ieref > 0) {
	wvnm1[iv] = gal->wnm1[vindex(connec(ie - ivel, pm * DEG),iv)];
	wvn[iv] = gal->wn[vindex(connec(ie - ivel, pm * DEG),iv)];
      } else {
	solexacte(gal->xnode[ino],tnow,wvnm1);
	solexacte(gal->xnode[ino],tnow+gal->dt,wvn);
	//printf("wR=(%f,%f)\n",creal(wv[iv]),cimag(wv[iv]));
      }

      dcmplx wloc[(DEG + 1) * M];
      for(int iloc = 0; iloc <= DEG; iloc++){
	for(int ivar = 0; ivar < M; ivar++){
	  wloc[ M * iloc + ivar] =
	    gal->wnm1[vindex(connec(ie, iloc),ivar)];
	}
      }


      local_dg(gal, wloc, wvnm1[iv], wvn[iv], ivel);
      
      // explicit part
      /* for(int iloc = 0; iloc <= DEG; iloc++)
	 res[iloc] = 0; 

	 res[DEG * (1- pm)] += (gal->wnm1[iw] - wvnm1[iv]) 
	 * VMAX * gal->dt * (THETA - 1) / gal->dx / wglop(DEG, 0);

	 for(int iloc = 0; iloc <= DEG; iloc++) {
	 for(int jloc = 0; jloc <= DEG; jloc++) {
	 res[iloc] += ivel * VMAX * gal->dt * (THETA - 1) / gal->dx *
	 dlag(DEG, jloc, iloc) * gal->wnm1[iw - DEG * (1 - pm) + jloc];
	 !	}
	 }

	 gal->wnm1[iw] -= wvn[iv] *  (-VMAX) *
	 gal->dt * THETA / gal->dx / wglop(DEG, DEG);
	 // local implicit scheme
	 iw -= DEG * (1 - pm); // return to first glop for matrix product
	 for(int iloc = 0; iloc <= DEG; iloc++){
	 gal->wn[iw + iloc] = 0;
	 for(int jloc = 0; jloc <= DEG; jloc++)
	 gal->wn[iw + iloc] += mat[n * iloc + jloc] *
	 (gal->wnm1[iw + jloc] + res[jloc]) ;
	 }*/

      for(int iloc = 0; iloc <= DEG; iloc++){
      	for(int ivar = 0; ivar < M; ivar++){
	  
	  gal->wn[vindex(connec(ie, iloc),ivar)] =
	    wloc[ M * iloc + ivar];
      	}
      }


    }
    // update
    for(int ie = 0; ie < NB_ELEMS ; ie++) {
      for(int iloc = 0; iloc <= DEG; iloc++){
  	for(int iw = 0; iw < WSIZE; iw++){
  	  int iw = vindex(connec(ie, iloc),iv);
  	  gal->wnm1[iw] = gal->wn[iw];
  	}
      }
    }
  

  }

  // update
  /*   for(int ie = 0; ie < NB_ELEMS ; ie++) { */
  /*     for(int iloc = 0; iloc <= DEG; iloc++){ */
  /* 	for(int iw = 0; iw < WSIZE; iw++){ */
  /* 	  int iw = vindex(connec(ie, iloc),iv); */
  /* 	  gal->wnm1[iw] = gal->wn[iw]; */
  /* 	} */
  /*     } */
  /*   } */
  
  /* } */
  


  
  for(int ino = 0; ino < NB_NODES; ino++){
    dcmplx wloc[M];
    for(int iv = 0; iv < M; iv++) wloc[iv] =gal->wnm1[vindex(ino,iv)]; 
    //bgk_relax(wloc, gal->dt);
    //bgk_project(wloc);
    for(int iv = 0; iv < M; iv++) gal->wn[vindex(ino,iv)]=wloc[iv]; 
  }


  // update
  for(int iw = 0; iw < WSIZE; iw++) {
    gal->wnm1[iw] = gal->wn[iw];
  }
}

void local_dg(galerkin *gal, dcmplx *wloc, dcmplx wvnm1, dcmplx wvn, int ivel){

  int n = M * (DEG + 1);
  dcmplx mat[n * n];
  loc_assembly_nonlin(gal, mat, ivel * VMAX);
  // explicit part
  dcmplx res[n];
  dcmplx res_imp[n];
  for(int i = 0; i < n; i++){
    res[i] = wloc[i];
    res_imp[i] = 0;
  }
  
  int iv = ivel + 1;
  int pm = (ivel + 1) / 2; // 0 if negative 1 if positive 

  
  res[M * DEG * (1- pm) + iv] +=  (wloc[M * DEG * (1- pm) + iv] - wvnm1)
    * VMAX * gal->dt * (THETA - 1) / gal->dx / wglop(DEG, 0);

  for(int iloc = 0; iloc <= DEG; iloc++) {
    int i = iloc * M + iv;
    for(int jloc = 0; jloc <= DEG; jloc++) {
      int j = jloc * M + iv;
      res[i] += ivel * VMAX * gal->dt * (THETA - 1) / gal->dx *
	dlag(DEG, jloc, iloc) * wloc[j];
    }
  }

  res_imp[M * DEG * (1- pm) + iv] -= wvn *  (-VMAX) *
    gal->dt * THETA / gal->dx / wglop(DEG, DEG);

  
  // local implicit scheme
  for(int iloc = 0; iloc <= DEG; iloc++){
    int i = M * iloc + iv;
    wloc[i] = 0;
    for(int jloc = 0; jloc <= DEG; jloc++){
      int j = M * jloc + iv;
      wloc[i] += mat[n * i + j] * (res_imp[j] + res[j]);	
    }
  }

}

 
  






void gal_plot(galerkin *gal)
{

  // commande de tracé
  // plot 'fe.dat'  using 1:9 w l, 'rho_u.dat' using 1:3 w l
  // plot 'fe.dat'  using 1:8 w l, 'rho_u.dat' using 1:2 w l
  FILE * gnufile;
  gnufile = fopen("fe.dat", "w" );

  FILE * ru_file;
  ru_file = fopen("fe.bin", "wb" );

  int temp = NB_ELEMS;
  fwrite(&temp, sizeof(int), 1, ru_file);
  temp = DEG;
  fwrite(&temp, sizeof(int), 1, ru_file);


  dcmplx wex[M];
  dcmplx wnum[M];
  for(int ino = 0; ino < NB_NODES; ino++){
    double x = gal->xnode[ino];
    fprintf(gnufile, "%f ", x);
    fwrite(&x, sizeof(double), 1, ru_file);
    double t=gal->t;
    solexacte(x, t, wex);
    for(int iv = 0; iv < M; iv++){
      int imem = vindex(ino,iv);
      wnum[iv] =  gal->wn[imem];
    }
    for(int iv = 0; iv < M; iv++){
      fprintf(gnufile, "%f ", creal(wex[iv]));
      fprintf(gnufile, "%f ", creal(wnum[iv]));
      fwrite(wnum + iv, sizeof(dcmplx), 1, ru_file);
    }
    dcmplx rho = wnum[0] + wnum[1] + wnum[2];
    dcmplx u = (wnum[2] - wnum[0]) / rho;
    dcmplx rhoex = wex[0] + wex[1] + wex[2];
    dcmplx uex = (wex[2] - wex[0]) / rho;
    fprintf(gnufile, "%f ", creal(rho));      // var 8 9
    fprintf(gnufile, "%f ", creal(u));        
    fprintf(gnufile, "%f ", creal(rhoex));      // var 10 11
    fprintf(gnufile, "%f ", creal(uex));        
    fwrite(&rho, sizeof(dcmplx), 1, ru_file);
    fwrite(&u, sizeof(dcmplx), 1, ru_file);
    fprintf(gnufile, "%f ", cimag(rho));
    fprintf(gnufile, "%f ", cimag(u));

    fprintf(gnufile, "\n");
  }
    
  fclose(gnufile);
  fclose(ru_file);


}

void solexacte(double x, dcmplx t, dcmplx* w){

  // smooth solution test case
  /* double xi = x + VMAX * t; */
  /* w[0] = 1 + cexp(- 30 * xi * xi); */
  /* xi = x; */
  /* w[1] = 1 + cexp(- 30 * xi * xi); */
  /* xi = x - VMAX * t; */
  /* w[2] = 1 + cexp(- 30 * xi * xi); */
  /* bgk_project(w); */
  double xi = x + VMAX * t;
  w[0] = 1 + cexp(- 30 * xi * xi);
  xi = x;
  w[1] = 2 * (1 - CSON * CSON) / CSON / CSON * (1 + cexp(- 30 * xi * xi));
  xi = x - VMAX * t;
  w[2] = 1 + cexp(- 30 * xi * xi);

  /* bgk_relax(w, 1.); */

  // debug test case
  /* double xi = x + VMAX * t; */
  /* xi = x; */
  /* w[0] = (xi-2)*(xi+2); */
  /* xi = x; */
  /* w[1] = (xi-2)*(xi+2); */
  /* xi = x; */
  /* w[2] = (xi-2)*(xi+2); */
  /* w[2] = 1; */
  
  
  // Riemann problem test case
  /* double xi; */
  /* xi = x + VMAX * t; */
  /* if (xi < 0) { */
  /*   w[0] = RL * UL * (UL - 1) / 2 + CSON * CSON * RL / 2; */
  /* } else { */
  /*   w[0] = RR * UR * (UR - 1) / 2 + CSON * CSON * RR / 2; */
  /* } */
  /* xi = x; */
  /* if (xi < 0) { */
  /*   w[1] = RL * (1 - UL * UL - CSON * CSON); */
  /* } else { */
  /*   w[1] = RR * (1 - UR * UR - CSON * CSON); */
  /* } */
  /* xi = x - VMAX * t; */
  /* if (xi < 0) { */
  /*   w[2] = RL * UL * (UL + 1) / 2 + CSON * CSON * RL / 2; */
  /* } else { */
  /*   w[2] = RR * UR * (UR + 1) / 2 + CSON * CSON * RR / 2; */
  /* } */

  /* dcmplx t0=2.; */
  /* dcmplx x0=1.; */
  /* dcmplx xi = (x - x0) / (t + t0); */
  /* dcmplx u = xi + CSON; */
  /* dcmplx r = RL * cexp((UL - u) / CSON); */
  /* w[0] = r * u * (u - 1) / 2 + CSON * CSON * r / 2; */
  /* w[1] = r * (1 - u * u - CSON * CSON); */
  /* w[2] = r * u * (u + 1) / 2 + CSON * CSON * r / 2; */
  //bgk_project(w);

}

void bgk_relax(dcmplx *w, dcmplx dt){

  dcmplx r = w[0] + w[1] + w[2];
  dcmplx q = w[2] - w[0];
  dcmplx u = q / r;

  dcmplx weq[3] = {q * (u - 1) / 2 + CSON * CSON * r / 2,
		   r * (1 - u * u - CSON * CSON),
		   q * (u + 1) / 2 + CSON * CSON * r / 2};

  dcmplx dw[3] = { w[0] - weq[0],
		   w[1] - weq[1],
		   w[2] - weq[2]};

  
  dcmplx z = dt / 2 * RELAX;
  w[0] = (1 - z) / (1 + z) * w[0] + 2 * z /(1 + z) * weq[0];
  w[1] = (1 - z) / (1 + z) * w[1] + 2 * z /(1 + z) * weq[1];
  w[2] = (1 - z) / (1 + z) * w[2] + 2 * z /(1 + z) * weq[2];

  /* w[0] = - w[0] + 2 * weq[0]; */
  /* w[1] = - w[1] + 2 * weq[1]; */
  /* w[2] = - w[2] + 2 * weq[2]; */
  /* dcmplx z = dt / (tau + dt / 2); */
  /* w[0] -= z * dw[0]; */
  /* w[1] -= z * dw[1]; */
  /* w[2] -= z * dw[2]; */

}


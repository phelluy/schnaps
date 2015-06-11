#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "skyline.h"


int sol_(real *vkgs, real *vkgd, real *
	vkgi, real *vfg, int *kld, real *vu, int neq, 
	 int ifac, int isol, int nsym, real *
	 energ, int *ier);

void InitSkyline(Skyline* sky, int n){

  sky->is_alloc=false;
  sky->is_sym=false;
  sky->is_lu=false;

  sky->neq=n;

  sky->nmem=0;

  sky->vkgs=NULL;

  sky->vkgd=malloc(n*sizeof(real));
  assert(sky->vkgd);
  for(int i=0;i<n;i++) sky->vkgd[i]=0;

  sky->vkgi=NULL;

  sky->prof=malloc(n*sizeof(int));
  assert(sky->prof);
  for(int i=0;i<n;i++) sky->prof[i]=0;

  sky->kld=malloc((n+1)*sizeof(int));
  assert(sky->kld);
  for(int i=0;i<n+1;i++) sky->kld[i]=0;

}


void SwitchOn(Skyline* sky,int i,int j){

  // update the profile
  sky->prof[j]= j-i > sky->prof[j] ? j-i : sky->prof[j] ;
  sky->prof[i]= i-j > sky->prof[i] ? i-j : sky->prof[i] ;

} 

void AllocateSkyline(Skyline* sky){

  assert(! sky->is_alloc);

  sky->kld[0]=0;
  for(int i=0;i<sky->neq;i++){
    sky->kld[i+1]=sky->kld[i]+sky->prof[i];
  }
  sky->nmem=sky->kld[sky->neq];

  sky->vkgs=malloc(sky->nmem * sizeof(real));
  assert(sky->vkgs);

  if (! sky->is_sym){
    sky->vkgi=malloc(sky->nmem * sizeof(real));
    assert(sky->vkgi);
  }
  else{
    sky->vkgi=sky->vkgs;
  }

  sky->is_alloc=true;

  // fill with zeros
  for(int k=0;k<sky->nmem;k++){
    sky->vkgs[k]=0;
     if (! sky->is_sym) sky->vkgi[k]=0;
  }
  

}

void SetSkyline(Skyline* sky,int i,int j,real val){

  assert(sky->is_alloc);


  if (i==j){
    sky->vkgd[i]+=val;
  }
  else if (j>i){
    int k=sky->kld[j+1]-j+i;
    sky->vkgs[k]+=val;
  }
  else {
    assert(!(sky->is_sym));
    int k=sky->kld[i+1]-i+j;
    sky->vkgi[k]+=val;
    //printf("i=%d j=%d k=%d nmem=%d\n",i,j,k,sky->nmem);
    //printf("i=%d j=%d k=%d v=%f\n",i,j,k,sky->vkgi[k]);
  }


} 

real GetSkyline(Skyline* sky,int i,int j){

  if (sky->is_sym && i>j){
    int temp=i;
    i=j;
    j=temp;
  }

  if (j-i > sky->prof[j] || i-j > sky->prof[i]){
    return 0;
  }
  else if (i==j){
    return sky->vkgd[i];
  }
  else if ( j>i){
    int k=sky->kld[j+1]-j+i;
    return sky->vkgs[k];
  }
  else {
    int k=sky->kld[i+1]-i+j;
    return sky->vkgi[k];
  }


}

void DisplaySkyline(Skyline* sky){

  int n=sky->neq;


  printf("profil=");
  for(int i=0;i<n;i++){
    printf("%d ",sky->prof[i]);
  }
  printf("\n");

  printf("kld=");
  for(int i=0;i<n+1;i++){
    printf("%d ",sky->kld[i]);
  }
  printf("\n");

  printf("vkgd=");
  for(int i=0;i<n;i++){
    printf("%f ",sky->vkgd[i]);
  }
  printf("\n");

  printf("vkgs=");
  for(int i=0;i<sky->nmem;i++){
    printf("%f ",sky->vkgs[i]);
  }
  printf("\n");

  printf("vkgi=");
  for(int i=0;i<sky->nmem;i++){
    printf("%f ",sky->vkgi[i]);
  }
  printf("\n");

  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf("%f ",GetSkyline(sky,i,j));
    }
    printf("\n");
  }

    
}



void FactoLU(Skyline* sky){

  real* vfg=NULL;
  real* vu=NULL;
  real energ;
  int ier;
  int ifac=1;
  int isol=0;
  int nsym=1;
  if (sky->is_sym) nsym=0;

  sol_(sky->vkgs,sky->vkgd, sky->vkgi,
       vfg, sky->kld, vu, sky->neq, 
        ifac, isol, nsym,
       &energ, &ier);

  sky->is_lu=true;

}

void SolveSkyline(Skyline* sky,real* vfg,real* vu){
  assert(sky->is_lu);

  real energ;
  int ier;
  int ifac=0;
  int isol=1;
  int nsym=1;
  if (sky->is_sym) nsym=0;

  sol_(sky->vkgs,sky->vkgd, sky->vkgi,
       vfg, sky->kld, vu, sky->neq, 
        ifac, isol, nsym,
       &energ, &ier);



}




void FreeSkyline(Skyline* sky){

  assert(sky->is_alloc);

  free(sky->vkgs);
  //printf("vkgd=%p\n",sky->vkgd);
  free(sky->vkgd);
  //assert(1==2);
  if (! sky->is_sym)  free(sky->vkgi);
  free(sky->prof);
  free(sky->kld);

  sky->is_alloc=false;
  //InitSkyline(sky,sky->neq);

}




/* sol.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

//#include "f2c.h"

/* Table of constant values */

static int c__1 = 1;

/* Subroutine */ int sol_(real *vkgs, real *vkgd, real *
	vkgi, real *vfg, int *kld, real *vu, int neq, 
	 int ifac, int isol, int nsym, real *
	energ, int *ier)
{
    /* Initialized data */

    static real vzero = 0.;

    /* Format strings */
    static char fmt_8000[] = "sol pivot nul equation";

    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Builtin functions */
    //int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static int i__;
    static real c1, c2;
    static int j1, j2, ic, ij, ik, jbk, jck, jhj, jhk, lhk, jhj1, jhk1, 
	    lhk1;
    extern real scal_(real *, real *, int *);
    static int imin, imax, imin1;
    static real cdiag;

    /* Fortran I/O blocks */
    //static cilist io___22 = { 0, 0, 0, fmt_8000, 0 };


/*   resolution d'un systeme lineaire symetrique ou non. la matrice est */
/*   stockee par ligne de ciel,en memoire dans les tables vkgs,vkgd,vkgi */

/*       entrees */
/*          vkgs,vkgd,vkgi    matrice du systeme : parties superieure, */
/*                            diagonale, inferieure (real precision) */
/*          vfg               second membre */
/*          kld               pointeurs vers les hauts de colonne */
/*          vu                vecteur solution (qui peut etre vfg) */
/*          neq               nombre d'equations */
/*          mp                unite logique d'impression */
/*          ifac              si ifac.eq.1 triangularisation de */
/*                            la matrice */
/*          isol              si isol.eq.1 calcul de la solution a */
/*                            partir de la matrice triangularisee */
/*          nsym              indice de probleme non symetrique */
/*       sorties */
/*          vkgs,vkgd,vkgi    matrice triangularisee (si ifac.eq.1) */
/*          vfg               solution (si isol.eq.1) */
/*          energ             energie du systeme (si nsym.eq.0) */
/*          ier               mis a 1 si pivot nul rencontre */

/* =========================== debut des declarations ==================== */
    /* Parameter adjustments */
    --vu;
    --kld;
    --vfg;
    --vkgi;
    --vkgd;
    --vkgs;


#define _Z 1

    /* Function Body */
/* =========================== debut du code executable ================== */

/* -------  traitement */

    ik = 1;
    if (vkgd[1] == vzero) {
	goto L800;
    }
    *energ = vzero;
    *ier = 0;
    if (isol == 1) {
	i__1 = neq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vu[i__] = vfg[i__];
	}
    }

/* -------  pour chaque colonne ik a modifier */

    jhk = 1;
    i__1 = neq;
    for (ik = 2; ik <= i__1; ++ik) {
      //printf("factolu %d/%d\n",ik,neq);

/* -------  pointeur du haut de la colonne suivante ik+1 */

	jhk1 = kld[ik + 1]+_Z;

/* -------  hauteur de la colonne ik (hors termes superieur et diagonal) */

	lhk = jhk1 - jhk;
	lhk1 = lhk - 1;

/* -------  ligne du premier terme a modifier dans la colonne ik */

	imin = ik - lhk1;
	imin1 = imin - 1;

/* -------  ligne du dernier terme a modifier dans la colonne ik */

	imax = ik - 1;
	if (lhk1 < 0) {
	    goto L100;
	}
	if (ifac != 1) {
	    goto L90;
	}
	if (nsym == 1) {
	    vkgi[jhk] /= vkgd[imin1];
	}
	if (lhk1 == 0) {
	    goto L40;
	}

/* -------  modifier les termes non diagonaux de la colonne ik */

	jck = jhk + 1;
	jhj = kld[imin]+_Z;

/* -------  pour chaque terme place en jck, correspondant a la colonne ij */

	i__2 = imax;
	for (ij = imin; ij <= i__2; ++ij) {
	    jhj1 = kld[ij + 1]+_Z;

/* -------  nombre de termes modificatifs du terme place en jck */

/* Computing MIN */
	    i__3 = jck - jhk, i__4 = jhj1 - jhj;
	    //ic = min(i__3,i__4);
	    ic = i__3 < i__4 ? i__3 : i__4;
	    if (ic <= 0 && nsym == 0) {
		goto L20;
	    }
	    c1 = vzero;
	    if (ic <= 0) {
		goto L17;
	    }
	    j1 = jhj1 - ic;
	    j2 = jck - ic;
	    if (nsym == 1) {
		goto L15;
	    }
	    vkgs[jck] -= scal_(&vkgs[j1], &vkgs[j2], &ic);
	    goto L20;
L15:
	    vkgs[jck] -= scal_(&vkgi[j1], &vkgs[j2], &ic);
	    c1 = scal_(&vkgs[j1], &vkgi[j2], &ic);
L17:
	    vkgi[jck] = (vkgi[jck] - c1) / vkgd[ij];
L20:
	    ++jck;
/* L30: */
	    jhj = jhj1;
	}

/* -------  modifier le terme diagonal */

L40:
	jck = jhk;
	cdiag = vzero;
	i__2 = imax;
	for (ij = imin1; ij <= i__2; ++ij) {
	    c1 = vkgs[jck];
	    if (nsym == 1) {
		goto L50;
	    }
	    c2 = c1 / vkgd[ij];
	    vkgs[jck] = c2;
	    goto L60;
L50:
	    c2 = vkgi[jck];
L60:
	    cdiag += c1 * c2;
/* L70: */
	    ++jck;
	}
	vkgd[ik] -= cdiag;
	if (vkgd[ik] == 0.f) {
	    goto L800;
	}

/* -------  resolution du systeme triangulaire inferieur */

L90:
	if (isol != 1) {
	    goto L100;
	}
	if (nsym != 1) {
	    vu[ik] = vfg[ik] - scal_(&vkgs[jhk], &vu[imin1], &lhk);
	}
	if (nsym == 1) {
	    vu[ik] = vfg[ik] - scal_(&vkgi[jhk], &vu[imin1], &lhk);
	}
L100:
	jhk = jhk1;
    }
    if (isol != 1) {
	goto L9999;
    }

/* -------  resolution du systeme diagonal : */

    if (nsym == 1) {
	goto L120;
    }
    i__1 = neq;
    for (ik = 1; ik <= i__1; ++ik) {
	c1 = vkgd[ik];
	if (c1 == vzero) {
	    goto L800;
	}
	c2 = vu[ik] / c1;
	vu[ik] = c2;
/* L110: */
	*energ += c1 * c2 * c2;
    }

/* -------  resolution du systeme triangulaire superieur */

L120:
    ik = neq + 1;
    jhk1 = kld[ik]+_Z;
L130:
    --ik;
    if (nsym == 1) {
	vu[ik] /= vkgd[ik];
    }
    if (ik == 1) {
	goto L9999;
    }
    c1 = vu[ik];
    jhk = kld[ik]+_Z;
    jbk = jhk1 - 1;
    if (jhk > jbk) {
	goto L150;
    }
    ij = ik - jbk + jhk - 1;
    i__1 = jbk;
    for (jck = jhk; jck <= i__1; ++jck) {
	vu[ij] -= vkgs[jck] * c1;
/* L140: */
	++ij;
    }
L150:
    jhk1 = jhk;
    goto L130;

/* -------  erreurs */

L800:
	/* io___22.ciunit = *mp; */
    printf("%s %d\n",fmt_8000,ik);
	/* s_wsfe(&io___22); */
	/* do_fio(&c__1, (char *)&ik, (ftnlen)sizeof(int)); */
	/* e_wsfe(); */
	
    *ier = 1;
    goto L9999;

/* -------  fin */

L9999:
    return 0;
/* ===========================   fin du module sol    ================== */
} /* sol_ */

real scal_(real *x, real *y, int *n)
{
    /* Initialized data */

    static real zero = 0.;

    /* System generated locals */
    int i__1;
    real ret_val;

    /* Local variables */
    static int i__;

/* ======================================================================= */
/* calcul du produit scalaire */
/* ======================================================================= */
    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
/* ----------------------------------------------------------------------- */
    ret_val = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += x[i__] * y[i__];
    }
    return ret_val;
} /* scal_ */


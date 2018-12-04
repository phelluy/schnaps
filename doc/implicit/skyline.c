#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "skyline.h"
#include "math.h"


int sol_(dcmplx *vkgs, dcmplx *vkgd, dcmplx *
	 vkgi, dcmplx *vfg, int *kld, dcmplx *vu, int neq, 
	 int ifac, int isol, int nsym, dcmplx *
	 energ, int *ier);

int mulku_(dcmplx *vkgs, dcmplx *vkgd, dcmplx *
	   vkgi, int *kld, dcmplx *vfg, int neq, int nsym, 
	   dcmplx *vres, int nsky);

dcmplx scal_(dcmplx *x, dcmplx *y, int *n);


void InitSkyline(Skyline* sky, int n){

  sky->is_alloc=false;
  sky->copy_is_alloc=false;
  sky->is_sym=false;
  sky->is_lu=false;

  sky->neq=n;

  sky->nmem=0;

  sky->vkgs=NULL;
  sky->copy_vkgs=NULL;

  sky->vkgd=calloc(n,sizeof(dcmplx));
  assert(sky->vkgd);
  for(int i=0;i<n;i++) sky->vkgd[i]=0;
  sky->copy_vkgd=NULL;

  sky->vkgi=NULL;
  sky->copy_vkgi=NULL;

  sky->prof=calloc(n,sizeof(int));
  assert(sky->prof);
  for(int i=0;i<n;i++) sky->prof[i]=0;

  sky->kld=calloc((n+1),sizeof(int));
  assert(sky->kld);
  for(int i=0;i<n+1;i++) sky->kld[i]=0;

}


void SwitchOn(Skyline* sky,int i,int j){

  // update the profile
  sky->prof[j]= j-i > sky->prof[j] ? j-i : sky->prof[j] ;
  sky->prof[i]= i-j > sky->prof[i] ? i-j : sky->prof[i] ;

} 

void AllocateSkyline(Skyline* sky){

  assert(!sky->is_alloc);

  sky->kld[0]=0;
  for(int i=0;i<sky->neq;i++){
    sky->kld[i+1]=sky->kld[i]+sky->prof[i];
  }
  sky->nmem=sky->kld[sky->neq];

  sky->vkgs=calloc(sky->nmem,sizeof(dcmplx));
  assert(sky->vkgs);

  if (! sky->is_sym){
    sky->vkgi=calloc(sky->nmem,sizeof(dcmplx));
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

void AllocateCopySkyline(Skyline* sky){

  assert(!sky->copy_is_alloc);

  sky->copy_vkgs=calloc(sky->nmem,sizeof(dcmplx));
  assert(sky->copy_vkgs);

  if (! sky->is_sym){
    sky->copy_vkgi=calloc(sky->nmem,sizeof(dcmplx));
    assert(sky->copy_vkgi);
  }
  else{
    sky->copy_vkgi=sky->copy_vkgs;
  }

  sky->copy_vkgd=calloc(sky->neq,sizeof(dcmplx));
  assert(sky->copy_vkgd);
  for(int i=0;i<sky->neq;i++) sky->copy_vkgd[i]=0;

  sky->copy_is_alloc=true;

  // fill with zeros
  for(int k=0;k<sky->nmem;k++){
    sky->copy_vkgs[k]=0;
    if (! sky->is_sym) sky->copy_vkgi[k]=0;
  }

  // Copying right after allocating.
  for (int i=0; i<sky->nmem; i++){
    sky->copy_vkgs[i] = sky->vkgs[i];
    if (!sky->is_sym) sky->copy_vkgi[i] = sky->vkgi[i];
  }
  for (int i=0; i<sky->neq; i++){
    sky->copy_vkgd[i] = sky->vkgd[i];
  }

}

void ZeroSkyline(Skyline* sky){

  assert(sky->is_alloc);

  for(int i=0;i<sky->neq;i++) sky->vkgd[i]=0;

  for(int k=0;k<sky->nmem;k++) sky->vkgs[k]=0;

  if (!sky->is_sym){
    for(int k=0;k<sky->nmem;k++) sky->vkgi[k]=0;
  }

}

void AddSkyline(Skyline* sky,int i,int j,dcmplx val){

  assert(sky->is_alloc);
  
  if ((j-i > sky->prof[j] || i-j > sky->prof[i]) && (val != 0)){
    printf("problem of profil with add %d %d\n",i,j);
  }
  
  if ((j-i > sky->prof[j] || i-j > sky->prof[i]) && val==0.0)
    {
      ;
    }
  else if (i==j){
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

void SetSkyline(Skyline* sky,int i,int j,dcmplx val){

  assert(sky->is_alloc);

  if ((j-i > sky->prof[j] || i-j > sky->prof[i]) && (val != 0)){
    printf("problem of profil with set %d %d\n",i,j);
  }
  
  if ((j-i > sky->prof[j] || i-j > sky->prof[i]) && val==0.0)
    {
      ;
    }
  else if (i==j){
    sky->vkgd[i]=val;
  }
  else if (j>i){
    int k=sky->kld[j+1]-j+i;
    sky->vkgs[k]=val;
  }
  else {
    assert(!(sky->is_sym));
    int k=sky->kld[i+1]-i+j;
    sky->vkgi[k]=val;
    //printf("i=%d j=%d k=%d nmem=%d\n",i,j,k,sky->nmem);
    //printf("i=%d j=%d k=%d v=%f\n",i,j,k,sky->vkgi[k]);
  }
  
  
} 

dcmplx GetSkyline(Skyline* sky,int i,int j){

  if (sky->is_sym && i>j){
    int temp=i;
    i=j;
    j=temp;
  }

  if (j-i > sky->prof[j] || i-j > sky->prof[i]){
    return 0.0;
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

  bool full = true;

  if (full){

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
      printf("(%.5e,%.5e) ",creal(sky->vkgd[i]),cimag(sky->vkgd[i]));
    }
    printf("\n");
  
    printf("vkgs=");
    for(int i=0;i<sky->nmem;i++){
      printf("(%.5e,%.5e) ",creal(sky->vkgs[i]),cimag(sky->vkgs[i]));
    }
    printf("\n");

    printf("vkgi=");
    for(int i=0;i<sky->nmem;i++){
      printf("(%.5e,%.5e) ",creal(sky->vkgi[i]),cimag(sky->vkgi[i]));
    }
    printf("\n");
    printf("\n");
    printf("\n");
    for(int i=0;i<n;i++){
      for(int j=0;j<n;j++){
	printf("(%.3e,%.3e) ",
	       creal(GetSkyline(sky,i,j)),
	       cimag(GetSkyline(sky,i,j)));
      }   
      printf("\n");
    }
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if (cabs(GetSkyline(sky,i,j)) < 1e-10) {
	printf("  ");
      } else {
	printf("><");
      }
    }   
    printf("\n");
  }

}



void FactoLU(Skyline* sky){

  dcmplx* vfg=NULL;
  dcmplx* vu=NULL;
  dcmplx energ;
  int ier;
  int ifac=1;
  int isol=0;
  int nsym=1;
  if (sky->is_sym) nsym=0;
  // Allocating and storing old matrix inside copies.
  //if (!sky->copy_is_alloc) AllocateCopySkyline(sky);
  //printf("LU factorization in progress...\n");
  //printf("address sky=%p\n",sky);

  sol_(sky->vkgs,sky->vkgd, sky->vkgi,
       vfg, sky->kld, vu, sky->neq, 
       ifac, isol, nsym,
       &energ, &ier);

  sky->is_lu=true;

}

void MatVectSkyline(Skyline * sky, dcmplx * x, dcmplx * prod) {

  assert(!sky->is_lu);

  int nsym=1;
  if (sky->is_sym) nsym=0;

  for(int i=0; i < sky->neq; i++) prod[i]=0;
  if (!sky->is_lu){
    mulku_(sky->vkgs, sky->vkgd, sky->vkgi,
	   sky->kld, x, sky->neq, nsym, 
	   prod, sky->nmem);
  }
  else
    {
      mulku_(sky->copy_vkgs, sky->copy_vkgd, sky->copy_vkgi,
	     sky->kld, x, sky->neq, nsym, 
	     prod, sky->nmem);
    }
}


// slow solve: better precision
void SolveSkyline(Skyline* sky,dcmplx* vfg,dcmplx* vu){
  assert(sky->is_lu);

  dcmplx energ;
  int ier,iter;
  int ifac=0;
  int isol=1;
  int nsym=1;
  dcmplx * vec_temp;
  dcmplx * sol_temp;
  dcmplx * sol_temp2;
  int nb_iterations=2;

  sol_temp=calloc(sky->neq,sizeof(dcmplx));
  sol_temp2=calloc(sky->neq,sizeof(dcmplx));  
  vec_temp=calloc(sky->neq,sizeof(dcmplx));
  
  if (sky->is_sym) nsym=0;

  sol_(sky->vkgs,sky->vkgd, sky->vkgi,
       vfg, sky->kld, sol_temp, sky->neq, 
       ifac, isol, nsym,&energ, &ier);

  /////////// Post treatment ////////////
  for(iter=0;iter<nb_iterations;iter++) {

    if(iter>0){
      for(int i=0; i < sky->neq; i++){
	sol_temp[i]=sol_temp2[i];
      }
    }
       
    MatVectSkyline(sky,sol_temp,vec_temp);


    for(int i=0; i < sky->neq; i++){
      vec_temp[i]=vec_temp[i]-vfg[i];
    }
      
    sol_(sky->vkgs,sky->vkgd, sky->vkgi,
	 vec_temp, sky->kld, sol_temp2, sky->neq, 
	 ifac, isol, nsym,&energ, &ier);
    
    double error=0.0;
    for(int i=0; i < sky->neq; i++){
      error=error+cabs(sol_temp2[i]);
    }
    printf("llllU %.13e %d \n",error,iter);
    
    for(int i=0; i < sky->neq; i++){
      sol_temp2[i]=sol_temp[i]-sol_temp2[i];
    }
   
    if(error <1.e-12 || iter==nb_iterations){
      break;
    }
  }

  if(iter>0){
    for(int i=0; i < sky->neq; i++){
      vu[i]=sol_temp2[i];
    }
  }
  else
    {
      for(int i=0; i < sky->neq; i++){
	vu[i]=sol_temp[i];
      }
    }
  free(sol_temp);
  free(sol_temp2); 
  free(vec_temp);


}

void FastSolveSkyline(Skyline* sky,dcmplx* vfg,dcmplx* vu){
  assert(sky->is_lu);

  dcmplx energ;
  int ier,iter;
  int ifac=0;
  int isol=1;
  int nsym=1;

  if (sky->is_sym) nsym=0;

  sol_(sky->vkgs,sky->vkgd, sky->vkgi,
       vfg, sky->kld, vu, sky->neq, 
       ifac, isol, nsym,&energ, &ier);


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

  if (sky->is_lu){
    free(sky->copy_vkgs);
    free(sky->copy_vkgd);
    if (! sky->is_sym)  free(sky->copy_vkgi);
  }
  sky->is_alloc=false;
  //InitSkyline(sky,sky->neq);

}

/* Table of constant values */

static int c__1 = 1;

/* Subroutine */ int sol_(dcmplx *vkgs, dcmplx *vkgd, dcmplx *
			  vkgi, dcmplx *vfg, int *kld, dcmplx *vu, int neq, 
			  int ifac, int isol, int nsym, dcmplx *
			  energ, int *ier)
{
  /* Initialized data */

  static dcmplx vzero = 0.0;

  /* Format strings */
  static char fmt_8000[] = "sol pivot nul equation";

  /* System generated locals */
  int i__1, i__2, i__3, i__4;

  /* Builtin functions */
  //int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);

  /* Local variables */
  static int i__;
  static dcmplx c1=0.0, c2=0.0;
  static int j1, j2, ic, ij, ik, jbk, jck, jhj, jhk, lhk, jhj1, jhk1, 
    lhk1;
  extern dcmplx scal_(dcmplx *, dcmplx *, int *);
  static int imin, imax, imin1;
  static dcmplx cdiag=0.0;

  /*   resolution d'un systeme lineaire symetrique ou non. la matrice est */
  /*   stockee par ligne de ciel,en memoire dans les tables vkgs,vkgd,vkgi */

  /*       entrees */
  /*          vkgs,vkgd,vkgi    matrice du systeme : parties superieure, */
  /*                            diagonale, inferieure (dcmplx precision) */
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

dcmplx scal_(dcmplx *x, dcmplx *y, int *n)
{
  /* Initialized data */

  static dcmplx zero = 0.0;

  /* System generated locals */
  int i__1;
  dcmplx ret_val=0.0;

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


/* muls.f -- translated by f2c (version 20100827).
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

/* Subroutine */ int mulku_(dcmplx *vkgs, dcmplx *vkgd, dcmplx *
			    vkgi, int *kld, dcmplx *vfg, int neq, int nsym, 
			    dcmplx *vres, int nsky)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  static dcmplx c__;
  static int j, i0, i1, ij, ik, jhk, lhk, jhk1;
  //extern dcmplx scal_(dcmplx *, dcmplx *, int *);

  /* =======================================================================MULK   2 */
  /*     CE SOUS-PROGRAMME AJOUTE AU VECTEUR RES LE PRODUIT DE LA          MULK   3 */
  /*     MATRICE KG PAR LE VECTEUR FG                                      MULK   4 */
  /*       ENTREES                                                         MULK   5 */
  /*          VKGS,VKGD,VKGI  MATRICE KG STOCKEE PAR LIGNE DE CIEL (SYM.   MULK   6 */
  /*                          OU NON SYM.)                                 MULK   7 */
  /*          KLD     TABLE DES POINTEURS DES HAUTS DE COLONNES DE KG      MULK   8 */
  /*          VFG     VECTEUR FG                                           MULK   9 */
  /*          NEQ     DIMENSION DES VECTEURS FG ET RES                     MULK  10 */
  /*          NSYM    .EQ.1 SI LE PROBLEME N'EST PAS SYMETRIQUE            MULK  11 */
  /*          VRES    VECTEUR RES                                          MULK  12 */
  /*       SORTIE                                                          MULK  13 */
  /*          VRES    VECTEUR RES                                          MULK  14 */
  /* =======================================================================MULK  15 */
  /* -----------------------------------------------------------------------MULK  18 */
  /* -------  POUR CHAQUE COLONNE DE LA MATRICE KG                          MULK  19 */
  /* Parameter adjustments */
  --vres;
  --vfg;
  --kld;
  --vkgd;
  --vkgi;
  --vkgs;

  /* Function Body */
  i__1 = neq;
  for (ik = 1; ik <= i__1; ++ik) {
    jhk = kld[ik]+_Z;
    jhk1 = kld[ik + 1]+_Z;
    lhk = jhk1 - jhk;
    /* -------  TERME DIAGONAL                                                MULK  24 */
    c__ = vkgd[ik] * vfg[ik];
    if (lhk <= 0) {
      goto L20;
    }
    i0 = ik - lhk;
    /* -------  TERMES DE LIGNE                                               MULK  28 */
    if (nsym != 1) {
      c__ += scal_(&vkgs[jhk], &vfg[i0], &lhk);
    }
    if (nsym == 1) {
      c__ += scal_(&vkgi[jhk], &vfg[i0], &lhk);
    }
    /* -------  TERMES DE COLONNE                                             MULK  31 */
    j = jhk;
    i1 = ik - 1;
    i__2 = i1;
    for (ij = i0; ij <= i__2; ++ij) {
      vres[ij] += vkgs[j] * vfg[ik];
      /* L10: */
      ++j; 
    }
  L20:
    vres[ik] += c__;
  }
  return 0;
} /* mulku_ */

void PLU_Square(int n, dcmplx *a, int *sigma){

  // init permutation to identity
  for(int p = 0; p < n; ++p) sigma[p] = p;
  
  // elimination in column p
  for(int p = 0; p < n - 1; ++p){

    // search max piv
    double sup = 0;
    int pmax = p;
    for(int i = p; i < n; i++) {
      if (sup <= cabs(a[i * n + p])){
	sup = cabs(a[i * n + p]);
	pmax = i;
      }	
    }  

    // swap two lines
    for(int j = 0; j < n; j++){
      dcmplx aux = a[p * n + j];
      a[p * n + j] = a[pmax * n + j];
      a[pmax * n + j] = aux;
    }

    // store permutation
    int temp = sigma[p];
    sigma[p] = sigma[pmax];
    sigma[pmax] = temp;

    for(int i = p + 1; i < n; i++){
      dcmplx c = a[i * n + p] / a[ p * n + p];
      for(int j = p + 1; j < n; j++)
	a[i * n + j]-= c * a[p * n + j];
      a[i * n + p] = c;
    }
  }

  //display
  /* for(int i = 0; i < n; ++i){ */
  /*   for(int j = 0; j < n; j++){ */
  /*     printf("(%f,%f) ", creal(a[i * n + j]), cimag(a[i * n + j])); */
  /*   } */
  /*   printf("\n"); */
  /* } */
  /* printf("permut="); */
  /* for(int i = 0; i < n; ++i) printf("%d ",sigma[i]); */
  /* printf("\n"); */
}

void PLU_Solve(int n, dcmplx *a, int *sigma, dcmplx *b,
	       dcmplx *x){

  x[0] = b[sigma[0]];
  for(int i = 1; i < n; ++i){
    x[i] = b[sigma[i]];
    for(int j = 0; j < i; j++){
      x[i] -= a[i * n + j] * x[j];
    }
  }

  x[n-1] /= a[n * n - 1];
  for(int i = n-2; i >= 0; --i){
    for(int j = i + 1; j < n; j++){
      x[i] -= a[i * n + j] * x[j];
    }
    x[i] /= a[i * n + i];
  }

  /* printf("sol="); */
  /* for(int i = 0; i < n; ++i) printf("(%f,%f) ",creal(x[i]),cimag(x[i])); */
  /* printf("\n"); */


}

void InvertSquare(int n, dcmplx *As, dcmplx *B){
  
  int sigma[n];

  dcmplx A[n*n];

  for(int k = 0; k < n * n; k++) A[k] = As[k];
  
  PLU_Square(n, A, sigma);
    
  
  dcmplx r[n], x[n];
  
  for(int j = 0; j < n; j++){
    
    for(int i = 0; i < n; i++){
      r[i] = (i==j);
    }
    

    PLU_Solve(n, A, sigma, r, x);
    
    for(int i = 0; i < n; i++) B[i * n + j] = x[i]; 
				 
    /* for(int i = 0; i < n; i++) printf("j=%d i=%d r=%f %f\n", */
    /* 				      j,i, */
    /* 				      creal(x[i]), */
    /* 				      cimag(x[i])); */


  }

  /* assert(1==2); */

  /* //display */
  /* for(int i = 0; i < n; ++i){ */
  /*   for(int j = 0; j < n; j++){ */
  /*     printf("(%f,%f) ", creal(B[i * n + j]), cimag(B[i * n + j])); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  
}

void TestPLU(void)
{
  dcmplx A[16] = {.2,-1,0,0,
		  -1,2,-1,0,
		  0,-1,2,-1,
		  0.1,0,-1,2};

  dcmplx As[16] = {.2,-1,0,0,
		  -1,2,-1,0,
		  0,-1,2,-1,
		  0.1,0,-1,2};

  int sigma[4] = {0,1,2,3};

  dcmplx b[4] = {-0.8,0,0,1};

  dcmplx x[4];



  //PLU_Square(4, A, sigma);
  //PLU_Solve(4, A, sigma, b, x);
  dcmplx B[16];
  InvertSquare(4, A, B);
  PLU_Square(4, As, sigma);
  PLU_Solve(4, As, sigma, b, x);
  int n = 4;
  double res = 0;
  for(int i = 0; i < n; i++){
    for(int k = 0; k < n; k++) x[i] -= B[i*n+k]*b[k];
    //printf("b=(%f,%f)\n",creal(x[i]),cimag(x[i]));
    res += cabs(x[i]);
  }
  
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      dcmplx v = 0;
      for(int k =  0; k < n; k++){
	v += A[n * i + k] * B[n * k + j];
      }
      printf("i=%d j=%d aij=(%f,%f)\n",i,j,
	     creal(v),cimag(v));
    }
  }


  printf("test inversion, erreur=%f\n",res);
  assert(res < 1e-14);
}


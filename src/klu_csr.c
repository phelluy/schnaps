#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "klu_csr.h"
#include "math.h"



void InitKLU(KLU * klu, int n)
{

  klu->is_alloc = false;
  klu->copy_is_alloc = false;
  klu->is_lu = false;

  klu->neq = n;

  klu_defaults(&klu->common);

  klu->T = cs_di_spalloc(0, 0, 1, 1, 1);

  klu->symbolic = NULL;

  klu->numeric = NULL;
  
  //klu->symb_Qinv = NULL;
  
  //klu->symb_Pinv = NULL;
  
  klu->display_orig = true;
  klu->display_orig_values = false;
  klu->display_klu_blocks = true;
}


void SwitchOnKLU(KLU * klu, int i, int j)
{

  assert(!klu->is_alloc);
  cs_di_entry(klu->T, i, j, 1);

}

void AllocateKLU(KLU * klu)
{

  assert(!klu->is_alloc);

  klu->A = cs_di_compress(klu->T);
  cs_di_spfree(klu->T);

  klu->is_alloc = true;

  cs_di_dupl(klu->A);
  // check that the matrix is square and its size conforms to the prescribed one
  assert(klu->A->m == klu->neq);
  assert(klu->A->n == klu->neq);
  // symbolic factorization
  klu->symbolic =
      klu_analyze(klu->neq, klu->A->p, klu->A->i, &(klu->common));
  assert(klu->symbolic);
  //
  printf(" KLU found %i blocks , %i nb non-zeroes in offdiagblocks , largest block size %i\n",
    klu->symbolic->nblocks,klu->symbolic->nzoff,klu->symbolic->maxblock);
  //
/*  klu->symb_Qinv = (int *) calloc(klu->neq,sizeof(int));*/
/*  assert(klu->symb_Qinv);*/
/*  klu->symb_Pinv = (int *) calloc(klu->neq,sizeof(int));*/
/*  assert(klu->symb_Pinv);*/
  //
/*  for (int i=0;i < klu->neq;i++){*/
/*    klu->symb_Qinv[klu->symbolic->Q[i]]=i;*/
/*  }*/
/*  for (int i=0;i < klu->neq;i++){*/
/*    klu->symb_Pinv[klu->symbolic->P[i]]=i;*/
/*  }*/
  //
}

void AllocateCopyKLU(KLU * klu)
{

  assert(!klu->copy_is_alloc);
  klu->Acopy = cs_di_add(klu->A, klu->A, 1, 0);

  klu->copy_is_alloc = true;
  //cs_di_print(klu->A, 0);


}

void AddKLU(KLU * klu, int i, int j, schnaps_real val)
{

  assert(klu->is_alloc);

  for (int iloc = klu->A->p[j]; iloc < klu->A->p[j + 1]; iloc++) {
    if (klu->A->i[iloc] == i) {
      klu->A->x[iloc] += val;
      return;
    }
  }
}

void SetKLU(KLU * klu, int i, int j, schnaps_real val)
{

  assert(klu->is_alloc);

  for (int iloc = klu->A->p[j]; iloc < klu->A->p[j + 1]; iloc++) {
    if (klu->A->i[iloc] == i) {
      klu->A->x[iloc] = val;
      return;
    }
  }
}

schnaps_real GetKLU(KLU * klu, int i, int j)
{

  assert(klu->is_alloc);
  for (int iloc = klu->A->p[j]; iloc < klu->A->p[j + 1]; iloc++) {
    if (klu->A->i[iloc] == i) {
      return klu->A->x[iloc];
    }
  }
  return 0.0;
}

void DisplayKLU(KLU * klu){
  int n = klu->neq;
  if (klu->display_orig){
    printf(" Original ordering\n");
    if (klu->display_orig_values){
      printf("\n");
      printf("\n");
      printf("\n");
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          printf("%.3e ", GetKLU(klu, i, j));
        }
        printf("\n");
      }
    }
    else{
      printf("\n");
      printf("\n");
      for (int i=0; i < n; i++){
        for (int j = 0; j < n; j++) {
          if (fabs(GetKLU(klu, i, j)) < 1e-8) {printf(" ");} else {printf("*");
          }
        }
        printf("\n");
      }
    }
  }
  if (klu->display_klu_blocks){
    printf(" KLU symbolic- pre-ordering\n");
    // Block structure
    int nb=klu->symbolic->nblocks;
    printf("KLU Block structure\n");
    for (int i=0;i<nb+1;i++){
      printf("\t%i",klu->symbolic->R[i]);
    }
    printf("\n");
    //
    for (int iblock = 0; iblock < nb; iblock++){
    int ibl_start = klu->symbolic->R[iblock];
    int ibl_end = klu->symbolic->R[iblock + 1];
    for (int i=ibl_start;i < ibl_end;i++){
      for (int jblock = 0; jblock < nb; jblock++){
        int jbl_start = klu->symbolic->R[jblock];
        int jbl_end = klu->symbolic->R[jblock + 1];
      for (int j=jbl_start;j < jbl_end;j++){
        if (fabs(GetKLU(klu, klu->symbolic->P[i], klu->symbolic->Q[j])) < 1e-8) {
          printf(" ");
        } else {
          printf("*");
        }
      } // j
      printf("|");
      } // jblock
      printf("\n");
      } // i
      for (int j=0; j < n;j++){
        printf("-");
      }
      printf("\n");
    } // iblock
  }
//
}



void FactoKLU(KLU * klu)
{


  klu->numeric =
      klu_factor(klu->A->p, klu->A->i, klu->A->x, klu->symbolic,
		 &(klu->common));
  klu->is_lu = true;
}

void ReFactoKLU(KLU * klu)
{

  assert(klu->is_lu);
  int res =
      klu_refactor(klu->A->p, klu->A->i, klu->A->x, klu->symbolic,
		   klu->numeric, &(klu->common));
  assert(res);

}

void MatVectKLU(KLU * klu, schnaps_real * x, schnaps_real * prod)
{

  for (int i = 0; i < klu->neq; i++) {
    prod[i] = 0.0;
  }
  for (int j = 0; j < klu->neq; j++) {
    for (int iloc = klu->A->p[j]; iloc < klu->A->p[j + 1]; iloc++) {
      prod[klu->A->i[iloc]] += klu->A->x[iloc] * x[j];
    }
  }
}


void SolveKLU(KLU * klu, schnaps_real * rhs, schnaps_real * sol)
{
  assert(klu->is_lu);
  // TODO? :evaluate  the need for an in place version avoiding the copy
  for (int i = 0; i < klu->neq; i++) {
    sol[i] = rhs[i];
  }
#ifndef _DOUBLE_PRECISION
  printf("warning: klu computations not available in single precision\n");
#endif
  klu_solve(klu->symbolic, klu->numeric, klu->neq, 1,
	    (double*) sol, &(klu->common));
}




void FreeKLU(KLU * klu)
{

  assert(klu->is_alloc);
  if (klu->symbolic) {
    klu_free_symbolic(&(klu->symbolic), &(klu->common));
    klu->symbolic = NULL;
  }
  if (klu->is_lu) {
    klu_free_numeric(&(klu->numeric), &(klu->common));
    klu->is_lu = false;
    klu->numeric = NULL;
  }
  //
  cs_di_spfree(klu->A);
  //
  if (klu->copy_is_alloc) {
    free(klu->Acopy);
  }
  //
/*  if (klu->symb_Pinv){*/
/*    free(klu->symb_Pinv);*/
/*    klu->symb_Pinv = NULL;*/
/*  };*/
/*  if (klu->symb_Qinv){*/
/*    free(klu->symb_Qinv);*/
/*    klu->symb_Qinv = NULL;*/
/*  }*/

}

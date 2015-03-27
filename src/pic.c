#include "pic.h"


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



#include "macromesh.h"

#include <stdio.h>
#include <assert.h>


void ReadMacroMesh(MacroMesh* m,char* filename){

  FILE* f=NULL;

  f=fopen(filename,"r");
  printf("%s\n",filename);

  assert(f != NULL);
  

}

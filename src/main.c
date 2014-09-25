#include "macromesh.h"
#include <stdio.h>
#include <assert.h>

int main(void) {

  MacroMesh m;

  ReadMacroMesh(&m,"../geo/disque.msh");
  PrintMacroMesh(&m);

  return 0;

} 

#ifndef _MACROMESH_H
#define _MACROMESH_H

#include "global.h"

// structure for managing the mesh obtained from
// gmsh. It is called a macromesh, because it will be refined
// afterward
typedef struct MacroMesh{
  //sizes
  int nbelems;
  int nbnodes;
  // connectivity
  int* elem2node;
  double* node;
} MacroMesh;

// get the mesh from a gmsh file
void ReadMacroMesh(MacroMesh* m,char* filename);
// list the mesh data
void PrintMacroMesh(MacroMesh* m);

#endif

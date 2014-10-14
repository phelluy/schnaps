#ifndef _MACROMESH_H
#define _MACROMESH_H

#include <stdbool.h>
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
  int* elem2elem;
  double* node;
} MacroMesh;

// a simple struct for modelling a four
// corner face with a left and a right element
typedef struct Face4Sort{
  int node[4];
  int left,right;
  int locfaceleft,locfaceright;
} Face4Sort;

// sort the node list of the face
void OrderFace4Sort(Face4Sort* f);

// compare two integers
int CompareInt(const void* a,const void* b);

// compare two ordered four-corner faces
int CompareFace4Sort(const void* a,const void* b);

// get the mesh from a gmsh file
void ReadMacroMesh(MacroMesh* m,char* filename);

// simple transformations of the mesh
void AffineMap(double* x);
void AffineMapMacroMesh(MacroMesh* m);

// detect if the mesh is 2D
// and then permut the nodes so that
// the z direction coincides in the reference
// or physical frame
bool Detect2DMacroMesh(MacroMesh* m);

// verify the validity and orientation of the mesh
void CheckMacroMesh(MacroMesh* m);
// list the mesh data
void PrintMacroMesh(MacroMesh* m);
void BuildConnectivity(MacroMesh* m);


#endif

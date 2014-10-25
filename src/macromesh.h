#ifndef _MACROMESH_H
#define _MACROMESH_H

#include <stdbool.h>
#include "global.h"

//! \brief structure for managing the mesh obtained from gmsh.
//! It is called a macromesh, because it will be refined
//! afterward
typedef struct MacroMesh{
  int nbelems; //!< number of macro elems
  int nbnodes; //!< number of nodes in the macromesh
  // connectivity
  int* elem2node; //!< elems to nodes connectivity (20 nodes/elem) 
  int* elem2elem; //!< elems to elems connectivity (along 6 faces)
  double* node; //!< nodes coordintes array
  bool is2d; //!< 2d computation detection
} MacroMesh;

//! \brief a simple struct for modelling a four
//! corner face with a left and a right element
//! used only by the connectivity builder
typedef struct Face4Sort{
  int node[4];
  int left,right;
  int locfaceleft,locfaceright;
} Face4Sort;

//!\brief  sort the node list of the face
//! used only by the connectivity builder
void OrderFace4Sort(Face4Sort* f);

//! \brief compare two integers
int CompareInt(const void* a,const void* b);

//! \brief compare two ordered four-corner faces
int CompareFace4Sort(const void* a,const void* b);

//! \brief get the mesh from a gmsh file
//! \param[inout] m pointer to a macromesh
//! \param[in] filename location of the gmsh file 
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
// for a given interpolation described in param[]
void CheckMacroMesh(MacroMesh* m,int param[7]);
// list the mesh data
void PrintMacroMesh(MacroMesh* m);
void BuildConnectivity(MacroMesh* m);


#endif

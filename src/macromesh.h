#ifndef _MACROMESH_H
#define _MACROMESH_H

#include <stdbool.h>
#include "global.h"

//! \brief structure for managing the mesh obtained from gmsh.
//! It is called a macromesh, because it will be refined
//! afterwards.
typedef struct MacroMesh{
  int nbelems; //!< number of macro elems
  int nbnodes; //!< number of nodes in the macromesh
  int nbfaces; //!< number of macrofaces
  // connectivity
  int* elem2node; //!< elems to nodes connectivity (20 nodes/elem)
  int* elem2elem; //!< elems to elems connectivity (along 6 faces)
  int* face2elem; //!< faces to elems connectivity (Left and Right)
  double* node; //!< nodes coordinates array
  bool is2d; //!< 2d computation detection
  bool periodic; //!< determines if the macromesh has periodic
		 //!boundary conditions
} MacroMesh;

//! \brief a simple struct for modelling a four
//! corner face with a left and a right element.
//! Used only by the connectivity builder.
typedef struct Face4Sort{
  //! 4 nodes
  int node[4];
  //! left elem index
  int left;
  //! right elem index
  //int right;  // FIXME: this is never referenced
  //! index of this face in the left elem
  int locfaceleft;
  //! index of this face in the right elem
  //int locfaceright;  // FIXME: this is never referenced
} Face4Sort;

//!\brief  sort the node list of the face.
//! Used only by the connectivity builder.
//! \param[in] f face whose nodes has to be sorted
void OrderFace4Sort(Face4Sort* f);

//! \brief compare two integers (used by quicksort).
//! \param[in] a first integer
//! \param[in] b second integer
//! \returns a value v, v<0 if a<b, v=0 if a==b, v>0 if a>b
int CompareInt(const void* a,const void* b);

//! \brief compare two ordered four-corner faces (used by quicksort).
//! Lexicographic order on the four nodes indices.
//! \param[in] a first face
//! \param[in] b second face
//! \returns a value v, v<0 if a<b, v=0 if a==b, v>0 if a>b
int CompareFace4Sort(const void* a,const void* b);

//! \brief get the mesh from a gmsh file.
//! \param[inout] m pointer to a macromesh
//! \param[in] filename location of the gmsh file
void ReadMacroMesh(MacroMesh* m,char* filename);

//! \brief compute additional connectivity arrays from
//! a basic connectivity given by gmsh.
//! \param[inout] m pointer to a macromesh
void BuildConnectivity(MacroMesh* m);

//! \brief affine transformation
//! \param[inout] x the transformed point
void AffineMap(double* x);
//! \brief simple transformations of the mesh
//! \param[inout] m the macromesh
void AffineMapMacroMesh(MacroMesh* m);

//! \brief detects if the mesh is 2D
//! and then permuts the nodes so that
//! the z direction coincides in the reference
//! or physical frame.
//! \param[inout] m a macromesh
//! \returns true if the mesh 2d
bool Detect2DMacroMesh(MacroMesh* m);

//! \brief verify the validity and orientation of the mesh
//! for  given interpolation parameters.
//! The function simply aborts if the mesh is bad because
//! going on with computations has no meaning.
//! \param[in] m a macromesh
//! \param[in] param interpolation parameters (m, degrees and refinements)
void CheckMacroMesh(MacroMesh* m,int param[7]);
//! \brief list the mesh data
//! \param[in] m a macromesh
void PrintMacroMesh(MacroMesh* m);


#endif

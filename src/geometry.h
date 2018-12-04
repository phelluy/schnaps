#ifndef _GEOMETRY_H
#define _GEOMETRY_H

// header of the geometry struct
// geometric transformations

#include "global.h"

//! \brief a struct for managing geometric mapping
typedef struct Geom{
  //! number of nodes in the reference element
  int nbrefnodes;
  /* //! nodes of the reference element */
  /* real physnode[_NB_REF_NODES][3]; */
  /* //! position of a point on the ref elem */
  /* real xref[3]; */
  /* //! gradient of a function on the ref elem */
  /* real dphiref[3]; */
  /* //! position of a point on the physical elem */
  /* real xphy[3]; */
  /* //! mapping jacobian, its comatrix and the determinant */
  /* real dtau[3][3]; */
  /* real codtau[3][3]; */
  /* real det; */
  /* //! gradient of a function on the physical elem */
  /* real dphi[3]; */
  /* //! normal vector in the physical elem on face ifa */
  /* int ifa; */
  /* real vnds[3]; */

  // some pointers to particular functions
  //! \brief pointer to particular Ref2Phy mapping
  void (*Ref2Phy)(schnaps_real physnode[][3],
                      schnaps_real xref[3],
                      schnaps_real dphiref[3],
                      int ifa,
                      schnaps_real xphy[3],
                      schnaps_real dtau[3][3],
                      schnaps_real codtau[3][3],
                      schnaps_real dphi[3],
                      schnaps_real vnds[3]);

  //! \brief pointer to a particular Phy2Ref inverse mapping
  void (*Phy2Ref)(schnaps_real physnode[][3],schnaps_real xphy[3],schnaps_real xref[3]);
} Geom;




//! \brief mapping tau from the reference point to the physical point.
//! If an optional variable is not used it HAS to be set to NULL
//! \param[in] physnode: coordinates of the physical nodes
//! \param[in] xref: coordinates of the mapped point in the reference frame
//! \param[in] dphiref: gradient of a function in the reference frame (optional)
//! \param[in] ifa: face index if computation of the normal vector (optional)
//! \param[out] xphy: coordinates of the mapped point in the physical frame
//! \param[out] dtau: jacobian of the mapping tau (optional)
//! \param[out] codtau: comatrix of dtau (optional)
//! \param[out] dphi: gradient of the function in the physical frame (optional)
//! \param[out] vnds: normal vector times the elementary surface ds (optional)
void schnaps_ref2phy(schnaps_real physnode[20][3],
             schnaps_real xref[3],
             schnaps_real dphiref[3],
             int ifa,
             schnaps_real xphy[3],
             schnaps_real dtau[3][3],
             schnaps_real codtau[3][3],
             schnaps_real dphi[3],
             schnaps_real vnds[3]);

//! \brief mapping tau from the reference point to the physical point.
//! Data encapsulation for more simplicity
//! \param[inout] g a Geom data structure
void GeomRef2Phy(Geom* g);


//! \brief inverse mapping tau from the physical point to the reference point.
//! solution by the Newton's method.
//! \param[in] physnode : coordinates of physical nodes
//! \param[in] xphy : coordinates of the mapped point in the physical frame
//! \param[out] xref: coordinates of the mapped point in the reference frame
void schnaps_phy2ref(schnaps_real physnode[20][3],schnaps_real xphy[3],schnaps_real xref[3]);

//! \brief inverse mapping tau from the physical point to the reference point.
//! more robust version of the function (Newton's method + homotopy)
//! \param[in] physnode : coordinates of physical nodes
//! \param[in] xphy : coordinates of the mapped point in the physical frame
//! \param[out] xref: coordinates of the mapped point in the reference frame
void RobustPhy2Ref(schnaps_real physnode[20][3],schnaps_real xphy[3],schnaps_real xref[3]);

//! \brief inverse mapping tau from the physical point to the reference point.
//! Function with encapsulation
//! \param[inout] g a Geom data structure
//void GeomPhy2Ref(Geom* g);

//! \brief dot product between two vectors
//! \param[in] a, b : the two points
//! \return the dot product
schnaps_real dot_product(schnaps_real a[3], schnaps_real b[3]);

//! \brief Length of a vector
//! \param[in] a: the vector
//! \return The length of the vector
schnaps_real norm(schnaps_real a[3]);

//! \brief normalize a vector
//! \param[inout] a: the vector
void Normalize(schnaps_real a[3]);

//! \brief periodic correction
//! \param[inout] xyz: the vector to be put inside the periodic box
//! \param[in] period:  sizes of the box in each direction
//! if period[dim]<0 -> non periodic in direction dim
//! the box is of the form [0,period[0]]x[0,period[1]]x[0,period[2]]
#pragma start_opencl
void PeriodicCorrection(schnaps_real xyz[3],schnaps_real period[3]);
#pragma end_opencl

//! \brief periodic correction
//! \param[inout] xyz the vector to be put inside the periodic box
//! \param[in] period  sizes of the box in each direction
//! \param[in] xmin  lower bounds of the box in each direction
//! \param[in] xmax  upper bounds of the box in each direction
//! if period[dim]<0 -> non periodic in direction dim
void PeriodicCorrectionB(schnaps_real xyz[3],schnaps_real period[3],schnaps_real xmin[3], schnaps_real xmax[3]);

//! \brief distance between two points
//! \param[in] a, b : the two points
//! \return the distance
schnaps_real Dist(schnaps_real a[3], schnaps_real b[3]);


//! \brief point coordinates on standard output
//! \param[in] x : the point
void PrintPoint(schnaps_real x[3]);

#endif

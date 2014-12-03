#ifndef _GEOMETRY_H
#define _GEOMETRY_H

// header of the geometry struct
// geometric transformations

//! \brief a struct for managing geometric mapping
typedef struct Geom{
  //! number of nodes in the reference element
#define _NB_REF_NODES 20
  //! nodes of the reference element
  double physnode[_NB_REF_NODES][3];
  //! position of a point on the ref elem
  double xref[3];
  //! gradient of a function on the ref elem
  double dphiref[3];
  //! position of a point on the physical elem
  double xphy[3];
  //! mapping jacobian, its comatrix and the determinant
  double dtau[3][3];
  double codtau[3][3];
  double det;
  //! gradient of a function on the physical elem
  double dphi[3];
  //! normal vector in the physical elem on face ifa
  int ifa;
  double vnds[3];

  // some pointers to particular functions
  //! \brief pointer to particular Ref2Phy mapping
  void (*Ref2Phy)(double physnode[_NB_REF_NODES][3],
                      double xref[3],
                      double dphiref[3],
                      int ifa,
                      double xphy[3],
                      double dtau[3][3],
                      double codtau[3][3],
                      double dphi[3],
                      double vnds[3]);

  //! \brief pointer to a particular Phy2Ref inverse mapping
  void (*Phy2Ref)(double physnode[_NB_REF_NODES][3],double xphy[3],double xref[3]);
} Geom;




//! \brief mapping tau from the reference point to the physical point.
//! If an optional variable is not used it HAS to be set to NULL
//! \param[in] physnode : coordinates of the physical nodes
//! \param[in] xref : coordinates of the mapped point in the reference frame
//! \param[in] dphiref : gradient of a function in the reference frame (optional)
//! \param[in] ifa : face index if computation of the normal vector (optional)
//! \param[out] xphy : coordinates of the mapped point in the physical frame
//! \param[out] dphi : gradient of the function in the physical frame (optional)
//! \param[out] dtau : jacobian of the mapping tau (optional)
//! \param[out] codtau : comatrix of dtau (optional)
//! \param[out] vnds : normal vector times the elementary surface ds (optional)
void Ref2Phy(double physnode[20][3],
             double xref[3],
             double dphiref[3],
             int ifa,
             double xphy[3],
             double dtau[3][3],
             double codtau[3][3],
             double dphi[3],
             double vnds[3]);

//! \brief mapping tau from the reference point to the physical point.
//! Data encapsulation for more simplicity
//! \param[inout] g a Geom data structure
void GeomRef2Phy(Geom* g);


//! \brief inverse mapping tau from the physical point to the reference point.
//!
//! \param[in] physnode : coordinates of physical nodes
//! \param[in] xphy : coordinates of the mapped point in the physical frame
//! \param[out] xref: coordinates of the mapped point in the reference frame
void Phy2Ref(double physnode[20][3],double xphy[3],double xref[3]);

//! \brief inverse mapping tau from the physical point to the reference point.
//! Function with encapsulation
//! \param[inout] g a Geom data structure
void GeomPhy2Ref(Geom* g);

//! \brief distance between two points
//! \param[in] x1,x2 : the two points
//! \return the distance
double Dist(double x1[3],double x2[3]);


//! \brief point coordinates on standard output
//! \param[in] x : the point
void PrintPoint(double x[3]);

#endif

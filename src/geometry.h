#ifndef _GEOMETRY_H
#define _GEOMETRY_H

// header of the geometry struct
// geometric transformations
//! \brief mapping tau from the reference point to the physical point.
//!
//! If an optional variable is not used it HAS to be set to NULL
//! \param[in] physnode[20*3] : coordinates of the physical nodes
//! \param[in] xref[3] : coordinates of the mapped point in the reference frame
//! \param[in] dphiref[3] : gradient of a function in the reference frame (optional)
//! \param[in] ifa=0..5 : face index if computation of the normal vector (optional)
//! \param[out] xphy[3]: coordinates of the mapped point in the physical frame
//! \param[out] dphi[3] : gradient of the function in the physical frame (optional)
//! \param[out] dtau[3*3] : jacobian of the mapping tau (optional)
//! \param[out] codtau[3*3] : comatrix of dtau (optional)
//! \param[out] vnds[3] : normal vector times the elementary surface ds (optional)
void Ref2Phy(double physnode[20][3],
             double xref[3],
             double dphiref[3],
             int ifa,
             double xphy[3],
             double dtau[3][3],
             double codtau[3][3],
             double dphi[3],
             double vnds[3]);

//! \brief inverse mapping tau from the physical point to the reference point.
//! \param[in] physnode[20*3] : coordinates of physical nodes
//! \param[in] xphy[3] : coordinates of the mapped point in the physical frame
//! \param[out] xref[3]: coordinates of the mapped point in the reference frame
void Phy2Ref(double physnode[20][3],double xphy[3],double xref[3]);

//! \brief distance between two points
//! \param[in] x1[3],x2[3] : the two points
//! \return the distance
double Dist(double x1[3],double x2[3]);


//! \brief point coordinates on standard output
//! \param[in] x[3] : the point
void PrintPoint(double x[3]);

#endif

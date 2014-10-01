#ifndef _GEOMETRY_H
#define _GEOMETRY_H

// header of the geometry struct
// geometric transformations

// mapping tau from the reference point to the physical point.
// If an optional variable is not used it HAS to be set to NULL
// input:
// physnode[20*3] : coordinates of physical nodes
// xref[3] : coordinates of the mapped point in the reference frame
// dphiref[3] : gradient of a function in the reference frame (optional)
// ifa=0..5 : face index if computation of the normal vector (optional)
// output:
// xphy[3]: coordinates of the mapped point in the physical frame
// dphi[3] : gradient of the function in the physical frame (optional)
// dtau[3*3] : jacobian of the mapping tau (optional)
// codtau[3*3] : comatrix of dtau (optional)
// vnds[3] : normal vector times the elementary surface ds (optional)
void Ref2Phy(double* physnode,
             double* xref,
             double* dphiref,
             int ifa,
             double* xphy,
             double* dtau,
             double* codtau,
             double* dphi,
             double* vnds);

// inverse mapping tau from the physical point to the reference point.
// input:
// physnode[20*3] : coordinates of physical nodes
// xphy[3] : coordinates of the mapped point in the physical frame
// output:
// xref[3]: coordinates of the mapped point in the reference frame
void Phy2Ref(double* physnode,double* xphy,double* xref);




#endif

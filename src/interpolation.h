#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H


// utilitie functions for computing interpolation
// on a macrocell
// in all these functions
// param is an input integer array
// param[0..2] : approximation degree in the 3 directions on 
// the reference cube
// param[3..5] : number of subcells in the 3 directions on 
// the reference cube
// param[6..nparam]: can be used for returning
// parameters (the user has to ensure that enough memory is reserved

//! \brief number of Gauss-LObatto Points (GLOPs) on the macro cell 
//! \param[in] param: the param list 
int NPG(int param[]);
// number of GLOPs on the face ifa of the macrocell
int NPGF(int* param,int ifa);

// return the reference coordinates xpg[3] and weight wpg of the GLOP ipg
void ref_pg_vol(int* param,int ipg,double* xpg,double* wpg);

// from a reference point find the nearest
// gauss point
int ref_ipg(int* param,double* xref);

// same function for the face 
// param[6] contains the volume GLOP index computed from face GLOP index.
void ref_pg_face(int* param,int ifa,int ipg,double* xpg,double* wpg);
// return the value psi  and the gradient dpsi[3] of the basis 
// function ib at point xref[3]. Warning: the value of the gradient is
// not reliable if xref is on the boundary 
//of a subcell (because the gradient is discontinuous)
void psi_ref(int* param, int ib, double* xref, double* psi, double* dpsiref);
// return the gradient dpsi[3] of the basis 
// function ib at GLOP ipg.
void grad_psi_pg(int* param,int ib,int ipg,double* dpsiref);







#endif

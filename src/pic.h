#ifndef _PIC_H
#define _PIC_H

#include "schnaps.h"
#include <math.h>

//! \brief struct for managing a set of particles
//! and the Particle-In-Cell (PIC) method
typedef struct PIC {

  //!  number of particles
  int nbparts;

  //! positions and velocity of particles (size=6*nbparts)
  real* xv;

  //! list of cell ids for each particles
  int* cell_id;
  int* old_cell_id;

  //! time step
  real dt; 

} PIC;


//! \brief init the PIC structure with 
//! space for n particles
//! \param[in] n number of particles
//! \param[inout] pic PIC object
void InitPIC(PIC* pic,int n);

//! \brief generate a 3d gaussian distribution of velocities
// ! Box-Muller algorithm
//! \param[in] k1 4 van der Corput parameters
//! \param[in] k2 4 van der Corput parameters
//! \param[out] v a pseudo-random gaussian vector 
void BoxMuller3d(real *v,int* k1, int* k2);



//! \brief free the allocated arrays
//! \param[inout] pic the PIC object
void FreePIC(PIC* pic);


//! \brief init the PIC structure with 
//! space for n particles
//! \param[in] n number of particles
//! \param[inout] pic PIC object
void CreateParticles(PIC* pic,MacroMesh *m);

//! \brief compute charge and current associated to particles
//! \param[in] pic a PIC struct containing the particles
//! \param[inout] f a maxwell field updated with charge and current sources
void AccumulateParticles(PIC* pic,field *f);


//! brief pseudo-random van der corput number generator
//! \param[in] n index of the number in the sequence
//! \param[in] k1 a prime number
//! \param[in] k2 a prime number k1 > k2 !!
real corput(int n,int k1,int k2);


//! brief create a gmsh file for plotting the particles
//! \param[in] pic a PIC struct
//! \param[in] m a macromesh
void PlotParticles(PIC* pic,MacroMesh *m);

//! brief push particles with a given field
//! \param[inout] pic a struct PIC describing the particles
//! \param[in] f a field
void PushParticles(field *f,PIC* pic);

#endif

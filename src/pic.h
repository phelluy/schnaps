#ifndef _PIC_H
#define _PIC_H

#include "schnaps.h"

//! \brief struct for managing a set of particles
//! and the Particle-In-Cell (PIC) method
typedef struct PIC {

  //!  number of particles
  int nbparts;

  //! positions and velocity of particles (size=6*nbparts)
  double* xv;

  //! list of cell ids for each particles
  int* cell_id;

} PIC;


//! \brief init the PIC structure with 
//! space for n particles
//! \param[in] n number of particles
//! \param[inout] pic PIC object
void InitPIC(PIC* pic,int n);


//! \brief free the allocated arrays
//! \param[inout] pic the PIC object
void FreePIC(PIC* pic);


//! \brief init the PIC structure with 
//! space for n particles
//! \param[in] n number of particles
//! \param[inout] pic PIC object
void CreateParticles(PIC* pic,Macromesh *m);


//! The prime number for the van der corput sequence
#define _K1 5   //!< k1 k2 pour van der corput dans la direction x
#define _K2 3

#define _K3 5    //!< k3 k4 pour van der corput dans la direction y
#define _K4 3

//! brief pseudo-random van der corput number generator
//! \param[in] n index of the number in the sequence
//! \param[in] k1 a prime number
//! \param[in] n index of the number in the sequence
double corput(int n,int k1,int k2);


#endif

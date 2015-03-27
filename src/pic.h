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



#endif

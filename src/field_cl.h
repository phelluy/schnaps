#ifndef _FIELD_CL_H
#define _FIELD_CL_H

#ifdef _WITH_OPENCL
#include "field.h"
#include "clinfo.h"

//! copy back the field to host memory
//! \param[inout] f a field
void CopyfieldtoCPU(field *f);

void dtfield_CL(field *f, cl_mem *dtwn_cl);
void DGVolume_CL(void *mcell, field *f, cl_mem *dtwn_cl);
void DGMacroCellInterface_CL(void *mface, field *f, cl_mem *wn_cl);
void DGMass_CL(void *mcell, field *f);

#endif // _WITH_OPENCL

#endif // _FIELD_CL_H

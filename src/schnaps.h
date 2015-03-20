//! \file schnaps.h Main header of the SCHNAPS library. Has to be included
//! by a program using the library.
#ifndef _SCHNAPS_H
#define _SCHNAPS_H

#include "global.h"
#include "geometry.h"
#include "field.h"

#ifdef _WITH_OPENCL
#include "field_cl.h"
#include "clinfo.h"
#endif

#include "maxwell.h"

#endif

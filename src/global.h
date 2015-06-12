#ifndef _GLOBAL_H
#define _GLOBAL_H
//! brief global variables and defs

#ifndef _OPENMP
// activate pthread if openmp is not here
//#define _WITH_PTHREAD
#endif

#ifdef _WITH_OPENCL
extern int nplatform_cl;
extern int ndevice_cl;
char numflux_cl_name[1024]; // FIXME: move to field struct.
char cl_buildoptions[1024];
#endif //_WITH_OPENCL

#define __constant
#define __local
#define __private

//#define _PERIOD 1

#endif // #ifndef _GLOBAL_H

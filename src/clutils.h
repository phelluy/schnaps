/*
  OpenCL Utility functions
*/

#ifndef CLUTILS_H_INC
#define CLUTILS_H_INC

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

const char* clErrorString(const cl_int err);

char* print_build_debug(cl_program* program, cl_device_id *device);
#endif /* CLUTILS_H_INC */


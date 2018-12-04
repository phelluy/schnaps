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

cl_ulong clv_starttime(cl_event clv);
cl_ulong clv_endtime(cl_event clv);
cl_ulong clv_duration(cl_event clv);

const char* clErrorString(const cl_int err);

char* print_build_debug(cl_program* program, cl_device_id *device);
#endif /* CLUTILS_H_INC */


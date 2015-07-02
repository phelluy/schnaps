#ifndef _CLINFO_H
#define _CLINFO_H
#include <stdbool.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>

//! \brief Data structure for managing the OpenCL
//! system informations
typedef struct CLInfo{
  cl_device_id device; //!< devices list
  cl_ulong devicememsize;  //!< GPU global memory size
  cl_ulong maxmembuffer;  //!< maximal size of a GPU buffer
  cl_ulong maxconstmem;  //!< maximal GPU constant memory
  cl_ulong cachesize;  //!< cache size of compute units

  cl_uint nbcomputeunits;  //!< number of compute units
  size_t maxworkgroupsize;  //!< maximal allowed size of a work group
  char platformname[1000]; //!< platform name
  char devicename[1000]; //!< accelerator name
  char clextensions[1000]; //!< list of OpenCL extensions
  cl_uint nbplatforms;  //!< number of platforms
  cl_uint platformnum; //!< platform number

  //cl_uint nbdevices;   //!< number of available devices
  cl_uint devicenum;   //!< chosen device number

  cl_context context;//!< OpenCL context
  cl_command_queue commandqueue;   //!< default command queue

  cl_program program;   //!< compute program

} CLInfo;


//! \brief Init the OpenCL framework
//! \param[inout] cli pointer to a CLInfo
//! \param[in] platform_id  Platform ID
//! \param[in] device_id  Compute device ID
void InitCLInfo(CLInfo *cli, int platform_id, int device_id);

//! \brief Display OpenCL informations
//! \param[in] cli pointer to a CLInfo
void PrintCLInfo(CLInfo *cli);

//! \brief Compile kernels source
//! \param[inout] cli pointer to a CLInfo
//! \param[in] program string containing the kernels sources
void BuildKernels(CLInfo *cli, char *program, char *buildoptions);

//! \brief scan all *.h and *.c in order to find the code
//! to be shared with opencl
//! such code is enclosed between \#pragma start_opencl and 
// #pragma end_opencl
//! get also field.cl
//! the result is put into schnaps.cl
void GetOpenCLCode(void);

//! \brief allocates and fills a string with a file content
void ReadFile(char filename[], char** s);

bool cldevice_is_acceptable(cl_uint ndevice, cl_uint nplatform);

double cl_dev_gflops(char *platform_name);
double cl_dev_bwidth(char *platform_name);
double kernel_min_time(double dev_flops, double bandwidth,
		       unsigned long int flop_count, 
		       unsigned long int io_count);
void print_kernel_perf(double dev_fflops, double bandwidth,
		       unsigned long int flop_count, unsigned long int io_count,
		       cl_ulong kernel_time_ns);

#endif

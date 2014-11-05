#ifndef _CLINFO_H
#define _CLINFO_H

#include <CL/cl.h>

//! \brief Data structure for managing the OpenCL
//! system informations
typedef struct CLInfo{

  cl_device_id* device; //!< devices list
  cl_ulong devicememsize;  //!< GPU global memory size
  cl_ulong maxmembuffer;  //!< maximal size of a GPU buffer
  cl_ulong cachesize;  //!< cache size of compute units

  cl_uint nbcomputeunits;  //!< number of compute units
  size_t  maxworkgroupsize;  //!< maximal allowed size of a work group
  char devicename[1000]; //!< accelerator name
  char clextensions[1000]; //!< list of OpenCL extensions
  cl_uint nbplatforms;  //!< number of platforms
  cl_uint platformid; //!< platform id

  cl_uint nbdevices;   //!< number of available devices
  cl_uint deviceid;   //!< chosen device ID 

  cl_context context;//!< OpenCL context

  cl_command_queue commandqueue;   //!< default command queue
  
  cl_program program;   //!< compute program

} CLInfo;


//! \brief Init the OpenCL framework
//! \param[inout] cli pointer to a CLInfo 
//! \param[in] platform_id  Platform ID
//! \param[in] device_id  Compute device ID
void InitCLInfo(CLInfo* cli,
                int platform_id,
                int device_id);

//! \brief Display OpenCL informations
//! \param[in] cli pointer to a CLInfo 
void PrintCLInfo(CLInfo* cli);

//! \brief Compile kernels source
//! \param[inout] cli pointer to a CLInfo 
//! \param[in] program string containing the kernels sources
void BuildKernels(CLInfo* cli,char* program);

#endif

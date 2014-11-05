#include "clinfo.h"
#include <stdio.h>
#include <assert.h>

void InitCLInfo(CLInfo* cli,int platform_id,int device_id){

  cl_int status;  // for checking OpenCL errors

  // numbers of platforms
  status = clGetPlatformIDs(0, NULL, &(cli->nbplatforms));
  assert (status == CL_SUCCESS);

  assert(cli->nbplatforms > 0);

  // platform array construction
  cl_platform_id* platforms =
    malloc(sizeof(cl_platform_id[cli->nbplatforms]));
  assert(platforms);

  status = clGetPlatformIDs(cli->nbplatforms, platforms, NULL);
  assert (status == CL_SUCCESS);

  printf("Available platforms:\n");

  char pbuf[2000];
  for (int i = 0; i < cli->nbplatforms; ++i) { 
    printf("Platform %d:\n",i);
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_NAME,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    assert (status == CL_SUCCESS);

    printf("%s\n",pbuf);

    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VENDOR,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    assert (status == CL_SUCCESS);

    printf("%s\n",pbuf);

    //  opencl version
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VERSION,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    assert (status == CL_SUCCESS);

    printf("%s\n",pbuf);

  }

  cli->platformid=platform_id;
  assert(cli->platformid < cli->nbplatforms);


  // devices count
  status = clGetDeviceIDs(platforms[platform_id],
			  CL_DEVICE_TYPE_ALL,
			  0,
			  NULL,
			  &(cli->nbdevices));
  assert (status == CL_SUCCESS);
  assert(cli->nbdevices > 0);

  // devices array construction
  cli->device = malloc(sizeof(cl_device_id) * cli->nbdevices);
  assert(cli->device);
  
  printf("DeviceID=%d NbDevices=%d\n",device_id,cli->nbdevices);
  assert(device_id < cli->nbdevices);

  printf("We choose device %d / %d ",device_id,cli->nbdevices);
  printf("of platform %d / %d\n",platform_id,cli->nbdevices);
  
  status = clGetDeviceIDs(platforms[platform_id],
			  CL_DEVICE_TYPE_ALL,
			  cli->nbdevices,
			  cli->device,
			  NULL);
  assert (status == CL_SUCCESS);

  // device name
  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_NAME,
			   sizeof(cli->devicename),
			   cli->devicename,
			   NULL);
  assert (status == CL_SUCCESS);
  printf("%s",cli->devicename);
  

  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_EXTENSIONS,
			   sizeof(cli->clextensions),
			   cli->clextensions,
			   NULL);
  assert (status == CL_SUCCESS);
  printf("CL extensions for that device:\n%s\n",cli->clextensions);


  // device memory
  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_GLOBAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->devicememsize),
			   NULL);
  assert (status == CL_SUCCESS);
  printf("Global memory: %f MB\n",cli->devicememsize/1024./1024.);
  
  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_MAX_MEM_ALLOC_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxmembuffer),
			   NULL);
  assert (status == CL_SUCCESS);
  printf("Max buffer size: %f MB\n",cli->maxmembuffer/1024./1024.);

  // compute unit size cache
  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_LOCAL_MEM_SIZE,
			   sizeof(cl_ulong),
                           &(cli->cachesize),
			   NULL);
  assert (status == CL_SUCCESS);  
  printf("Cache size: %f KB\n",cli->cachesize/1024.);


  // nb of compute units
  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_MAX_COMPUTE_UNITS,
                             sizeof(cl_uint),
                             (void*)&(cli->nbcomputeunits),
			   NULL);
  assert (status == CL_SUCCESS);
  printf("Nb of compute units: %d\n",cli->nbcomputeunits);
  


  // max workgroup size
  status = clGetDeviceInfo(
			   cli->device[device_id],
			   CL_DEVICE_MAX_WORK_GROUP_SIZE,
			   sizeof(size_t),
			   (void*)&(cli->maxworkgroupsize),
			   NULL);
  assert (status == CL_SUCCESS);
  printf("Max workgroup size: %zu\n",cli->maxworkgroupsize);

  // first opencl context
  cli->context = clCreateContext(NULL, // no context properties
    1,         // only one device in the list
    &cli->device[device_id], // device list
    NULL, // callback function 
    NULL, // function arguments
    &status);

  assert (status == CL_SUCCESS);

  // command queue
  cli->commandqueue = clCreateCommandQueue(
				      cli->context,
				      cli->device[device_id],
				      CL_QUEUE_PROFILING_ENABLE,
				      &status);
  assert (status == CL_SUCCESS);

  

}


void PrintCLInfo(CLInfo* cli){

   cl_int status;  // for checking OpenCL errors

   char pbuf[2000];

  printf("%s\n",cli->devicename);

  // device memory
  printf("Global memory: %f MB\n",cli->devicememsize/1024./1024.);
  printf("Max buffer size: %f MB\n",cli->maxmembuffer/1024./1024.);

  printf("Cache size: %f KB\n",cli->cachesize/1024.);


  printf("Nb of compute units: %d\n",cli->nbcomputeunits);
  


  printf("Max workgroup size: %zu\n",cli->maxworkgroupsize);


   printf("CL extensions:\n%s\n",cli->clextensions);



}


void BuildKernels(CLInfo* cli,char* strprog){

  cl_int err;

  // kernels creation
  cli->program = 
    clCreateProgramWithSource(
                              cli->context,
                              1,
                              (const char **) &strprog,
                              NULL,
                              &err);
  if (!(cli->program)) {
    printf("Failed to create program.\n");
  }

  // compilation
  err = clBuildProgram(cli->program, 0, NULL, NULL, NULL, NULL);
  // if not successfull: display the errors
  if (err != CL_SUCCESS) {
    size_t len;
    char buffer[204800];
    printf("Failed to build program.\n");
    clGetProgramBuildInfo(cli->program, cli->device[cli->deviceid], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    printf("%s\n",buffer);
  }
    assert( err == CL_SUCCESS);


}

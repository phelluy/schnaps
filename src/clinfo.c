#include "clinfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "clutils.h"

void InitCLInfo(CLInfo* cli, int platform_id, int device_id)
{
  cl_int status;
  
  /* numbers of platforms */
  status = clGetPlatformIDs(0, NULL, &(cli->nbplatforms));
  assert(status == CL_SUCCESS);
  
  assert(cli->nbplatforms > 0);
 
  /* platform array construction */
  cl_platform_id* platforms = malloc(sizeof(cl_platform_id[cli->nbplatforms]));
  assert(platforms);
  
  status = clGetPlatformIDs(cli->nbplatforms, platforms, NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  
  printf("\nAvailable OpenCL platforms:\n");
  
  char pbuf[2000];
  for(int i = 0; i < cli->nbplatforms; ++i) {
    printf("\nPlatform %d:\n", i);
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_NAME,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);
    printf("\t%s\n", pbuf);

    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VENDOR,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);
    printf("\t%s\n",pbuf);

    //  opencl version
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VERSION,
			       sizeof(cli->platformname),
			       cli->platformname,
			       NULL);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status == CL_SUCCESS);
    printf("\t%s\n",cli->platformname);
  }

  cli->platformid=platform_id;
  assert(cli->platformid < cli->nbplatforms);

  // devices count
  status = clGetDeviceIDs(platforms[platform_id],
			  CL_DEVICE_TYPE_ALL,
			  0,
			  NULL,
			  &(cli->nbdevices));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  assert(cli->nbdevices > 0);

  // devices array construction
  cli->device = malloc(sizeof(cl_device_id) * cli->nbdevices);
  assert(cli->device);

  printf("\nWe choose device %d/%d ", device_id, cli->nbdevices-1);
  assert(device_id < cli->nbdevices);

  printf("of platform %d/%d\n",platform_id,cli->nbplatforms-1);
  status = clGetDeviceIDs(platforms[platform_id],
			  CL_DEVICE_TYPE_ALL,
			  cli->nbdevices,
			  cli->device,
			  NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  // device name
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_NAME,
			   sizeof(cli->devicename),
			   cli->devicename,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("%s\n", cli->devicename);

  cl_device_type dtype;
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_TYPE,
			   sizeof(dtype), 
			   &dtype, 
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tDevice type: ");
  switch(dtype) {
  case CL_DEVICE_TYPE_CPU:
    printf("CL_DEVICE_TYPE_CPU (%d)\n", (int)dtype);
    break;
  case CL_DEVICE_TYPE_GPU:
    printf("CL_DEVICE_TYPE_GPU (%d)\n", (int)dtype);
    break;
  case CL_DEVICE_TYPE_ACCELERATOR:
    printf("CL_DEVICE_TYPE_ACCELERATOR (%d)\n", (int)dtype);
    break;
  case CL_DEVICE_TYPE_DEFAULT:
    printf("CL_DEVICE_TYPE_DEFAULT (%d)\n", (int)dtype);
    break;
  default:
    printf("ERROR: unknown OpenCL device type %d\n", (int)dtype);
  }

  // device memory
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_GLOBAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->devicememsize),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tGlobal memory: %f MB\n",cli->devicememsize/1024./1024.);

  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_MEM_ALLOC_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxmembuffer),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tMax buffer size: %f MB\n",cli->maxmembuffer/1024./1024.);

  // compute unit size cache
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_LOCAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->cachesize),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tCache size: %f KB\n",cli->cachesize/1024.);

  // get maxconstmem
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxconstmem),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tConst mem: %f KB\n",cli->maxconstmem/1024.);

  // get maxconst args
  int maxcstargs;
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_CONSTANT_ARGS,
			   sizeof(cl_ulong),
			   &maxcstargs,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tMax Const args: %d \n",maxcstargs);

  // nb of compute units
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_COMPUTE_UNITS,
			   sizeof(cl_uint),
			   (void*)&(cli->nbcomputeunits),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tNb of compute units: %d\n",cli->nbcomputeunits);

  // max workgroup size
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_WORK_GROUP_SIZE,
			   sizeof(size_t),
			   (void*)&(cli->maxworkgroupsize),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tMax workgroup size: %zu\n",cli->maxworkgroupsize);

  // OpenCL extensions
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_EXTENSIONS,
			   sizeof(cli->clextensions),
			   cli->clextensions,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  printf("\tOpenCL extensions: %s\n", cli->clextensions);

  // First opencl context
  cli->context = clCreateContext(NULL, // no context properties
				 1,         // only one device in the list
				 &cli->device[device_id], // device list
				 NULL, // callback function
				 NULL, // function arguments
				 &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

#if 0
  //#ifdef CL_VERSION_2_0
  cl_queue_properties queue_properties = CL_QUEUE_PROFILING_ENABLE;
  cli->commandqueue 
    = clCreateCommandQueueWithProperties(cli->context,
					 cli->device[device_id],
					 &queue_properties,
					 &status);
#else
    cli->commandqueue = clCreateCommandQueue(cli->context,
					   cli->device[device_id],
					   CL_QUEUE_PROFILING_ENABLE,
					   &status);
#endif
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);

  printf("\tOpenCL Init OK\n\n");
}

void PrintCLInfo(CLInfo *cli){
  cl_int status;  // for checking OpenCL errors
  char pbuf[2000];

  printf("%sPlatform: \n",cli->platformname);
  printf("\tDevice: %s\n",cli->devicename);

  // device memory
  printf("\tGlobal memory: %f MB\n", cli->devicememsize/1024./1024.);
  printf("\tMax buffer size: %f MB\n", cli->maxmembuffer/1024./1024.);
  printf("\tCache size: %f KB\n", cli->cachesize/1024.);
  printf("\tNb of compute units: %d\n", cli->nbcomputeunits);
  printf("\tMax workgroup size: %zu\n", cli->maxworkgroupsize);
  printf("\tOpenCL extensions:\n%s\n", cli->clextensions);
}

void BuildKernels(CLInfo *cli, char *strprog, char *buildoptions)
{
  cl_int status;

  cli->program = clCreateProgramWithSource(cli->context,
					   1,
					   (const char **) &strprog,
					   NULL,
					   &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status == CL_SUCCESS);
  
  if(!(cli->program)) 
    printf("Failed to create program.\n");

  status = clBuildProgram(cli->program,
			  0,               // one device
			  NULL,
			  buildoptions,//NULL, 
			  NULL, NULL);

  /* cl_int clBuildProgram(	cl_program program, */
  /* 				cl_uint num_devices, */
  /* 				const cl_device_id *device_list, */
  /* 				const char *options, */
  /* 				void(CL_CALLBACK *pfn_notify)(cl_program program, void *user_data), */
  /* 				void *user_data) */

  if(status != CL_SUCCESS) {
    /* printf("%s\n", strprog); */
    printf("%s\n", clErrorString(status));
    printf("Compilation output:\n%s\n", 
	   print_build_debug(&(cli->program), &cli->device[cli->deviceid]));
  }
  assert(status == CL_SUCCESS);
}

void ReadFile(char filename[], char **s){
  FILE *f = fopen(filename , "r");
  assert(f);

  fseek(f , 0L, SEEK_END);
  int size = ftell(f);
  rewind(f);

  /* allocate memory for entire content */
  *s = calloc(size + 2, sizeof(char));
  assert(*s);
  /* copy the file into the buffer */
  int status = fread(*s , sizeof(char), size , f);
  // assert(status==0);
  //assert(s[size]=='\0');
  fclose(f);
}

//! \brief scan all *.h and *.c in order to find the code
//! to be shared with opencl
//! such code is enclosed between #pragma start_opencl and 
//! #pragma end_opencl
void GetOpenCLCode(void){
  int status;
  status = system("sh get_opencl_code.sh");
  assert(!status);
}


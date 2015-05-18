#include "clinfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "clutils.h"
#include <stdbool.h>

bool cldevice_exists(cl_platform_id *platform, cl_uint ndevice)
{
  cl_int status;
  cl_uint ndevices;
  status = clGetDeviceIDs(*platform,
			  CL_DEVICE_TYPE_ALL,
			  0,
			  NULL,
			  &ndevices);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("ndevices:%d\n",ndevices);
  printf("ndevice:%d\n",ndevice);
  if(ndevice < ndevices)
    return true;
  return false;
}

cl_uint get_nbplatforms()
{
  cl_int status;
  cl_uint nbplatforms;
  status = clGetPlatformIDs(0, NULL, &(nbplatforms));
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  return nbplatforms;
}

bool clplatform_exists(cl_uint nplatform)
{ 
  cl_uint nbplatforms = get_nbplatforms();
  if(nplatform < nbplatforms)
    return true;
  else
    return false;
}

cl_platform_id get_platform_id(cl_uint nplatform)
{
  cl_uint nbplatforms = get_nbplatforms();
  cl_platform_id *platforms = malloc(sizeof(cl_platform_id[nbplatforms]));
  
  cl_int status;
  status = clGetPlatformIDs(nbplatforms, platforms, NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  cl_platform_id the_platform = platforms[nplatform];

  free(platforms);

  return the_platform;
}

cl_uint get_nbdevices(cl_platform_id *platform_id)
{
  cl_uint nbdevices;
  cl_int status;
  status = clGetDeviceIDs(*platform_id,
			  CL_DEVICE_TYPE_ALL,
			  0,
			  NULL,
			  &nbdevices);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  return nbdevices;
}

cl_device_id get_device_id(cl_platform_id platform, cl_uint ndevice)
{
  cl_uint nbdevices = get_nbdevices(&platform);
  cl_device_id *devices;
  devices = malloc(sizeof(cl_device_id) * nbdevices);

  cl_int status;
  status = clGetDeviceIDs(platform,
			  CL_DEVICE_TYPE_ALL,
			  nbdevices,
			  devices,
			  NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  cl_device_id device = devices[ndevice];

  free(devices);

  return device;
}

void get_cldevice_extensions(cl_device_id device, char * buf, size_t bufsize)
{
  cl_int status;
  status = clGetDeviceInfo(device,
			   CL_DEVICE_EXTENSIONS,
			   bufsize,
			   buf,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

bool cldevice_supports_double(cl_device_id *device)
{
  char clextensions[1000];
  get_cldevice_extensions(*device, clextensions, sizeof(clextensions));
  printf("\tOpenCL extensions: %s\n", clextensions);
  if(strstr(clextensions, "cl_khr_fp64") != NULL) {
    printf("\t\tDouble-precision enabled (cl_khr_fp64).\n");
    return true;
  } else {
    printf("\t\tDouble-precision not enabled (cl_khr_fp64): ABORT!\n");
    return false;
  }
}

bool cldevice_is_acceptable(cl_uint nplatform, cl_uint ndevice)
{
  if(!clplatform_exists(nplatform)) {
    printf("Invalid clplatform\n");
    return false;
  }

  cl_platform_id platform = get_platform_id(nplatform);
  if(!cldevice_exists(&platform, ndevice)) {
    printf("Invalid cldevice\n");
    return false;
  }

#if real == double
  cl_device_id device = get_device_id(platform, ndevice);
  if(!cldevice_supports_double(&device)) {
    printf("cldevice does not support double\n");
    return false;
  }
#endif

  return true;
}

void InitCLInfo(CLInfo *cli, int platform_id, int device_id)
{
  cl_int status;
  
  /* numbers of platforms */
  cli->nbplatforms = get_nbplatforms();
  assert(cli->nbplatforms > 0);
 
  /* platform array construction */
  cl_platform_id *platforms = malloc(sizeof(cl_platform_id[cli->nbplatforms]));
  assert(platforms);

  // Store the platform_ids in platforms
  status = clGetPlatformIDs(cli->nbplatforms, platforms, NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
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
    assert(status >= CL_SUCCESS);
    printf("\t%s\n", pbuf);

    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VENDOR,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    printf("\t%s\n",pbuf);

    //  opencl version
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VERSION,
			       sizeof(cli->platformname),
			       cli->platformname,
			       NULL);
    if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    printf("\t%s\n",cli->platformname);
  }

  cli->platformid = platform_id;
  assert(cli->platformid < cli->nbplatforms);

  cl_platform_id platform = get_platform_id(platform_id);

  // devices count
  cli->nbdevices = get_nbdevices(&platform);
  assert(cli->nbdevices > 0);

  // devices array construction
  cli->device = malloc(sizeof(cl_device_id) * cli->nbdevices);
  assert(cli->device);

  printf("\nWe choose device %d/%d ", device_id, cli->nbdevices-1);
  assert(device_id < cli->nbdevices);
  printf("of platform %d/%d\n",platform_id,cli->nbplatforms-1);

  // Get all of the device_ids and store them in cli->device 
  status = clGetDeviceIDs(platform,
			  CL_DEVICE_TYPE_ALL,
			  cli->nbdevices,
			  cli->device,
			  NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  // device name
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_NAME,
			   sizeof(cli->devicename),
			   cli->devicename,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("%s\n", cli->devicename);

  cl_device_type dtype;
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_TYPE,
			   sizeof(dtype), 
			   &dtype, 
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
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
  assert(status >= CL_SUCCESS);
  printf("\tGlobal memory: %f MB\n",cli->devicememsize/1024./1024.);

  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_MEM_ALLOC_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxmembuffer),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tMax buffer size: %f MB\n",cli->maxmembuffer/1024./1024.);

  // compute unit size cache
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_LOCAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->cachesize),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tLocal memory size: %f KB\n",cli->cachesize/1024.);

  // get maxconstmem
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxconstmem),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tConst mem: %f KB\n",cli->maxconstmem/1024.);

  // get global cache size
  cl_ulong global_mem_cache_size;
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,
			   sizeof(cl_ulong),
			   &global_mem_cache_size,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tGlobal memory cache size: %d bytes\n", (int)global_mem_cache_size);

  // get maxconst args
  int maxcstargs;
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_CONSTANT_ARGS,
			   sizeof(cl_ulong),
			   &maxcstargs,
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tMax Const args: %d \n",maxcstargs);

  // nb of compute units
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_COMPUTE_UNITS,
			   sizeof(cl_uint),
			   (void*)&(cli->nbcomputeunits),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tNb of compute units: %d\n",cli->nbcomputeunits);

  // max workgroup size
  status = clGetDeviceInfo(cli->device[device_id],
			   CL_DEVICE_MAX_WORK_GROUP_SIZE,
			   sizeof(size_t),
			   (void*)&(cli->maxworkgroupsize),
			   NULL);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tMax workgroup size: %zu\n",cli->maxworkgroupsize);

  // OpenCL extensions
  get_cldevice_extensions(cli->device[device_id], cli->clextensions,
			  sizeof(cli->clextensions));
  printf("\tOpenCL extensions: %s\n", cli->clextensions);
  
  // First opencl context
  cli->context = clCreateContext(NULL, // no context properties
				 1,         // only one device in the list
				 &cli->device[device_id], // device list
				 NULL, // callback function
				 NULL, // function arguments
				 &status);
  if(status != CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

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
  assert(status >= CL_SUCCESS);

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
  printf("\tLocal memory size: %f KB\n", cli->cachesize/1024.);
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
  assert(status >= CL_SUCCESS);
  
  if(!(cli->program)) 
    printf("Failed to create program.\n");

  int deflen = 20;
  char *buildoptions0 
    = (buildoptions == NULL) ? 
    malloc(deflen + 1) :  
    malloc(deflen + strlen(buildoptions) + 1);

#if real == double
  sprintf(buildoptions0, "-D real=double ");
#else
  sprintf(buildoptions0, "-D real=float ");
#endif
  if(buildoptions != NULL)
    strcat(buildoptions0, buildoptions);

  printf("OpenCL compilation arguments: %s\n", buildoptions0);
  status = clBuildProgram(cli->program,
			  0,               // one device
			  NULL,
			  buildoptions0,
			  NULL, NULL);

  free(buildoptions0);

  if(status < CL_SUCCESS) {
    //printf("%s\n", strprog); // Print the OpenCL code
    printf("%s\n", clErrorString(status));
    printf("Compilation output:\n%s\n",
  	   print_build_debug(&(cli->program), &cli->device[cli->deviceid]));
  }
  assert(status >= CL_SUCCESS);
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

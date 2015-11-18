#include <stdlib.h>
#include "clinfo.h"
#include <stdio.h>
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
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
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
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  cl_device_id device = devices[ndevice];
  free(devices);

#if 0
  // Split the device into subdevices.
  // NB: requires OpenCL >= 1.2
  unsigned int nsubdev = 8;
  cl_device_partition_property *partprop =
    malloc(sizeof(cl_device_partition_property) * 3);
  partprop[0] = CL_DEVICE_PARTITION_EQUALLY;
  //partprop[0] = CL_DEVICE_PARTITION_BY_COUNTS;
  partprop[1] = 1;
  partprop[2] = 0;
  cl_device_id *subdevice =  malloc(sizeof(cl_device_id) * nsubdev);
  status = clCreateSubDevices(device,
			      partprop,
			      nsubdev,
			      subdevice,
			      NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  return subdevice[0];
#else
  return device;
#endif
}

void get_cldevice_extensions(cl_device_id device, char * buf, size_t bufsize)
{
  cl_int status;
  status = clGetDeviceInfo(device,
			   CL_DEVICE_EXTENSIONS,
			   bufsize,
			   buf,
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
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

  if(sizeof(real) == sizeof(double)) {
      cl_device_id device = get_device_id(platform, ndevice);
      if(!cldevice_supports_double(&device)) {
	printf("cldevice does not support double\n");
	return false;
      }
    }

  return true;
}

void print_platforms(CLInfo *cli)
{
  cl_int status;
  
  int nplatforms = get_nbplatforms();
  
  /* platform array construction */
  cl_platform_id *platforms = malloc(nplatforms * sizeof(cl_platform_id));
    
  status = clGetPlatformIDs(cli->nbplatforms, platforms, NULL);
  
  char pbuf[2000];
  for(int i = 0; i < nplatforms; ++i) {
    printf("\nPlatform %d:\n", i);
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_NAME,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    printf("\t%s\n", pbuf);

    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VENDOR,
			       sizeof(pbuf),
			       pbuf,
			       NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    printf("\t%s\n",pbuf);

    //  opencl version
    status = clGetPlatformInfo(platforms[i],
			       CL_PLATFORM_VERSION,
			       sizeof(cli->platformname),
			       cli->platformname,
			       NULL);
    if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
    assert(status >= CL_SUCCESS);
    printf("\t%s\n",cli->platformname);
  }

  free(platforms);
}

void set_device_name(CLInfo *cli)
{
  cl_int status;
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_NAME,
			   sizeof(cli->devicename),
			   cli->devicename,
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void print_device_type(CLInfo *cli)
{
  cl_int status;
  
  cl_device_type dtype;
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_TYPE,
			   sizeof(dtype), 
			   &dtype, 
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
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
}

void print_opencl_version(CLInfo *cli)
{
  cl_int status;
  
  char buf[128];
  status = clGetDeviceInfo(cli->device, CL_DEVICE_VERSION,
			   128, buf, NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tOpenCL version: %s\n", buf);
}

void set_device_memory(CLInfo *cli)
{
  cl_int status;
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_GLOBAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->devicememsize),
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_max_buffer_size(CLInfo *cli)
{
  cl_int status;
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_MAX_MEM_ALLOC_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxmembuffer),
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_local_memory_size(CLInfo *cli)
{
  cl_int status;

  // compute unit size cache
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_LOCAL_MEM_SIZE,
			   sizeof(cl_ulong),
			   &(cli->cachesize),
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_max_constant_memory_size(CLInfo *cli)
{
  cl_int status;

  // get maxconstmem
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE,
			   sizeof(cl_ulong),
			   &(cli->maxconstmem),
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void print_global_mem_cache_size(CLInfo *cli)
{
  cl_int status;

  cl_ulong global_mem_cache_size;
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,
			   sizeof(cl_ulong),
			   &global_mem_cache_size,
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tGlobal memory cache size: %d bytes\n", (int)global_mem_cache_size);
}
  
void print_max_constant_arguments(CLInfo *cli)
{
  cl_int status;
  
  // get maxconst args
  int maxcstargs;
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_MAX_CONSTANT_ARGS,
			   sizeof(cl_ulong),
			   &maxcstargs,
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  printf("\tMax Const args: %d \n",maxcstargs);
}

void set_max_compute_units(CLInfo *cli)
{
  cl_int status;

  // nb of compute units
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_MAX_COMPUTE_UNITS,
			   sizeof(cl_uint),
			   (void*)&(cli->nbcomputeunits),
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_max_workgroup_size(CLInfo *cli)
{
  cl_int status;

  // max workgroup size
  status = clGetDeviceInfo(cli->device,
			   CL_DEVICE_MAX_WORK_GROUP_SIZE,
			   sizeof(size_t),
			   (void*)&(cli->maxworkgroupsize),
			   NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
}

void set_clinfo_data(CLInfo *cli)
{
  set_clinfo_data(cli);
  set_device_name(cli);
  set_device_memory(cli);
  set_max_buffer_size(cli);
  set_local_memory_size(cli);
  set_max_constant_memory_size(cli);
  set_max_compute_units(cli);
  set_max_workgroup_size(cli);
}

void PrintCLInfo(CLInfo *cli)
{
  cl_int status;  // for checking OpenCL errors

  printf("%sPlatform: \n",cli->platformname);
  printf("\tDevice: %s\n",cli->devicename);

  // device memory
  printf("\tGlobal memory: %f MB\n", cli->devicememsize/1024./1024.);
  printf("\tMax buffer size: %f MB\n", cli->maxmembuffer/1024./1024.);
  printf("\tLocal memory size: %f KB\n", cli->cachesize/1024.);
  printf("\tNb of compute units: %d\n", cli->nbcomputeunits);
  printf("\tMax workgroup size: %zu\n", cli->maxworkgroupsize);
  print_global_mem_cache_size(cli);
  print_max_constant_arguments(cli);
  print_device_type(cli);
  print_opencl_version(cli);
  get_cldevice_extensions(cli->device, cli->clextensions,
			  sizeof(cli->clextensions));
  printf("\tOpenCL extensions:\n%s\n", cli->clextensions);
}

void InitCLInfo(CLInfo *cli, int platform_num, int device_num)
{
  cl_int status;
  
  cli->nbplatforms = get_nbplatforms();
  cl_platform_id *platforms = malloc(cli->nbplatforms * sizeof(cl_platform_id));
 
  // Store the platform_ids in platforms
  status = clGetPlatformIDs(cli->nbplatforms, platforms, NULL);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  printf("\nAvailable OpenCL platforms:\n");
  print_platforms(cli);
  
  cli->platformnum = platform_num;
  cli->platform = get_platform_id(platform_num);

  int ndevices = get_nbdevices(&cli->platform);
  assert(ndevices > 0);

  printf("\nWe choose device %d/%d ", device_num, ndevices);
  assert(device_num < ndevices);
  printf("of platform %d/%d\n", platform_num, cli->nbplatforms - 1);

  cli->device = get_device_id(cli->platform, device_num);

  PrintCLInfo(cli);
  
  // First opencl context
  cli->context = clCreateContext(NULL, // no context properties
				 1,         // only one device in the list
				 &cli->device, // device list
				 NULL, // callback function
				 NULL, // function arguments
				 &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  cli->commandqueue 
    = clCreateCommandQueue(cli->context,
			   cli->device,
			   CL_QUEUE_PROFILING_ENABLE
			   | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE,
			   &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);

  printf("\tOpenCL Init OK\n\n");
}

void BuildKernels(CLInfo *cli, char *strprog, char *buildoptions)
{
  cl_int status;

  cli->program = clCreateProgramWithSource(cli->context,
					   1,
					   (const char **) &strprog,
					   NULL,
					   &status);
  if(status < CL_SUCCESS) printf("%s\n", clErrorString(status));
  assert(status >= CL_SUCCESS);
  
  if(!(cli->program)) 
    printf("Failed to create program.\n");

  int deflen = 20;
  char *buildoptions0
    = (buildoptions == NULL)
    ?  malloc(deflen + 1)
    : malloc(deflen + strlen(buildoptions) + 1);
  
  if(sizeof(real) == sizeof(double))
    sprintf(buildoptions0, "-D real=double ");
  else
    sprintf(buildoptions0, "-D real=float ");

  if(buildoptions != NULL) {
    strcat(buildoptions0, buildoptions);
  }
    
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
    printf("Compilation output:\n");
    print_build_debug(&cli->program, &cli->device);
  } else {
    printf("Compilation output:\n");
    print_build_debug(&cli->program, &cli->device);
  }
    
  assert(status >= CL_SUCCESS);
}

void ReadFile(char filename[], char **s) {
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
//! such code is enclosed between \#pragma start_opencl and 
// \#pragma end_opencl
void GetOpenCLCode(void){
  int status;
  status = system("sh ../get_opencl_code.sh");
  assert(status == 0);
}

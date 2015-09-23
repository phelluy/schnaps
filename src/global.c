#include "global.h"

// Set default platform choices via cmake, for example:
// cmake -D_CL_PLATFORM=1 -D_CL_DEVICE=1 ../
#ifndef _CL_PLATFORM
int nplatform_cl = 0;
#else
int nplatform_cl = _CL_PLATFORM;
#endif

#ifndef _CL_DEVICE
int ndevice_cl = 0;
#else
int ndevice_cl = _CL_DEVICE;
#endif

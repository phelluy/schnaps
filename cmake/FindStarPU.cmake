# Find StarPU
#
# Once done this will define
# STARPU_FOUND - system has StarPU
# STARPU_INCLUDE_DIRS - the StarPU include directory
# STARPU_LIB - link these to use StarPU
# HWLOC_LIB - link to hwloc libs


FIND_PACKAGE(PackageHandleStandardArgs)
find_path(STARPU_INCLUDE_DIR
  NAMES starpu.h
  HINTS /usr/include/starpu/1.1/ )
find_library(STARPU_LIB
  NAMES starpu-1.1)
find_library(HWLOC_LIB
  NAMES hwloc)
set(STARPU_INCLUDE_DIR ${STARPU_INCLUDE_DIR})
set(STARPU_LIB ${STARPU_LIB} ${HWLOC_LIB})
message(${STARPU_LIB})
message(${STARPU_INCLUDE_DIR})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(StarPU DEFAULT_MSG
  STARPU_LIB
  STARPU_INCLUDE_DIR)
mark_as_advanced(
  STARPU_LIB
  STARPU_INCLUDE_DIR
  )

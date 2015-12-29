set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required(VERSION 2.8)

#set(CMAKE_VERBOSE_MAKEFILE ON)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g -D_POSIX_C_SOURCE=200112L")
set(CMAKE_C_FLAGS "-std=gnu99 -D_REENTRANT -O3 -g")
#set(CMAKE_C_FLAGS " -std=c99 -D_REENTRANT -O3")


#set(CMAKE_C_FLAGS " -std=c99 -D_REENTRANT -O3")

IF(DEFINED _CL_PLATFORM)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_CL_PLATFORM=${_CL_PLATFORM}")
ENDIF()
IF(DEFINED _CL_DEVICE)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_CL_DEVICE=${_CL_DEVICE}")
ENDIF()


OPTION(SINGLE_PRECISION
  "Compile using floats"
  OFF)

IF(SINGLE_PRECISION)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dreal=float")
ELSE()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dreal=double")
ENDIF()

OPTION(USE_OPENCL
  "Compile with the OpenCL library"
  ON)

IF(USE_OPENCL)# set with: cmake . -DUSE_OPENCL:BOOL=ON
  MESSAGE(STATUS "Compiling with OpenCL Library")
  find_package (OpenCL)
  if(OPENCL_FOUND)
    MESSAGE(STATUS "OpenCL found")
  else()
    MESSAGE(STATUS "OpenCL not found")
  endif(OPENCL_FOUND)
ENDIF(USE_OPENCL)

if(OPENCL_FOUND)
  INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_OPENCL")
else()
  INCLUDE_DIRECTORIES(".")
endif(OPENCL_FOUND)

if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wno-deprecated-declarations")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

SET(BUILD_PARALUTION      OFF CACHE BOOL   "Build paralution library")

IF(BUILD_PARALUTION)

  ADD_DEFINITIONS(-DPARALUTION)
  INCLUDE(ExternalProject)
  SET(PARALUTION_VERSION 1.0.0 CACHE STRING "Paralution version number")
  EXTERNALPROJECT_ADD( paralution-${PARALUTION_VERSION}
     URL  http://www.paralution.com/downloads/paralution-${PARALUTION_VERSION}.tar.gz
     CMAKE_ARGS -DSUPPORT_CUDA=OFF -DSUPPORT_OMP=OFF
   -DBUILD_EXAMPLES=OFF -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -O3
     INSTALL_COMMAND ""
  )

  INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution-${PARALUTION_VERSION}-prefix/src/paralution-${PARALUTION_VERSION}/src/)
  LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/paralution-${PARALUTION_VERSION}-prefix/src/paralution-${PARALUTION_VERSION}-build/lib)
  ADD_LIBRARY(paralution_c   src/paralution_c.cpp )
  ADD_DEPENDENCIES(paralution_c  paralution-${PARALUTION_VERSION})
  TARGET_LINK_LIBRARIES(paralution_c paralution)

ENDIF()


include_directories(src)

enable_testing()

list(APPEND
  source_files
  src/schnaps.h
  src/field.c
  src/field.h
  src/geometry.c
  src/geometry.h
  src/global.c
  src/global.h
  src/model.c
  src/model.h
  src/interpolation.c
  src/interpolation.h
  src/macromesh.c
  src/macromesh.h
  src/h20phi.h
  src/maxwell.h
  src/maxwell.c
  src/mhd.h
  src/mhd.c
  src/gyro.h
  src/gyro.c
  src/pic.h
  src/pic.c
  src/skyline.h
  src/dpackfgmres.h
  src/dpackfgmres.c
  src/linear_solver.h
  src/linear_solver.c
  src/solverpoisson.h
  src/quantities_vp.h
  src/diagnostics_vp.h
  src/skyline.c
  src/solverpoisson.c
  src/quantities_vp.c
  src/diagnostics_vp.c
  src/collision.h
  src/collision.c
  )


if (OPENCL_FOUND)
  list(APPEND
    source_files
    src/field_cl.c
    src/field_cl.h
    src/clutils.c
    src/clutils.h
    src/clinfo.c
    src/clinfo.h
    )
endif (OPENCL_FOUND)

#schnaps library
add_library (_schnaps ${source_files})

IF(BUILD_PARALUTION)
  target_link_libraries(_schnaps paralution_c paralution)
ENDIF(BUILD_PARALUTION)

set(LIBS "")

if (OPENCL_FOUND)
  list(APPEND LIBS _schnaps pthread ${OPENCL_LIBRARIES} m)
else()
  list(APPEND LIBS _schnaps pthread m)
endif (OPENCL_FOUND)

# nearest neighbour fast search library
#find_package(flann)
if(FLANN_FOUND)
  MESSAGE(STATUS "Delicious!")
  list(APPEND LIBS ${FLANN_LIBRARY})
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_WITH_FLANN")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -D_WITH_FLANN")
endif()

# main program schnaps
add_executable(
  schnaps
  src/schnaps.c
  )
target_link_libraries(schnaps ${LIBS})
add_test(schnaps schnaps)


# Find all of the files ending in .c in the test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "test/")
file(GLOB globtest test/*.c)
MESSAGE(STATUS "globtest")
foreach(test IN LISTS globtest)
  string(REPLACE  ".c"  ""  test  ${test})
  string(REGEX REPLACE  "([^/]+$)"  ""  dir  ${test})
  string(REPLACE  ${dir}  ""  test  ${test})
  #MESSAGE(STATUS ${test})
  list(APPEND testlist ${test})
endforeach(test)

foreach(test IN LISTS testlist)
  MESSAGE(STATUS "Adding test " ${test})
  add_executable(
    ${test}
    test/test.h
    test/${test}.c
    )
  target_link_libraries(${test} ${LIBS})
  add_test(test/${test} test/${test})
endforeach(test)


# Find all of the files ending in .c in the example directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "example/")
file(GLOB globexample example/*.c)
MESSAGE(STATUS "globexample")
foreach(example IN LISTS globexample)
  string(REPLACE  ".c"  ""  example  ${example})
  string(REGEX REPLACE  "([^/]+$)"  ""  dir  ${example})
  string(REPLACE  ${dir}  ""  example  ${example})
  #MESSAGE(STATUS ${example})
  list(APPEND examplelist ${example})
endforeach(example)

foreach(example IN LISTS examplelist)
  MESSAGE(STATUS "Adding example " ${example})
  add_executable(
    ${example}
    test/test.h
    example/${example}.c
    )
  target_link_libraries(${example} ${LIBS})
endforeach(example)

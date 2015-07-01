set(CMAKE_LEGACY_CYGWIN_WIN32 0) # Remove when CMake >= 2.8.4 is required
cmake_minimum_required(VERSION 2.8)

project(schnaps)

#set(CMAKE_C_COMPILER "tau_cc.sh")
#set(CMAKE_C_FLAGS "-optRevert -optTauSelectFile=select.tau -std=c99 -fopenmp -Ofast")
#set(CMAKE_C_COMPILER "/opt/local/bin/gcc-mp-4.8")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#set(CMAKE_C_FLAGS "-std=c99 -D_REENTRANT -O3 -g -D_POSIX_C_SOURCE=200112L")
#set(CMAKE_C_FLAGS "-std=gnu99 -D_REENTRANT -O3 -g")
#set(CMAKE_C_FLAGS " -std=c99 -D_REENTRANT -O3 -g")
set(CMAKE_C_FLAGS " -std=c99 -D_REENTRANT -g")


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Dreal=double")

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
find_package(flann)
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

# list of tests
add_executable(
  testmacromesh
  test/test.h
  test/testmacromesh.c
  )
target_link_libraries(testmacromesh ${LIBS})
add_test(testmacromesh testmacromesh)

add_executable(
  testgeometry
  test/test.h
  test/testgeometry.c
  )
target_link_libraries(testgeometry ${LIBS})
add_test(testgeometry testgeometry)

add_executable(
  testinterpolation
  test/test.h
  test/testinterpolation.c
  )
target_link_libraries(testinterpolation ${LIBS})
add_test(testinterpolation testinterpolation)

add_executable(
  testmodel
  test/test.h
  test/testmodel.c
  )
target_link_libraries(testmodel ${LIBS})

add_test(testfield testfield)
add_executable(
  testfield
  test/test.h
  test/testfield.c
  )
target_link_libraries(testfield ${LIBS})
add_test(testfield testfield)

add_test(testfielddg testfielddg)
add_executable(
  testfielddg
  test/test.h
  test/testfielddg.c
  )
target_link_libraries(testfielddg ${LIBS})
add_test(testfielddg testfielddg)

add_test(testfieldrk2 testfieldrk2)
add_executable(
  testfieldrk2
  test/test.h
  test/testfieldrk2.c
  )
target_link_libraries(testfieldrk2 ${LIBS})
add_test(testfieldrk2 testfieldrk2)

add_test(testfieldrk4 testfieldrk4)
add_executable(
  testfieldrk4
  test/test.h
  test/testfieldrk4.c
  )
target_link_libraries(testfieldrk4 ${LIBS})
add_test(testfieldrk4 testfieldrk4)

add_test(test2dmeshdetection test2dmeshdetection)
add_executable(
  test2dmeshdetection
  test/test.h
  test/test2dmeshdetection.c
  )
target_link_libraries(test2dmeshdetection ${LIBS})
add_test(test2dmeshdetection test2dmeshdetection)

add_test(testfieldrk2-2d testfieldrk2-2d)
add_executable(
  testfieldrk2-2d
  test/test.h
  test/testfieldrk2-2d.c
  )
target_link_libraries(testfieldrk2-2d ${LIBS})
add_test(testfieldrk2-2d testfieldrk2-2d)

add_test(testsubcelldgvol testsubcelldgvol)
add_executable(
  testsubcelldgvol
  test/test.h
  test/testsubcelldgvol.c
  )
target_link_libraries(testsubcelldgvol ${LIBS})
add_test(testsubcelldgvol testsubcelldgvol)

add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)

add_executable(
  testfieldrk2-2d-subcell
  test/test.h
  test/testfieldrk2-2d-subcell.c
  )
target_link_libraries(testfieldrk2-2d-subcell ${LIBS})
add_test(testfieldrk2-2d-subcell testfieldrk2-2d-subcell)

add_test(testmeq2 testmeq2)
add_executable(
  testmeq2
  test/test.h
  test/testmeq2.c
  )
target_link_libraries(testmeq2 ${LIBS})
add_test(testmeq2 testmeq2)

add_test(testmaxwell2d testmaxwell2d)
add_executable(
  testmaxwell2d
  test/test.h
  test/testmaxwell2d.c
  )
target_link_libraries(testmaxwell2d ${LIBS})
add_test(testmaxwell2d testmaxwell2d)

add_executable(
  testgyro
  test/test.h
  test/testgyro.c
  )
target_link_libraries(testgyro ${LIBS})
add_test(testgyro testgyro)

add_executable(
  testpic
  test/test.h
  test/testpic.c
  )
target_link_libraries(testpic ${LIBS})
add_test(testpic testpic)

add_executable(
  testpic_accumulate
  test/test.h
  test/testpic_accumulate.c
  )
target_link_libraries(testpic_accumulate ${LIBS})
add_test(testpic_accumulate testpic_accumulate)

add_executable(
  testperiodic
  test/test.h
  test/testperiodic.c
)
target_link_libraries(testperiodic ${LIBS})
# This test is too slow.  Please re-enable after decreasing execution time.
#add_test(testperiodic testperiodic)


add_test(testskyline testskyline)
add_executable(
  testskyline
  test/test.h
  test/testskyline.c
  )
target_link_libraries(testskyline ${LIBS})
add_test(testskyline testskyline)


add_executable(
  testpoisson
  test/test.h
  test/testpoisson.c
)
target_link_libraries(testpoisson ${LIBS})
add_test(testpoisson testpoisson)

add_executable(
  testpoisson2d
  test/test.h
  test/testpoisson2d.c
)
target_link_libraries(testpoisson2d ${LIBS})
add_test(testpoisson2d testpoisson2d)

add_executable(
  test_transportvp
  test/test.h
  test/test_transportvp.c
)
target_link_libraries(test_transportvp ${LIBS})
add_test(test_transportvp test_transportvp)

add_executable(
  testlandau_damping_1D
  test/test.h
  test/testlandau_damping_1D.c
)
target_link_libraries(testlandau_damping_1D ${LIBS})
# This test is too slow.  Please re-enable after decreasing execution time.
#add_test(testlandau_damping_1D testlandau_damping_1D)


 add_executable(
    testwave2d
    test/test.h
    test/testwave2d.c
    )
  target_link_libraries(testwave2d ${LIBS})

##########  OpenCL unit tests  ##########
if(OPENCL_FOUND)
  add_test(testclinfo testclinfo)
  add_executable(
    testclinfo
    test/test.h
    test/testclinfo.c
    )
  target_link_libraries(testclinfo ${LIBS})

  add_test(testmanyv testmanyv)
  add_executable(
    testmanyv
    test/test.h
    test/testmanyv.c
    )
  target_link_libraries(testmanyv ${LIBS})
  add_test(testmanyv testmanyv)

  add_test(testdtfield_cl testdtfield_cl)
  add_executable(
    testdtfield_cl
    test/test.h
    test/testdtfield_cl.c
    )
  target_link_libraries(testdtfield_cl ${LIBS})

  add_test(testcache_cl testcache_cl)
  add_executable(
    testcache_cl
    test/test.h
    test/testcache_cl.c
    )
  target_link_libraries(testcache_cl ${LIBS})

  add_test(testmacrocellinterface testmacrocellinterface)
  add_executable(
    testmacrocellinterface
    test/test.h
    test/testmacrocellinterface.c
    )
  target_link_libraries(testmacrocellinterface ${LIBS})

  add_test(testfieldrk4_cl testfieldrk4_cl)
  add_executable(
    testfieldrk4_cl
    test/test.h
    test/testfieldrk4_cl.c
    )
  target_link_libraries(testfieldrk4_cl ${LIBS})

  add_test(testkernelinterface testkernelinterface)
  add_executable(
    testkernelinterface
    test/test.h
    test/testkernelinterface.c
    )
  target_link_libraries(testkernelinterface ${LIBS})

  add_test(testfieldrk2_cl testfieldrk2_cl)
  add_executable(
    testfieldrk2_cl
    test/test.h
    test/testfieldrk2_cl.c
    )
  target_link_libraries(testfieldrk2_cl ${LIBS})

  add_test(testkernel testkernel)
  add_executable(
    testkernel
    test/test.h
    test/testkernel.c
    )
  target_link_libraries(testkernel ${LIBS})

  add_test(testkernelvolume testkernelvolume)
  add_executable(
    testkernelvolume
    test/test.h
    test/testkernelvolume.c
    )
  target_link_libraries(testkernelvolume ${LIBS})

  add_test(testkernelflux testkernelflux)
  add_executable(
    testkernelflux
    test/test.h
    test/testkernelflux.c
    )
  target_link_libraries(testkernelflux ${LIBS})

  add_test(testmeq2_cl testmeq2_cl)
  add_executable(
    testmeq2_cl
    test/test.h
    test/testmeq2_cl.c
    )
  target_link_libraries(testmeq2_cl ${LIBS})


  add_test(testmhd testmhd)
  add_executable(
    testmhd
    test/test.h
    test/testmhd.c
    )
  target_link_libraries(testmhd ${LIBS})
  
  add_test(testmhd1d testmhd1d)
  add_executable(
    testmhd1d
    test/test.h
    test/testmhd1d.c
    )
  target_link_libraries(testmhd1d ${LIBS})

endif(OPENCL_FOUND)


IF(BUILD_PARALUTION)


  SET(MATRIX_FILE gr_30_30.mtx)

  IF(EXISTS ${MATRIX_FILE})
  ELSE()
    FILE(DOWNLOAD "ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/laplace/${MATRIX_FILE}.gz" 
                  ${CMAKE_BINARY_DIR}/${MATRIX_FILE}.gz)

    FIND_PROGRAM(GZIP_TOOL
               NAMES gzip
               PATHS /bin /usr/bin /usr/local/bin /opt/local/bin)

    EXECUTE_PROCESS( COMMAND ${GZIP_TOOL} -d ${MATRIX_FILE}.gz
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  ENDIF()

  ADD_EXECUTABLE(paralution_cg test/cg.cpp)
  TARGET_LINK_LIBRARIES(paralution_cg paralution)

  ADD_TEST(NAME testparalution COMMAND paralution_cg ${MATRIX_FILE} 4)

 add_executable(
  testlinearsolver
  test/test.h
  test/testlinearsolver.c
  )
target_link_libraries(testlinearsolver paralution_c ${LIBS})
add_test(testlinearsolver testlinearsolver)

ELSE()

 add_executable(
  testlinearsolver
  test/test.h
  test/testlinearsolver.c
  )
target_link_libraries(testlinearsolver ${LIBS})
add_test(testlinearsolver testlinearsolver)


ENDIF()
